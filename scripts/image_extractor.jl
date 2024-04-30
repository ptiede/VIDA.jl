using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using Distributed
@everywhere begin
    using Pkg; Pkg.activate(@__DIR__)
end

@everywhere using VIDA

using ArgParse
using CSV
using DataFrames
using Random
@everywhere begin
    using OptimizationBBO
    using OptimizationMetaheuristics
end



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "arg1"
            help = "file list of fits images to read"
            arg_type = String
            required = true
        "--stride"
             help = "Checkpointing stride, i.e. number of steps."
             arg_type = Int
             default = 200
        "--blur"
             help = "Blur images before extracting fitted model parameters"
             arg_type = Float64
             default = 0.0
        "--out"
            help = "name of output files with extractions"
            default = "fit_summaries.csv"
        "--restart"
            help = "Tells the sampler to read in the old files and restart the run."
            action = :store_true
        "--template"
            help = "Parses sting with models separated by spaces. For instance to\n"*
                   "run a model with 1 m=(1,4) m-ring, 2 gaussian and a stretched disk\n"*
                   "one could do `--template mring_1_4 gauss_2 disk_1`. The current model\n"*
                   "options are: \n"*
                   "  - `[stretch]mring_n_m`: adds a [stretched] m-ring of order `n` `m` thickenss and azimuth resp.\n"*
                   "  - `gauss_n`           : add `n` asymmetric gaussian components to the template\n"*
                   "  - `[stretch]disk_n    : add `n` [stretched] flat top disks to the template\n"
            action = :store_arg
            nargs = '*'
            arg_type = String
            required = true
        "--seed"
            help = "Random seed for initial positions in extract"
            arg_type = Int
            default = 42
    end
    return parse_args(s)
end

function main()
    #Parse the command line arguments
    parsed_args = parse_commandline()
    #Assign the arguments to specific variables
    fitsfiles = parsed_args["arg1"]
    out_name = parsed_args["out"]
    seed = parsed_args["seed"]
    stride = parsed_args["stride"]
    blur = μas2rad(parsed_args["blur"])
    templates = parsed_args["template"]
    @info "Template types $templates"
    restart = parsed_args["restart"]
    println("Using options: ")
    println("list of files: $fitsfiles, ")
    println("output name: $out_name, ")
    println("random seed $seed")
    println("Checkpoint stride $stride")
    println("Blurring Gaussian kernel width in µas: $(rad2μas(blur))")

    #Read in a file and create list of images to template
    #the last line is the termination of the file
    files = split(read(fitsfiles,String),"\n")

    #check if the last entry of files is an empty string
    if files[end] == ""
        files = files[1:end-1]
    end

    println("Starting fit")


    #Now run on the files for real
    main_sub(files, out_name,
             templates,
             seed,
             restart, stride, blur)
    println("Done! Check $out_name for summary")
    return 0
end



function make_initial_templates(templates...)
    res = map(make_initial_template, templates)
    syms = ntuple(i->Symbol(:model_, i), length(res))
    templates = NamedTuple{syms}(getindex.(res, 1))
    lowers    = NamedTuple{syms}(getindex.(res, 2))
    uppers    = NamedTuple{syms}(getindex.(res, 3))

    templates_const = merge(templates, (c = x->(x.floor*Constant(μas2rad(150.0))),))
    lowers_const = merge(lowers, (c = (floor = 1e-6,),))
    uppers_const = merge(uppers, (c = (floor = 100.0,),))

    syms_const = (syms..., :c)

    temp = let t = templates_const, s=syms_const
        x->begin
            sum(s) do si
                getproperty(t, si)(getproperty(x, si))
            end
        end
    end
    return temp, lowers_const, uppers_const
end



function make_initial_template(template)
    if occursin("mring", template)
        stretch = occursin("stretch", template)
        type = parse.(Int, split(template, "_")[2:3])
        @info "Template includes a m-ring of order $(type) with stretch $stretch"
        return make_initial_template_mring(type[1], type[2], stretch)
    elseif occursin("gauss", template)
        ngauss = parse.(Int, split(template, "_")[end])
        @info "Template includes $ngauss gaussians"
        return make_initial_template_gauss(ngauss)
    elseif occursin("disk", template)
        ndisk = parse.(Int, split(template, "_")[end])
        stretch = occursin("stretch", template)
        @info "Template includes $ndisk disks with stretch $stretch"
        return make_initial_template_disk(ndisk, stretch)
    else
        @error "Template $template not available"
    end
end






function make_initial_template_gauss(n::Int)
    gauss = @inline function (θ)
        mapreduce(+, 1:n) do i
            modify(Gaussian(),
                Stretch(θ.σ[i]*sqrt(1-θ.τ[i]), θ.σ[i]/sqrt(1-θ.τ[i])),
                Rotate(θ.ξ[i]),
                Shift(θ.x0[i], θ.y0[i])
            )
        end
    end
    lower = (σ=ntuple(_-> μas2rad(0.5), n),
             τ=ntuple(_->0.001, n),
             ξ=ntuple(_->-π/2, n),
             x0=ntuple(_->-μas2rad(60.0), n),
             y0=ntuple(_->-μas2rad(60.0), n),
             )
    upper = (σ=ntuple(_-> μas2rad(35.0), n),
             τ=ntuple(_->0.5, n),
             ξ=ntuple(_->π/2, n),
             x0=ntuple(_->μas2rad(60.0), n),
             y0=ntuple(_->μas2rad(60.0), n),
             )

    return gauss, lower, upper
end


function make_initial_template_disk(n::Int, stretch)
    if stretch
        disk = @inline function (θ)
            mapreduce(+, 1:n) do i
                mod = sqrt(1-θ.τ[i])
                modify(GaussDisk(θ.σ[i]/θ.r0[i]),
                        Stretch(θ.r0[i]*mod, θ.r0[i]/mod),
                          Rotate(θ.ξ[i]),
                          Shift(θ.x0[i], θ.y0[i])
                         )
            end
        end
        lower = (
                 r0=ntuple(_-> μas2rad(0.5), n),
                 σ =ntuple(_->μas2rad(0.1), n),
                 τ=ntuple(_->0.001, n),
                 ξ=ntuple(_->-π/2, n),
                 x0=ntuple(_->-μas2rad(60.0), n),
                 y0=ntuple(_->-μas2rad(60.0), n),
                )
        upper = (
                 r0=ntuple(_->μas2rad(35.0), n),
                 σ =ntuple(_->μas2rad(20.0), n),
                 τ=ntuple(_->0.5, n),
                 ξ=ntuple(_->π/2, n),
                 x0=ntuple(_->μas2rad(60.0), n),
                 y0=ntuple(_->μas2rad(60.0), n),
                )
    else
        disk = @inline function (θ)
            mapreduce(+, 1:n) do i
                modify(GaussDisk(θ.σ[i]/θ.r0[i]),
                          Stretch(θ.r0[i], θ.r0[i]),
                          Shift(θ.x0[i], θ.y0[i])
                         )
            end
        end
        lower = (
                 r0=ntuple(_-> μas2rad(0.5), n),
                 σ =ntuple(_->μas2rad(0.1), n),
                 x0=ntuple(_->-μas2rad(60.0), n),
                 y0=ntuple(_->-μas2rad(60.0), n),
                )
        upper = (
                 r0=ntuple(_-> μas2rad(35.0), n),
                 σ =ntuple(_->μas2rad(20.0), n),
                 x0=ntuple(_->μas2rad(60.0), n),
                 y0=ntuple(_->μas2rad(60.0), n),
                )
    end
    return disk, lower, upper
end


function make_initial_template_mring(N::Int, M::Int, stretch)

    lower = (
             r0 = μas2rad(10.0),
             σ0 = μas2rad(1.0),
             σ  = ntuple(_->μas2rad(0.0), N),
             ξσ = ntuple(_->-1π, N),
             s  = ntuple(_->0.0, M),
             ξs = ntuple(_->-1π, M),
             x0 = -μas2rad(60.0),
             y0 = -μas2rad(60.0)
             )

    upper = (
             r0 = μas2rad(35.0),
             σ0 = μas2rad(20.0),
             σ  = ntuple(_->μas2rad(5.0), N),
             ξσ = ntuple(_->1π, N),
             s  = ntuple(_->0.999, M),
             ξs = ntuple(_->1π, M),
             x0 = μas2rad(60.0),
             y0 = μas2rad(60.0)
             )


    if stretch
        mring = @inline function (θ)
            mod = sqrt(1-θ.τ)
            return modify(CosineRing(θ.σ0/θ.r0, θ.σ./θ.r0, θ.ξσ .- θ.ξτ, θ.s, θ.ξs .- θ.ξτ),
                            Stretch(θ.r0*mod, θ.r0/mod),
                            Rotate(θ.ξτ),
                            Shift(θ.x0, θ.y0)
                         )
        end
        lower = merge(lower, (τ = 0.001, ξτ = -π/2))
        upper = merge(upper, (τ = 0.5, ξτ = π/2))
    else
        mring = @inline function (θ)
            return modify(CosineRing(θ.σ0/θ.r0, θ.σ./θ.r0, θ.ξσ, θ.s, θ.ξs),
                            Stretch(θ.r0),
                            Shift(θ.x0, θ.y0)
                         )
        end
    end
    return mring, lower, upper
end

flattenvals(x::NamedTuple) = flattenvals(values(x)...)
flattenvals(x, y...) = (flattenvals(x)..., flattenvals(y...)...)
flattenvals(x::Tuple) = flattenvals(x...)
flattenvals(x) = (x,)
flattenvals() = ()

function flattenkeys(nt::NamedTuple{N, T}) where {N, T<:Tuple}
    k = keys(nt)
    v = values(nt)
    t  = ntuple(i->combine_symbol(k[i], v[i]), length(k))
    flattenvals(t)
end

combine_symbol(a::Symbol, ::Number) = a
combine_symbol(a::Symbol, b::Symbol) = Symbol(a, "_", b)
combine_symbol(a::Symbol, ::NTuple{N, <:Real}) where {N} = ntuple(i->Symbol(a, "_", i), N)
function combine_symbol(a::Symbol, b::NamedTuple{N, T}) where {N, T}
    k = keys(b)
    v = values(b)
    ntuple(i->combine_symbol(Symbol(a, "_", k[i]), v[i]), length(k))
end


function create_initial_df(fitsfiles, names, restart, outname)
    start_indx = 1
    nfiles = length(fitsfiles)
    df = DataFrame()
    if !restart
      #we want the keynames to match the model parameters
      for i in eachindex(names)
        insertcols!(df, ncol(df)+1, names[i] => zeros(nfiles))
      end
      #fill the data frame with some likely pertinent information
      df[:, :divmin] = zeros(nfiles)
      df[:, :fitsfiles] =  fitsfiles
    else
      df = DataFrame(CSV.File(outname))
      start_indx = findfirst(isequal(0.0), df[:,1])
      println("Restarting run for $outname at index $start_indx")
    end
    return df, start_indx
end

@everywhere function fit_template(file, template, lower, upper, blur)
    println("Extracting $file")
    image = load_image(string(file))
    if blur > 0.0
        image = VIDA.blur(image, blur) # INI: blur the image with Gaussian kernel of given fwhm
    end
    rimage = VIDA.regrid(image, μas2rad(200.0), μas2rad(200.0), 64, 64)
    cimage = VIDA.clipimage(0.0,rimage)
    div = VIDA.LeastSquares(cimage)
    t = @elapsed begin
        prob = VIDAProblem(div, template, lower, upper)
        xopt,  _, _ = vida(prob, BBO_adaptive_de_rand_1_bin(); maxiters=20_000)
        # xopt2, θ, divmin = vida(prob, CMAEvolutionStrategyOpt(); init_params=xopt, maxiters=1_000)
        xopt2, θ, divmin = vida(prob, ECA(;options=Options(f_calls_limit = 1000, f_tol = 1e-5)); init_params=xopt, use_initial=true)
    end
    println("This took $t seconds")
    println("The minimum divergence is $divmin")
    return xopt2, θ, divmin
end


function main_sub(fitsfiles, out_name,
                  templates,
                  seed,
                  restart, stride, blur)

    #"Define the template I want to use and the var bounds"
    model, lower, upper = make_initial_templates(templates...)

    # get the flat keys that we will save everything to
    k = flattenkeys(lower)

    #Set up the data frame to hold the optimizer output that
    #will be saved
    df, start_indx = create_initial_df(fitsfiles, k, restart, out_name)

    #Now fit the files!
    indexpart = Iterators.partition(start_indx:length(fitsfiles), stride)
    for ii in indexpart
      results = pmap(fitsfiles[ii]) do f
                fit_template(f, model, lower, upper, blur)
      end
      df[ii,1:length(k)] .= reduce(hcat, collect.(flattenvals.(first.(results))))'
      df[ii,end-1] = last.(results)
      df[ii,end] = fitsfiles[ii]
      #save the file
      println("Checkpointing $(ii)")
      CSV.write(out_name, df)
    end
    CSV.write(out_name, df)

    return df
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
