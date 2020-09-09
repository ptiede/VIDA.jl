using Distributed
using ArgParse
Distributed.@everywhere using VIDA
using CSV
using DataFrames
using Random


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
        "--out"
            help = "name of output files with extractions"
            default = "fit_summaries.txt"
        "--restart"
            help = "Tells the sampler to read in the old files and restart the run."
            action = :store_true
        "-c"
            help = "lower percentage of flux to clip from image"
            arg_type = Float64
            default = 0.0
        "--div"
            help = "Divergence type to be used in image reconstruction"
            arg_type = String
            default = "Bh"
        "--filter"
            help = "Cosine filter to use for feature extraction. Expects two numbers"
            action = :store_arg
            nargs = 2
            required=true
        "--seed"
            help = "Random seed for initial positions in extract"
            arg_type = Int
            default = 42
        "--istart"
            help = "start index of fits file"
            arg_type = Int
            default = 1
        "--iend"
            help = "end index of fits file"
            arg_type = Int
            default = -1
    end
    return parse_args(s)
end

function main()
    #Parse the command line arguments
    parsed_args = parse_commandline()
    #Assign the arguments to specific variables
    fitsfiles = parsed_args["arg1"]
    out_name = parsed_args["out"]
    clip_percent = parsed_args["c"]
    startindx = parsed_args["istart"]
    endindx = parsed_args["iend"]
    seed = parsed_args["seed"]
    div_type = parsed_args["div"]
    stride = parsed_args["stride"]
    if (parsed_args["div"]=="KL")
      div_type = "KL"
    elseif (parsed_args["div"]=="Bh")
      div_type = "Bh"
    else
      error("$div_type not found! Must be Bh or KL")
    end
    filter_type = parse.(Int, parsed_args["filter"])
    println("Filter type $filter_type")
    restart = parsed_args["restart"]
    println("Using options: ")
    println("list of files: $fitsfiles, ")
    println("output name: $out_name, ")
    println("Clipping $(clip_percent*100) percent")
    println("random seed $seed")
    println("divergence type $div_type")
    println("Checkpoint stride $stride")

    #Read in a file and create list of images to filter
    #the last line is the termination of the file
    files = split(read(fitsfiles,String),"\n")

    #check if the last entry of files is an empty string
    if files[end] == ""
        files = files[1:end-1]
    end

    if startindx == 0
      startindx = 1
    end

    if endindx == -1 || endindx > length(files)
      endindx = length(files)
    end

    println("starting index $startindx")
    println("ending index $endindx")
    println("Starting fit")
    files = files[startindx:endindx]


    #Now run on the files for real
    main_sub(files, out_name,
             div_type,filter_type,
             clip_percent, seed,
             restart, stride)
    println("Done! Check $out_name for summary")
    return 0
end

function make_initial_filter(filter_type)
    @show filter_type
    lower_σ = [0.01, [ -2.0 for i in 1:filter_type[1]]... ]
    upper_σ = [15.0, [ 2.0 for i in 1:filter_type[1]]... ]
    lower_ξσ = Float64[]
    upper_ξσ = Float64[]
    if filter_type[1] > 0
        lower_ξσ = Float64[-π for i in 1:filter_type[1] ]
        upper_ξσ = Float64[ π for i in 1:filter_type[1] ]
    end
    lower_s = [0.001, [-0.99 for i in 2:filter_type[2]]...]
    upper_s = [0.999, [0.99 for i in 2:filter_type[2]]...]
    lower_ξs = [-π for i in 1:filter_type[2] ]
    upper_ξs = [π for i in 1:filter_type[2] ]

    lower = [5.0 ,
             lower_σ..., lower_ξσ...,
             0.001, 0.0,
             lower_s..., lower_ξs...,
             -60.0, -60.0
            ]
    upper = [ 30.0 ,
              upper_σ..., upper_ξσ...,
              0.999, π,
              upper_s..., upper_ξs...,
              60.0, 60.0
            ]
    filter = CosineRing{filter_type[1],filter_type[2]}(lower)
    return (filter, lower, upper)
end

function create_initial_df!(start_indx, fitsfiles, filter, restart, out_name)
    start_indx = 1
    df = DataFrame()
    nfiles = length(fitsfiles)
    if !restart
      #we want the keynames to match the model parameters
      key_names = fieldnames(typeof(filter))
      n = length(filter.θ1.σ)
      m = length(filter.θ1.s)
      key_names = [:r0,
                    [:σ for i in 1:n]...,
                    [:ξσ for i in 2:n]...,
                    :τ,
                    :ξτ,
                    [:s for i in 1:m]...,
                    [:ξs for i in 1:m]...,
                    :x0,
                    :y0,
                    :Irel
                ]
      for i in 1:length(key_names)
        insertcols!(df, ncol(df)+1, Symbol(key_names[i]) => zeros(nfiles); makeunique=true)
      end
      @show key_names
      #fill the data frame with some likely pertinent information
      setproperty!(df, :divmin, zeros(nfiles))
      setproperty!(df, :fitsfiles,  fitsfiles)
    else
      df = DataFrame!(CSV.File(out_name, delim=";"))
      start_indx = findfirst(isequal(0.0), df[:,1])
      println("Restarting run for $out_name at index $start_indx")
    end
    return df, start_indx
end

@everywhere function fit_func(filter,lower, upper, div_type, clip_percent)
    function (file)
        d_type = :none
        if (div_type == "KL")
            d_type = :KL
        elseif (div_type == "Bh")
            d_type = :Bh
        else
            error("$(div_type) not found! Must be KL or Bh")
        end
        println("Extracting $file using $d_type divergence")
        image = load_fits(string(file))
        cimage = VIDA.clipimage(clip_percent,image)
        div = VIDA.Bhattacharyya(cimage)
        if (d_type == :KL)
            div = VIDA.KullbackLeibler(cimage)
        end
        θ,divmin,_,_ = bbextract(div, filter, lower, upper;
                                 TraceMode=:silent, MaxFuncEvals=50*10^3)
        θ,divmin,_,_ = extract(div, θ, lower, upper)
        return VIDA.unpack(θ),divmin
    end
end

function main_sub(fitsfiles, out_name,
                  div_type, filter_type,
                  clip_percent, seed,
                  restart, stride)

    #"Define the filter I want to use and the var bounds"
    matom = make_initial_filter(filter_type)
    model = matom[1] + 1.0*Constant()
    lower = [matom[2]...,1e-6]
    upper = [matom[3]...,1]
    @show lower
    @show upper
    @show model

    @assert length(lower) == length(upper) "Bounds must have equal size"
    @assert VIDA.size(typeof(model)) == length(lower) "Number of filter params $(length(unpack(model))) != $(length(lower))"

    #Need to make sure all the procs know this information
    @everywhere model = $(model)
    @everywhere lower = $(lower)
    @everywhere upper = $(upper)
    @everywhere div_type = $(div_type)
    @everywhere clip_percent = $(clip_percent)

    #Set up the data frame to hold the optimizer output that
    #will be saved
    start_indx = 1
    df,start_indx = create_initial_df!(start_indx,fitsfiles, model, restart, out_name)
    rng = MersenneTwister(seed) #set the rng

    #Now fit the files!
    @everywhere fit = fit_func(model,lower,upper, div_type, clip_percent)
    indexpart = Iterators.partition(start_indx:length(fitsfiles), stride)
    for ii in indexpart
      results = pmap(fit, fitsfiles[ii])
      println(results)
      df[ii,1:length(lower)] = hcat(first.(results)...)'
      df[ii,length(lower)+1] = last.(results)
      df[ii,end] = fitsfiles[ii]
      #save the file
      println("Checkpointing $(ii)")
      CSV.write(out_name, df, delim=';')
    end

    return df
end




main()
