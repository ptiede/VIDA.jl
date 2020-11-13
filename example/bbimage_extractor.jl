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
        "-b"
            help = "Likelihood power L^b, so high b increases surface contrast"
            arg_type = Float64
            default = 1.0
        "--div"
            help = "Divergence type to be used in image reconstruction"
            arg_type = String
            default = "Bh"
        "--filter"
            help = "Filter to use for the image extraction.\nOptions are, Gen, TIDA, Slash, Ellip, Circ, AsymG"
            action = :append_arg
            nargs = '*'
            arg_type = String
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
    breg = parsed_args["b"]
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
    filter_type = vcat(parsed_args["filter"]...)
    println("Filter type $filter_type")
    restart = parsed_args["restart"]
    println("Using $(Threads.nthreads()) threads")
    println("Using options: ")
    println("list of files: $fitsfiles, ")
    println("output name: $out_name, ")
    println("Clipping $(clip_percent*100) percent")
    println("likelihood multiplication $breg")
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
             clip_percent, breg, seed,
             restart, stride)
    println("Done! Check $out_name for summary")
    return 0
end

function make_initial_filter(filter_type)
    if filter_type == "TIDA"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [5.0, 0.01, 1e-3, -0.999, -π, -50.0, -50.0]
      upper = [35.0, 20.0, 0.999, 0.999, π, 50.0, 50.0]
      filter = TIDAGaussianRing(20.0, 5.0,
                                0.5, 0.5,
                                0.0, 0.0, 0.0)
      return (filter, lower, upper)
    elseif filter_type == "Gen"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [5.0, 0.01, 1e-3, 0.0, 1e-3,-π, -50.0, -50.0]
      upper = [35.0, 20.0, 0.999, π,  0.99, π, 50.0, 50.0]
      filter = GeneralGaussianRing(20.0, 5.0,
                                   0.5, 0.0,
                                   0.5, 0.0,
                                   0.0, 0.0)
      return (filter, lower, upper)
    elseif filter_type == "Slash"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [5.0, 0.01, 1e-3, -π, -50.0, -50.0]
      upper = [35.0, 20.0, 0.999, π,  50.0, 50.0]
      filter = SlashedGaussianRing(20.0, 5.0,
                                   0.5, 0.0,
                                   0.0, 0.0)
      return (filter, lower, upper)
    elseif filter_type == "Ellip"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [5.0, 0.01, 1e-3, 0.0, -50.0, -50.0]
      upper = [35.0, 20.0, 0.999, π,  50.0, 50.0]
      filter = EllipticalGaussianRing(20.0, 5.0,
                                      0.5, 0.0,
                                      0.0, 0.0)
      return (filter, lower, upper)
    elseif filter_type == "Circ"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [5.0, 0.01, -40.0, -50.0]
      upper = [35.0, 20.0, 50.0, 50.0]
      filter = GaussianRing(20.0, 5.0,0.0, 0.0)
      return (filter, lower, upper)
    elseif filter_type == "AsymG"
        lower = [0.5, 0.001, 0.0, -60.0,-60.0]
        upper = [35.0, 0.999, π, 60.0, 60.0 ]
        filter = AsymGaussian(5.0,0.001, 0.001, 0.0,0.0)
        return (filter, lower, upper)
    elseif filter_type == "Disk"
        lower = [0.5 , 0.001, -60.0, -60.0]
        upper = [35.0, 20.0 ,  60.0,  60.0 ]
        filter = Disk(5.0, 1.0, 0.0, 0.0)
        return (filter, lower, upper)
    else
      error("$filter_type not found must be Circ, Ellip, Slash, TIDA, Gen, Disk, AsymG")
    end
end

function create_initial_df!(start_indx, fitsfiles, filter, restart, out_name)
    start_indx = 1
    df = DataFrame()
    nfiles = length(fitsfiles)
    if !restart
      #we want the keynames to match the model parameters
      key_names = fieldnames(typeof(filter))
      for i in 1:length(key_names)
        insertcols!(df, ncol(df)+1, Symbol(key_names[i]) => zeros(nfiles); makeunique=true)
      end
      #fill the data frame with some likely pertinent information
      setproperty!(df, :divmin, zeros(nfiles))
      setproperty!(df, :fitsfiles,  fitsfiles)
    else
      df = CSV.read(out_name, delim=";") |> DataFrame
      start_indx = findfirst(isequal(0.0), df[:,1])
      println("Restarting run for $out_name at index $start_indx")
    end
    return df, start_indx
end

@everywhere function fit_func(filter,lower, upper, div_type, clip_percent, breg)
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
        prob = ExtractProblem(div, filter, lower, upper)
        θ,divmin = extractor(prob, BBO(tracemode=:silent, maxevals=50_000))
        prob_new = ExtractProblem(div, θ, lower, upper)
        θ,divmin = extractor(prob_new, CMAES(cov_scale=0.01, verbosity=0))
        return VIDA.unpack(θ),divmin
    end
end

function main_sub(fitsfiles, out_name,
                  div_type, filter_type,
                  clip_percent, breg, seed,
                  restart, stride)

    #"Define the filter I want to use and the var bounds"
    matom = make_initial_filter.(filter_type)
    model = 1.0*matom[1][1]
    lower = [matom[1][2]...,1e-10]
    upper = [matom[1][3]...,1e2]
    if length(matom) > 1
      for i in 2:length(matom)
        model = model + 1.0*matom[i][1]
        lower = [lower..., matom[i][2]..., 1e-10]
        upper = [upper..., matom[i][3]..., 1e2]
      end
    end
    model = model + 1.0*Constant()
    lower = typeof(model)([lower..., 1e-10])
    upper = typeof(model)([upper..., 1])

    #Need to make sure all the procs know this information
    @everywhere model = $(model)
    @everywhere lower = $(lower)
    @everywhere upper = $(upper)
    @everywhere div_type = $(div_type)
    @everywhere clip_percent = $(clip_percent)
    @everywhere breg = $(breg)

    #Set up the data frame to hold the optimizer output that
    #will be saved
    start_indx = 1
    df,start_indx = create_initial_df!(start_indx,fitsfiles, model, restart, out_name)

    #Now fit the files!
    @everywhere fit = fit_func(model,lower,upper, div_type, clip_percent, breg)
    #Partition the list for checkpointing purposes
    indexpart = Iterators.partition(start_indx:length(fitsfiles), stride)
    for ii in indexpart
      results = pmap(fit, fitsfiles[ii])
      df[ii,1:VIDA.size(typeof(lower))] = hcat(first.(results)...)'
      df[ii,VIDA.size(typeof(lower))+1] = last.(results)
      df[ii,end] = fitsfiles[ii]
      #save the file
      println("Checkpointing $(ii)")
      CSV.write(out_name, df, delim=';')
    end

    return df
end




main()
