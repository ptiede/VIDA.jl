using Distributed
using ArgParse
Distributed.@everywhere using VIDA
using CSV
using DataFrames
using PyPlot
using Optim
using Random
using SharedArrays


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "arg1"
            help = "file list of fits images to read"
            arg_type = String
            required = true
        "--ninits", "-n"
            help = "number of initial starting locations"
            arg_type = Int
            default = 10
        "--out"
            help = "name of output files with extractions"
            default = "fit_summaries.txt"
        "--plot"
            help = "Switch that turns on plotting triptics as fitting"
            action = :store_true
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
            help = "Filter to use for the image extraction.\nOptions are, Gen, TIDA, Slash, Ellip, Circ"
            arg_type = String
            default = "TIDA"
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
    nstart = parsed_args["ninits"]
    out_name = parsed_args["out"]
    plotbool = parsed_args["plot"]
    clip_percent = parsed_args["c"]
    startindx = parsed_args["istart"]
    endindx = parsed_args["iend"]
    breg = parsed_args["b"]
    seed = parsed_args["seed"]
    div_type = parsed_args["div"]
    if (parsed_args["div"]=="KL")
      div_type = "KL"
    elseif (parsed_args["div"]=="Bh")
      div_type = "Bh"
    else
      error("$div_type not found! Must be Bh or KL")
    end
    filter_type = parsed_args["filter"]
    restart = parsed_args["restart"]
    println("Using $(Threads.nthreads()) threads")
    println("Using options: ")
    println("list of files: $fitsfiles, ")
    println("nstart: $nstart, ")
    println("output name: $out_name, ")
    println("plotting: $plotbool")
    println("Clipping $(clip_percent*100) percent")
    println("likelihood multiplication $breg")
    println("random seed $seed")
    println("divergence type $div_type")

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
    main_sub(files, nstart, out_name,
             div_type,filter_type,
             clip_percent, breg, seed,
             plotbool, restart)
    println("Done! Check $out_name for summary")
    return 0
end

function make_initial_filter(filter_type)
    if filter_type == "TIDA"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [10.0, 0.01, 1e-3, -0.999, -π, -30.0, -30.0, 1e-6]
      upper = [35.0, 15.0, 0.999, 0.999, π, 30.0, 30.0, 1e3]
      filter = TIDAGaussianRing(20.0, 5.0,
                                0.5, 0.5,
                                0.0, 0.0, 0.0) + 1.0*Constant()
      return (filter, lower, upper)
    elseif filter_type == "Gen"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [10.0, 0.01, 1e-3, 0.0, 1e-3,-π, -30.0, -30.0, 1e-6]
      upper = [35.0, 15.0, 0.999, π,  0.99, π, 30.0, 30.0, 1e3]
      filter = GeneralGaussianRing(20.0, 5.0,
                                   0.5, 0.0,
                                   0.5, 0.0,
                                   0.0, 0.0) + 1.0*Constant()
      return (filter, lower, upper)
    elseif filter_type == "Slash"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [10.0, 0.01, 1e-3, -π, -30.0, -30.0, 1e-6]
      upper = [35.0, 15.0, 0.999, π,  30.0, 30.0, 1e3]
      filter = SlashedGaussianRing(20.0, 5.0,
                                   0.5, 0.0,
                                   0.0, 0.0) + 1.0*Constant()
      return (filter, lower, upper)
    elseif filter_type == "Ellip"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [10.0, 0.01, 1e-3, 0.0, -30.0, -30.0, 1e-6]
      upper = [35.0, 15.0, 0.999, π,  30.0, 30.0, 1e3]
      filter = EllipticalGaussianRing(20.0, 5.0,
                                      0.5, 0.0,
                                      0.0, 0.0) + 1.0*Constant()
      return (filter, lower, upper)
    elseif filter_type == "Circ"
      #parameter bounds, be aggressive with these. If they
      #are too small the optimizer can struggle
      lower = [10.0, 0.01, -30.0, -30.0, 1e-6]
      upper = [35.0, 15.0, 30.0, 30.0, 1e3]
      filter = GaussianRing(20.0, 5.0,0.0, 0.0) + 1.0*Constant()
      return (filter+1.0*Constant(), lower, upper)
    else
      error("$filter_type not found must be Circ, Ellip, Slash, TIDA, Gen")
    end
end

function create_initial_df!(start_indx, fitsfiles, filter, restart)
    if !restart
      start_indx = 1
      df = DataFrame()
      nfiles = length(fitsfiles)
      #we want the keynames to match the model parameters
      key_names = fieldnames(typeof(filter))
      for i in 1:length(key_names)
          setproperty!(df, key_names[i], zeros(nfiles))
      end
      #fill the data frame with some likely pertinent information
      setproperty!(df, :divmin, zeros(nfiles))
      setproperty!(df, :fitsfiles,  fitsfiles)
    else
      df = CSV.read(out_name, delim=";")
      start_indx = findfirst(isequal(0.0), df[!,:r0])
      println("Restarting run for $out_name at index $start_indx")
    end
    return df
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
        image = VIDA.load_ehtimfits(string(file))
        cimage = VIDA.clipimage(clip_percent,image)
        div = VIDA.make_div(cimage, d_type, breg)
        θ,divmin,_,_ = bbextract(div, filter, lower, upper;
                                 TraceMode=:silent, MaxFuncEvals=2*10^4)
        return VIDA.unpack(θ),divmin
    end
end

function main_sub(fitsfiles, nstart, out_name,
                  div_type, filter_type,
                  clip_percent, breg, seed,
                  plotbool, restart)

    #"Define the filter I want to use and the var bounds"
    model,lower,upper = make_initial_filter(filter_type)
    @everywhere model = $(model)
    @everywhere lower = $(lower)
    @everywhere upper = $(upper)
    @everywhere div_type = $(div_type)
    @everywhere clip_percent = $(clip_percent)
    @everywhere breg = $(breg)
    #prixt(d_type)
    #Set up the data frame to hold the optimizer output that
    #will be saved
    start_indx = 1
    df = create_initial_df!(start_indx,fitsfiles, model, restart)

    rng = MersenneTwister(seed) #set the rng
    #Now fit the files!
    @everywhere fit = fit_func(model,lower,upper, div_type, clip_percent, breg)
    results = pmap(fit, fitsfiles)
    println(first.(results))
    df[:,1:length(lower)] = hcat(first.(results)...)'
    df[:,length(lower)+1] = last.(results)
    df[:,end] = fitsfiles
    #save the file
    CSV.write(out_name, df, delim=';')

    return df
end




main()
