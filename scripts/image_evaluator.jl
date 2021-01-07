using ArgParse
using CSV
using DataFrames
using Random
using VIDA


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "arg1"
            help = "file list of fits images to read"
            arg_type = String
            required = true
        "--ntries"
             help = "Number of instances of optimizer to run."
             arg_type = Int
             default = 8
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
        "--truth"
            help = "File name for the ground truth hdf5 movie"
            arg_type = String
            required=true
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
    recon_file = parsed_args["arg1"]
    out_name = parsed_args["out"]
    clip_percent = parsed_args["c"]
    seed = parsed_args["seed"]
    div_type = parsed_args["div"]
    ntries = parsed_args["ntries"]
    if (parsed_args["div"]=="KL")
      div_type = "KL"
    elseif (parsed_args["div"]=="Bh")
      div_type = "Bh"
    else
      error("$div_type not found! Must be Bh or KL")
    end
    truth_file = parsed_args["truth"]
    println("Truth image filename $truth_file")
    restart = parsed_args["restart"]
    println("Using $(Threads.nthreads()) threads")
    println("Using options: ")
    println("Recon movie: $recon_file, ")
    println("output name: $out_name, ")
    println("Clipping $(clip_percent*100) percent")
    println("random seed $seed")
    println("divergence type $div_type")
    println("Number of tries $ntries")

    #Read recon_movie
    mov_recon = load_hdf5(recon_file)

    #Read in truth movie
    mov_truth = load_hdf5(truth_file)


    #Now run on the files for real
    main_sub(mov_recon, out_name,
             div_type, mov_truth,
             clip_percent,
             restart, ntries)
    println("Done! Check $out_name for summary")
    return 0
end

function make_problem(div, truthimg::EHTImage)
    lower = ImageFilter(-80.0, -80.0, truthimg)
    upper = ImageFilter(80.0, 80.0, truthimg)
    filter = ImageFilter(0.0, 0.0, truthimg)
    return ExtractProblem(div, filter, lower, upper)
end

function create_initial_df!(start_indx, times, restart, out_name)
    start_indx = 1
    df = DataFrame()
    nfiles = length(times)
    if !restart
      #we want the keynames to match the model parameters
      key_names = fieldnames(ImageFilter)
      for i in 1:length(key_names)
        insertcols!(df, ncol(df)+1, Symbol(key_names[i]) => zeros(nfiles); makeunique=true)
      end
      #fill the data frame with some likely pertinent information
      setproperty!(df, :divmin, zeros(nfiles))
      setproperty!(df, :time,  times)
    else
      df = CSV.File(out_name) |> DataFrame
      start_indx = findfirst(isequal(0.0), df[:,1])
      println("Restarting run for $out_name at index $start_indx")
    end
    return df, start_indx
end


function main_sub(mov_recon, out_name,
                  div_type, mov_truth,
                  clip_percent,
                  restart, ntries)

    #Get reconstruction times and images
    times = get_times(mov_recon)
    images = VIDA.clipimage.(clip_percent, get_frames(mov_recon))
    println("I am going to evaluate $(length(times)) images")

    # Create initial data frames
    start_indx = 1
    df,start_indx = create_initial_df!(start_indx, times, restart, out_name)

    for (i,t) in enumerate(times)
        imtruth = get_image(mov_truth, t)
        bh = Bhattacharyya(images[i])
        prob = make_problem(bh, imtruth)
        res = threaded_extractor(ntries, prob, CMAES(cov_scale=20.0, verbosity=0))
        df[i,1:2] = unpack(res[1])
        df[i,3]   = res[2]
        println("Done image $i/$(length(times))")
    end

    CSV.write(out_name, df)

    return df
end

main()
