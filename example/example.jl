using ArgParse
using LaVIDA
using CSV
using DataFrames
using Optim
using Random


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
        "-c"
            help = "lower percentage of flux to clip from image"
            arg_type = Float64
            default = 0.0
        "-b"
            help = "Likelihood power L^b, so high b increases surface contrast"
            arg_type = Float64
            default = 1e6
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
    nstart = parsed_args["ninits"]
    out_name = parsed_args["out"]
    plotbool = parsed_args["plot"]
    clip_percent = parsed_args["c"]
    breg = parsed_args["b"]
    seed = parsed_args["seed"]
    println("Using options: ")
    println("list of files: $fitsfiles, ")
    println("nstart: $nstart, ")
    println("output name: $out_name, ")
    println("plotting: $plotbool")
    println("Clipping $(clip_percent*100) percent")
    println("likelihood multiplication $breg")
    println("random seed $seed")
    println("Starting fit")

    #Read in a file and create list of images to filter
    #the last line is the termination of the file
    files = string.(Base.split(read(fitsfiles,String),"\n"))

    #check if the last entry of files is an empty string
    if files[end] == ""
        files = files[1:end-1]
    end


    #Now run on the files for real
    main_sub(files, nstart, out_name,
             clip_percent, breg, seed,
             plotbool)
    println("Done! Check $out_name for summary")
    return 0
end

function main_sub(fitsfiles, nstart, out_name,
                  clip_percent, breg, seed, plotbool)
   "Define the filter I want to use and the starting values"
   filter = TIDAGaussianRing(20.0, 5.0,
                              0.5, 0.5,
                              0.0, 0.0, 0.0)
    #parameter bounds, be aggressive with these. If they
    #are too small the optimizer can struggle
    lower = [10.0, 0.01, 1e-3, -0.999, -π, -50.0, -50.0]
    upper = [35.0, 20.0, 0.999, 0.999, π, 50.0, 50.0]

    #Set up the data frame to hold the optimizer output that
    #will be saved
    df = DataFrame()
    nfiles = length(fitsfiles)
    #we want the keynames to match the model parameters
    key_names = fieldnames(typeof(filter))
    for i in 1:length(key_names)
        setproperty!(df, key_names[i], zeros(2*nfiles))
    end
    #fill the data frame with some likely pertinent information
    setproperty!(df, :ℓmax, zeros(2*nfiles))
    setproperty!(df, :threadid, zeros(2*nfiles))
    setproperty!(df, :converged, Vector{Bool}(undef, 2*nfiles))
    setproperty!(df, :iterations, zeros(2*nfiles))
    setproperty!(df, :fitsfiles,  fill("sdf",2*nfiles))


    rng = MersenneTwister(seed) #set the rng
    #Now fits the files!
    for (i,f) in enumerate(fitsfiles)
        println("Fitting image $f")
        #load the image and clip it to remove low cost contaminants
        image = load_ehtimfits(string(f))
        image = clipimage(clip_percent,image,:relative)

        #Create the measure and fitting function using a closure
        bh = make_div(image, :Bh)
        kl = make_div(image, :KL)

        #Extract! and fill the dataframe
        results = extract(rng, nstart, bh, filter, lower, upper,
                          SAMIN(),
                          Optim.Options(iterations = 10^6))
        df[i,1:end-1] = results[1,:]
        df[i,end] = f

        #If true we will plot the triptic
        if plotbool
            pname = replace(f, ".fits"=>"_bh_triptic.png")
            θfit = TIDAGaussianRing(results[1,1:length(key_names)])
            triptic(image, θfit)
            savefig(pname)
            closeall()
        end


        results = extract(rng, nstart, kl, filter, lower, upper,
                          Fminbox(Optim.LBFGS()),
                          Optim.Options(iterations = 20))
        df[i+1,1:end-1] = results[1,:]
        df[i+1,end] = f

        #If true we will plot the triptic
        if plotbool
            pname = replace(f, ".fits"=>"_kl_triptic.png")
            θfit = TIDAGaussianRing(results[1,1:length(key_names)])
            triptic(image, θfit)
            savefig(pname)
            closeall()
        end
    end
    df[!,:ℓmax] /=  breg #divide out breg to give real max

    #save the file
    CSV.write(out_name, df, delim=';')

    return df
end

main()
