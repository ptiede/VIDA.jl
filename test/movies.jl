using Test,VIDA
include("common.jl")

@testset "MovieRead" begin
    times = collect(1.0:0.1:2.0)
    phi = times*π
    images = EHTImage{Float64, Array{Float64,2}}[]
    for p in phi
        f = SlashedGaussianRing(r0=r0, σ=σ, s=s, ξ = p, x0=0.0, y0=0.0)
        im = VIDA.make_image(f, 64, [-60.0, 60.0], [-60.0, 60.0])
        push!(images, im)
    end
    mov = join_frames(times, images)
    #Test the save and read hdf5
    save_hdf5("test.hdf5", mov)
    mov_read = load_hdf5("test.hdf5")
    @test mov_read.frames.itp.coefs ≈ mov.frames.itp.coefs
    # Kill the temp file
    rm("test.hdf5")

    #Get an image and make sure everything matches
    img = get_image(mov_read, times[1])
    frames = get_frames(mov)
    #@inferred get_frames(mov)
    @inferred getindex(frames, 1)
    @test typeof(frames[1]) === typeof(images[1])
    @test frames[1].img == images[1].img
    @test img.img == images[1].img

    flux(mov, 10.0)

    #Now lets blur the movie
    bmov = blur(mov, 10.0)
    #regrid the movie
    rmov = regrid(mov, 128, [-50.0, 50.0], [-50.0, 50.0])

end
