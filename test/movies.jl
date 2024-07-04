using Test,VIDA
include("common.jl")

@testset "MovieRead" begin
    times = collect(1.0:0.1:2.0)
    phi = times*π
    frames = map(phi) do p
        f = SlashedGaussianRing(r0, σ, s, p, 0.0, 0.0)
        return intensitymap(f, imagepixels(120.0, 120.0, 64, 64))
    end
    mov = join_frames(times, frames)
    VIDAMovie(times, frames)
    @test times == get_times(mov)
    show(mov)
    println(mov)
    #Test the save and read hdf5
    save_hdf5("test.hdf5", mov)
    mov_read = load_hdf5("test.hdf5")
    @test mov_read.frames ≈ mov.frames
    # Kill the temp file
    rm("test.hdf5")

    #Get an image and make sure everything matches
    img = get_image(mov_read, times[1])
    images = get_frames(mov)
    #@inferred get_frames(mov)
    @inferred getindex(frames, 1)
    # @test typeof(images[T=1]) === typeof(frames[1])
    @test frames[1] == images[T=1]
    @test img == images[T=1]

    flux(mov, 10.0)

    #Now lets blur the movie
    bmov = VIDA.blur(mov, 10.0)

    regrid(mov, imagepixels(100.0, 100.0, 32, 32))

end
