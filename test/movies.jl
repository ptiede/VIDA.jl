using Test,VIDA
include("common.jl")

@testset "MovieRead" begin
    times = 1:0.1:2
    phi = times*phi
    images = EHTImage[]
    for p in phi
        im = SlashedGaussianRing(r0=r0, σ=σ, s=s, ξs = p, x0=0.0, y0=0.0)
        push!(images, im)
    end
    mov = join_frames(times, images)
    #Test the save and read hdf5
    save_hdf5("test.jl", mov)
    mov_read = load_hdf5("test.jl")
    @test mov_read.frames.coef .== mov.frames.coef
    # Kill the temp file
    rm("test.hdf5")

    #Get an image and make sure everything matches
    img = get_image(times[1], img)
    frames = get_frames(mov)
    @test frames[1].img .== images[1].img
    @test img.img .== images[1].img

end
