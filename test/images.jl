using Test,VIDA
include("common.jl")


@testset "ReadWriteImage" begin
    img = VIDA.load_image(joinpath(@__DIR__, "../example/data/example_image.fits"))
    xcent,ycent = centroid(img)
    @test xcent ≈ μas2rad(-1.7307380406643236)
    @test ycent ≈ μas2rad(0.5182951434801205)
    θ = GaussianRing(r0,σ,x0,y0)
    img_fake = intensitymap(θ, axiskeys(img))
    ComradeBase.save("tmp", img)
    rm("tmp")
end

@testset "ImageModifiers" begin
    img = VIDA.load_image(joinpath(@__DIR__, "../example/data/example_image.fits"))
    cimg = clipimage(0.1, img)
    cimg = clipimage(0.0, img, :absolute)
end

@testset "ImageModifiers Asymmetric" begin
    img = VIDA.load_image("../example/data/asymetric_image.fits")
    cimg = clipimage(0.1, img)
    cimg = clipimage(0.0, img, :absolute)
end

@testset "BlurImages" begin
    filt = GaussianRing(1.0)
    img = intensitymap(filt, 10.0, 10.0, 128, 128)
    bimg = VIDA.blur(img, σ*2*sqrt(2*log(2)))
    @test size(img) == (128,128)
end
