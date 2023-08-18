using Test,VIDA
include("common.jl")


@testset "ReadWriteImage" begin
    img = VIDA.load_fits("../example/data/example_image.fits")
    _ = VIDA.load_image("../example/data/example_image.fits")
    xcent,ycent = centroid(img)
    @test xcent == -1.7307380406643236
    @test ycent == 0.5182951434801205
    θ = GaussianRing(r0,σ,x0,y0)
    img_fake = VIDA.make_image(θ, 64, xlim, ylim, img)
    save_fits(img, "tmp")
    rm("tmp")
end

@testset "ImageModifiers" begin
    img = VIDA.load_fits("../example/data/example_image.fits")
    cimg = clipimage(0.1, img)
    cimg = clipimage(0.0, img, :absolute)
    regrid(img, npix, xlim, ylim)
    ra,dec = pixelloc(img)
end

@testset "ImageModifiers Asymmetric" begin
    img = VIDA.load_fits("../example/data/asymetric_image.fits")
    cimg = clipimage(0.1, img)
    cimg = clipimage(0.0, img, :absolute)
    rescale(img, npix, xlim, ylim)
    ra,dec = pixelloc(img)
end

@testset "BlurImages" begin
    filt = AsymGaussian(σ=σ, τ=0.0, ξ=0.0, x0=0.0, y0=0.0)
    img = VIDA.make_image(filt, 1024, [-100.0, 100.0], [-100.0, 100.0])
    bimg = blur(img, σ*2*sqrt(2*log(2)))
    s2 = sqrt(inertia(bimg)[1])
    @test isapprox(s2,sqrt(2)*σ, rtol=1e-3)
    @test size(img) == (1024,1024)
end
