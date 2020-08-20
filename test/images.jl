using Test,VIDA
include("common.jl")


@testset "ReadWriteImage" begin
    img = load_ehtimfits("../example/elliptical_gaussian_rot-45.00m87Scale_seed_23_simobs_netcal_scanavg-z0.6-s100-t1-v100-l10-p40-e0.000.fits")
    xcent,ycent = centroid(img)
    @test xcent == -2.0442682583272647
    @test ycent == 1.2458998319861223
    θ = GaussianRing(r0,σ,x0,y0)
    img_fake = VIDA.make_ehtimage(θ, 64, xlim, ylim, img)
    save_ehtimfits(img, "tmp")
    rm("tmp")
end

@testset "ImageModifiers" begin
    img = VIDA.load_ehtimfits("../example/elliptical_gaussian_rot-45.00m87Scale_seed_23_simobs_netcal_scanavg-z0.6-s100-t1-v100-l10-p40-e0.000.fits")
    cimg = VIDA.clipimage(0.1, img)
    dcimg = VIDA.downsample(2, cimg)
    wdcimg = VIDA.window_image([-40,50], [-30,20], img)

end
