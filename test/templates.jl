using Test,VIDA
using LinearAlgebra
using Plots
include(joinpath(@__DIR__, "common.jl"))

function test_template(θ::ComradeBase.AbstractModel; npix=128, fov=120.0, )
    @inferred ComradeBase.intensity_point(θ, (X=0.0, Y=0.0))
    @test ComradeBase.intensity_point(θ, (X=0.0, Y=0.0)) isa AbstractFloat
    @inferred ComradeBase.radialextent(θ)
    @test ComradeBase.radialextent(θ) isa AbstractFloat

    img = intensitymap(θ, fov, fov, npix, npix)
    triptic(img, θ)
end


@testset "GaussianRing" begin
    t = RingTemplate(RadialGaussian(0.1), AzimuthalUniform())
    test_template(t)
    test_template(GaussianRing(0.1))
    test_template(GaussianRing(5.0, 0.1, 0.1, 0.1))
    g = imagepixels(fovx, fovy, npix, npix)

    @test intensitymap(modify(GaussianRing(0.1/5), Stretch(5.0), Shift(0.1, 0.1)), g) ≈
          intensitymap(GaussianRing(5.0, 0.1, 0.1, 0.1), g)
end

@testset "EllipticalGaussianRing" begin
    t = modify(RingTemplate(RadialGaussian(0.1), AzimuthalUniform()), Stretch(0.5, 2.0))
    test_template(t)
    test_template(EllipticalGaussianRing(5.0, 0.1, 0.5, 0.0, 0.1, 0.1))
    g = imagepixels(fovx, fovy, npix, npix)
    # @test intensitymap(modify(GaussianRing(0.1/5), Stretch(5.0), Shift(0.1, 0.1)), g) ==
    #       intensitymap(GaussianRing(5.0, 0.1, 0.1, 0.1), g)
end

@testset "RingTemplate" begin
    gr = RadialGaussian(0.1)
    dr = RadialDblPower(3.00, 5.0)
    tr = RadialTruncExp(1.0)
    rads = (gr, dr, tr)

    au = AzimuthalUniform()
    ac1 = AzimuthalCosine(0.5, π/2)
    ac2 = AzimuthalCosine((0.5, 0.1), (0.0, π/3))
    azs = (au, ac1, ac2)

    for r in rads, a in azs
        test_template(RingTemplate(gr, a) + 0.1*Constant(1.0))
    end

end

@testset "CosineRing" begin
    t = CosineRing(0.1, (0.1, 0.2, 0.3), (0.0, 0.25, 1.5), (0.5, 0.1, 0.2), (0.0, -0.25, -2.0))
    test_template(t)
    t1 = CosineRing(0.1, (), (), (0.5,), (0.0,))
    test_template(t1)
    t2 = SlashedGaussianRing(0.1, 0.5)
    t3 = CosineRing(0.1, (0.1,), (0.0,), (), (), )
    test_template(t3)

    g = imagepixels(10.0, 10.0, 64, 64)
    img1 = intensitymap(t1, g)
    img2 = intensitymap(t2, g)
    @test img1 ≈ img2

    tr = CosineRing(5.0, 1.0, (0.1, 0.2, 0.3), (0.0, 0.25, 1.5), (0.5, 0.1, 0.2), (0.0, -0.25, -2.0), 1.0, 2.0)
    plot(tr)
end

@testset "CosineRing With floors" begin
    tr1 = CosineRingwFloor(5.0, 1.0, (0.1, 0.2, 0.3), (0.0, 0.25, 1.5), (0.5, 0.1, 0.2), (0.0, -0.25, -2.0), 0.0, 1.0, 2.0)
    CosineRingwFloor(5.0, 1.0, (), (), (0.5, 0.1, 0.2), (0.0, -0.25, -2.0), 0.0, 1.0, 2.0)
    CosineRingwFloor(5.0, 1.0, (0.1, 0.2, 0.3), (0.0, 0.25, 1.5), (), (), 0.0, 1.0, 2.0)
    test_template(tr1)

    tr2 = CosineRingwGFloor(5.0, 1.0, (0.1, 0.2, 0.3), (0.0, 0.25, 1.5), (0.5, 0.1, 0.2), (0.0, -0.25, -2.0), 0.0, 1.0, 1.0, 2.0)
    test_template(tr2)

    g = imagepixels(10.0, 10.0, 64, 64)
    img1 = intensitymap(tr1, g)
    img2 = intensitymap(tr2, g)
    @test img1 ≈ img2
end


@testset "TemplateLogSpiral" begin
    θ = LogSpiral(r0, τ, σ, δϕ, ξs, x0, y0)
    test_template(θ)
end

@testset "TemplateConstant" begin
    θ = Constant(1.0)
    @test ComradeBase.intensity_point(θ, (1.0, 1.0)) == 1.0
    test_template(θ)
end
