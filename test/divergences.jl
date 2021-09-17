using Test,VIDA
include("common.jl")


@testset "DivergenceKL" begin
    θ = GaussianRing(r0,σ,x0,y0)
    img = VIDA.make_image(θ, 64, xlim, ylim)
    kl = KullbackLeibler(img)
    @test kl(θ) < 1e-8
    @inferred kl(θ)
end

@testset "DivergenceBh" begin
    θ = GaussianRing(r0,σ,x0,y0)
    img = VIDA.make_image(θ, 64, xlim, ylim)
    bh = Bhattacharyya(img)
    @test bh(θ) < 1e-8
    @inferred bh(θ)
end

@testset "DivergenceLS" begin
    θ = GaussianRing(r0,σ,x0,y0)
    img = VIDA.make_image(θ, 64, xlim, ylim)
    bh = LeastSquares(img)
    @test bh(θ) < 1e-10
    @inferred bh(θ)
end

@testset "DivergenceRy" begin
    θ = GaussianRing(r0,σ,x0,y0)
    θa = GaussianRing(r0*1.2, σ/1.2, x0, y0)
    img = VIDA.make_image(θ, 64, xlim, ylim)
    ry = Renyi(img, 0.5)
    bh = Bhattacharyya(img)
    @test ry(θ) < 1e-10
    @test isapprox(bh(θa)*2, ry(θa), rtol=1e-6)
    @inferred ry(θ)
end
