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
