using Test,VIDA
include("common.jl")


@testset "DivergenceKL" begin
    θ = GaussianRing(r0,σ,x0,y0)
    img = VIDA.make_ehtimage(θ, 64, xlim, ylim)
    kl = KullbackLeibler(img)
    @test kl(θ) < 1e-8
end

@testset "DivergenceBh" begin
    θ = GaussianRing(r0,σ,x0,y0)
    img = VIDA.make_ehtimage(θ, 64, xlim, ylim)
    kl = Bhattacharyya(img)
    @test kl(θ) < 1e-8
end
