using Test,VIDA
include("common.jl")


@testset "DivergenceKL" begin
    θ = GaussianRing(r0,σ,x0,y0)
    img = intensitymap(θ, fovx, fovy, npix, npix)
    div = KullbackLeibler(img)
    @test divergence(div, θ) < 1e-8
    @inferred divergence(div, θ)
end

@testset "DivergenceBh" begin
    θ = GaussianRing(r0,σ,x0,y0)
    img = intensitymap(θ, fovx, fovy, npix, npix)
    div = Bhattacharyya(img)
    @test divergence(div, θ) < 1e-8
    @inferred divergence(div, θ)
end

@testset "DivergenceLS" begin
    θ = GaussianRing(r0,σ,x0,y0)
    img = intensitymap(θ, fovx, fovy, npix, npix)
    div = LeastSquares(img)
    @test divergence(div, θ) < 1e-8
    @inferred divergence(div, θ)
end

@testset "DivergenceRy" begin
    θ = GaussianRing(r0,σ,x0,y0)
    θa = GaussianRing(r0*1.2, σ/1.2, x0, y0)
    img = intensitymap(θ, fovx, fovy, npix, npix)
    ry = Renyi(img, 0.5)
    bh = Bhattacharyya(img)
    @test divergence(ry, θ) < 1e-10
    @test isapprox(divergence(bh, θa)*2, divergence(ry,θa), rtol=1e-6)
    @inferred divergence(ry, θ)
end
