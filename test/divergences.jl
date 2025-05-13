using Test, VIDA
include("common.jl")


@testset "DivergenceKL" begin
    θ = GaussianRing(r0, σ, x0, y0)
    g = imagepixels(fovx, fovy, npix, npix)
    img = intensitymap(θ, g)
    div = KullbackLeibler(img)
    @test divergence(div, θ) < 1.0e-8
    @inferred divergence(div, θ)
end

@testset "DivergenceBh" begin
    θ = GaussianRing(r0, σ, x0, y0)
    g = imagepixels(fovx, fovy, npix, npix)
    img = intensitymap(θ, g)
    div = Bhattacharyya(img)
    @test divergence(div, θ) < 1.0e-8
    @inferred divergence(div, θ)
end

@testset "DivergenceLS" begin
    θ = GaussianRing(r0, σ, x0, y0)
    g = imagepixels(fovx, fovy, npix, npix)
    img = intensitymap(θ, g)
    div = LeastSquares(img)
    @test divergence(div, θ) < 1.0e-8
    @inferred divergence(div, θ)
end

@testset "DivergenceRy" begin
    θ = GaussianRing(r0, σ, x0, y0)
    θa = GaussianRing(r0 * 1.2, σ / 1.2, x0, y0)
    g = imagepixels(fovx, fovy, npix, npix)
    img = intensitymap(θ, g)
    ry = Renyi(img, 0.5)
    bh = Bhattacharyya(img)
    @test divergence(ry, θ) < 1.0e-10
    @test isapprox(divergence(bh, θa), divergence(ry, θa), rtol = 1.0e-6)
    @inferred divergence(ry, θ)
end


@testset "DivergenceNxCorr" begin
    θ = GaussianRing(r0, σ, x0, y0)
    g = imagepixels(fovx, fovy, npix, npix)
    img = intensitymap(θ, g)
    div = NxCorr(img)
    @test divergence(div, θ) < 1.0e-10
    @inferred divergence(div, θ)
end
