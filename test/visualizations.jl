using Test,VIDA,Plots
include("common.jl")

@testset "VisualizationsTest" begin
    θ = GaussianRing(r0, σ, x0, y0)
    img = intensitymap(θ, fovx, fovxy, npix, npix)
    plot(θ)
    plot(img)
    triptic(img, θ)
end
