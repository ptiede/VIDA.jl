using Test,VIDA,Plots
include("common.jl")

@testset "VisualizationsTest" begin
    θ = GaussianRing(r0, σ, x0, y0)
    g = imagepixels(fovx, fovy, npix, npix)
    img = intensitymap(θ, g)
    plot(θ)
    plot(img)
    triptic(img, θ)
end
