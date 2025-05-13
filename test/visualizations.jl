using Test, VIDA, CairoMakie
include("common.jl")

@testset "VisualizationsTest" begin
    θ = GaussianRing(r0, σ, x0, y0)
    g = imagepixels(fovx, fovy, npix, npix)
    img = intensitymap(θ, g)
    image(g, θ)
    imageviz(img)
    triptic(img, θ)
end
