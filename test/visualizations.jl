using Test,VIDA
include("common.jl")

@testset "VisualizationsTest" begin
    θ = GaussianRing(r0,σ,x0,y0)
    img = VIDA.make_ehtimage(θ, 64, xlim, ylim)
    plot(θ)
    plot(img)
    triptic(img, θ)
end
