using Test,VIDA
using OptimizationBBO
using OptimizationCMAEvolutionStrategy
include("common.jl")

@testset "VIDA" begin
    temp(θ) = SlashedGaussianRing(θ.r0, θ.σ, θ.s, θ.ξs, θ.x0, θ.y0)+
               modify(Gaussian(),
                     Stretch(θ.σG*sqrt(1-θ.τG), θ.σG/sqrt(1- θ.τG)),
                     Rotate(θ.ξG),
                     Shift(θ.xG, θ.yG),
                     Renormalize(θ.fG)
                     )

    lower = (r0 = 5.0, σ = 0.1, s = 0.001, ξs = -1π, x0 = -60.0, y0 = -60.0,
             σG = 0.1, τG = 0.001, ξG = -π/2, xG = -60.0, yG = -60.0, fG = 0.0)
    upper = (r0 = 30.0, σ = 10.0, s = 0.999, ξs = 1π, x0 = 60.0, y0 = 60.0,
             σG = 30.0, τG = 0.999, ξG = π/2, xG = 60.0, yG = 60.0, fG = 50.0)

    p0 = (;r0, σ, s, ξs, x0, y0,
         σG=σ*2, τG=τ, ξG=π/4, xG = 10.0, yG = 10.0, fG=2.0)
    θ = temp(p0)
    g = imagepixels(fovx, fovy, 64, 64)
    fimg = intensitymap(θ, g)

    bh = Bhattacharyya(fimg)

    prob = VIDAProblem(bh, temp, lower, upper)
    xopt, θopt, dmin = vida(prob, BBO_adaptive_de_rand_1_bin(); maxiters=100_000)
    map((x,y)->@test(isapprox(x, y, atol=1e-2)), xopt, p0)

    xopt, θopt, dmin = vida(prob, CMAEvolutionStrategyOpt(); init_params=p0, unit_cube=false, maxiters=100_000)
    map((x,y)->@test(isapprox(x, y, atol=1e-2)), xopt, p0)
end
