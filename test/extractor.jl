using Test,VIDA, Optim
include("common.jl")

@testset "OptimExtractor" begin
    θ = SlashedGaussianRing(r0, σ, s, ξs, x0, y0) +
        0.1*AsymGaussian(σ*2, τ, ξτ, x0-2.0, y0+10.0)
    lower = typeof(θ)([1.0, 0.1, 0.001, -π, -60.0, -60.0,
             0.1, 0.001, 0.0, -60.0, -60.0, 1e-6])
    upper = typeof(θ)([40.0, 20.0, 0.999, π, 60.0, 60.0,
             30.0, 0.999, π, 60.0, 60.0, 1.0])
    fimg = VIDA.make_image(θ, 64, xlim, ylim)
    bh = Bhattacharyya(fimg)
    θ0 = SlashedGaussianRing(r0*1.05, σ*1.05, s*0.95, ξs, x0, y0) +
                            0.11*AsymGaussian(σ*1.5, τ*0.95, ξτ*1.1, x0, y0)
    prob = ExtractProblem(bh, θ0, lower, upper)
    rθ,divmin = threaded_extractor(1, prob, Opt(Fminbox(LBFGS())))
    rθ,divmin = extractor(prob, Opt(Fminbox(LBFGS())))
    @test isapprox(unpack(θ), unpack(θ), rtol=ϵ)
end


@testset "BBExtractor" begin
    θ = SlashedGaussianRing(r0, σ, s, ξs, x0, y0) +
        0.1*AsymGaussian(σ*2, τ, ξτ, x0-2.0, y0+10.0)
    lower = typeof(θ)([1.0, 0.1, 0.001, -π, -60.0, -60.0,
             0.1, 0.001, 0.0, -60.0, -60.0, 1e-6])
    upper = typeof(θ)([40.0, 20.0, 0.999, π, 60.0, 60.0,
             30.0, 0.999, π, 60.0, 60.0, 1.0])
    fimg = VIDA.make_image(θ, 64, xlim, ylim)
    bh = Bhattacharyya(fimg)
    θ0 = SlashedGaussianRing(r0*2, σ*1.5, s*1.1, 0.5, 0.0, 0.0) +
         0.5*AsymGaussian(5.0, 0.5, 0.25, 0.0, 0.0)
    prob = ExtractProblem(bh, θ0, lower, upper)
    rθ,divmin = threaded_extractor(4, prob, BBO(tracemode=:silent))
    prob2 = ExtractProblem(bh, rθ, lower, upper)
    rθ,divmin = threaded_extractor(4, prob2, CMAES(verbosity=0,ftol=1e-20, cov_scale=0.01))
    @test isapprox(unpack(rθ), unpack(θ), rtol=1e-1)
end
