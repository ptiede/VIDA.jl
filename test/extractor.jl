using Test,VIDA
include("common.jl")

@testset "Extractor" begin
    θ = SlashedGaussianRing(r0, σ, s, ξs, x0, y0) +
        0.1*AsymGaussian(σ*2, τ, ξτ, x0-2.0, y0+10.0)
    lower = [1.0, 0.1, 0.001, -π, -60.0, -60.0,
             0.1, 0.001, 0.0, -60.0, -60.0, 1e-6]
    upper = [40.0, 20.0, 0.999, π, 60.0, 60.0,
             30.0, 0.999, π, 60.0, 60.0, 1.0]
    fimg = VIDA.make_ehtimage(θ, 64, xlim, ylim)
    bh = Bhattacharyya(fimg)
    θ0 = SlashedGaussianRing(r0*1.05, σ*1.05, s*0.95, ξs, x0, y0) +
                            0.11*AsymGaussian(σ*1.5, τ*0.95, ξτ*1.1, x0, y0)
    df = extract(1,bh, θ0, lower, upper)
    rθ,divmin,_,_ = extract(bh, θ0, lower, upper)
    @test isapprox(unpack(θ), unpack(θ), rtol=ϵ)
end


@testset "BBExtractor" begin
    θ = SlashedGaussianRing(r0, σ, s, ξs, x0, y0) +
        0.1*AsymGaussian(σ*2, τ, ξτ, x0-2.0, y0+10.0)
    lower = [1.0, 0.1, 0.001, -π, -60.0, -60.0,
             0.1, 0.001, 0.0, -60.0, -60.0, 1e-6]
    upper = [40.0, 20.0, 0.999, π, 60.0, 60.0,
             30.0, 0.999, π, 60.0, 60.0, 1.0]
    fimg = VIDA.make_ehtimage(θ, 64, xlim, ylim)
    bh = Bhattacharyya(fimg)
    θ0 = SlashedGaussianRing(r0*2, σ*1.5, s*1.1, 0.5, 0.0, 0.0) +
         0.5*AsymGaussian(5.0, 0.5, 0.25, 0.0, 0.0)
    rθ,divmin,_,_ = bbextract(bh, θ0, lower, upper, MaxFuncEvals=40000)
    rθ,divmin,_,_ = extract(bh, rθ, lower, upper)
    @test isapprox(unpack(rθ), unpack(θ), rtol=ϵ)
end
