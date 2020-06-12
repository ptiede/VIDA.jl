using Test,VIDA

include("common.jl")

@testset "FilterAsymGaussian" begin
    θ = AsymGaussian(σ, τ, ξτ, x0,y0)
    θ1 = AsymGaussian(x0=x0,σ=σ,y0=y0, τ=τ, ξ=ξτ)
    θ2 = AsymGaussian([σ,τ,ξτ,x0,y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)

    @test imagefilter(x0,y0,θ) == 1.0
    @test length(fieldnames(AsymGaussian)) == VIDA.size(AsymGaussian)
end

@testset "FilterGaussianRing" begin
    θ = GaussianRing(r0,σ,x0,y0)
    θ1 = GaussianRing(r0=r0,x0=x0,σ=σ,y0=y0)
    θ2 = GaussianRing([r0,σ,x0,y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)

    @test imagefilter(r0+x0,y0,θ) == 1.0
    @test length(fieldnames(GaussianRing)) == VIDA.size(GaussianRing)
end

@testset "FilterSlashedGaussianRing" begin
    θ = SlashedGaussianRing(r0,σ,s, ξs, x0, y0)
    θ1 = SlashedGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0, s=s, ξ=ξs)
    θ2 = SlashedGaussianRing([r0, σ, s, ξs, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    xrot,yrot = VIDA.rotate(r0,0, π-ξs)
    @test imagefilter(xrot+x0,yrot+y0, θ) == 1.0
    @test length(fieldnames(SlashedGaussianRing)) == VIDA.size(SlashedGaussianRing)
end


@testset "FilterEllipticalGaussianRing" begin
    θ = EllipticalGaussianRing(r0,σ,τ, ξτ, x0, y0)
    θ1 = EllipticalGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0, τ=τ, ξ=ξτ)
    θ2 = EllipticalGaussianRing([r0, σ, τ, ξτ, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    xrot,yrot = VIDA.rotate(r0/sqrt(1-τ),0, π-ξτ)
    @test imagefilter(xrot+x0,yrot+y0, θ) == 1.0
    @test length(fieldnames(EllipticalGaussianRing)) == VIDA.size(EllipticalGaussianRing)
end


@testset "FilterTIDAGaussianRing" begin
    θ = TIDAGaussianRing(r0,σ,τ, s, ξτ, x0, y0)
    θ1 = TIDAGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0,s=s, τ=τ, ξ=ξτ)
    θ2 = TIDAGaussianRing([r0, σ, τ, s, ξτ, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    xrot,yrot = VIDA.rotate(r0/sqrt(1-τ),0, π-ξτ)
    @test imagefilter(xrot+x0,yrot+y0, θ) == 1.0
    @test length(fieldnames(TIDAGaussianRing)) == VIDA.size(TIDAGaussianRing)
end

@testset "FilterGeneralGaussianRing" begin
    θ = GeneralGaussianRing(r0,σ,τ,ξτ, s, ξs, x0, y0)
    θ1 = GeneralGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0,s=s,ξs=ξs, τ=τ, ξτ=ξτ)
    θ2 = GeneralGaussianRing([r0, σ, τ, ξτ, s, ξs, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    @test (sum(VIDA.filter_image(θ,npix,xlim,ylim)[3]) - 193.31484083354985) < ϵ
    @test length(fieldnames(GeneralGaussianRing)) == VIDA.size(GeneralGaussianRing)
end

@testset "FilterAddMul" begin
    θ1 = GaussianRing(r0,σ,x0,y0)
    θ2 = AsymGaussian(σ,τ,ξτ,x0,y0)
    θ = θ1+1.0*θ2
    θs = cat(θ1,θ2)
    @test θ == θs
    θarr = split(θ)
    @test θ1 == θarr[1]
    @test θ2 == θarr[2].θ
    @test unpack(θ) == unpack(θs)
    @test imagefilter(x0,y0,θ1)+imagefilter(x0,y0,θ2) == imagefilter(x0,y0,θ)
    @test imagefilter(x0,y0,θ1)+imagefilter(x0,y0,θ2) == imagefilter(x0,y0,θs)
end
