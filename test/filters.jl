using Test,VIDA
using LinearAlgebra
include("common.jl")

@testset "FilterAsymGaussian" begin
    θ = AsymGaussian(σ, τ, ξτ, x0,y0)
    θ1 = AsymGaussian(x0=x0,σ=σ,y0=y0, τ=τ, ξ=ξτ)
    θ2 = AsymGaussian([σ,τ,ξτ,x0,y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)

    @test θ(x0,y0) == 1.0
    @test length(fieldnames(AsymGaussian)) == VIDA.size(AsymGaussian)
    img = VIDA.make_ehtimage(θ, npix, xlim, ylim)
    @test isapprox(flux(img), 1.0, atol=ϵ)
    xcent,ycent = centroid(img)
    @test isapprox(xcent, x0; rtol=ϵ)
    @test isapprox(ycent, y0; rtol=ϵ)
    img_inert = VIDA.inertia(img, true)
    σx2 = σ^2/(1-τ)
    σy2 = σ^2*(1-τ)
    ed = eigen(Symmetric(img_inert))
    @test isapprox(σx2, maximum(eigvals(ed)); rtol=ϵ)
    @test isapprox(σy2, minimum(eigvals(ed)); rtol=ϵ)

end

@testset "FilterDisk" begin
    θ = Disk(r0, α, x0,y0)
    θ1 = Disk(r0=r0,α=α,y0=y0, x0=x0)
    θ2 = Disk([r0,α,x0,y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)

    @test θ(x0,y0) == 1.0
    @test length(fieldnames(Disk)) == VIDA.size(Disk)
    img = VIDA.make_ehtimage(θ, npix, xlim, ylim)
    @test isapprox(flux(img), 1.0, atol=ϵ)
    xcent,ycent = centroid(img)
    @test isapprox(xcent, x0; rtol=ϵ)
    @test isapprox(ycent, y0; rtol=ϵ)
end


@testset "FilterGaussianRing" begin
    θ = GaussianRing(r0,σ,x0,y0)
    θ1 = GaussianRing(r0=r0,x0=x0,σ=σ,y0=y0)
    θ2 = GaussianRing([r0,σ,x0,y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)

    @test θ(r0+x0,y0) == 1.0
    @test length(fieldnames(GaussianRing)) == VIDA.size(GaussianRing)
end

@testset "FilterSlashedGaussianRing" begin
    θ = SlashedGaussianRing(r0,σ,s, ξs, x0, y0)
    θ1 = SlashedGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0, s=s, ξ=ξs)
    θ2 = SlashedGaussianRing([r0, σ, s, ξs, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    xrot,yrot = VIDA.rotate(r0,0, π-ξs)
    @test θ(xrot+x0,yrot+y0) == 1.0
    @test length(fieldnames(SlashedGaussianRing)) == VIDA.size(SlashedGaussianRing)
end


@testset "FilterEllipticalGaussianRing" begin
    θ = EllipticalGaussianRing(r0,σ,τ, ξτ, x0, y0)
    θ1 = EllipticalGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0, τ=τ, ξ=ξτ)
    θ2 = EllipticalGaussianRing([r0, σ, τ, ξτ, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    xrot,yrot = VIDA.rotate(r0/sqrt(1-τ),0, π-ξτ)
    @test θ(xrot+x0,yrot+y0) == 1.0
    @test length(fieldnames(EllipticalGaussianRing)) == VIDA.size(EllipticalGaussianRing)
end


@testset "FilterTIDAGaussianRing" begin
    θ = TIDAGaussianRing(r0,σ,τ, s, ξτ, x0, y0)
    θ1 = TIDAGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0,s=s, τ=τ, ξ=ξτ)
    θ2 = TIDAGaussianRing([r0, σ, τ, s, ξτ, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    xrot,yrot = VIDA.rotate(r0/sqrt(1-τ),0, π-ξτ)
    @test θ(xrot+x0,yrot+y0) == 1.0
    @test length(fieldnames(TIDAGaussianRing)) == VIDA.size(TIDAGaussianRing)
end

@testset "FilterGeneralGaussianRing" begin
    θ = GeneralGaussianRing(r0,σ,τ,ξτ, s, ξs, x0, y0)
    θ1 = GeneralGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0,s=s,ξs=ξs, τ=τ, ξτ=ξτ)
    θ2 = GeneralGaussianRing([r0, σ, τ, ξτ, s, ξs, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    fimg = VIDA.filter_image(θ,npix,xlim,ylim)
    @test (abs(sum(fimg[3])*step(fimg[1])*step(fimg[2]))-730.1726987113645) < ϵ
    @test length(fieldnames(GeneralGaussianRing)) == VIDA.size(GeneralGaussianRing)
end

@testset "FilterConstant" begin
    θ = Constant()
    @test θ(1,1) == 1.0
    @test length(unpack(θ)) == VIDA.size(typeof(θ))
end

@testset "FilterCosineRing" begin
    θg = GeneralGaussianRing(r0,σ,τ, ξτ, s, ξs, x0, y0)
    θsg = CosineRing{1,1}([r0, σ, τ, ξτ, s, ξs, x0, y0])
    _,_,fimg_g = VIDA.filter_image(θg, 64, xlim, ylim)
    _,_,fimg_sg = VIDA.filter_image(θsg, 64, xlim, ylim)
    @test sum(fimg_g/sum(fimg_g) - fimg_sg/sum(fimg_sg)) < 1e-8
    θ = CosineRing{N,M}(r0,
                        [σ, σ_1],
                        [ξσ],
                        τ,
                        ξτ,
                        [s,s_1,s_2],
                        [ξs, ξs_1, ξs_2],
                        x0, y0)
    θ1 = CosineRing{N,M}(r0=r0,
                        σ=[σ, σ_1],
                        ξσ=[ξσ],
                        τ=τ,
                        ξτ=ξτ,
                        s=[s,s_1,s_2],
                        ξs = [ξs, ξs_1, ξs_2],
                        x0=x0, y0=y0)
    θ2 = CosineRing{N,M}([r0,
                         σ, σ_1,
                         ξσ,
                         τ,
                         ξτ,
                         s,s_1,s_2,
                         ξs, ξs_1, ξs_2,
                         x0, y0])
    @test VIDA.size(typeof(θ)) == length(unpack(θ))
    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    #fimg = VIDA.filter_image(θ,npix,xlim,ylim)
    #@test (abs(sum(fimg[3])*step(fimg[1])*step(fimg[2]))-730.1726987113645) < ϵ
    #@test length(fieldnames(GeneralGaussianRing)) == VIDA.size(GeneralGaussianRing)
end

@testset "FilterMulAdd" begin
    θ1 = GaussianRing(r0,σ,x0,y0)
    θ2 = AsymGaussian(σ,τ,ξτ,x0,y0)
    θ = θ1+1.0*θ2
    @test 1.0*θ2 == θ2*1.0
    θs = stack(θ1,θ2)
    println(1.0*θ2)

    @test fieldnames(typeof(1.0*θ2)) == [fieldnames(typeof(θ2))..., :Irel]
    @test fieldnames(typeof(θ1+θ2)) == [fieldnames(typeof(θ1))...,fieldnames(typeof(θ2))...]

    @test VIDA.size(typeof(1.0*θ2)) == length(unpack(1.0*θ2))
    @test VIDA.size(typeof(θ1+θ2)) == length(unpack(θ1+θ2))
    @test θ == θs
    θarr = split(θ)
    @test θ1 == θarr[1]
    @test θ2 == θarr[2].θ
    @test θarr[2].σ == σ
    @test unpack(θ) == unpack(θs)
    @test θ1(x0,y0) + θ2(x0,y0) == θ(x0,y0)
    @test θ1(x0,y0) + θ2(x0,y0) == θs(x0,y0)
end
