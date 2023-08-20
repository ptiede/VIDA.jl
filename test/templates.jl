using Test,VIDA
using LinearAlgebra
using Interpolations: Lanczos
include(joinpath(@__DIR__, "common.jl"))

function test_template(θ::AbstractTemplate; npix=128, fov=120.0, )
    @inferred ComradeBase.intensity_point(θ, (X=0.0, Y=0.0))
    @test ComradeBase.intensity_point(θ, (X=0.0, Y=0.0)) isa AbstractFloat
    @inferred ComradeBase.radialextent(θ)
    @test ComradeBase.radialextent(θ) isa AbstractFloat
end


@testset "GaussianRing" begin
    θ = AsymGaussian(σ, τ, ξτ, x0,y0)
    θ1 = AsymGaussian(x0=x0,σ=σ,y0=y0, τ=τ, ξ=ξτ)
    θ2 = AsymGaussian([σ,τ,ξτ,x0,y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)

    @test θ(x0,y0) == 1.0
    @test length(propertynames(θ)) == VIDA.size(AsymGaussian)
    @test size(θ) == VIDA.size(AsymGaussian)
    img = VIDA.make_image(θ, npix, xlim, ylim)
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

@testset "TemplateDisk" begin
    θ = Disk(r0, α, x0,y0)
    θ1 = Disk(r0=r0,α=α,y0=y0, x0=x0)
    θ2 = Disk([r0,α,x0,y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)

    @test θ(x0,y0) == 1.0
    @test length(propertynames(θ)) == VIDA.size(Disk)
    img = VIDA.make_image(θ, npix, xlim, ylim)
    @test isapprox(flux(img), 1.0, atol=ϵ)
    xcent,ycent = centroid(img)
    @test isapprox(xcent, x0; rtol=ϵ)
    @test isapprox(ycent, y0; rtol=ϵ)
end


@testset "TemplateGaussianRing" begin
    θ = GaussianRing(r0,σ,x0,y0)
    θ1 = GaussianRing(r0=r0,x0=x0,σ=σ,y0=y0)
    θ2 = GaussianRing([r0,σ,x0,y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)

    @test θ(r0+x0,y0) == 1.0
    @test length(propertynames(θ)) == VIDA.size(GaussianRing)
end

@testset "TemplateSlashedGaussianRing" begin
    θ = SlashedGaussianRing(r0,σ,s, ξs, x0, y0)
    θ1 = SlashedGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0, s=s, ξ=ξs)
    θ2 = SlashedGaussianRing([r0, σ, s, ξs, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    xrot,yrot = VIDA.rotate(r0,0, π-ξs)
    @test θ(xrot+x0,yrot+y0) == 1.0
    @test length(propertynames(θ)) == VIDA.size(SlashedGaussianRing)
end


@testset "TemplateEllipticalGaussianRing" begin
    θ = EllipticalGaussianRing(r0,σ,τ, ξτ, x0, y0)
    θ1 = EllipticalGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0, τ=τ, ξ=ξτ)
    θ2 = EllipticalGaussianRing([r0, σ, τ, ξτ, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    xrot,yrot = VIDA.rotate(r0/sqrt(1-τ),0, π-ξτ)
    @test θ(xrot+x0,yrot+y0) == 1.0
    @test length(propertynames(θ)) == VIDA.size(EllipticalGaussianRing)
end


@testset "TemplateTIDAGaussianRing" begin
    θ = TIDAGaussianRing(r0,σ,τ, s, ξτ, x0, y0)
    θ1 = TIDAGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0,s=s, τ=τ, ξ=ξτ)
    θ2 = TIDAGaussianRing([r0, σ, τ, s, ξτ, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    xrot,yrot = VIDA.rotate(r0/sqrt(1-τ),0, π-ξτ)
    @test θ(xrot+x0,yrot+y0) == 1.0
    @test length(propertynames(θ)) == VIDA.size(TIDAGaussianRing)
end

@testset "TemplateGeneralGaussianRing" begin
    θ = GeneralGaussianRing(r0,σ,τ,ξτ, s, ξs, x0, y0)
    θ1 = GeneralGaussianRing(r0=r0,x0=x0,σ=σ,y0=y0,s=s,ξs=ξs, τ=τ, ξτ=ξτ)
    θ2 = GeneralGaussianRing([r0, σ, τ, ξτ, s, ξs, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    fimg = VIDA.template_image(θ,npix,xlim,ylim)
    @test (abs(sum(fimg[3])*step(fimg[1])*step(fimg[2]))-730.1726987113645) < ϵ
    @test length(propertynames(θ)) == VIDA.size(GeneralGaussianRing)
end

@testset "TemplateLogSpiral" begin
    θ = LogSpiral(r0, τ, σ, δϕ, ξs, x0, y0)
    θ1 = LogSpiral(r0=r0,x0=x0,σ=σ,y0=y0,κ=τ,ξ=ξs, δϕ=δϕ)
    θ2 = LogSpiral([r0, τ, σ, δϕ, ξs, x0, y0])

    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    fimg = VIDA.template_image(θ,npix,xlim,ylim)
    @test length(propertynames(θ)) == VIDA.size(typeof(θ))
end


@testset "TemplateImageTemplate" begin
    tmp = GaussianRing(r0, σ, x0, y0)
    img = VIDA.make_image(tmp, 60, (-60.0, 60.0), (-60.0, 60.0))
    θ = ImageTemplate(0.0, 0.0, img)
    sitp = VIDA.image_interpolate(img, Lanczos())
    θ2 = ImageTemplate(0.0, 0.0, sitp)
    bh = Bhattacharyya(img)
    @test bh(θ) < 1e-6
end

@testset "TemplateConstant" begin
    θ = Constant()
    @test θ(1,1) == 1.0
    @test length(unpack(θ)) == VIDA.size(typeof(θ))
end

@testset "TemplateSymCosineRing" begin
    θg = SlashedGaussianRing(r0,σ, s, ξs, x0, y0)
    θsg = SymCosineRing{0,1}([r0, σ, s, ξs, x0, y0])
    _,_,fimg_g = VIDA.template_image(θg, 64, xlim, ylim)
    _,_,fimg_sg = VIDA.template_image(θsg, 64, xlim, ylim)
    @test sum(fimg_g/sum(fimg_g) - fimg_sg/sum(fimg_sg)) < 1e-8
    θ = SymCosineRing{N-1,M}(r0,
                        [σ, σ_1],
                        [ξσ],
                        [s,s_1,s_2],
                        [ξs, ξs_1, ξs_2],
                        x0, y0)
    θ1 = SymCosineRing{N-1,M}(r0=r0,
                        σ=[σ, σ_1],
                        ξσ=[ξσ],
                        s=[s,s_1,s_2],
                        ξs = [ξs, ξs_1, ξs_2],
                        x0=x0, y0=y0)
    θ2 = SymCosineRing{N-1,M}([r0,
                         σ, σ_1,
                         ξσ,
                         s,s_1,s_2,
                         ξs, ξs_1, ξs_2,
                         x0, y0])
    @inferred θ2(1.0, 1.0)
    @test VIDA.size(typeof(θ)) == length(unpack(θ))
    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    #fimg = VIDA.template_image(θ,npix,xlim,ylim)
    #@test (abs(sum(fimg[3])*step(fimg[1])*step(fimg[2]))-730.1726987113645) < ϵ
    #@test length(propertynames(GeneralGaussianRing)) == VIDA.size(GeneralGaussianRing)
end

@testset "TemplateSymCosineRingwFloor" begin
    θg = SlashedGaussianRing(r0,σ, s, ξs, x0, y0)
    θsg = SymCosineRingwFloor{0,1}([r0, σ, s, ξs, 0.0, x0, y0])
    _,_,fimg_g = VIDA.template_image(θg, 64, xlim, ylim)
    _,_,fimg_sg = VIDA.template_image(θsg, 64, xlim, ylim)
    @test sum(fimg_g/sum(fimg_g) - fimg_sg/sum(fimg_sg)) < 1e-8
    θ = SymCosineRingwFloor{N-1,M}(r0,
                        [σ, σ_1],
                        [ξσ],
                        [s,s_1,s_2],
                        [ξs, ξs_1, ξs_2],
                        0.1,
                        x0, y0)
    θ1 = SymCosineRingwFloor{N-1,M}(r0=r0,
                        σ=[σ, σ_1],
                        ξσ=[ξσ],
                        s=[s,s_1,s_2],
                        ξs = [ξs, ξs_1, ξs_2],
                        floor = 0.1,
                        x0=x0, y0=y0)
    θ2 = SymCosineRingwFloor{N-1,M}([r0,
                         σ, σ_1,
                         ξσ,
                         s,s_1,s_2,
                         ξs, ξs_1, ξs_2,
                         0.1,
                         x0, y0])
    @inferred θ2(1.0, 1.0)
    @test VIDA.size(typeof(θ)) == length(unpack(θ))
    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    #fimg = VIDA.template_image(θ,npix,xlim,ylim)
    #@test (abs(sum(fimg[3])*step(fimg[1])*step(fimg[2]))-730.1726987113645) < ϵ
    #@test length(propertynames(GeneralGaussianRing)) == VIDA.size(GeneralGaussianRing)
end


@testset "TemplateSymCosineRingwGFloor" begin
    θg = SlashedGaussianRing(r0,σ, s, ξs, x0, y0)
    θsg = SymCosineRingwGFloor{0,1}([r0, σ, s, ξs, 0.0, 30.0, x0, y0])
    _,_,fimg_g = VIDA.template_image(θg, 64, xlim, ylim)
    _,_,fimg_sg = VIDA.template_image(θsg, 64, xlim, ylim)
    @test sum(fimg_g/sum(fimg_g) - fimg_sg/sum(fimg_sg)) < 1e-8
    θ = SymCosineRingwGFloor{N-1,M}(r0,
                        [σ, σ_1],
                        [ξσ],
                        [s,s_1,s_2],
                        [ξs, ξs_1, ξs_2],
                        0.1, 30.0,
                        x0, y0)
    θ1 = SymCosineRingwGFloor{N-1,M}(r0=r0,
                        σ=[σ, σ_1],
                        ξσ=[ξσ],
                        s=[s,s_1,s_2],
                        ξs = [ξs, ξs_1, ξs_2],
                        floor = 0.1, σg=30.0,
                        x0=x0, y0=y0)
    θ2 = SymCosineRingwGFloor{N-1,M}([r0,
                         σ, σ_1,
                         ξσ,
                         s,s_1,s_2,
                         ξs, ξs_1, ξs_2,
                         0.1, 30.0,
                         x0, y0])
    @inferred θ2(1.0, 1.0)
    @test VIDA.size(typeof(θ)) == length(unpack(θ))
    @test unpack(θ) == unpack(θ1)
    @test unpack(θ1) == unpack(θ2)
    #fimg = VIDA.template_image(θ,npix,xlim,ylim)
    #@test (abs(sum(fimg[3])*step(fimg[1])*step(fimg[2]))-730.1726987113645) < ϵ
    #@test length(propertynames(GeneralGaussianRing)) == VIDA.size(GeneralGaussianRing)
end



@testset "TemplateCosineRing" begin
    θg = GeneralGaussianRing(r0,σ,τ, ξτ, s, ξs, x0, y0)
    θsg = CosineRing{0,1}([r0, σ, τ, ξτ, s, ξs, x0, y0])
    _,_,fimg_g = VIDA.template_image(θg, 64, xlim, ylim)
    _,_,fimg_sg = VIDA.template_image(θsg, 64, xlim, ylim)
    @test sum(fimg_g/sum(fimg_g) - fimg_sg/sum(fimg_sg)) < 1e-8
    θ = CosineRing{N-1,M}(r0,
                        [σ, σ_1],
                        [ξσ],
                        τ,
                        ξτ,
                        [s,s_1,s_2],
                        [ξs, ξs_1, ξs_2],
                        x0, y0)
    θ1 = CosineRing{N-1,M}(r0=r0,
                        σ=[σ, σ_1],
                        ξσ=[ξσ],
                        τ=τ,
                        ξτ=ξτ,
                        s=[s,s_1,s_2],
                        ξs = [ξs, ξs_1, ξs_2],
                        x0=x0, y0=y0)
    θ2 = CosineRing{N-1,M}([r0,
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
    #fimg = VIDA.template_image(θ,npix,xlim,ylim)
    #@test (abs(sum(fimg[3])*step(fimg[1])*step(fimg[2]))-730.1726987113645) < ϵ
    #@test length(propertynames(GeneralGaussianRing)) == VIDA.size(GeneralGaussianRing)
end


@testset "TemplateMods" begin
    θ = SlashedGaussianRing(r0, σ, s, ξs, -20.0, 10.0)
    θ0 = SlashedGaussianRing(r0, σ, s, 0.0,-20.0, 10.0)
    rθ = VIDA.rotate(θ, -ξs)
    rimg = intensitymap(rθ, 320.0, 320.0, 512, 512)
    img0 = intensitymap(θ0, 320.0, 320.0, 512, 512)
    @test rimg.img ≈ img0.img
end
