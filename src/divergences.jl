
"""
    $(TYPEDEF)
An abstract class for a divergence of a function. This expects
that a subtype has a field with an EHTImage object and a flux type.
The struct is then assumed to be a **functor** and have a function
that computes the divergence of the image and a template.

For example
```julia
    struct MyDiv{T,F,S} <: AbstractDivergence
        img::EHTImage{T,F}
        flux::S
    end
    function (bh::MyDiv)(θ::AbstractTemplate)
        ...
    end
```

"""
abstract type AbstractDivergence end

"""
    $(TYPEDEF)
Type for the Bhattacharyya divergence. It constructed from an `EHTImage` i.e. data.
Additionally to evaluate the divergence we use a functor approach where if θ
is your
### Details
This computes the Bhattacharyya divergence which is related to Hellinger distance between
two distributions. In fact, they are both minimized at the same point. The Bhattacharyya
divergence is defined as

```math
Bh(f_\\theta||\\hat{I}) = -\\log\\int \\sqrt{f_\\theta(x,y)\\hat{I}(x,y)}dxdy,
```
where ``\\hat{I}`` is defined as the image normalized to unit flux.

"""
struct Bhattacharyya{T<:EHTImage,S} <: AbstractDivergence
    """
    Abstract image class
    """
    img::T
    flux::S
end
function Bhattacharyya(img::T) where {T<:EHTImage}
    Bhattacharyya(img, flux(img))
end
function (bh::Bhattacharyya)(θ::T) where {T<:AbstractTemplate}
    @unpack img, flux = bh
    bsum = zero(eltype(img.img))
    template_norm = zero(eltype(img.img))
    xstart = (-img.nx*img.psize_x + img.psize_x)/2.0
    ystart = (-img.ny*img.psize_y + img.psize_y)/2.0
    for i in 1:img.nx
        @simd for j in 1:img.ny
            x = xstart + img.psize_x*(i-1)
            y = ystart + img.psize_y*(j-1)
            template_value = abs(θ(x,y))
            @inbounds bsum += sqrt(template_value*img.img[j,i])
            template_norm += template_value
        end
    end
    return -log(bsum/sqrt(template_norm*flux))
end






"""
    $(TYPEDEF)
Type for the KL divergence. It constructed from an `EHTImage` i.e. data.
Additionally to evaluate the divergence we use a functor approach where if θ
is your
### Details
This computes the KL divergence which is related to Hellinger distance between
two distributions. In fact, they are both minimized at the same point. The Bhattacharyya
divergence is defined as

```math
KL(f_\\theta||\\hat{I}) = -\\log\\int f_{\\theta}(x,y)\\log
        \\left(\\frac{f_{\\theta}(x,y)}{\\hat{I}(x,y)}\\right)dxdy,
```
where ``\\hat{I}`` is defined as the image normalized to unit flux.

This struct is also a functor.
"""
struct KullbackLeibler{T,S} <: AbstractDivergence
    """
    Abstract image class
    """
    img::T
    flux::S
end
function KullbackLeibler(img::T) where {T<:EHTImage}
    KullbackLeibler(img, flux(img))
end


function (kl::KullbackLeibler)(θ::T) where {T<:AbstractTemplate}
    @unpack img, flux = kl
    klsum = zero(eltype(img.img))
    template_norm = zero(eltype(img.img))
    xstart = (-img.nx*img.psize_x + img.psize_x)/2.0
    ystart = (-img.ny*img.psize_y + img.psize_y)/2.0

    @inbounds for i in 1:img.nx
        @simd for j in 1:img.ny
            x = xstart + img.psize_x*(i-1)
            y = ystart + img.psize_y*(j-1)
            template_value = θ(x,y)+1e-12
            klsum += template_value*log(template_value/(img.img[j,i]+1e-12))
            template_norm += template_value
        end
    end
    return (klsum/template_norm - log(template_norm/flux))
end


"""
    $(TYPEDEF)
Type for the least squares divergence. It constructed from an `EHTImage` i.e. data.
Additionally to evaluate the divergence we use a functor approach where if θ
is your
### Details
This computes the squared 2 norm between your image and template, both of which
are normalized to unit flux.

To construct this just pass it an image object
```julia
ls = LeastSquares(img::EHTImage)
```

# Notes
We have a template internal matrix the prevents internal allocations during computation
This is a internal feature that the user does not need to worry about.
"""
struct LeastSquares{T,S,V} <: AbstractDivergence
    img::T
    flux::S
    template::V
end

function LeastSquares(img::T) where {T<:EHTImage}
    LeastSquares(img, flux(img), zeros(size(img)))
end


function (ls::LeastSquares)(θ::AbstractTemplate)
    @unpack img, flux, template = ls
    template_norm = zero(eltype(img.img))
    xstart = (-img.nx*img.psize_x + img.psize_x)/2.0
    ystart = (-img.ny*img.psize_y + img.psize_y)/2.0

    @inbounds for i in 1:img.nx
        @simd for j in 1:img.ny
            x = xstart + img.psize_x*(i-1)
            y = ystart + img.psize_y*(j-1)
            template_value = θ(x,y)+1e-12
            template[j,i] = template_value
            template_norm += template_value
        end
    end
    template .= template./template_norm .- img./flux
    return sum(abs2, template)
end
