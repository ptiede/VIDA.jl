export divergence, nxcorr

"""
    $(TYPEDEF)

Defines an `AbstractDivergence` type. This is the basic cost function of `VIDA`.

A subtype of `AbstractDivergence` must be a `struct` with at least two fields
  - `img::IntensityMap` which holds the image you wish to run feature extraction on
  - `mimg::IntensityMap` which is an internal buffer that stores the model image.

The divergence is evaluated roughly as

```julia
    normalize_div(div, sum(divergence_point(div, image, model)))
```


Therefore a user must implement the following methods

  - [`divergence_point`](@ref): Which evaluates the divergence at a single pixel
  - [`normalize_div`](@ref): Which normalizes and modifies the

"""
abstract type AbstractDivergence end

struct InplaceDivergence{D<:AbstractDivergence, T} <: AbstractDivergence
    div::D
    img::T
    mimg::T
    function InplaceDivergence(div::D, img::IntensityMap) where {D<:AbstractDivergence}
        all(>(0), img) || throw(ArgumentError("All intensities must be positive.")) 
        nimg = img./flux(img)
        T = typeof(nimg)
        return new{D,T}(div, nimg, zero(nimg))
    end
end

"""
    divergence_point(div::AbstractDivergence, p, q)

Evaluate the divergence `div` at a single pixel with value `p` and `q` for image and
template. The total divergence
"""
function divergence_point end

"""
    normalize_div(div, val)

Returns the proper divergence whose values have been normalized so that the value is always
postive semi-definite and is only zero then the template and the image are equal.
"""
function normalize_div end


"""
    $(SIGNATURES)

Computes the divergence of `d` with respect to the template `m`. This measure how similar
the divergence with respect to the image stored in `d` is compared to the template.

# Arguments
  - `d`: An divergence or a subtype of [`VIDA.AbstractDivergence`](@ref)
  - `m`: A template which `ComradeBase.AbstractModel` or more commonly
         a [`VIDA.AbstractImageTemplate`](@ref)
"""
@inline function divergence(d::AbstractDivergence, m::ComradeBase.AbstractModel)
    _divergence(d, m)
end

function divergence(d::AbstractDivergence, m::IntensityMap{<:Real})
    _divergence(d, m)
end

# @inline function _divergence(::IsAnalytic, d::AbstractDivergence, m::AbstractModel)
#     return divergence_analytic(d, m)
# end

# @inline function _divergence(::NotAnalytic, d::AbstractDivergence, m::AbstractModel)
#     return divergence_numeric(d, m)
# end


# function _divergence(d::AbstractDivergence, m::AbstractModel)
#     (;img) = d
#     g = CB.imagegrid(img)
#     div, fm = mapreduce(.+, zip(g, img)) do (p, I)
#         Imod = CB.intensity_point(m, p)
#         return divergence_point(d, I, Imod), Imod
#     end
#     return normalize_div(d, div, fm)
# end

function _divergence(d::InplaceDivergence, m::ComradeBase.AbstractModel)
    (;img, mimg) = d
    CB.intensitymap!(mimg, m)
    fm = sum(x->max(x, 0), mimg)
    return __divergence(d.div, img, mimg, fm)
end

function _divergence(d::InplaceDivergence, mimg::IntensityMap{<:Real})
    (;img) = d
    fm = sum(x->max(x, 0), mimg)
    return __divergence(d.div, img, mimg, fm)
end


function __divergence(d, img, mimg, fm)
    div  = mapreduce(+, parent(img), parent(mimg)) do ii, im
        return divergence_point(d, ii, max(im/fm, 0))
    end
    # div = sum(divergence_point(d, ii, max(im/fm, 0))))
    return normalize_div(d, div)
end




"""
    $(TYPEDEF)
Type for the Bhattacharyya divergence. It constructed from an `IntensityMap` i.e. data.
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
struct Bhattacharyya <: AbstractDivergence end

function Bhattacharyya(img::T) where {T<:IntensityMap}
    InplaceDivergence(Bhattacharyya(), img)
end


@inline function divergence_point(::Bhattacharyya, p, q)
    return sqrt(p*q)
end
@inline normalize_div(::Bhattacharyya, div) = -log(div)


"""
    $(TYPEDEF)
Type for the KL divergence. It constructed from an `IntensityMap` i.e. data.
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
struct KullbackLeibler <: AbstractDivergence end
function KullbackLeibler(img::IntensityMap)
    return InplaceDivergence(KullbackLeibler(), img)
end

@inline divergence_point(::KullbackLeibler, p, q) = q*log(q/(p+eps(typeof(p))))
@inline normalize_div(::KullbackLeibler, div) = div


struct Renyi{S} <: AbstractDivergence
    α::S
end

"""
    Renyi(img::IntensityMap, α)
Construct the Renyi divergence with parameter α. It constructed from an `IntensityMap` i.e. data.

### Details
This computes the KL divergence which is related to Hellinger distance between
two distributions. In fact, they are both minimized at the same point. The Bhattacharyya
divergence is defined as

```math
Ry(f_\\theta||\\hat{I}) = \\frac{1}{α-1}\\log\\int
        \\left(\\frac{f_{\\theta}(x,y)^\\alpha}{\\hat{I}(x,y)^{\\alpha-1}}\\right)dxdy,
```
where ``\\hat{I}`` is defined as the image normalized to unit flux.

This is a very flexible divergence that reduces to many of the other divergences implemented.
 - `α = 1` corresponds to the KL divergence
 - `α = 1/2` corresponds to the Bhattacharyya divergence up to a multiplicative factor of 2

Typically we find that `α=1.5` works well, as it focusses on the bright regions of the images
moreso than the Bh and KL divergence. For `α>2` the measure tends to devolve in something
akin the to sup norm and fails to match the image structure.
"""
function Renyi(img::T, α) where {T<:IntensityMap}
    @assert !(α-1 ≈ 0) "α=1 is the KL divergence use that instead"
    InplaceDivergence(Renyi(α), img)
end

@inline divergence_point(d::Renyi, p, q) = q*(q/(p + eps()))^(d.α-1)
@inline function normalize_div(d::Renyi, div)
    α = d.α
    return inv(α-1)*log(div)/2
end



"""
    $(TYPEDEF)
Type for the least squares divergence. It constructed from an `IntensityMap` i.e. data.
Additionally to evaluate the divergence we use a functor approach where if θ
is your
### Details
This computes the squared 2 norm between your image and template, both of which
are normalized to unit flux.

To construct this just pass it an image object
```julia-repl
julia> ls = LeastSquares(img::IntensityMap)
```

# Notes
We have a template internal matrix the prevents internal allocations during computation
This is a internal feature that the user does not need to worry about.
"""
struct LeastSquares <: AbstractDivergence end

function LeastSquares(img::SpatialIntensityMap)
    return InplaceDivergence(LeastSquares(), img)
end

function divergence_point(::LeastSquares, p, q)
    return abs2(p - q)
end

@inline normalize_div(::LeastSquares, div) = div

struct NxCorr <: AbstractDivergence end

"""
    NxCorr(img::IntensityMap)

Construct the normalized cross correlation (NXCORR) divergence with respect to the image `img`.
To maximize the NXCorr we instead compute the -log(|NxCorr|) as the divergence

NxCorr is defined as:
    NXCORR(n, m) = (Nσₙσₘ) Σᵢ (nᵢ - μₙ)(mᵢ - μₘ)
where `n` and `m` are the two images `μ` is the mean and `σ` is the pixelwise standard deviation
"""
function NxCorr(img::T) where {T<:IntensityMap}
    InplaceDivergence(NxCorr(), img)
end

function __divergence(d::NxCorr, img, mimg, fm)
    return -log(abs(nxcorr(img, mimg)))
end

"""
    $(SIGNATURES)

Returns the normalized cross correlation (NXCORR) between `img1` and `img2`.
NXCORR is defined as

    NXCORR(n, m) = (Nσₙσₘ) Σᵢ (nᵢ - μₙ)(mᵢ - μₘ)

where `n` and `m` are the two images `μ` is the mean and `σ` is the pixelwise standard deviation
"""
function nxcorr(img1::IntensityMap{T}, img2::IntensityMap{T}) where {T<:Real}
    m1, s1 = mean_and_std(parent(img1))
    m2, s2 = mean_and_std(parent(img2))
    xcorr =  mapreduce(+, parent(img1), parent(img2); init=zero(T)) do I, J
        (I - m1)*(J - m2)
    end
    # dm1 = img1 .- m1
    # dm2 = img2 .- m2
    # return sum(dm1.*dm2)/(length(img1)*s1*s2)
    return xcorr/((length(img1)-1)*s1*s2)
end
