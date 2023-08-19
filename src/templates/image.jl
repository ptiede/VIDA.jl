


"""
    $(TYPEDEF)
Extrememly flexible symmetric ring model. The thickness is modeled as a cosine
expansion with `N` terms and the slash by a expansion with `M` terms.


# Details
The ring is forced to be symmetric for a significant speed boost over CosineRing.
The thickness of the ring is modeled by a cosine expansion in azimuthal
angle. `N` specifies the number of cosine modes to fit, where the first
mode is the constant thickness portion and so has no corresponding angle.
The slash is modeled as a separate cosine expansion, with `M` terms.
Here the zero order term is forced to be unity, so `M` defines the `M`
additional terms.

"""
Base.@kwdef struct SymCosineRing{N,M,T} <: AbstractImageTemplate
    σ0::T
    """Standard deviation expansion (length N) of the width of the Gaussian ring"""
    σ::NTuple{N, T}
    """Orientations of the cosine expansion width (length N)"""
    ξσ::NTuple{N, T}
    """Slash of Gaussian ring (length M)."""
    s::NTuple{M, T}
    """Slash orientations (length M) in radians measured north of east"""
    ξs::NTuple{M, T}
end
CB.radialextent(d::SymCosineRing{N,M,T}) where {N,M,T} = one(T) + 3*d.σ[1]

function SymCosineRing(
    r0,
    σ0,
    σ::NTuple{N}, ξσ::NTuple{N},
    s::NTuple{M}, ξs::NTuple{M}, x0, y0) where {N, M}
    return modify(SymCosineRing(σ0/r0, σ/r0, ξσ, s, ξs), Stretch(r0), Shift(x0, y0))
end


function CB.intensity_point(θ::SymCosineRing{N,M,T}, p) where {N,M,T}
    (;X, Y) = p
    ex = X
    ey = Y
    ϕ = atan(-ey,-ex)
    r = sqrt(ex^2 + ey^2)
    dr2 = (r-1)^2
    #construct the slash
    n = one(T)
    @inbounds for i in 1:M
        n -= θ.s[i]*cos(i*(ϕ - θ.ξs[i]))
    end

    σ = θ.σ[1]
    @inbounds for i in 1:N
        σ += θ.σ[i+1]*cos(i*(ϕ - θ.ξσ[i]))
    end

    return abs(n)*exp(-dr2/(2*σ^2))
end




"""
    $(SIGNATURES)
Symmetric ring model of radius `r0`  whose azimuthal brightness and thickness is described using a
`M` and `N` order cosine expansion. In addition a disk with radius matching the diameter
of the ring with flux `floor` is added to the center of the ring.

# Details
The thickness of the ring is modeled by a cosine expansion in azimuthal
angle. `N` specifies the number of cosine modes to fit, where the first
mode is the constant thickness portion and so has no corresponding angle.
The slash is modeled as a separate cosine expansion, with `M` terms.
Here the zero order term is forced to be unity, so `M` defines the `M`
additional terms.

"""
function SymCosineRingwFloor(
    r0,
    σ0,
    σ::NTuple{N}, ξσ::NTuple{N},
    s::NTuple{M}, ξs::NTuple{M},
    floor, x0, y0) where {N, M}
    return modify(SymCosineRing(σ0/r0, σ/r0, ξσ, s, ξs), Stretch(r0), Shift(x0, y0)) +
           floor*modify(GaussDisk(), Stretch(r0), Shift(x0, y0))
end

"""
    $(SIGNATURES)
Symmetric ring model of radius `r0` whose azimuthal brightness and thickness is described using a
`M` and `N` order cosine expansion. In addition a Gaussian component with thickness
`σG` and flux fraction `floor` is added to the center of the ring.

# Details
The thickness of the ring is modeled by a cosine expansion in azimuthal
angle. `N` specifies the number of cosine modes to fit, where the first
mode is the constant thickness portion and so has no corresponding angle.
The slash is modeled as a separate cosine expansion, with `M` terms.
Here the zero order term is forced to be unity, so `M` defines the `M`
additional terms.

"""
function SymCosineRingwGFloor(
        r0,
        σ0,
        σ::NTuple{N}, ξσ::NTuple{N},
        s::NTuple{M}, ξs::NTuple{M}, floor, σG, x0, y0) where {N, M}
    return modify(SymCosineRing(σ0/r0, σ/r0, ξσ, s, ξs), Stretch(r0), Shift(x0, y0)) +
           floor*modify(Gaussian(), Stretch(σG/r0), Shift(x0, y0))
end


"""
    $(TYPEDEF)
Template type for a logarithmic spiral segment

## Fields
$(FIELDS)
"""
Base.@kwdef struct LogSpiral{T<:Real} <: AbstractImageTemplate
    """ Unit curvature of the logarithmic spiral """
    κ::T
    """ thickness of the Gaussian spiral arm """
    σ::T
    """ Azimuthal extent of the spiral arm """
    δϕ::T
end
CB.radialextent(d::LogSpiral{T}) where {T} = d.r0*exp(d.κ*d.δϕ)

function LogSpiral(r0, κ, σ, δϕ, ξ, x0, y0)
    return modify(LogSpiral(κ, σ/r0, δϕ), Stretch(r0), Rotate(ξ), Shift(x0, y0))
end


@inline function CB.intensity_point(θ::LogSpiral, p)
    (;X, Y) = p
    (;κ, σ, δϕ) = θ
    #Set up the spiral
    k = sqrt(1-κ*κ)/κ
    rc = exp(k*10π) #This finds where we should start our spiral arm from
    a = inv(rc) #Get on the correct logspiral

    r = hypot(X, Y)
    α = atan(Y,X)

    #Now I need to find the distance from the closest spiral arm
    n = (log(r/a)/k - α)/(2π)
    nc = ceil(n)
    nf = floor(n)
    rc = a*exp(k*(α + nc*2π))
    rf = a*exp(k*(α + nf*2π))
    r1,r2 = abs(rc-r),abs(rf-r)
    if r1 < r2
        nn = nc
        dist = r1
    else
        nn = nf
        dist = r2
    end
    #Get the angular extent
    dtheta = (10π - (α + nn*2π))
    return exp(-dist^2/(2*σ^2) -dtheta^2/(2*(δϕ/2)^2))
end

#Load the templates
"""
    $(TYPEDEF)
An constant template.

### Details
Defines an image that just has constant flux. This is very useful for soaking up
low levels of flux in image reconstructions that can bias the results.

Since images or normalized to unity, this means the `Constant` template has no
additional parameters.
"""
struct Constant{T} <: AbstractImageTemplate
    scale::T
end
@inline CB.intensity_point(c::Constant{T}, p) where {T} = inv(c.scale)^2
CB.radialextent(::Constant{T}) where {T} = one(T)

"""
    $(TYPEDEF)
A smoothed disk model

### Details
Defines a template for an image that has a smoothed disk model.


"""
Base.@kwdef struct GaussDisk{T} <: AbstractImageTemplate
    """
    Disk edge standard deviation
    """
    α::T
end

@inline function CB.intensity_point(θ::GaussDisk{T}, p) where {T}
    (;α) = θ
    r = hypot(p.X, p.Y)
    if ( r < 1)
        return one(T)
    else
        return exp(-(r-1)^2/(2*α^2))
    end
end
CB.radialextent(d::GaussDisk{T}) where {T} = one(T) + 3*d.α

GaussDisk(r0, σ, x0, y0) = modify(GaussDisk(σ/r0), Stretch(r0), Shift(x0, y0))


"""
    $(TYPEDEF)
Symmetric gaussian ring template. This is the most basic template and just attempts
to recover a location `x0`,`y0`, radius `r0` and thickness `σ` from some image.
### Fields
    $(FIELDS)
### Example
```julia
GaussianRing(5.0)
```
"""
Base.@kwdef struct GaussianRing{T} <: AbstractImageTemplate
    """Standard deviation of Gaussian ring"""
    σ::T #standard deviation of the Gaussian ring
    function GaussianRing(σ)
        @assert σ>0 "Ring thickness must be positive"
        new{typeof(σ)}(σ)
    end
end

 @inline function CB.intensity_point(θ::GaussianRing, p)
    r = hypot(p.X, p.Y)
    return exp(-(r - 1)^2/(2*θ.σ^2))
end
CB.radialextent(θ::GaussianRing{T}) where {T} = one(T) + 3*θ.σ

GaussianRing(r0, σ, x0, y0) = modify(GaussianRing(σ/r0), Stretch(r0), Shift(x0, y0))


"""
    $(TYPEDEF)
Implements the slashed gaussian ring template, that uses a cosine
to symmetrically implement the slash. While this is marginally more
complicated that a linear slash, it has a number of benefits such as
mainting the azimuthal and smooth structure of the image.

### Fields
    $(FIELDS)

"""
Base.@kwdef struct SlashedGaussianRing{T} <: AbstractImageTemplate
    """Standard deviation of the unit Gaussian ring"""
    σ::T
    """Slash strength of Gaussiang ring. 0 means no slash"""
    s::T
    function SlashedGaussianRing(σ,s)
        @assert σ>0 "SlashedGaussianRing: Ring thickness must be positive"
        @assert s<=1 && s>=0 "SlashedGaussianRing: Slash strength, $s, is bounded in [0,1]"
        T = promote_type(typeof(σ), typeof(s))
        new{T}(T(σ), T(s))
    end
end

CB.radialextent(d::SlashedGaussianRing{T}) where {T} = one(T) + 3*d.σ
SlashedGaussianRing(r0, σ, s, ξ, x0, y0) = modify(SlashedGaussianRing(σ/r0, s), Stretch(r0), Rotate(ξ), Shift(x0, y0))

#Template function
@inline function CB.intensity_point(θ::SlashedGaussianRing, p)
    (;X, Y) = p
    r = hypot(X, Y)

    ϕ = atan(X,-Y)
    n = (1-θ.s*cos(ϕ))/(θ.s + 1)

    return n*exp(-(r-1)^2/(2*θ.σ^2))
end

"""
    $(SIGNATURES)
Implements the elliptical gaussian ring template. Where the ellipticity `τ` is defined
as one minus ratio between the semi-minor and semi-major axis.

## Details
Adds ellipticity to the ring. The radius `r0` of the ring is now defined as the
geometric mean of the semi-major (a) and semi-minor (b) axis lengths i.e.

r0 = √(a*b).

The ellipticity `τ` is given by τ = 1-b/a.

## Arguments
 - `r0`: geometric mean radius (√ab) of the Gaussian ring
 - `σ` : standard deviation of the Gaussian ring
 - `τ` : asymmetry of the Gaussian ring τ = 1-b/a
 - `ξτ`: semi-major axis measured east of north
 - `x0`: location of the ring center horizontally
 - `y0`: location of the ring center vertically
"""
function EllipticalGaussianRing(r0, σ, τ, ξτ, x0, y0)
    return modify(GaussianRing(σ/r0),
        Stretch(r0*sqrt(1-τ), r0/sqrt(1-τ)),
        Rotate(ξτ),
        Shift(x0, y0))
end

"""
    $(SIGNATURES)

Combination of the elliptical and slashed gaussian ring.

### Details
Adds ellipticity to the ring. The radius `r0` of the ring is now defined as the
geometric mean of the semi-major (a) and semi-minor (b) axis lengths i.e.

r0 = √(a*b).

The ellipticity `τ` is given by τ = 1-b/a.

## Arguments
 - `r0`: geometric mean radius (√ab) of the Gaussian ring
 - `σ` : standard deviation of the Gaussian ring
 - `τ` : asymmetry of the Gaussian ring τ = 1-b/a
 - `ξτ`: semi-major axis measured east of north
 - `s` : Slash amplitude. 0 means no slash, and 1 is maximal.
 - `ξs`: azimuthal peak brightness P.A. measured east of north
 - `x0`: location of the ring center horizontally
 - `y0`: location of the ring center vertically

"""
function EllipticalSlashedGaussianRing(r0, σ, τ, ξτ, s, ξs, x0, y0)
    return modify(SlashedGaussianRing(σ/r0, s),
        Rotate(ξs-ξτ),
        Stretch(r0*sqrt(1-τ), r0/sqrt(1-τ)),
        Rotate(ξτ),
        Shift(x0, y0)
        )
end

function EllipticalCosineRing(
    r0,
    σ0,
    σ::NTuple{N}, ξσ::NTuple{N},
    τ, ξτ,
    s::NTuple{M}, ξs::NTuple{M}, x0, y0) where {N, M}
    return modify(
            SymCosineRing(σ0/r0, σ/r0, ξσ, s, ξs .- ξτ),
            Stretch(r0*sqrt(1-τ), r0/sqrt(1-τ)),
            Rotate(ξτ),
            Shift(x0, y0))
end
