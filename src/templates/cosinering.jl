export CosineRing, CosineRingwFloor, CosineRingwGFloor, GaussDisk

"""
    $(TYPEDEF)
A flexible symmetric ring model. The thickness is modeled as a cosine
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
struct CosineRing{N,M,T} <: AbstractImageTemplate
    """
    0th order ring thickness of the Gaussian ring
    """
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
CB.radialextent(d::CosineRing{N,M,T}) where {N,M,T} = one(T) + 3*(d.σ0 + maximum(d.σ))
CB.radialextent(d::CosineRing{0,M,T}) where {M,T} = one(T) + 3*(d.σ0)

function CosineRing(
    r0,
    σ0,
    σ::NTuple{N}, ξσ::NTuple{N},
    s::NTuple{M}, ξs::NTuple{M}, x0, y0) where {N, M}
    return modify(CosineRing(σ0/r0, σ./r0, ξσ, s, ξs), Stretch(r0), Shift(x0, y0))
end


function CB.intensity_point(θ::CosineRing{N,M,T}, p) where {N,M,T}
    (;X, Y) = p
    ex = X
    ey = Y
    ϕ = ringphase(X,Y)
    r = sqrt(ex^2 + ey^2)
    dr2 = (r-1)^2
    #construct the slash
    n = one(T)
    @inbounds for i in 1:M
        n -= θ.s[i]*cos(i*(ϕ - θ.ξs[i]))
    end

    σ = θ.σ0
    @inbounds for i in 1:N
        σ += θ.σ[i]*cos(i*(ϕ - θ.ξσ[i]))
    end

    return abs(n)*exp(-dr2/(2*σ^2))
end




"""
    $(SIGNATURES)

A Cosine ring with a [`GaussDisk`](@ref) at the center of the ring that matches the
`r0` and `σ0` of the CosineRing.


This is a convienence constructor of the more basic VLBISkyModel image constructors.
It is equivalent to

```julia-repl
julia> CosineRing(r0, σ0, σ, ξσ, s, ξs, x0, y0) + floor*modify(GaussDisk(σ0/r0), Stretch(r0), Shift(x0, y0))
```

### Details
Adds ellipticity to the ring. The radius `r0` of the ring is now defined as the
geometric mean of the semi-major (a) and semi-minor (b) axis lengths i.e.

r0 = √(a*b).

The ellipticity `τ` is given by τ = 1-b/a.

## Arguments

 - `r0`: geometric mean radius (√ab) of the ring
 - `σ0`: standard deviation of the Gaussian ring
 - `σ` : amplitudes of the `N` order cosine expansion of the ring thickness
 - `ξσ`: phase of the `N` order cosine expansion of the ring thickness
 - `s` : amplitudes of the `M` order cosine expansion of the azimuthal brightness
 - `ξs`: phase of the `M` order cosine expansion of the azimuthal brightness
 - `floor`: The flux of the center Gaussian. This is relative to the CosineRing.
 - `x0`: location of the ring center horizontally
 - `y0`: location of the ring center vertically


"""
function CosineRingwFloor(
    r0,
    σ0,
    σ::NTuple{N}, ξσ::NTuple{N},
    s::NTuple{M}, ξs::NTuple{M},
    floor, x0, y0) where {N, M}
    return CosineRing(r0, σ0, σ, ξσ, s, ξs, x0, y0) + floor*modify(GaussDisk(σ0/r0), Stretch(r0), Shift(x0, y0))
end

"""
    $(SIGNATURES)

A Cosine ring with a Gaussian blob at the center of the ring.
This is a convienence constructor of the more basic VLBISkyModel image constructors.
It is equivalent to

```julia-repl
julia> CosineRing(r0, σ0, σ, ξσ, s, ξs, x0, y0) + floor*modify(Gaussian(), Stretch(σG), Shift(x0, y0))
```

### Details
Adds ellipticity to the ring. The radius `r0` of the ring is now defined as the
geometric mean of the semi-major (a) and semi-minor (b) axis lengths i.e.

r0 = √(a*b).

The ellipticity `τ` is given by τ = 1-b/a.

## Arguments
 - `r0`: geometric mean radius (√ab) of the ring
 - `σ0`: standard deviation of the Gaussian ring
 - `σ` : amplitudes of the `N` order cosine expansion of the ring thickness
 - `ξσ`: phase of the `N` order cosine expansion of the ring thickness
 - `s` : amplitudes of the `M` order cosine expansion of the azimuthal brightness
 - `ξs`: phase of the `M` order cosine expansion of the azimuthal brightness
 - `floor`: The flux of the center Gaussian. This is relative to the CosineRing.
 - `σG`: The size of the central Gaussian.
 - `x0`: location of the ring center horizontally
 - `y0`: location of the ring center vertically

"""
function CosineRingwGFloor(
        r0,
        σ0,
        σ::NTuple{N}, ξσ::NTuple{N},
        s::NTuple{M}, ξs::NTuple{M}, floor, σG, x0, y0) where {N, M}
    return CosineRing(r0, σ0, σ, ξσ, s, ξs, x0, y0) +
           floor*modify(Gaussian(), Stretch(σG), Shift(x0, y0))
end


"""
    $(SIGNATURES)

An Elliptical Cosine ring. This is a convienence constructor of the more basic
VLBISkyModel image constructors. It is equivalent to

```julia-repl
julia> modify(CosineRing(σ0/r0, σ/r0, ξσ, s, ξs .- ξτ),
            Stretch(r0*sqrt(1-τ), r0/sqrt(1-τ)),
            Rotate(ξτ),
            Shift(x0, y0)
        )
```

### Details
Adds ellipticity to the ring. The radius `r0` of the ring is now defined as the
geometric mean of the semi-major (a) and semi-minor (b) axis lengths i.e.

r0 = √(a*b).

The ellipticity `τ` is given by τ = 1-b/a.

## Arguments
 - `r0`: geometric mean radius (√ab) of the ring
 - `σ0`: standard deviation of the Gaussian ring
 - `σ` : amplitudes of the `N` order cosine expansion of the ring thickness
 - `ξσ`: phase of the `N` order cosine expansion of the ring thickness
 - `τ` : asymmetry of the ring τ = 1-b/a
 - `ξτ`: semi-major axis measured east of north
 - `s` : amplitudes of the `M` order cosine expansion of the azimuthal brightness
 - `ξs`: phase of the `M` order cosine expansion of the azimuthal brightness
 - `x0`: horizontal location of the ring center
 - `y0`: vertical location of the ring center vertically

"""
function EllipticalCosineRing(
    r0,
    σ0,
    σ::NTuple{N}, ξσ::NTuple{N},
    τ, ξτ,
    s::NTuple{M}, ξs::NTuple{M}, x0, y0) where {N, M}
    return modify(
            CosineRing(σ0/r0, σ./r0, ξσ .- ξτ, s, ξs .- ξτ),
            Stretch(r0*sqrt(1-τ), r0/sqrt(1-τ)),
            Rotate(ξτ),
            Shift(x0, y0))
end
