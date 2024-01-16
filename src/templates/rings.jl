export RingTemplate,
       RadialGaussian, RadialDblPower, RadialTruncExp,
       AzimuthalUniform, AzimuthalCosine


# Convience function that creates a phase that matches the EHT convention
ringphase(x, y) = atan(-x,-y)

abstract type AbstractRadial end
abstract type AbstractAzimuthal end

"""
    $(TYPEDEF)

A flexible ring template that forms a ring by taking the product of
a radial and azimuthal brightness profile.

A list of radial profiles is given by `subtypes(AbstractRadial)`

A list of azimuthal profiles is given by `subtypes(AbstractAzimuthal)`

## Examples
```julia-repl
julia> rad = RadialGaussian(0.1)
julia> azi = AzimuthalUniform()
julia> ring = modify(RingTemplate(rad, azi), Stretch(10.0), Shift(1.0, 2.0))
```

## Fields
$(FIELDS)
"""
struct RingTemplate{R<:AbstractRadial, A<:AbstractAzimuthal} <: AbstractImageTemplate
    """
    Radial profile of the ring
    """
    radial::R
    """
    Azimuthal profile of the ring
    """
    azimuthal::A
end

@inline function CB.intensity_point(d::RingTemplate, p)
    (;X, Y) = p
    r = hypot(X, Y)
    ϕ = ringphase(X, Y)
    fr = radial_profile(d.radial, r)
    fϕ = azimuthal_profile(d.azimuthal, ϕ)
    return fr*fϕ
end

@inline CB.radialextent(d::RingTemplate) = CB.radialextent(d.radial)


"""
    RadialGaussian(σ)

Create a profile that is a radial gaussian ring with radius unity
and standard deviation `σ`.

## Notes
This is usually couple with a azimuthal profile to create a general ring template

```julia-repl
julia> rad = RadialGaussian(0.1)
julia> azi = AzimuthalUniform()
julia> t = RingTemplate(rad, azi)
```

## Arguments
  - `σ`: The standard deviation for the Gaussian ring.

"""
struct RadialGaussian{T} <: AbstractRadial
    σ::T
end

@inline function radial_profile(d::RadialGaussian, r)
    return exp(-(r - 1)^2/(2*d.σ^2))
end

CB.radialextent(d::RadialGaussian) = 1 + 3*d.σ


"""
    RadialDblPower(αinner, αouter)

Radial profile that given by a double power-law with the function form

```julia
    r^αinner*inv(1+ r^(αinner + αouter + 1))
```

## Notes
This is usually couple with a azimuthal profile to create a general ring template

```julia-repl
julia> rad = RadialDblPower(3.0, 3.0)
julia> azi = AzimuthalUniform()
julia> t = RingTemplate(rad, azi)
```

## Arguments
  - `αinner` the power law index for `r<1`.
  - `αouter` the power law index for `r≥1`.

"""
struct RadialDblPower{T} <: AbstractRadial
    αinner::T
    αouter::T
end

CB.radialextent(d::RadialDblPower{T}) where {T} = 1 + convert(T, 50^(1/d.αouter))


@inline function radial_profile(d::RadialDblPower, r)
    (;αinner, αouter,) = d
    return r^αinner/(1 + r^(αouter + αinner + 1))
end


"""
    RadialTruncExp(σ)

Radial profile that has a exponential profile up to a cutoff radius of 1 where for `r<1`
the profile is identically zero.

This evalutes to
```julia
    exp(-(r-1)/σ)
```
## Notes
This is usually couple with a azimuthal profile to create a general ring template

```julia-repl
julia> rad = RadialTruncExp(2.0)
julia> azi = AzimuthalUniform()
julia> t = RingTemplate(rad, azi)
```

## Arguments
  - `σ`: Exponential inverse fall off parameter.

"""
struct RadialTruncExp{T} <: AbstractRadial
    σ::T
end

function radial_profile(d::RadialTruncExp{T}, r) where {T}
    r < 1 && return zero(T)
    return exp(-(r - 1)/(d.σ))
end

CB.radialextent(d::RadialTruncExp{T}) where {T} = 1 + log(convert(T, 50))*d.σ

struct AzimuthalUniform{T} <: AbstractAzimuthal end

"""
    AzimuthalUniform()
    AzimuthalUniform{T}()

A azimuthal profile that is uniform for all angles.

## Notes
This is usually couple with a radial profile to create a general ring template

```julia-repl
julia> rad = RadialDblPower(3.0, 3.0)
julia> azi = AzimuthalUniform() # Defaults to Float64
julia> t = RingTemplate(rad, azi)
```
"""
AzimuthalUniform() = AzimuthalUniform{Float64}()

@inline azimuthal_profile(::AzimuthalUniform{T}, ϕ) where {T} = one(T)/2π


"""
    AzimuthalCosine(s::NTuple{N,T}, ξ::NTuple{N, T}) where {N, T}

A azimuthal profile that is  given by a cosine expansion of order `N`.
The expansion is given by

```
    1 - ∑ₙ sₙcos(nϕ - ξₙ)
```

## Notes
This is usually couple with a radial profile to create a general ring template

```julia-repl
julia> rad = RadialDblPower(3.0, 3.0)
julia> azi = AzimuthalCosine((0.5, 0.2), (0.0, π/4)) # Defaults to Float64
julia> t = RingTemplate(rad, azi)
```

## Arguments
  - `s` : amplitudes of the `N` order cosine expansion of the azimuthal brightness
  - `ξs`: phase of the `N` order cosine expansion of the azimuthal brightness


"""
struct AzimuthalCosine{T, N} <: AbstractAzimuthal
    s::NTuple{N,T}
    ξs::NTuple{N,T}
end

AzimuthalCosine(s::Real, ξs::Real) = AzimuthalCosine((s, ), (ξs,))

@inline function azimuthal_profile(d::AzimuthalCosine{T, N}, ϕ) where {T,N}
    (;s, ξs) = d
    mapreduce(+, 1:N; init=one(T)) do n
        return -s[n]*cos(n*ϕ - ξs[n])
    end
end

"""
    $(SIGNATURES)

Implements the gaussian ring template.

## Notes
This is a convienence constructor that uses [`RingTemplate`](@ref) under
the hood. To create this function your self do

```julia
RingTemplate(RadialGaussian(σ), AzimuthalUniform())
```

## Arguments
- `σ` : standard deviation of the Gaussian ring

"""
@inline GaussianRing(σ) = RingTemplate(RadialGaussian(σ), AzimuthalUniform())

"""
    $(SIGNATURES)

Implements the gaussian ring template.

## Notes
This is a convienence constructor that uses [`RingTemplate`](@ref) under
the hood. To create this function your self do

```julia
modify(GaussianRing(σ/r0), Stretch(r0), Shift(x0, y0))
```

## Arguments
- `r0`: radius of the ring
- `σ` : standard deviation of the Gaussian ring
- `x0`: location of the ring center horizontally
- `y0`: location of the ring center vertically


"""
@inline GaussianRing(r0, σ, x0, y0) = modify(GaussianRing(σ/r0), Stretch(r0), Shift(x0, y0))


"""
    $(SIGNATURES)

Implements the slashed gaussian ring template, that uses a cosine
to implement the slash and has a brightness position angle of 0 degrees
east of north.

## Notes
This is a convienence constructor that uses [`RingTemplate`](@ref) under
the hood. To create this function your self do

## Arguments
- `σ` : standard deviation of the Gaussian ring
- `s` : Slash amplitude. 0 means no slash, and 1 is maximal.


```julia
RingTemplate(RadialGaussian(σ), AzimuthalCosine((s,), (zero(s),)))
```

"""
@inline SlashedGaussianRing(σ, s) = RingTemplate(RadialGaussian(σ), AzimuthalCosine((s,), (zero(s),)))


"""
    $(SIGNATURES)

Implements the slashed gaussian ring template, that uses a cosine
to implement the slash.

## Arguments
- `r0`: radius of the ring
- `σ` : standard deviation of the Gaussian ring
- `s` : Slash amplitude. 0 means no slash, and 1 is maximal.
- `ξs`: azimuthal peak brightness P.A. measured east of north
- `x0`: location of the ring center horizontally
- `y0`: location of the ring center vertically

## Notes
This is a convienence constructor that uses [`RingTemplate`](@ref) under
the hood. To create this function your self do

```julia
modify(SlashedGaussianRing(σ/r0, s), Stretch(r0), Rotate(ξ), Shift(x0, y0))
```

"""
SlashedGaussianRing(r0, σ, s, ξ, x0, y0) = modify(SlashedGaussianRing(σ/r0, s), Stretch(r0), Rotate(ξ), Shift(x0, y0))



"""
    $(SIGNATURES)
Implements the elliptical gaussian ring template with
semi-minor axis `a` and semi-major axis `b`.


## Arguments
 - `r0`: geometric mean radius (√ab) of the ring
 - `σ` : standard deviation of the Gaussian ring
 - `τ` : asymmetry of the ring τ = 1-b/a
 - `ξτ`: semi-major axis measured east of north
 - `x0`: location of the ring center horizontally
 - `y0`: location of the ring center vertically

 ## Notes
 This is a convienence constructor that uses [`RingTemplate`](@ref) under
 the hood. To create this function your self do

 ```julia-repl
 julia> modify(GaussianRing(σ/r0), Stretch(r0*sqrt(1-τ), r0/sqrt(1-τ)), Rotate(ξτ), Shift(x0, y0))
 ```

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

The brightness asymmetry uses a cosine to implement the slash.

## Arguments
 - `r0`: geometric mean radius (√ab) of the Gaussian ring
 - `σ` : standard deviation of the Gaussian ring
 - `τ` : asymmetry of the Gaussian ring τ = 1-b/a
 - `ξτ`: semi-major axis measured east of north
 - `s` : Slash amplitude. 0 means no slash, and 1 is maximal.
 - `ξs`: azimuthal peak brightness P.A. measured east of north
 - `x0`: location of the ring center horizontally
 - `y0`: location of the ring center vertically

 ## Notes
 This is a convienence constructor that uses [`RingTemplate`](@ref) under
 the hood. To create this function your self do

 ```julia
 modify(SlashedGaussianRing(σ/r0, s),
            Rotate(ξs-ξτ),
            Stretch(r0*sqrt(1-τ), r0/sqrt(1-τ)),
            Rotate(ξτ),
            Shift(x0, y0)
        )
 ```


"""
function EllipticalSlashedGaussianRing(r0, σ, τ, ξτ, s, ξs, x0, y0)
    return modify(SlashedGaussianRing(σ/r0, s),
        Rotate(ξs-ξτ),
        Stretch(r0*sqrt(1-τ), r0/sqrt(1-τ)),
        Rotate(ξτ),
        Shift(x0, y0)
        )
end
