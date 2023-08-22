"""
    $(TYPEDEF)
Template type for a logarithmic spiral segment

## Fields
$(FIELDS)
"""
struct LogSpiral{T<:Real} <: AbstractImageTemplate
    """ Unit curvature of the logarithmic spiral """
    κ::T
    """ thickness of the Gaussian spiral arm """
    σ::T
    """ Azimuthal extent of the spiral arm """
    δϕ::T
end
CB.radialextent(d::LogSpiral{T}) where {T} = exp(d.κ*d.δϕ)

"""
    $(SIGNATURES)

Constructs a `LogSpiral` template.

## Notes
This is a convienence constructor that uses `VLBISkyModels` `modify`
under the hood. To create this function yourself do

```julia
modify(LogSpiral(κ, σ/r0, δϕ),
       Stretch(r0), Rotate(ξ), Shift(x0, y0))
```

## Arguments
  - `r0`: Radius of the start of the logarithmic spiral
  - `κ`:  Unit curvature of the spiral
  - `σ`:  Thickness of the logarithmic spiral
  - `δϕ`: The angular extent of the spiral
  - `ξ`:  The position angle of the start of the spiral
  - `x0`: The horizontal location of the start of the spiral
  - `y0`: The vertical location of the start of the spiral

"""
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
    α = ringphase(X,Y)

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

## Fields
$(FIELDS)

"""
struct Constant{T} <: AbstractImageTemplate
    """
    Sets the angular scale of the image. This is usually equal to the image FOV
    """
    scale::T
end
@inline CB.intensity_point(c::Constant{T}, p) where {T} = inv(c.scale)^2
CB.radialextent(::Constant{T}) where {T} = one(T)

"""
    $(TYPEDEF)
A smoothed disk model

## Details
Defines a template for an image that has a smoothed disk model.

## Fields
$(FIELDS)

"""
struct GaussDisk{T} <: AbstractImageTemplate
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

"""
    $(SIGNATURES)

Constructs a disk template with a Gaussian fall off.

## Notes
This is a convienence constructor that uses `VLBISkyModels` `modify`
under the hood. To create this function yourself do

```julia
modify(GaussDisk(σ/r0), Stretch(r0), Shift(x0, y0))
```

## Arguments
  - `r0`: Radius of flat top of the disk.
  - `σ`:  Standard deviation of the Gaussian of the disk edge.
  - `x0`: The horizontal location of the disk center
  - `y0`: The vertical location of the disk center

"""
GaussDisk(r0, σ, x0, y0) = modify(GaussDisk(σ/r0), Stretch(r0), Shift(x0, y0))
