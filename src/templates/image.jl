@doc """
    $(SIGNATURES)
Unpacks the parameters of the template `θ`

Returns the parameters in a vector.
"""
function unpack(θinit::T) where {T<:AbstractImageTemplate}
    n = Base.size(T)
    fields = propertynames(θinit)
    p = zeros(n)
    for i in 1:n
        p[i] = getproperty(θinit,fields[i])
    end
    return p
end


@inline function imagecenter(θ::AbstractImageTemplate)
    return θ.x0, θ.y0
end


@doc """
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
@with_kw struct SymCosineRing{N,M} <: AbstractImageTemplate
    """Radius of the Gaussian ring"""
    r0::Float64
    """Standard deviations (length N+1) of the width of the Gaussian ring"""
    σ::Vector{Float64}
    """Orientations of the cosine expansion width (length N)"""
    ξσ::Vector{Float64}
    """Slash of Gaussian ring (length M)."""
    s::Vector{Float64}
    """Slash orientations (length M) in radians measured north of east"""
    ξs::Vector{Float64}
    """x position of the center of the ring in μas"""
    x0::Float64
    """y position of the center of the ring in μas"""
    y0::Float64
    function SymCosineRing{N,M}(r0,σ, ξσ, s, ξs,x0,y0) where {N, M}
        #@assert N isa Integer
        #@assert M isa Integer
        new{N,M}(float(r0),σ, ξσ, s, ξs,float(x0),float(y0))
    end
end
@doc """
    CosineRing{N,M}(p::AbstractArray) where {N,M}
Takes in a vector of paramters describing the template.
# Details
The order of the vector must be
 - p[1] = `r0`
 - p[2:(N+1)] = `σ`
 - p[(N+2):(2N)] = `ξσ`
 - p[2N+1] = `τ`
 - p[2N+2] = `ξτ`
 - p[2N+3:2N+M+2] = `s`
 - p[2N+3+M:2N+2+2M] = `ξs`
 - p[2N+3+2M] = `x0`
 - p[2N+4+2M] = `y0`
"""
function SymCosineRing{N,M}(p::AbstractArray) where {N,M}
    #@assert length(p) == size(CosineRing{N,M})
    SymCosineRing{N,M}(p[1], p[2:(N+2)],
                    p[(N+3):(2N+2)],
                    p[2N+3:2N+2+M], p[2N+3+M:2N+2+2M],
                    p[2N+3+2M], p[2N+4+2M]
    )
end
Base.size(::Type{SymCosineRing{N,M}}) where {N, M} = 3 + N+1 + N + 2*M

@inline function (θ::SymCosineRing{N,M})(x,y) where {N, M}
    ex = x-θ.x0
    ey = y-θ.y0
    ϕ = atan(-ey,-ex)
    r = sqrt(ex^2 + ey^2)
    dr2 = (r-θ.r0)^2
    #construct the slash
    n = one(θ.r0)
    @inbounds for i in 1:M
        n -= θ.s[i]*cos(i*(ϕ - θ.ξs[i]))
    end

    σ = θ.σ[1]
    @inbounds for i in 1:N
        σ += θ.σ[i+1]*cos(i*(ϕ - θ.ξσ[i]))
    end

    return abs(n)*exp(-dr2/(2.0*σ^2+1e-2))
end

function unpack(θ::SymCosineRing{N,M}) where {N,M}
    n = Base.size(typeof(θ))
    p = zeros(n)
    p[1] = θ.r0
    p[2:(N+2)] = θ.σ
    if N>0
        p[(N+3):(2N+2)] = θ.ξσ
    end
    p[2N+3:2N+2+M] = θ.s
    p[2N+3+M:2N+2+2M] = θ.ξs
    p[2N+3+2M] = θ.x0
    p[2N+4+2M] = θ.y0

    return p
end

@doc """
    $(SIGNATURES)
Template type for a logarithmic spiral segment

## Fields
$(FIELDS)
"""
@with_kw struct LogSpiral{T<:Real} <: AbstractImageTemplate
    """ Radius of the spiral peak brightness """
    r0::T
    """ Unit curvature of the logarithmic spiral """
    κ::T
    """ thickness of the Gaussian spiral arm """
    σ::T
    """ Azimuthal extent of the spiral arm """
    δϕ::T
    """ peak brightness location """
    ξ::T
    """ x location of disk center in μas """
    x0::T
    """ y location of disk center in μas """
    y0::T
end
function LogSpiral(p::Vector{T}) where {T<:Real}
    @assert length(p) == 7
    LogSpiral{T}(p[1],p[2],p[3],p[4],p[5],p[6],p[7])
end
Base.size(::Type{LogSpiral{T}}) where {T} = 7


@inline function (θ::LogSpiral)(x,y)
    @unpack κ, σ, r0, δϕ, ξ, x0, y0 = θ
    x′,y′ = x-x0,y-y0
    #Set up the spiral
    k = sqrt(1-κ*κ)/κ
    rc = exp(k*10π) #This finds where we should start our spiral arm from
    a = r0/rc #Get on the correct logspiral

    r = hypot(x′,y′)
    α = (atan(y′,x′)) - ξ

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

"""
    $(SIGNATURES)
Type for an image template. This takes an EHTImage or any image object
and creates a template out of it. The parameters of the image are
the center of the image `x0`, `y0`.

# Fields
$(FIELDS)
"""
@with_kw struct ImageTemplate{T} <: AbstractImageTemplate
    x0::Float64
    y0::Float64
    itp::T
end
function ImageTemplate(x0, y0, img::EHTImage; interp=BSpline(Cubic(Line(OnGrid()))))
    itp = image_interpolate(img, interp)
    return ImageTemplate(float(x0), float(y0), itp)
end

@inline function Base.size(::Type{ImageTemplate{T}}) where {T<:Interpolations.AbstractInterpolation}
    return 2
end
function (θ::ImageTemplate{T})(x,y) where {T<:Interpolations.AbstractInterpolation}
    @unpack x0, y0 = θ
    itp = getfield(θ, :itp)
    return itp(y-y0, x-x0)
end

function _update(θ::ImageTemplate, p)
    return typeof(θ)(p[1], p[2], getfield(θ, :itp))
end

function Base.getproperty(θ::ImageTemplate, symbol::Symbol)
    names = propertynames(θ)
    if symbol ∈ names
        return getfield(θ, symbol)
    end
    throw("type ImageTemplate has no field $symbol")
end


function Base.propertynames(::ImageTemplate)
    return (:x0, :y0)
end

function unpack(θ::ImageTemplate)
    p = zeros(2)
    p[1] = θ.x0
    p[2] = θ.y0
    return p
end

#Load the templates
@doc """
    $(TYPEDEF)
An constant template.

### Details
Defines an image that just has constant flux. This is very useful for soaking up
low levels of flux in image reconstructions that can bias the results.

Since images or normalized to unity, this means the `Constant` template has no
additional parameters.
"""
struct Constant <: AbstractImageTemplate end
Constant(p) = Constant()
Base.size(::Type{Constant}) = 0
@inline (θ::Constant)(x,y) = 1

@doc """
    $(TYPEDEF)
A smoothed disk model

### Details
Defines a template for an image that has a smoothed disk model.


"""
@with_kw struct Disk <: AbstractImageTemplate
    """
    Radius of the disk
    """
    r0::Float64
    """
    Disk edge standard deviation in μas
    """
    α::Float64
    """
    x location of disk center in μas
    """
    x0::Float64
    """
    y location of disk center in μas
    """
    y0::Float64
end
function Disk(p)
    @assert length(p) == 4 "There are 4 parameters for the Disk template."
    Disk(p[1],p[2],p[3],p[4])
end
 @inline function (θ::Disk)(x,y)
    @unpack r0, α, x0, y0 = θ
    r = sqrt((x-x0)^2 + (y-y0)^2)
    if ( r < r0 )
        return one(eltype(r))
    else
        return exp(-(r-r0)^2/(2.0*α^2)) + 1e-50
    end
end
Base.size(::Type{Disk}) = 4


@doc """
    $(TYPEDEF)
An asymmetric Gaussian blob.

### Details
Defines a asymmetric Gaussian image. This is useful if the image has some
non-ring emission in the and you need to soak up some of the flux.

The parameters of the model follow very closely to those used in Themis.
The Gaussian size `σ` is given by
    σ = √(σxσy)
    τ = 1-σy/σx,
where σx,σy are the semi-major,minor axis lenght respectively. This is similar to
how the asymmetry for the `EllipticalGaussianRing`.

### Fields
$(FIELDS)
"""
@with_kw struct AsymGaussian <: AbstractImageTemplate
    """Gaussian size in μas"""
    σ::Float64
    """Gaussian asymmetry"""
    τ::Float64
    """Gaussian orientation in radians measured north of east"""
    ξ::Float64
    """x position of Gaussian center in μas"""
    x0::Float64
    """y position of Gaussian center in μas"""
    y0::Float64
    function AsymGaussian(σ,τ,ξ,x0,y0)
        @assert τ>=0 && τ < 1 "τ must be in [0,1)"
        new(σ, τ, ξ, x0, y0)
    end
end
function AsymGaussian(p)
    @assert length(p) == 5 "There are 5 parameters for the AsymGaussian template."
    AsymGaussian(p[1],p[2],p[3],p[4],p[5])
end

 @inline function (θ::AsymGaussian)(x,y)
    x′,y′ = rotate(x-θ.x0,y-θ.y0,θ.ξ)
    σx2 = θ.σ *θ.σ/(1.0-θ.τ)
    σy2 = θ.σ*θ.σ*(1.0-θ.τ)
    d2 = x′*x′/σx2 + y′*y′/σy2
    return exp(-0.5*d2) + 1e-50
end
Base.size(::Type{AsymGaussian}) = 5


@doc """
    $(TYPEDEF)
Symmetric gaussian ring template. This is the most basic template and just attempts
to recover a location `x0`,`y0`, radius `r0` and thickness `σ` from some image.
### Fields
    $(FIELDS)
### Example
```julia
GaussianRing(r0=20.0,σ=5.0,x0=0.0,y0=-10.0)
```
"""
@with_kw struct GaussianRing <: AbstractImageTemplate
    """Radius of Gaussian ring in μas"""
    r0::Float64
    """Standard deviation of Gaussian ring in μas"""
    σ::Float64 #standard deviation of the Gaussian ring
    """x location of the center of the Gaussian ring in μas"""
    x0::Float64 #center position of the gaussian ring x componenent
    """y location of the center of the Gaussian ring in μas"""
    y0::Float64 #center position of the gaussian ring y compoenent
    function GaussianRing(r0,σ,x0,y0)
        @assert r0>=0 "Ring radius must be positive"
        @assert σ>0 "Ring thickness must be positive"
        new(float(r0),float(σ),float(x0),float(y0))
    end
end
function GaussianRing(p)
    @assert length(p)==4 "GaussianRing: Template requires 4 parameters"
    GaussianRing(p[1],p[2],p[3],p[4])
end

Base.size(::Type{GaussianRing}) = 4

 @inline function (θ::GaussianRing)(x, y)
    r = sqrt((x - θ.x0)^2 + (y - θ.y0)^2)
    return exp( -(r - θ.r0)^2/(2*θ.σ^2)) + 1e-50
end


"""
    ($TYPEDEF)
Diffuse background template. This consists of a number of Gaussians at
fixed grid locations based on the field of view passed during construction.

There are two ways to construct a DiffuseBack template
```julia
fovx,fovy = 80.0, 90.0
θ = DiffuseBack(10.0, zeros(6,6), fovx, fovy)
```
which construct a 6x6 Gaussian grid over a field of view given by fovx and fovy.

The other way build the fovx based off of an EHTImage object
```julia
θ = DiffuseBack(10.0, zeros(6,6), img)
```
This construct the background to match the field of view of the image object.

# Internals
The Px and Py parameters in the struct are the pixel sizes in the x and y direction,
and passed through constant propogation thanks to Julia.

"""
@with_kw struct DiffuseBack{N,M,Px,Py} <: AbstractImageTemplate
    width::Float64
    intensities::Matrix{Float64}
end

function DiffuseBack(width, intensity::Matrix{Float64}, fovx, fovy)
    N,M = size(intensity)
    dx = abs(fovx/max(M-1,1))
    dy = fovy/max(N-1,1)
    return DiffuseBack{N,M,dx,dy}(width, intensity)
end

function DiffuseBack(width, intensity::Matrix{Float64}, img::EHTImage)
    fovx,fovy = field_of_view(img)
    N,M = size(intensity)
    dx = abs(fovx/max(M-1,1))
    dy = fovy/max(N-1,1)
    return DiffuseBack{N,M,dx,dy}(width, intensity)
end
Base.size(::Type{DiffuseBack{N,M,dx,dy}}) where {N,M,dx,dy} = N*M+1


function DiffuseBack{N,M,Px,Py}(p) where {N,M,Px,Py}
    return DiffuseBack{N,M,Px,Py}(p[1],reshape(@view(p[2:end]),N,M))
end

function _update(θ::DiffuseBack{N,M,dx,dy}, p) where {N,M,dx,dy}
    w = first(p)
    impix = @view p[2:end]
    DiffuseBack{N,M,dx,dy}(w, reshape(impix, N, M))
end

@inline function (θ::DiffuseBack{N,M,Px,Py})(x,y) where {N,M,Px,Py}
    sum = zero(eltype(θ.intensities))
    width = θ.width
    intensities = θ.intensities
    xstart = (-M*Px + Px)/2.0
    ystart = (-N*Py + Py)/2.0

    for i in 1:M
        for j in 1:N
            x0 = xstart + Px*(i-1)
            y0 = ystart + Py*(j-1)
            dr = (x-x0)^2 + (y-y0)^2
            sum += exp(-dr/(2*width^2))*intensities[j,i]
        end
    end
    return sum
end

function unpack(θ::DiffuseBack{N,M,Px,Py}) where {N,M,Px,Py}
    p = zeros(N*M+1)
    p[1] = θ.width
    p[2:end] = θ.intensities
    return p
end






@doc """
    $(TYPEDEF)
Implements the slashed gaussian ring template, that uses a cosine
to symmetrically implement the slash. While this is marginally more
complicated that a linear slash, it has a number of benefits such as
mainting the azimuthal and smooth structure of the image.

### Fields
    $(FIELDS)

"""
@with_kw struct SlashedGaussianRing <: AbstractImageTemplate
    """Radius of the ring in μas"""
    r0::Float64
    """Standard deviation of Gaussian ring in μas"""
    σ::Float64
    """Slash strength of Gaussiang ring. 0 means no slash"""
    s::Float64
    """Rotation angle in radians of slash direction, measured north of west"""
    ξ::Float64
    """x position of the center of the ring in μas"""
    x0::Float64
    """y position of the center of the ring in μas"""
    y0::Float64
    function SlashedGaussianRing(r0,σ,s,ξ,x0,y0)
        @assert r0>0 "SlashedGaussianRing: Ring radius must be positive"
        @assert σ>0 "SlashedGaussianRing: Ring thickness must be positive"
        @assert s<=1 && s>=0 "SlashedGaussianRing: Slash strength, $s, is bounded in [0,1]"
        new(float(r0),float(σ),float(s),float(ξ),float(x0),float(y0))
    end
end

function SlashedGaussianRing(p)
        @assert length(p)==6 "SlashedGaussianRing: Template requires 6 parameters"
        SlashedGaussianRing(p[1],p[2],p[3],p[4],p[5],p[6])
end
Base.size(::Type{SlashedGaussianRing}) = 6
#Template function
 @inline function (θ::SlashedGaussianRing)(x,y)
    ex = x - θ.x0
    ey = y - θ.y0
    r = sqrt((ex)^2 + (ey)^2)

    #rotate the image so slash is on the x axis
    #xrot,yrot = rotate(x-θ.x0,y-θ.y0,θ.ξ)
    #construct the slash
    ϕ = atan(-ey,-ex)
    n = (1-θ.s*cos(θ.ξ-ϕ))/(θ.s + 1)

    return n*exp(-(r-θ.r0)^2/(2*θ.σ^2)) + 1e-50
end


@doc """
    $(TYPEDEF)
Implements the elliptical gaussian ring template. Where the ellipticity `tau` is defined
as one minus ratio between the semi-minor and semi-major axis.

### Details
Adds ellipticity to the ring. The radius `r0` of the ring is now defined as the
geometric mean of the semi-major (a) and semi-minor (b) axis lengths i.e.

r0 = √(a*b).

The ellipticity `τ` is given by τ = 1-b/a.

### Fields
$(FIELDS)

### Notes
There is no normalization since the ellipticity makes it impossible to
normalize analytically. In fact the distance from the ellipse is implemented
numerically using an algorithm adapted from
[git](https://github.com/0xfaded/ellipse_demo/issues/1#issuecomment-405078823)
"""
@with_kw struct EllipticalGaussianRing <: AbstractImageTemplate
    """Radius of the Gaussian ring"""
    r0::Float64 #geometric mean of the semi-major, a, and semi-minor axis, b, r0=√ab
    """Standard deviation of the width of the Gaussian ring"""
    σ::Float64 #standard deviation of the Gaussian ring
    """Asymmetry of the Gaussian ring defined as ``1-b/a``"""
    τ::Float64 #asymmetry of the Gaussian ring τ = 1-b/a
    """Asymmetry orientation in radians measured north of east"""
    ξ::Float64 #slash orientation measured from north of east
    """x position of the center of the ring in μas"""
    x0::Float64
    """y position of the center of the ring in μas"""
    y0::Float64
    function EllipticalGaussianRing(r0,σ,τ,ξ,x0,y0)
        @assert r0>0 "EllipticalGaussianRing: r0 must be positive"
        @assert σ>0 "EllipticalGaussianRing: σ must be positive"
        @assert τ<=1 && τ>=0 "EllipticalGaussianRing: τ must be in [0,1)"
        new(float(r0),float(σ),float(τ),float(ξ),float(x0),float(y0))
    end
end
function EllipticalGaussianRing(p)
    @assert length(p)==6 "EllipticalGaussianRing: Template requires 6 parameters"
    EllipticalGaussianRing(p[1],p[2],p[3],p[4],p[5],p[6])
end
Base.size(::Type{EllipticalGaussianRing}) = 6

 @inline function (θ::EllipticalGaussianRing)(x,y)
    ex = x-θ.x0
    ey = y-θ.y0
    ex′,ey′ = rotate(ex,ey,θ.ξ)
    a = θ.r0/sqrt(1.0-θ.τ)
    b = θ.r0*sqrt(1.0-θ.τ)
    distance = ellipse_sqdist(ex′,ey′,a, b)
    return exp(-distance/(2.0*θ.σ^2)) + 1e-50
end

@doc """
    $(TYPEDEF)
Creates the template from the Paper I am writing. It is a combination of
the elliptical and slashed gaussian ring. The slash and the semi-major axis
are either aligned if the slash parameter `s`>0 or antialigned if `s`<0.

### Details
Adds ellipticity to the ring. The radius `r0` of the ring is now defined as the
geometric mean of the semi-major (a) and semi-minor (b) axis lengths i.e.

r0 = √(a*b).

The ellipticity `τ` is given by τ = 1-b/a.

### Fields
$(FIELDS)

"""
@with_kw struct TIDAGaussianRing <: AbstractImageTemplate
    """Radius of the Gaussian ring"""
    r0::Float64
    """Standard deviation of the width of the Gaussian ring"""
    σ::Float64
    """Asymmetry of the Gaussian ring defined as ``1-b/a``"""
    τ::Float64
    """Slash of Gaussian ring."""
    s::Float64
    """Slash/Asymmetry orientation in radians measured north of east"""
    ξ::Float64
    """x position of the center of the ring in μas"""
    x0::Float64
    """y position of the center of the ring in μas"""
    y0::Float64
    function TIDAGaussianRing(r0,σ,τ,s,ξ,x0,y0)
        @assert r0>0 "TIDAGaussianRing: r0 must be positive"
        @assert σ>0 "TIDAGaussianRing: σ must be positive"
        @assert τ<=1 && τ>=0 "TIDAGaussianRing: τ must be in [0,1)"
        @assert s>=-1 && s <= 1 "TIDAGaussianRing: s must be in [-1,1]"
        new(float(r0),float(σ),float(τ),float(s),float(ξ),float(x0),float(y0))
    end
end
function TIDAGaussianRing(p)
    @assert length(p)==7 "TIDAGaussianRing: Template requires 7 parameters"
    TIDAGaussianRing(p[1],p[2],p[3],p[4],p[5],p[6],p[7])
end

Base.size(::Type{TIDAGaussianRing}) = 7

#Template function for the TIDAGaussianRing
 @inline function (θ::TIDAGaussianRing)(x,y)
    ex = x-θ.x0
    ey = y-θ.y0
    ex′,ey′ = rotate(ex,ey,θ.ξ)
    a = θ.r0/sqrt(1.0-θ.τ)
    b = θ.r0*sqrt(1.0-θ.τ)
    distance = ellipse_sqdist(ex′,ey′,a, b)


    #construct the slash
    ϕ = atan(ey′,ex′)
    if θ.s >= 0
        n = (1-θ.s*cos(ϕ))/(θ.s + 1)
    else
        ϕ = mod(ϕ+π/2,2π)-π
        n = (1+θ.s*cos(ϕ))/(-θ.s + 1)
    end

    return n*exp(-distance/(2.0*θ.σ^2)) + 1e-50
end


@doc """
    $(TYPEDEF)
Creates the most general elliptical slashed gaussian ring model. It is a combination of
the elliptical and slashed gaussian ring. The direction of the slash and the ellipticity are
not aligned or anti-aligned like with the TIDAGaussianRing type.

### Details
Adds ellipticity to the ring. The radius `r0` of the ring is now defined as the
geometric mean of the semi-major (a) and semi-minor (b) axis lengths i.e.

r0 = √(a*b).

The ellipticity `τ` is given by τ = 1-b/a.

### Fields
$(FIELDS)
"""
@with_kw struct GeneralGaussianRing <: AbstractImageTemplate
    """Radius of the Gaussian ring"""
    r0::Float64
    """Standard deviation of the width of the Gaussian ring"""
    σ::Float64
    """Asymmetry of the Gaussian ring defined as ``1-b/a``"""
    τ::Float64
    """Asymmetry orientation in radians, measured north of east"""
    ξτ::Float64
    """Slash of Gaussian ring."""
    s::Float64
    """Slash orientation in radians measured north of east"""
    ξs::Float64
    """x position of the center of the ring in μas"""
    x0::Float64
    """y position of the center of the ring in μas"""
    y0::Float64
    function GeneralGaussianRing(r0,σ,τ,ξτ,s,ξs,x0,y0)
        @assert r0>0 "GeneralGaussianRing: r0 must be positive"
        @assert σ>0 "GeneralGaussianRing: σ must be positive"
        @assert τ<1 && τ>=0 "GeneralGaussianRing: τ must be in [0,1)"
        @assert s>=0 && s <= 1 "GeneralGaussianRing: s must be in [0,1]"
        new(float(r0),float(σ),float(τ),float(ξτ),float(s),float(ξs),float(x0),float(y0))
    end
end
function GeneralGaussianRing(p)
    @assert length(p)==8 "GeneralGaussianRing: Template requires 8 parameters"
    GeneralGaussianRing(p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8])
end
Base.size(::Type{GeneralGaussianRing}) = 8

 @inline function (θ::GeneralGaussianRing)(x,y)
    ex = x-θ.x0
    ey = y-θ.y0
    ex′,ey′ = rotate(ex,ey,θ.ξτ)
    a = θ.r0/sqrt(1.0-θ.τ)
    b = θ.r0*sqrt(1.0-θ.τ)
    distance = ellipse_sqdist(ex′,ey′,a, b)

    #construct the slash
    ex′,ey′ = rotate(ex,ey,θ.ξs)
    ϕ = atan(ey′,ex′)
    n = (1-θ.s*cos(ϕ))/(θ.s + 1)

    return n*exp(-distance/(2.0*θ.σ^2)) + 1e-50
end



@doc """
    $(TYPEDEF)
Extrememly flexible ring model. The thickness is modeled as a cosine
expansion with `N` terms and the slash by a expansion with `M` terms.


# Details
The ring is allowed to be elliptical.
The thickness of the ring is modeled by a cosine expansion in azimuthal
angle. `N` specifies the number of cosine modes to fit, where the first
mode is the constant thickness portion and so has no corresponding angle.
The slash is modeled as a separate cosine expansion, with `M` terms.
Here the zero order term is forced to be unity, so `M` defines the `M`
additional terms.

"""
@with_kw struct CosineRing{N,M} <: AbstractImageTemplate
    """Radius of the Gaussian ring"""
    r0::Float64
    """Standard deviations (length N+1) of the width of the Gaussian ring"""
    σ::Vector{Float64}
    """Orientations of the cosine expansion width (length N)"""
    ξσ::Vector{Float64}
    """Asymmetry of the Gaussian ring defined as ``1-b/a``"""
    τ::Float64
    """Asymmetry orientation in radians, measured north of east"""
    ξτ::Float64
    """Slash of Gaussian ring (length M)."""
    s::Vector{Float64}
    """Slash orientations (length M) in radians measured north of east"""
    ξs::Vector{Float64}
    """x position of the center of the ring in μas"""
    x0::Float64
    """y position of the center of the ring in μas"""
    y0::Float64
    function CosineRing{N,M}(r0,σ, ξσ, τ, ξτ,s, ξs,x0,y0) where {N, M}
        #@assert N isa Integer
        #@assert M isa Integer
        new{N,M}(float(r0),σ, ξσ, float(τ),float(ξτ), s, ξs,float(x0),float(y0))
    end
end

@doc """
    CosineRing{N,M}(p::AbstractArray) where {N,M}
Takes in a vector of paramters describing the template.
# Details
The order of the vector must be
 - p[1] = `r0`
 - p[2:(N+1)] = `σ`
 - p[(N+2):(2N)] = `ξσ`
 - p[2N+1] = `τ`
 - p[2N+2] = `ξτ`
 - p[2N+3:2N+M+2] = `s`
 - p[2N+3+M:2N+2+2M] = `ξs`
 - p[2N+3+2M] = `x0`
 - p[2N+4+2M] = `y0`
"""
function CosineRing{N,M}(p::AbstractArray) where {N,M}
    #@assert length(p) == size(CosineRing{N,M})
    CosineRing{N,M}(p[1], p[2:(N+2)],
                    p[(N+3):(2N+2)],
                    p[2N+3], p[2N+4],
                    p[2N+5:2N+4+M], p[2N+5+M:2N+4+2M],
                    p[2N+5+2M], p[2N+6+2M]
    )
end


Base.size(::Type{CosineRing{N,M}}) where {N, M} = 5 + N+1 + N + 2*M

 @inline function (θ::CosineRing{N,M})(x,y) where {N, M}
    ex = x-θ.x0
    ey = y-θ.y0
    ϕ = atan(-ey,-ex)
    ex′,ey′ = rotate(ex,ey,θ.ξτ)
    a = θ.r0/sqrt(1.0-θ.τ)
    b = θ.r0*sqrt(1.0-θ.τ)
    d2 = ellipse_sqdist(ex′,ey′,a, b)

    #construct the slash
    n = one(θ.r0)
    for i in 1:M
        n -= θ.s[i]*cos(i*(ϕ - θ.ξs[i]))
    end

    σ = θ.σ[1]
    for i in 1:N
        σ += θ.σ[i+1]*cos(i*(ϕ - θ.ξσ[i]))
    end

    return abs(n)*exp(-d2/(2.0*σ^2+1e-2))
end

function unpack(θ::CosineRing{N,M}) where {N,M}
    n = Base.size(typeof(θ))
    p = zeros(n)
    p[1] = θ.r0
    p[2:(N+2)] = θ.σ
    if N>0
        p[(N+3):(2N+2)] = θ.ξσ
    end
    p[2N+3] = θ.τ
    p[2N+4] = θ.ξτ
    p[2N+5:2N+4+M] = θ.s
    p[2N+5+M:2N+4+2M] = θ.ξs
    p[2N+5+2M] = θ.x0
    p[2N+6+2M] = θ.y0

    return p
end
