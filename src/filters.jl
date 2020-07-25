
#Filter information and documentation
"""
    $(TYPEDEF)
An absract type that will contain the filter information, such as the parameters.
Specific instanstantiations will need to be defined for you to use this.

### Details
    This defined the highest function type. If you wish to implement your own filter you
    need to define a a couple of things
    1. The filter type <: AbstractFilter
    2. an functor of the type that computes the filter function
    3. an `size` function that defines the number of parameters of the filter.

An example is given by:
```julia
#All of our composite type are defined using the Paramters.jl package to you
can directly refer to the struct parameters when creating it, although this isn't
actually used anywhere in the code.
@with_kw mutable struct Gaussian <: AbstractFilter
    σ::Float64
    x0::Float64
    y0::Float64
end

#Typically we inline and force the function to use fastmath
@fastmath @inline function (θ::Gaussian)(x,y)
    return 1.0/(2π*σ^2)*exp(-0.5*( (x-x0)^2+(y-y0)^2 ) )
end
size(::Type{Gaussian}) = 3
```
"""
abstract type AbstractFilter end


#Load the filters
"""
    $(TYPEDEF)
An constant filter.

### Details
Defines an image that just has constant flux. This is very useful for soaking up
low levels of flux in image reconstructions that can bias the results.

Since images or normalized to unity, this means the `Constant` filter has no
additional parameters.
"""
struct Constant <: AbstractFilter end
Constant(p) = Constant()
size(::Type{Constant}) = 0
@inline (θ::Constant)(x,y) = 1

"""
    $(TYPEDEF)
A smoothed disk model

### Details
Defines a filter for an image that has a smoothed disk model.


"""
@with_kw mutable struct Disk <: AbstractFilter
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
    @assert length(p) == 4 "There are 4 parameters for the Disk filter."
    Disk(p[1],p[2],p[3],p[4])
end
@fastmath @inline function (θ::Disk)(x,y)
    @unpack r0, α, x0, y0 = θ
    r = sqrt((x-x0)^2 + (y-y0)^2)
    if ( r < r0 )
        return one(eltype(r))
    else
        return exp(-(r-r0)/(2.0*α))
    end
end
size(::Type{Disk}) = 4


"""
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
@with_kw mutable struct AsymGaussian <: AbstractFilter
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
    @assert length(p) == 5 "There are 5 parameters for the AsymGaussian filter."
    AsymGaussian(p[1],p[2],p[3],p[4],p[5])
end

@fastmath @inline function (θ::AsymGaussian)(x,y)
    x′,y′ = rotate(x-θ.x0,y-θ.y0,θ.ξ)
    σx2 = θ.σ *θ.σ/(1.0-θ.τ)
    σy2 = θ.σ*θ.σ*(1.0-θ.τ)
    d2 = x′*x′/σx2 + y′*y′/σy2
    return exp(-0.5*d2)
end
size(::Type{AsymGaussian}) = 5


"""
    $(TYPEDEF)
Symmetric gaussian ring filter. This is the most basic filter and just attempts
to recover a location `x0`,`y0`, radius `r0` and thickness `σ` from some image.
### Fields
    $(FIELDS)
### Example
```julia
GaussianRing(r0=20.0,σ=5.0,x0=0.0,y0=-10.0)
```
"""
@with_kw mutable struct GaussianRing <: AbstractFilter
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
    @assert length(p)==4 "GaussianRing: Filter requires 4 parameters"
    GaussianRing(p[1],p[2],p[3],p[4])
end

size(::Type{GaussianRing}) = 4

@fastmath @inline function (θ::GaussianRing)(x, y)
    r = sqrt((x - θ.x0)^2 + (y - θ.y0)^2)
    return exp( -(r - θ.r0)^2/(2*θ.σ^2))
end


"""
    $(TYPEDEF)
Implements the slashed gaussian ring filter, that uses a cosine
to symmetrically implement the slash. While this is marginally more
complicated that a linear slash, it has a number of benefits such as
mainting the azimuthal and smooth structure of the image.

### Fields
    $(FIELDS)

"""
@with_kw mutable struct SlashedGaussianRing <: AbstractFilter
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
        @assert length(p)==6 "SlashedGaussianRing: Filter requires 6 parameters"
        SlashedGaussianRing(p[1],p[2],p[3],p[4],p[5],p[6])
end
size(::Type{SlashedGaussianRing}) = 6
#Filter function
@fastmath @inline function (θ::SlashedGaussianRing)(x,y)
    r = sqrt((x - θ.x0)^2 + (y - θ.y0)^2)

    #rotate the image so slash is on the x axis
    xrot,yrot = rotate(x-θ.x0,y-θ.y0,θ.ξ)
    #construct the slash
    ϕ = atan(yrot,xrot)
    n = (1-θ.s*cos(ϕ))/(θ.s + 1)

    return n*exp(-(r-θ.r0)^2/(2*θ.σ^2))
end


"""
    $(TYPEDEF)
Implements the elliptical gaussian ring filter. Where the ellipticity `tau` is defined
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
@with_kw mutable struct EllipticalGaussianRing <: AbstractFilter
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
    @assert length(p)==6 "EllipticalGaussianRing: Filter requires 6 parameters"
    EllipticalGaussianRing(p[1],p[2],p[3],p[4],p[5],p[6])
end
size(::Type{EllipticalGaussianRing}) = 6

@fastmath @inline function (θ::EllipticalGaussianRing)(x,y)
    ex = x-θ.x0
    ey = y-θ.y0
    ex′,ey′ = rotate(ex,ey,θ.ξ)
    a = θ.r0/sqrt(1.0-θ.τ)
    b = θ.r0*sqrt(1.0-θ.τ)
    distance = ellipse_sqdist(ex′,ey′,a, b)
    return exp(-distance/(2.0*θ.σ^2))
end

"""
    $(TYPEDEF)
Creates the filter from the Paper I am writing. It is a combination of
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
@with_kw mutable struct TIDAGaussianRing <: AbstractFilter
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
    @assert length(p)==7 "TIDAGaussianRing: Filter requires 7 parameters"
    TIDAGaussianRing(p[1],p[2],p[3],p[4],p[5],p[6],p[7])
end

size(::Type{TIDAGaussianRing}) = 7

#Filter function for the TIDAGaussianRing
@fastmath @inline function (θ::TIDAGaussianRing)(x,y)
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

    return n*exp(-distance/(2.0*θ.σ^2))
end


"""
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
@with_kw mutable struct GeneralGaussianRing <: AbstractFilter
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
    @assert length(p)==8 "GeneralGaussianRing: Filter requires 8 parameters"
    GeneralGaussianRing(p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8])
end
size(::Type{GeneralGaussianRing}) = 8

@fastmath @inline function (θ::GeneralGaussianRing)(x,y)
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

    return n*exp(-distance/(2.0*θ.σ^2))
end



#Rotates our points. Note that we use astronomer conventions
@fastmath @inline function rotate(x,y,ξ)
    s,c = sincos(π-ξ)
    x′ = c*x - s*y
    y′ = s*x + c*y

    return (x′ ,y′)
end

"""
    $(SIGNATURES)
Find the minimum square distance between an ellipse centered at (0,0) with semi-major
axis `a` and semi-minor axis `b` and the point (`x`, `y`).
Uses an iterative method with accuracy `ϵ` which defaults to 1e-6.
# Credit:
Algorithm taken from
https://github.com/0xfaded/ellipse_demo/issues/1#issuecomment-405078823
except written in Julia and using an adaptive termination condition for accuracy
"""
#@fastmath
@fastmath @inline function ellipse_sqdist(x,y, a, b, ϵ=1e-6)
    #For simplicity we will only look at the positive quadrant
    px = abs(x)
    py = abs(y)

    #initial guess
    tx = 0.707
    ty = 0.707
    err = 1.0
    n = 0
    while err > ϵ && n < 10
        tx′, ty′ = dist_ellipse_unit(px,py,tx,ty,a,b)
        #err = hypot((tx′-tx), ty′-ty)
        err = sqrt((tx-tx′)*(tx-tx′) + (ty-ty′)*(ty-ty′))
        tx = tx′
        ty = ty′
        n+=1
    end
    return (a*tx-px)*(a*tx-px) + (b*ty-py)*(b*ty-py)
end

#Finds the closest point on the ellipse. This is an internal function
#used for ellipse_distance.
#@fastmath
@fastmath @inline function dist_ellipse_unit(px, py, tx, ty, a, b)
    x′ = a*tx
    y′ = b*ty

    ex = (a*a - b*b)*tx*tx*tx/a
    ey = (b*b - a*a)*ty*ty*ty/b

    rx = x′ - ex
    ry = y′ - ey

    qx = px - ex
    qy = py - ey

    r = sqrt(rx*rx + ry*ry)
    q = sqrt(qx*qx + qy*qy)

    xx = min(1.0, max(0.0, (qx*r/q + ex)/a) )
    yy = min(1.0, max(0.0, (qy*r/q + ey)/b) )
    #t = hypot(xx,yy)
    t = sqrt(xx*xx + yy*yy)
    xx /= t
    yy /= t

    return xx,yy
end

"""
    $(SIGNATURES)
Unpacks the parameters of the filter `θ`, except the last which is
the boolean for if the filter is analytically normalized.


Returns the parameters in a vector.
"""
function unpack(θinit::T) where {T<:AbstractFilter}
    n = length(fieldnames(T))
    fields = fieldnames(T)
    p = zeros(n)
    for i in 1:n
        p[i] = getfield(θinit,fields[i])
    end
    return p
end


"""
    $(TYPEDEF)
Combines two filters together into one object. Since addition is
assoiciative this can actually we used to hold multiple different filters.

### Details
Overloads the Base.:+ function so you can easily add two filters together.

### Example
```julia
θ1 = GaussianRing(10,5,0,0)
θ2 = SlashedGaussianRing(15,5,0.5,π/4,0,0)
θ12 = θ1+θ2
```
"""
struct AddFilter{T1<:AbstractFilter,T2<:AbstractFilter} <: AbstractFilter
    θ1::T1
    θ2::T2
end

function AddFilter{T1,T2}(p) where {T1<:AbstractFilter,T2<:AbstractFilter}
    p1 = @view p[1:size(T1)]
    θ1 = T1(p1)
    p2 = @view p[(size(T1)+1):end]
    θ2 = T2(p2)
    AddFilter{T1,T2}(θ1,θ2)
end

function size(::Type{AddFilter{T1,T2}}) where {T1<:AbstractFilter, T2<:AbstractFilter}
    return size(T1) + size(T2)
end

function Base.fieldnames(::Type{AddFilter{T1,T2}}) where {T1<:AbstractFilter, T2<:AbstractFilter}
    return [fieldnames(T1)...,fieldnames(T2)...]
end

Base.:+(x1::T1,x2::T2) where {T1<:AbstractFilter,T2<:AbstractFilter} = AddFilter(x1,x2)
function (θ::AddFilter)(x,y)
    return θ.θ1(x,y) + θ.θ2(x,y)
end

function unpack(θinit::AddFilter)
    p1 = unpack(θinit.θ1)
    p2 = unpack(θinit.θ2)

    return append!(p1,p2)
end

"""
    $(TYPEDEF)
Multiplies filter by a constant. This is useful when combining with
AddFilter since it will change the relative weights of each filter.

### Details
Overloads the Base.:* function so you can easily multiple a filter by a number.

### Example
```julia
θ = GaussianRing(15,5,0.0,0.0)
2*θ
```
"""
struct MulFilter{S<:Number,T<:AbstractFilter} <: AbstractFilter
    a::S
    θ::T
end

function MulFilter{S,T}(p) where {S<:Number,T<:AbstractFilter}
    MulFilter{S,T}(p[end],T(@view p[1:end-1]))
end

function Base.fieldnames(::Type{MulFilter{S,T}}) where {S<:Number, T<:AbstractFilter}
    return [fieldnames(T)...,:Irel]
end

function size(::Type{MulFilter{S,T}}) where {S<:Number,T<:AbstractFilter}
    return 1 + size(T)
end

Base.:*(a,x::T) where {T<:AbstractFilter} = MulFilter(a,x)
Base.:*(x::T,a) where {T<:AbstractFilter} = MulFilter(a,x)
function (θ::MulFilter)(x,y)
    return θ.a*θ.θ(x,y)
end
function unpack(θinit::MulFilter)
    return append!(unpack(θinit.θ), θinit.a)
end

"""
    $(SIGNATURES)
Stacks filters together so you can easily combine multiple filters.
It does this by calling the :+ and :* method. Every filter added will
include an additional parameter that controls the relative weight of each filter.
"""
function Base.cat(θ::T, θ1...) where {T<:AbstractFilter}
    return θ+mapreduce(x->1.0*x, + , θ1)
end

"""
    $(SIGNATURES)
Splits the filter into an array with its subcomponents so you can easily access them.
"""
function Base.split(θ::AbstractFilter)
    return [θ]
end

function Base.split(θ::AddFilter)
    return [split(θ.θ1)..., split(θ.θ2)...]
end


"""
    $(SIGNATURES)
Creates an npix×npix rasterized image of the filter `θ` with
limits `xlim` and `ylim`

Returns the tuple (xitr,yitr,image) where xitr,yitr are the iterators
defining the pixel locations (which are centered) and the rasterized image,
 in Jy/μas^2.

# Note
I use the pixel size definition field_of_view/npix, but the image is evaluated
at the pixel centers.

We also use the astronomer orientation and ordering.

"""
function filter_image(θ::AbstractFilter,
                      npix::Int, xlim, ylim)
    fovx = xlim[2]-xlim[1]
    fovy = ylim[2]-ylim[1]

    px = fovx/(npix)
    py = fovy/(npix)

    xitr = (fovx/2-px/2):-px:(-fovx/2)
    yitr = (-fovy/2+py/2):py:(fovy/2)
    img = Matrix{Float64}(undef,npix,npix)
    for (i,x) in enumerate(xitr)
        for (j,y) in enumerate(yitr)
            img[j,i] = θ(x,y)
        end
    end
    return (xitr,yitr,img)
end
