@inline basetemplate(θ::AbstractModifierTemplate) = θ.θ
function unpack(θ::AbstractModifierTemplate)
    return push!(unpack(basetemplate(θ)), modifiers(θ)...)
end

"""
    StretchMod(θ::AbstractTemplate, τ)
Modifies a filter by stretching it in the x-direction by 1/√(1-τ) and
y-direction by sqrt(1-τ). Typically call the stretch function instead.
"""
@with_kw struct StretchMod{T<:AbstractTemplate, S<:Number} <: AbstractModifierTemplate
    """
    unmodified template
    """
    θ::T
    """
    Stretch in the x-direction
    """
    scx::S
    """
    Stretch in the y-direction
    """
    scy::S


end


function imagecenter(θ::AbstractModifierTemplate)
    return imagecenter(basetemplate(θ))
end

function StretchMod(θ::T, τ::S) where {T<:AbstractTemplate,S<:Number}
    scx = 1/sqrt(1-τ)
    scy = 1/scx
    return StretchMod(θ, scx, scy)
end


function StretchMod{T,S}(p::AbstractVector) where {T<:AbstractTemplate, S<:Number}
    StretchMod(T(@view p[1:end-1]), p[end])
end

"""
    $(SIGNATURES)
Stretch the template m, using the τ factor. This will stretch
the x-axis by 1/√(1-τ) and y-axis by √(1-τ).
"""
function stretch(m::AbstractTemplate, τ)
    scx = 1/sqrt(1-τ)
    scy = sqrt(1-τ)
    StretchMod(m, τ)
end
modifiers(θ::StretchMod) = 1-θ.scy/θ.scx
function Base.size(::Type{StretchMod{T,S}}) where {T<:AbstractTemplate, S<:Number}
    return 1 + size(T)
end

@inline function (θ::StretchMod)(x,y)
    return basetemplate(θ)(x/θ.scx,y/θ.scy)
end

function Base.propertynames(θ::StretchMod{T,S}) where {T<:AbstractTemplate, S<:Number}
    return (propertynames(basetemplate(θ))...,:τ)
end


@doc """
    RotateMod(θ::AbstractTemplate, ξ)
Modifies a filter by rotating it using the rotation matrix with angle ξ.
Typically do not call this function directly but instread the rotate
function.
"""
@with_kw struct RotateMod{T<:AbstractTemplate,S} <: AbstractModifierTemplate
    """
    unmodified template
    """
    θ::T
    """
    cosine of rotation angle
    """
    c::S
    """
    sin of rotation angle
    """
    s::S

end
function RotateMod(θ::T, ξ::S) where {T<:AbstractTemplate,S<:Number}
    s,c = sincos(ξ)
    return RotateMod(θ, c, s)
end

function RotateMod{T,S}(p::AbstractVector) where {T<:AbstractTemplate, S<:Number}
    RotateMod(T(@view p[1:end-1]), p[end])
end

"""
    $(SIGNATURES)
Rotate the template m by ξ. This rotates north of east.
"""
rotate(m::AbstractTemplate, ξ) = RotateMod(m, ξ)
modifiers(θ::RotateMod) = atan(θ.s, θ.c)
function Base.size(::Type{RotateMod{T,S}}) where {T<:AbstractTemplate, S<:Number}
    return 1 + size(T)
end

function Base.propertynames(θ::RotateMod{T,S}) where {T<:AbstractTemplate, S<:Number}
    return (propertynames(basetemplate(θ))...,:ξ)
end

@inline function (θ::RotateMod)(x,y)
    x0,y0 = imagecenter(θ)
    ex,ey = x-x0, y-y0
    x′ = θ.c*ex + θ.s*ey
    y′ = -θ.s*ex + θ.c*ey
    return basetemplate(θ)(x′+x0,y′+y0)
end

"""
    $(SIGNATURES)
Stretch and then rotate a template producing a modified template.
This is equivalent to rotate(stretch(m, τ), ξ).
"""
@inline function stretchrotate(m, τ, ξ)
    return rotate(stretch(m, τ), ξ)
end



@doc """
    $(TYPEDEF)
Multiplies template by a constant. This is useful when combining with
AddTemplate since it will change the relative weights of each template.

### Details
Overloads the Base.:* function so you can easily multiple a template by a number.

### Example
```julia
θ = GaussianRing(15,5,0.0,0.0)
2*θ
```
"""
struct MulTemplate{T<:AbstractTemplate,S<:Number} <: AbstractModifierTemplate
    θ::T
    Irel::S
end
modifiers(θ::MulTemplate) = θ.Irel

function Base.show(io::IO,θ::MulTemplate{T,S}) where {T<:AbstractTemplate, S<:Number}
    println(io,"VIDA.MulTemplate{$T,$S}")
    print(io,"θ: ")
    show(io,θ.θ)
    println(io,"Irel: $S $(θ.Irel)")
end

function MulTemplate{T,S}(p) where {T<:AbstractTemplate, S<:Number}
    MulTemplate{T,S}(T(@view p[1:end-1]), p[end])
end

function Base.propertynames(θ::MulTemplate{T,S}) where {T<:AbstractTemplate, S<:Number}
    return (propertynames(basetemplate(θ))...,:Irel)
end

function Base.size(::Type{MulTemplate{T,S}}) where {T<:AbstractTemplate, S<:Number}
    return 1 + size(T)
end

Base.:*(a,x::T) where {T<:AbstractTemplate} = MulTemplate(x,a)
Base.:*(x::T,a) where {T<:AbstractTemplate} = MulTemplate(x,a)

@inline function (θ::MulTemplate)(x,y)
    return θ.Irel*basetemplate(θ)(x,y)
end



function Base.getproperty(θmul::MulTemplate{T,S}, field::Symbol) where {T<:AbstractTemplate, S<:Number}
    if field == :θ
        return getfield(θmul,:θ)
    elseif field == :Irel
        return getfield(θmul, :Irel)
    elseif field in propertynames(basetemplate(θmul))
        return getfield(θmul.θ, field)
    else
        throw(KeyError(field))
    end
end
