


@doc """
    $(TYPEDEF)
Combines two templates together into one object. Since addition is
assoiciative this can actually we used to hold multiple different templates.

### Details
Overloads the Base.:+ function so you can easily add two templates together.

### Example
```julia
θ1 = GaussianRing(10,5,0,0)
θ2 = SlashedGaussianRing(15,5,0.5,π/4,0,0)
θ12 = θ1+θ2
```
"""
struct AddTemplate{T1<:AbstractTemplate,T2<:AbstractTemplate} <: AbstractCompositeTemplate
    θ1::T1
    θ2::T2
end

function AddTemplate{T1,T2}(p) where {T1<:AbstractTemplate,T2<:AbstractTemplate}
    p1 = @view p[1:Base.size(T1)]
    θ1 = T1(p1)
    p2 = @view p[(Base.size(T1)+1):end]
    θ2 = T2(p2)
    AddTemplate{T1,T2}(θ1,θ2)
end

function Base.size(::Type{AddTemplate{T1,T2}}) where {T1<:AbstractTemplate, T2<:AbstractTemplate}
    return Base.size(T1) + Base.size(T2)
end

function Base.propertynames(θ::AddTemplate{T1,T2}) where {T1<:AbstractTemplate, T2<:AbstractTemplate}
    return (propertynames(θ.θ1)...,propertynames(θ.θ2)...)
end

Base.:+(x1::T1,x2::T2) where {T1<:AbstractTemplate,T2<:AbstractTemplate} = AddTemplate(x1,x2)
function (θ::AddTemplate)(x,y)
    return θ.θ1(x,y) + θ.θ2(x,y)
end

function unpack(θinit::AddTemplate)
    p1 = unpack(θinit.θ1)
    p2 = unpack(θinit.θ2)

    return append!(p1,p2)
end

@doc """
    $(SIGNATURES)
Stacks templates together so you can easily combine multiple templates.
It does this by calling the :+ and :* method. Every template added will
include an additional parameter that controls the relative weight of each template.
"""
function stack(θ::T, θ1...) where {T<:AbstractTemplate}
    return θ+mapreduce(x->1.0*x, + , θ1)
end

@doc """
    $(SIGNATURES)
Splits the template into an array with its subcomponents so you can easily access them.
"""
function Base.split(θ::AbstractTemplate)
    return [θ]
end

function Base.split(θ::AddTemplate)
    return [split(θ.θ1)..., split(θ.θ2)...]
end
