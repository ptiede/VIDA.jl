
#Template information and documentation

"""
    $(TYPEDEF)
An abstract type that defines super template type.
"""
abstract type AbstractTemplate <: CB.AbstractModel end

# Hook into ComradeBase interface
CB.visanalytic(::Type{<:AbstractTemplate}) = CB.NotAnalytic()
CB.imanalytic(::Type{<:AbstractTemplate}) = CB.IsAnalytic()
CB.isprimitive(::Type{<:AbstractTemplate}) = CB.IsPrimitive()

CB.flux(::AbstractTemplate) = 1.0


"""
    size(f::AbstractTemplate)
Get the number of parameters for the template f
"""
Base.size(f::T) where {T<:AbstractTemplate} = Base.size(T)


# Updates the paramters for a given template even though they are immutable
function _update(θ::AbstractTemplate, p)
    return typeof(θ)(p)
end

#Load the variety of utils needed
include(joinpath(@__DIR__, "utils.jl"))


"""
    $(TYPEDEF)

An abstract type that will contain the template information, such as the parameters.
Specific instanstantiations will need to be defined for you to use this.

### Details
    This defined the highest function type. If you wish to implement your own template you
    need to define a a couple of things
    1. The template type <: AbstractTemplate
    2. an functor of the type that computes the template function
    3. an `Base.size` function that defines the number of parameters of the template.

An example is given by:
```julia
#All of our composite type are defined using the Paramters.jl package to you
can directly refer to the struct parameters when creating it, although this isn't
actually used anywhere in the code.
@with_kw struct Gaussian <: AbstractImageTemplate
    σ::Float64
    x0::Float64
    y0::Float64
end

#Typically we inline and force the function to use fastmath
 @inline function (θ::Gaussian)(x,y)
    return 1.0/(2π*σ^2)*exp(-0.5*( (x-x0)^2+(y-y0)^2 ) )
end
Base.size(::Type{Gaussian}) = 3
```
"""
abstract type AbstractImageTemplate <: AbstractTemplate end
include(joinpath(@__DIR__, "image.jl"))

"""
    $(TYPEDEF)
Meta-template that accepts a abstract template and modifies it through
some analytic transformation. This is useful for model composability.
"""
abstract type AbstractModifierTemplate <: AbstractTemplate end
include(joinpath(@__DIR__, "modifiers.jl"))

"""
    $(TYPEDEF)
Abstract type for taking two templates composing them in some manner
and returning a new template.
"""
abstract type AbstractCompositeTemplate <: AbstractTemplate end
include(joinpath(@__DIR__,"composite.jl"))
