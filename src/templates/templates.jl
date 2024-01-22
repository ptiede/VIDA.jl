"""
    $(TYPEDEF)

An abstract `ComradeBase.AbstractModel` that serves as the parent type of all
the `VIDA` templates implemented in this repo.

By default this model assumes the following `ComradeBase` traits
```julia
ComradeBase.visanalytic(::Type{<:AbstractImageTemplate}) = NoAnalytic()
ComradeBase.imanalytic(::Type{<:AbstractImageTemplate})  = IsAnalytic()
ComradeBase.ispolarized(::Type{<:AbstractImageTemplate}) = NotPolarized()
```

As a result if a user wishes the implement their own subtype (e.g., `MyTemplate`) of `AbstratImageTemplate`
they will need to implement the following methods
  - `ComradeBase.intensity_point(m::MyTemplate, p)`: which computes the potentially unormalized brightness of the template at the point `p`.
  - `ComradeBase.radialextent(m::MyTemplate)`: which computes the rough radial extent of the model `m`.

For more information about the total interface see [VLBISkyModels.jl](https://ehtjulia.github.io/VLBISkyModels.jl/stable/interface/)
"""
abstract type AbstractImageTemplate <: ComradeBase.AbstractModel end

# Hook into ComradeBase interface
CB.visanalytic(::Type{<:AbstractImageTemplate}) = CB.NotAnalytic()
CB.imanalytic(::Type{<:AbstractImageTemplate})  = CB.IsAnalytic()
CB.isprimitive(::Type{<:AbstractImageTemplate}) = CB.IsPrimitive()

function VLBISkyModels.__extract_tangent(m::AbstractImageTemplate)
    nt = ntfromstruct(m)
    if ntm isa NamedTuple{(), Tuple{}}
        return ZeroTangent()
    else
        return Tangent{typeof(dm)}(;nt...)
    end
end


include(joinpath(@__DIR__, "image.jl"))
include(joinpath(@__DIR__, "rings.jl"))
include(joinpath(@__DIR__, "cosinering.jl"))
