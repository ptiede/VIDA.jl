using HypercubeTransform, Random
using Distributions: Uniform, product_distribution

export VIDAProblem, vida

"""
    $(TYPEDEF)

A composite type that holds various properties that define the optimization process
to extract the optimal filter.

## Fields
$(FIELDS)
"""
Base.@kwdef struct VIDAProblem{D<:AbstractDivergence, F, N, B}
    """
    Divergence you will fit
    """
    div::D
    """
    Function that defines the parametric family of templates
    """
    f::F
    """
    Type of autodiff to use when optimizing if any
    """
    autodiff::N = nothing
    """
    The lower bounds of the parameter ranges to search over
    """
    lb::B
    """
    The upper bounds of parameter ranges to search over
    """
    ub::B
end

"""
    $(SIGNATURES)

Defines a `VIDAProblem` for optimization.

## Arguments
  - `div`: The divergence you wish to fit. This is an instantiation of [`VIDA.AbstractDivergence`](@ref)
  - `f`:   The function that defines the parametric family of templates you will fit. The function
           must accept a named tuple as its first argument, whose names define the parameters.
  - `lb`:  A NamedTuple whose names match the argument of `f` and whose values define the lower
           bounds of the parameters range you want to search over
   - `lb`: A NamedTuple whose names match the argument of `f` and whose values define the upper
           bounds of the parameters range you want to search over

## Example
```julia-repl
julia> f(x) = SlashedGaussianRing(x.σ, x.s)
julia> div = Renyi(img, 0.75)
julia> lb = (x = 0.1, s = 0.001)
julia> ub = (x = 10.0, s = 0.999)
julia> prob = VIDAProblem(div, f, lb, ub)
```
"""
function VIDAProblem(div, f, lb, ub)
    return VIDAProblem(div, f, nothing, lb, ub)
end

_distize(x::Real, y::Real) = Uniform(x, y)
_distize(x::Tuple, y::Tuple) = ntuple(i->_distize(x[i], y[i]), length(x))
_distize(x::NamedTuple{N}, y::NamedTuple{N}) where {N} = NamedTuple{N}(_distize(values(x), values(y)))

function _distize(x::AbstractArray, y::AbstractArray)
    dists = map((x,y)->Uniform(x, y), x, y)
    return product_distribution(dists)
end


function build_opt(prob, unit_cube)
    T = eltype(prob.div.img)
    dist = map((x,y)->_distize(x, y), prob.lb, prob.ub)
    if unit_cube
        t = ascube(dist)
        bounds_lower = fill(zero(T), dimension(t))
        bounds_upper = fill(one(T), dimension(t))
    else
        t = asflat(dist)
        bounds_lower = fill(-convert(T, 20), dimension(t))
        bounds_upper = fill(convert(T, 20), dimension(t))

    end
    f = let t=t, div = prob.div, mod = prob.f, lb=bounds_lower, ub=bounds_upper
        x->begin
            for i in eachindex(x)
                (lb[i] > x[i] || ub[i] < x[i]) && return convert(T, Inf)
            end
            return divergence(div, mod(transform(t, x)))
        end
    end
    return f, t, (bounds_lower, bounds_upper)
end

function initial_point(rng, T::Type, t::HypercubeTransform.AbstractHypercubeTransform, init_params)
    !isnothing(init_params) && return inverse(t, init_params)
    return rand(rng, T, dimension(t))
end

function initial_point(rng, T::Type, t::HypercubeTransform.TransformVariables.AbstractTransform, init_params)
    !isnothing(init_params) && return inverse(t, init_params)
    return randn(rng, T,  dimension(t))
end

"""
    $(SIGNATURES)

Runs the `VIDA` algorithm to find the optimal template defined by the problem `prob`.
The optimization will be run using the `Optimization.jl` `optimizer`. You can optionally
pass a random number generator used to specify the initial points or can explicitly pass
the initial location of the optimization using `init_params`. By default `vida` first
transforms the parameter space to the unit hypercube. If you do not wish to do this
set `unit_cube = false`.

The remaining `kwargs...` are forwarded to the `Optimization.solve` function.

!!! warn
    To use this function a user must first load Optimization, i.e. `using Optimization`

## Arguments
   - `prob`: Defines the problem you wish to solve.
   - `optimizer`: Specifies the optimizer you want to use. Any optimizer from `Optimization.jl` works.

## Optional Keyword arguments
   - `rng`: Specifies the random number generator you want to use to select the initial points.
            Note that is not forwarded to the optimizers since not all can use a specific rng.
   - `init_params`: Specify the initial point of the optimization routine. The default `nothing`
                    will randomly draw a starting point from the parameter space
   - `unit_cube`: If true the parameters are first transformed to the unit hypercube. If false
                  they are transformed to ℝⁿ using `TransformVariables`
   - `kwargs...`: Additional options to be passed to the `solve` function from `Optimization.jl`
"""
function vida end

