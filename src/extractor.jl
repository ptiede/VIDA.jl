using HypercubeTransform, Random
using Distributions: Uniform, product_distribution
using Optimization

export VIDAProblem, vida, threaded_vida

Base.@kwdef struct VIDAProblem{D<:AbstractDivergence, F, N, B}
    div::D
    f::F
    autodiff::N = SciMLBase.NoAD()
    lb::B
    ub::B
end

function VIDAProblem(div, f, lb, ub)
    return VIDAProblem(div, f, SciMLBase.NoAD(), lb, ub)
end

_distize(x::Real, y::Real) = Uniform(x, y)
_distize(x::Tuple, y::Tuple) = ntuple(i->_distize(x[i], y[i]), length(x))
_distize(x::NamedTuple{N}, y::NamedTuple{N}) where {N} = NamedTuple{N}(_distize(values(x), values(y)))

function _distize(x::AbstractArray, y::AbstractArray)
    dists = map((x,y)->Uniform(x, y), x, y)
    return product_distribution(dists)
end


function build_opt(prob, unit_cube)
    dist = map((x,y)->_distize(x, y), prob.lb, prob.ub)
    if unit_cube
        t = ascube(dist)
        bounds_lower = fill(1e-3, dimension(t))
        bounds_upper = fill(0.999, dimension(t))
    else
        t = asflat(dist)
        bounds_lower = fill(-20.0, dimension(t))
        bounds_upper = fill(20.0, dimension(t))

    end
    f = let t=t, div = prob.div, mod = prob.f, lb=bounds_lower, ub=bounds_upper
        x->begin
            for i in eachindex(x)
                (lb[i] > x[i] || ub[i] < x[i]) && return Inf
            end
            return divergence(div, mod(transform(t, x)))
        end
    end
    return f, t, (bounds_lower, bounds_upper)
end

function initial_point(rng, t::HypercubeTransform.AbstractHypercubeTransform, init_params)
    !isnothing(init_params) && return inverse(t, init_params)
    return rand(rng, dimension(t))
end

function initial_point(rng, t::HypercubeTransform.TransformVariables.AbstractTransform, init_params)
    !isnothing(init_params) && return inverse(t, init_params)
    return randn(rng, dimension(t))
end


function vida(prob::VIDAProblem, optimizer; rng=Random.default_rng(), init_params=nothing, unit_cube=true, kwargs...)
    f, t, (lb, ub) = build_opt(prob, unit_cube)
    x0 = initial_point(rng, t, init_params)
    fopt = OptimizationFunction((x,p)->f(x), prob.autodiff)
    optprob = OptimizationProblem(fopt, x0, nothing; lb=lb, ub=ub)
    xopt, min =  _vida(fopt, t, optprob, optimizer; kwargs...)
    return xopt, prob.f(xopt), min
end


function _vida(fopt, t, optprob, optimizer; kwargs...)
    sol = solve(optprob, optimizer; kwargs...)
    xopt = transform(t, sol.u)
    return xopt, sol.objective
end
