using HypercubeTransform, Random, Distributions
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


function build_opt(prob, unit_cube)
    dist = map((x,y)->Uniform(x, y), prob.lb, prob.ub)
    if unit_cube
        t = ascube(dist)
        bounds_lower = fill(1e-3, dimension(t))
        bounds_upper = fill(0.999, dimension(t))
    else
        t = asflat(dist)
        bounds_lower = nothing
        bounds_upper = nothing

    end
    f = let t=t, div = prob.div, mod = prob.f
        x->divergence(div, mod(transform(t, x)))
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
    optprob = OptimizationProblem(fopt, x0; lb, ub)
    xopt, min =  _vida(fopt, t, optprob, optimizer; kwargs...)
    return xopt, prob.f(xopt), min
end


function _vida(fopt, t, optprob, optimizer; kwargs...)
    sol = solve(optprob, optimizer; kwargs...)
    xopt = transform(t, sol.u)
    return xopt, sol.minimum
end

function threaded_vida(nstart::Int, prob::VIDAProblem, optimizer; rng=Random.default_rng(), init_params=nothing, unit_cube=true, kwargs...)
    f, t, (lb, ub) = build_opt(prob, unit_cube)
    x0 = initial_point(rng, t, init_params)
    fopt = OptimizationFunction((x,p)->f(x), autodiff)

    tasks = map(1:nstart) do i
        x0 = init_params(rng, t, init_params[i])
        prob = OptimizationProblem(fopt, x0, lb, ub)
        xopt, lp = _vida(fopt, prob, optimizer; kwargs...)
        return xopt, lp
    end

    sols = fetch.(tasks)
    xopts = first.(sols)
    lopts = last.(sols)
    inds = sortperm(lopts)
    return xopts[inds], lops[inds]
end
