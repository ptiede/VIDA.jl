module VIDAOptimizationExt

using VIDA
using Optimization: OptimizationFunction, OptimizationProblem, solve, SciMLBase
using Random

function VIDA.vida(prob::VIDAProblem, optimizer; rng=Random.default_rng(), init_params=nothing, unit_cube=true, kwargs...)
    T = eltype(prob.div.img)
    f, t, (lb, ub) = VIDA.build_opt(prob, unit_cube)
    x0 = VIDA.initial_point(rng, T, t, init_params)
    ad = isnothing(prob.autodiff) ? SciMLBase.NoAD() : prob.autodiff
    fopt = OptimizationFunction((x,p)->f(x), ad)
    optprob = OptimizationProblem(fopt, x0, nothing; lb=lb, ub=ub)
    xopt, min =  _vida(fopt, t, optprob, optimizer; kwargs...)
    return xopt, prob.f(xopt), min
end


function _vida(fopt, t, optprob, optimizer; kwargs...)
    sol = solve(optprob, optimizer; kwargs...)
    xopt = VIDA.HypercubeTransform.transform(t, sol.u)
    return xopt, sol.objective
end


end