@doc """
    $(TYPEDEF)
Abstract optimizer class that defines the optimizer to use.
This is the interface I'll need to define all the different optimizer
classes below.
"""
abstract type Optimizer end

"""
    $(TYPEDEF)
Defines the interface for the BlackBoxOptim interface.
This requires that the user to have imported the BlackBoxOptim package.

# Fields
 - popsize: The population size DEFAULT 64
 - maxevals: The maximum number of times to evaluate the divergence before terminiation.
 - tracemode: The output option. Default is silent, i.e. no output.
 Other options are `:compact` and `:verbose`


# Notes
This uses the default BlackBoxOptim optimizer i.e. adaptive_de_rand_1_bin_radiuslimited.
Currently other options aren't implemented since I found that this version
tended to work the best.
"""
@with_kw struct BBO <: Optimizer
    method::Symbol =:adaptive_de_rand_1_bin_radiuslimited
    popsize::Int = 64
    maxevals::Int = 25_000
    tracemode::Symbol = :compact
end


@doc """
 $(TYPEDEF)
Defines the interface for the LBFGS-B interface from the Optim.jl package.
This requires that the user to have imported the **Optim package**.

 # Fields
  $(FIELDS)

"""
@with_kw struct Opt{O<:Optim.AbstractConstrainedOptimizer} <: Optimizer
    opt::O
    options::Optim.Options = Optim.Options()
end

function Opt(opt::O) where {O<:Optim.AbstractConstrainedOptimizer}
    return Opt(opt=opt, options=Optim.Options())
end

@doc """
    $(TYPEDEF)
Defines the interface for the LBFGS-B interface from the Optim.jl package.
This requires that the user to have imported the **CMAEvolutionStrategy package**.

Typically I have found that this works very well. This the usual first optimizer to try.

# Fields
$(FIELDS)
"""
@with_kw struct CMAES <: Optimizer
    """ PopulationSize """
    popsize::Int = 64
    """ Initial covariance scale """
    cov_scale::Float64 = 1.0
    """ target divergence values. `nothing` means there is no target """
    ftarget::Union{Nothing,Float64} = nothing
    """xtol"""
    xtol::Union{Nothing,Float64} = nothing
    """ftol"""
    ftol::Union{Nothing,Float64} = 1e-11
    """maximum number of divergence evals, nothing means run until termination"""
    maxevals::Union{Nothing,Float64} = nothing
    """verbosity of output"""
    verbosity=1
end



@doc """
    ExtractProblem{T<:AbstractDivergence, S<:AbstractTemplate}
Defines a feature extraction problem to minimize, with an abstract filte
and an abstract divergence. This is needed to interface with the extractor
minimizer, which will minimize the divergence to find the optimal template.

# Fields
$(FIELDS)
"""
struct ExtractProblem{T<:AbstractDivergence, S<:AbstractTemplate}
    """ Divergence function to minimize """
    div::T
    """ Initial location of the optimizer """
    θinit::S
    """ Lower bound of the search region """
    θlower::S
    """ Upper bound of the search region """
    θupper::S
end

@doc """
    threaded_extractor(nstart::Int, prob::ExtractProblem, optimizer::Optimizer)
A threaded multi-start version of the extractor method. This will run `nstart`
instances of extractor, where the initial location of chosen uniformly within the
bounds defined in `prob`.

# Outputs
This outputs the best template and minimum divergence of all the extractors run.
"""
function threaded_extractor(nstart::Int, prob::ExtractProblem{S1,S2}, optimizer::T) where {S1,S2,T<:Optimizer}
    lbounds, ubounds = _bounds(prob)
    θmin = prob.θinit
    divmin = prob.div(θmin)
    Threads.@threads for i in 1:nstart
        p0 = lbounds + (ubounds-lbounds)*rand()
        θ0 = _update(prob.θinit, p0)
        prob = ExtractProblem(prob.div, θ0, prob.θlower, prob.θupper)
        res = extractor(prob, optimizer)
        if divmin > res[2]
            θmin = res[1]
            divmin = res[2]
        end
    end
    return θmin, divmin
end


@doc """
    extractor(prob::ExtractProblem, optimizer::Optimizer)
This extracts the optimal template defined by the `prob` problem.
This will minimize the divergence in prob and return the optimal template
and minimum divergence in a tuple.

`optimizer` is one of VIDA's optimizer types. Typically I would recommend the BBO()
optimizer

# Examples
```julia
    θopt, divmin = extractor(prob, BBO())
```
"""
function extractor(prob, optimizer) end

function extractor(prob::ExtractProblem{S,T}, optimizer::BBO) where {S,T}
    lbounds, ubounds, f, pinit = _create_opts(prob)
    search_range = [ (lbounds[i],ubounds[i])  for i in 1:length(lbounds)]
    ndim = length(lbounds)
    resbb =  bboptimize(f; NumDimensions=ndim,
                        Method=optimizer.method,
                        MaxFuncEvals=optimizer.maxevals, TraceMode=optimizer.tracemode,
                        PopulationSize=optimizer.popsize,
                        SearchRange=search_range)
    θinit2 = _update(prob.θinit, best_candidate(resbb))
    return (θinit2, best_fitness(resbb))

end




function extractor(prob::ExtractProblem{S,T}, optimizer::Opt) where {S,T}
    lbounds, ubounds, f, pinit = _create_opts(prob)
    results =  optimize(f, lbounds, ubounds, pinit, optimizer.opt, optimizer.options)
    θ = _update(prob.θinit, Optim.minimizer(results))
    ℓmax = Optim.minimum(results)
    return (θ, ℓmax)
end


function extractor(prob::ExtractProblem{S,T}, optimizer::CMAES) where {S,T}
    lbounds, ubounds, f, pinit = _create_opts(prob)
    results =  minimize(f, pinit, optimizer.cov_scale;
                        lower=lbounds,
                        upper=ubounds,
                        maxfevals=optimizer.maxevals,
                        xtol=optimizer.xtol,
                        ftol=optimizer.ftol,
                        ftarget=optimizer.ftarget,
                        verbosity=optimizer.verbosity
                    )
    θ = _update(prob.θinit, xbest(results))
    ℓmax = fbest(results)
    return (θ, ℓmax)
end


@inline function _bounds(prob::ExtractProblem)
    return unpack(prob.θlower), unpack(prob.θupper)
end

@inline function _create_opts(prob::ExtractProblem{S,T}) where {S,T}
    lbounds, ubounds = _bounds(prob)
    f(p) = prob.div(_update(prob.θinit, p))
    #unpack the starting location
    pinit = unpack(prob.θinit)
    return lbounds, ubounds, f, pinit
end
