@doc """
    $(SIGNATURES)
Function that uses Optim.jl to minimize our divergence to extract the image features.
`divergence`. `θinit`,
is the initial filter to use and must be a subtype of AbstractFilter. `lbounds` and
`ubounds` are the upper and lower bounds of the problem. `method` is the maximization
algorithm to use. For a list of availible methos see the Optim.jl package, as well
as for the other various args, and kwargs that can be passed.

# Notes
`θlower` and `θupper` form the lower and upper bounds of the parameter space.
For example, if the filter is `GaussianRing` then `r0` will be bounded in
`[θlower.r0, θupper.r0]`.


# Returns
Returns a tuple (θmin, divmin, converged, iterations), where θmin is the
optimal filter, divmin is the divergence evaluated as the optimal filter,
converged is a boolean that says whether the algorithm converged, and
iterations is the number of iterations it ran for.

"""
function extract(divergence::S, θinit::T, θlower::T, θupper::T,
                 args...; method=Fminbox(BFGS()),
                 kwargs...) where {S<:AbstractDivergence, T<:AbstractFilter}
    lower = unpack(θlower)
    upper = unpack(θupper)
    extract(divergence, θinit, lower, upper, args...; method=method, kwargs...)
end

function extract(divergence, θinit::T, lbounds, ubounds,
                 args...; method=Fminbox(BFGS()),
                 kwargs...) where {T<:AbstractFilter}
    #Construct the function that takes in vector of params
    f(p) = divergence(T(p))
    #unpack the starting location
    pinit = unpack(θinit)
    try
        results =  optimize(f, lbounds, ubounds, pinit,
                            method,
                            args...; kwargs...)
        θ = T(Optim.minimizer(results))
        ℓmax = Optim.minimum(results)
        converged = Optim.converged(results)
        iterations = Optim.iterations(results)
        return (θ, ℓmax, converged, iterations)
    catch e
        println("Extractor failed for initializer:")
        println(pinit)
        println(e)
        return (θinit, f(pinit), false, 1)
    end
end




@doc """
    $(SIGNATURES)
Function that uses Optim.jl to minimize our divergence to extract the image features.
`divergence`. `θinit`,
is the initial filter to use and must be a subtype of AbstractFilter. `lbounds` and
`ubounds` are the upper and lower bounds of the problem. `method` is the maximization
algorithm to use. For a list of availible methos see the Optim.jl package, as well
as for the other various args, and kwargs that can be passed.

# Returns
This returns a data frame with the optimal parameters


# Notes
`θlower` and `θupper` form the lower and upper bounds of the parameter space.
For example, if the filter is `GaussianRing` then `r0` will be bounded in
`[θlower.r0, θupper.r0]`.

This is also the threaded version of the code. If you want to just run a single
case don't pass nstart.

Also most of the filters aren't autodiffable right now so be careful with the autodiff feature.
"""
function extract(nstart::Int, divergence::S, θinit::T, θlower::T, θupper::T,
                 args...; kwargs...) where {S<:AbstractDivergence,T<:AbstractFilter}
    return extract(GLOBAL_RNG, nstart, divergence, θinit, θlower, θupper,
                    args...; kwargs...)
end


function extract(nstart::Int, divergence, θinit::T, lbounds, ubounds,
                 args...; kwargs...) where {T<:AbstractFilter}
  return extract(GLOBAL_RNG, nstart, divergence, θinit, lbounds, ubounds,
                 args...; kwargs...)
end

function extract(rng::AbstractRNG, nstart::Int, divergence::S, θinit::T,
                 θlower::T, θupper::T, args...; kwargs...) where {S<:AbstractDivergence, T<:AbstractFilter}
    lower = unpack(θlower)
    upper = unpack(θupper)
    extract(rng, nstart, divergence, θinit, lower, upper, args...; kwargs...)
end

function extract(rng::AbstractRNG, nstart::Int, divergence, θinit::T,
                 lbounds, ubounds,
                 args...; kwargs...) where {T<:AbstractFilter}
    pinit = unpack(θinit)
    start_states = rand(rng, length(pinit), nstart)

    #Create data frame that holds the results
    df = DataFrame()
    key_names = fieldnames(typeof(θinit))
    for i in 1:length(key_names)
        insertcols!(df,i, Symbol(key_names[i])=>zeros(nstart), makeunique=true)
    end
    setproperty!(df, :ℓmax, zeros(nstart))
    setproperty!(df, :threadid, zeros(nstart))
    setproperty!(df, :converged, Vector{Bool}(undef, nstart))
    setproperty!(df, :iterations, Vector{Int}(undef, nstart))
    @assert length(pinit)==length(ubounds) "Number of parameters doesn't equal number of bounds"
    @assert length(pinit)==length(lbounds) "Number of parameters doesn't equal number of bounds"
    Threads.@threads for i in 1:nstart
        p = zeros(length(pinit))
        for j in 1:length(pinit)
            p[j] = (ubounds[j]-lbounds[j])*start_states[j,i] + lbounds[j]
        end
        θinit = T(p)
        θ, ℓmax, converged, iterations = extract(divergence, θinit,
                             lbounds, ubounds,
                             args..., kwargs...)
        df[i,1:length(p)] = unpack(θ)
        df[i,length(p)+1] = ℓmax
        df[i,length(p)+2] = Threads.threadid()
        df[i,length(p)+3] = converged
        df[i,length(p)+4] = iterations
    end
    return sort!(df, length(pinit)+1)
end

@doc """
    $(SIGNATURES)

Function uses the BlackBoxOptim package to minimize the `divergence` function.
The output from this is then passed to extract to use a deterministic minimizer
to find the true minimum.

# Notes
`θlower` and `θupper` form the lower and upper bounds of the parameter space.
For example, if the filter is `GaussianRing` then `r0` will be bounded in
`[θlower.r0, θupper.r0]`.

# Returns
Returns a tuple (θmin, divmin, converged, iterations), where θmin is the
optimal filter, divmin is the divergence evaluated as the optimal filter,
converged is a boolean that says whether the algorithm converged
(currently always returns false), and
iterations is the number of iterations it ran for.

"""
function bbextract(divergence::S, θ::T, θlower::T, θupper::T,
                   args...; kwargs...) where {S<:AbstractDivergence, T<:AbstractFilter}
    lower = unpack(θlower)
    upper = unpack(θupper)
    return bbextract(divergence, θ, lower, upper, args...; kwargs...)
end

function bbextract(divergence::S, θ::T, lbounds::B, ubounds::B,
                   args...; kwargs...) where {S<:AbstractDivergence, T<:AbstractFilter,B<:AbstractArray}
    @assert length(lbounds)==length(ubounds) "lbounds and ubounds must have the same length"
    @assert length(lbounds)==length(unpack(θ)) "lbounds must have same number of params as θ"
    search_range = [ (lbounds[i],ubounds[i])  for i in 1:length(lbounds)]
    #Create function to optimize
    #println("Using $θ")
    ndim = length(lbounds)
    f(p) = divergence(T(p))
    resbb =  bboptimize(f, args...; NumDimensions=ndim,
                        Method=:adaptive_de_rand_1_bin_radiuslimited,
                        MaxFuncEvals=20000, TraceMode=:silent,
                        SearchRange=search_range,kwargs...)
    θinit2 = T(best_candidate(resbb))
    return (θinit2, best_fitness(resbb), false, 1)
end
