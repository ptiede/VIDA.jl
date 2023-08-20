# Interfaces


## Templates

The template describe the image features we wish to extract. 
As of 0.11 VIDA uses the [ComradeBase.jl](https://ptiede.github.io/Comrade.jl/stable/base_api/) 
and [VLBISkyModels.jl](https://ehtjulia.github.io/VLBISkyModels.jl/stable/) interface. This means that any model that obeys that interface can be used within VIDA.
Additionally VIDA defines a number of additional templates that are useful. For a complete list see the [API](@ref) page.

For specifically `VIDA` we have created a `AbstractImageTemplate` subtype of the `ComradeBase.AbstractModel` type and partially implemented some functions. 
For instance, we assume that `ComradeBase.imanalytic(::Type{<:AbstractImageTemplate}) = IsAnalytic()`. As such if an end user wants to implement a new
feature they just need to implement

```julia
ComradeBase.intensity_point(m::MyNewTemplate, p)
ComradeBase.radialextent(m::MyNewTemplate)
```

Note we do not implement the flux of the templates since they are often difficult to calculate.

## Ring Templates

As of VIDA 0.11 we also include a composite image template class called [`RingTemplate`](@ref).

## Divergence `AbstractDivergence`

In order to extract a feature you need to create a probability divergence function. Currently the divergences are defined using a [`AbstractDivergence`](@ref) type. The general user-facing interface is

```julia
bh = Bhattacharyya(image) #make the Bh divergence
kl = KullbackLeibler(image) #makes the KL divergence
```

to initialize the divergence. To evaluate the divergence on a template you use the [`divergence`](@ref) function

```julia
divergence(bh, θ::AbstractTemplate)
divergence(kl, θ::AbstractTemplate)
```


## Extract features `vida` and `VIDAProblem`

The main goal of `VIDA` is to extract image features. To do this we need to define the template and parameterization we want to use. 
The first step is to create a template function that takes in a `NamedTuple` and returns an `<:ComradeBase.AbstractModel`.
For example
```julia
temp(θ) = GaussianRing(θ.r0, θ.σ, θ.x0, θ.y0)
```

For our search we also need to provide the domain over which we want to search
```julia
lower = (r0 = 5.0, σ = 0.1, x0 = -60.0, y0 = 60.0)
upper = (r0 = 30.0, σ = 5.0, x0 = -60.0, y0 = 60.0)
```

We can then form our [`VIDAProblem`](@ref) using the divergence defined above
```julia
prob = VIDAProblem(bh, temp, lower, upper)
```

Finally, to get the optimal parameters and template we can call the [`vida`](@ref) function

```julia
using OptimizationMetaheuristics
xopt, opt_template, divmin = vida(prob, ECA())
```

!!! warn
   Older versions of `VIDA` also included a threaded version of `vida`. This no longer exists
   and would give a race condition if someone tried to use it at this point due to a difference
   in how `divergence` is calculated. If you want to run multiple copies of VIDA at once please
   use Julia's `Distributed` functionality.