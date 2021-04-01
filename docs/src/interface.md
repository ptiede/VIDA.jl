# Interface

The interface is based off of using probability divergences to extract features from images. The main functions a user needs to be aware of are:

``` julia
#load the image
image = load_fits(file::String)
#create the divergence we want (Bhattacharyya divergence)
bh = Bhattacharyya(image)
#Create the template we want to use to extract for example an slashed elliptical Gaussian
template = GeneralGaussianRing(p::Array) <: AbstractTemplate
#To extract the features we define the ExtractProblem and run extractor
prob = ExtractProblem(divergence, template::T, template_lower::T, template_upper::T) where {T<:AbstractTemplate}
opt_template, divmin = extractor(prob, opt::Optimizer)
```

Let's dive into what each piece means

## Image type: `EHTImage`

- This is the image you want to fit. Currently we only have support for fits images that are similar to the fits image objects created by [ehtim](https://github.com/achael/eht-imaging). These images are loaded with the function

 ```julia
  #load fits image
  image = load_fits("fitsname.fits")
```

- There are an additional number of tools available for image processing, such as clipping flux and rescaling the image. For a list of such functions please refer to the [API](@ref) page.

## Template type `AbstractTemplate`

- This is the template that will be used to extract an image feature. Basically it will find the template that is the closest to. Currently there are 4 types of templates implemented, but each are a subset of the other:

```julia
    GaussianRing(p) #Gaussian circular ring
    SlashedGaussianRing(p) #Gaussian circular ring with flux asymmetry
    EllipticalGaussianRing(p) #Gaussian elliptical ring
    TIDAGaussianRing(p)  #Gaussian elliptical ring with flux asymmetry where the orientations relative to each other are fixed
    GeneralGaussianRing(p) #Gaussian elliptical slashed ring where orientation of the slash and asymmetry are independent.
```

- For the parameters that it takes please use the julia ? mode.

The plotting is done through the recipes macros in Plots.jl. So it should
just work! In addition to the plot function there is a new recipe called
`triptic(img,θ)` that will produce a comparison between the template and
the true image. This can be useful when comparing the best template to the
image.

Additionally any other function that dispatches on the template type should just work! One thing to note is that the weight between the two templates is relative. Namely, total intensity will always be normalized, so the above code says that θ2 has 5 times the relative flux compared to the first.

## Divergence `AbstractDivergence`

In order to extract a feature you need to create a probability divergence function. Currently the divergences are defined using a AbstractDivergence type. Currently we have two divergences implemented [Bhattacharyya divergence (Bh)](https://en.wikipedia.org/wiki/Bhattacharyya_distance) and the [Kullback-Leibler divergence (KL)](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence). In order to construct the divergence we first need to specify the `image` that we are trying to fit.

```julia
bh = Bhattacharyya(image) #make the Bh divergence
kl = KullbackLeibler(image) #makes the KL divergence
```

which creates a functor that depends on the image. The functor itself take a template, e.g.,

```julia
bh(θ::AbstractTemplate)
```

and `bh` will use multiple dispatch to figure out which template function to use.

## Extract features `extractor` and `ExtractProblem` and `Optimizer`

We then use the extractor function to extract the image feature. To call extractor you need to define a `ExtractProblem` type with the call

```julia
prob = ExtractProblem(divergence, template, template_lower, template_upper)
```

where `divergence` is your probability divergence, `template` is the template function and `template_lower` and `template_upper` are the templates that represent the
lower and upper bounds you want to search over. Once you have your problem defined you can find the optimal template and extract image features using the `extractor`
function

```julia
opt_template, divmin = extractor(prob, BBO())
```

Here `BBO()` used the [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) package to optimize. Currently, I have 3 packages that are incorporated with VIDA

- `BBO()`, [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) which uses genetic and evolutionary algorithms to find global maximum.
- `CMAES()`, [CMAESEvolutionStrategy.jl](https://github.com/jbrea/CMAEvolutionStrategy.jl) which uses the CMA-ES genetic algorithm to find the global maximum.
- `Opt()`,  [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl), which allows you to use and box-constrained optimizer from the Optim.jl package, e.g. SAMIN().

In the future I may add interfaces. In terms of which optimizer to select I would suggest to use BBO() as your default, maybe with the `maxevals` increased from the default `25_000`. For other options you can see the documentation and help.

Note there is a threaded version of the extractor `threaded_extractor(nstart, prob, BBO())` that run `nstart` runs of extractor, with different initial conditions randomly sampled in the box-constraints.
