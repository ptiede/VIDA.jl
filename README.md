# VIDA

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ptiede.github.io/VIDA.jl/dev)
[![Build Status](https://github.com/ptiede/VIDA.jl/workflows/CI/badge.svg)](https://github.com/ptiede/VIDA.jl/actions)
[![Coverage](https://codecov.io/gh/ptiede/VIDA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ptiede/VIDA.jl)


**Warning v0.6 features a breaking change to the optimizer funtionality. If upgrading from a previous version your scripts WILL break**

VIDA.jl or the *Variational Image Domain Analysis* provides a interface to extracting features from fits images created for the EHT, using the notion of probabilty divergences similar to variational inference, hence the name. The currently implemented divergences are the Bhattacharyya distance/divergence as well as the Kullback-Leiber divergence. These are used to extract ring-like features from image reconstructions of black holes such as from M87. A paper applying this to various images ring-like images is in preparation.

## Interface
The inteface is based off of using probability divergences to extract features from images. The main functions a user needs to be aware of are:

``` julia
#load the image
image = load_fits(file::String)
#create the divergence we want (Bhattacharyya divergence)
bh = Bhattacharyya(image)
#Create the filter we want to use to extract for example an slashed elliptical Gaussian
filter = GeneralGaussianRing(p::Array) <: AbstractFilter
#To extract the features we define the ExtractProblem and run extractor
prob = ExtractProblem(divergence, filter::T, filter_lower::T, filter_upper::T) where {T<:AbstractFilter}
opt_filter, divmin = extractor(prob, opt::Optimizer)
```
Let's dive into what each piece means

### Image type: `EHTImage`
 - This is the image you want to fit. Currently we only have support for fits images that are similar to the fits image objects created by [ehtim](https://github.com/achael/eht-imaging). These images are loaded with the function
 ```julia
  #load fits image
  image = load_fits("fitsname.fits")
```
 - There are an additional number of tools available for image processing, such as clipping flux and increasing constrast in the image. That will be included in the documentation when I get around to making it.

### Filter type `AbstractFilter`
 - This is the filter that will be used to extract an image feature. Basically it will find the filter that is the closest to. Currently there are 4 types of filters implemented, but each are a subset of the other:
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
`triptic(img,θ)` that will produce a comparison between the filter and
the true image. This can be useful when comparing the best filter to the 
image.


Additionally any other function that dispatches on the filter type should just work! One thing to note is that the weight between the two filters is relative. Namely, total intensity will always be normalized, so the above code says that θ2 has 5 times the relative flux compared to the first.


### Divergence `AbstractDivergence`
In order to extract a feature you need to create a probability divergence function. Currently the divergences are defined using a AbstractDivergence type. Currently we have two divergences implemented [Bhattacharyya divergence (Bh)](https://en.wikipedia.org/wiki/Bhattacharyya_distance) and the [Kullback-Leiber divergence (KL)](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence). In order to construct the divergence we first need to specify the `image` that we are trying to fit. 
```julia
bh = Bhattacharyya(image) #make the Bh divergence
kl = KullbackLeibler(image) #makes the KL divergence
```
which creates a functor that depends on the image. The functor itself take a filter, e.g.,
```julia
bh(θ::AbstractFilter)
```
and `bh` will use multiple dispatch to figure out which filter function to use.


### Extract features `extractor` and `ExtractProblem` and `Optimizer`
We then use the extractor function to extract the image feature. To call extractor you need to define a `ExtractProblem` type with the call
```julia
prob = ExtractProblem(divergece, filter, filter_lower, filter_upper)
```
where `divergence` is your probability divergence, `filter` is the filter function and `filter_lower` and `filter_upper` are the filters that represent the
lower and upper bounds you want to search over. Once you have your problem defined you can find the optimal filter and extract image features using the `extractor`
function
```julia
opt_filter, divmin = extractor(prob, BBO())
```
Here `BBO()` used the [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) package to optimize. Currently, I have 3 packages that are incorporated with VIDA
 - `BBO()`, [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) which uses genetic and evolutionary algorithms to find global maximum.
 - `CMAES()`, [CMAESEvolutionStrategy.jl](https://github.com/jbrea/CMAEvolutionStrategy.jl) which uses the CMA-ES genetic algorithm to find the global maximum.
 - `Opt()`,  [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl), which allows you to use and box-constrained optimizer from the Optim.jl package, e.g. SAMIN().
In the future I may add interfaces. In terms of which optimizer to select I would suggest to use BBO() as your default, maybe with the maxevals increased from the default `25_000`. For other options you can see the documentation and help.

Note I also have a threaded version of the extractor `threaded_extractor(nstart, prob, BBO())` that run `nstart` runs of extractor, with different initial conditions randomly sampled in the box-constraints.







## A minimal example of extracting ring features
We have provided a minimal example of how to run the filter in examples using command line arguments.
A simpler example is
```julia
  using VIDA
  #load the image and plot it
  image = load_fits("examples/data/elliptical_gaussian_rot-0.00.fits")
  plot(image)

  #Create the filter to use
  filter = GaussianRing(20.0, #20 μas Gaussian ring
                            5.0,  #std dev is 5.0 μas
                            0.0,  #RA (x) location of ring center in μas
                            0.0   #DEC (y) locatin of ring center in μas
                           )
  #Plot the filter
  plot(filter)

  #make the measure you can choose from :KL or :Bh currently.
  bh = Bhattacharyya(image)
  
  #To call the function bh
  bh(filter)

  #Now lets define our extraction problem, i.e. the filter, bounds, and divergence
  lower = GaussianRing(0.1, 0.01, -60.0, -60.0)
  upper = GaussianRing(40.0, 20.0, 60.0, 60.0)
  prob = ExtractProblem(bh, filter, lower, upper)
  #Now we run the extractor with the BlackBoxOptim optimizer
  opt_filter, divmin = extractor(prob, BBO())

  #Now lets run 8 extractors using threads and the CMAES optimizer
  opt_filter, divmin = threaded_extractor(prob, CMAES(ftol=1e-20, cov_scale=10))
  
  #plot the results
  triptic(image, opt_filter)
```

### Image Evaluation
We now also have an ImageFilter type that allows for divergences to be used for image evaluation.

### Distributed computing
In the examples folder we have a complete script that shows how to use VIDA on a cluster to extract image features from multiple images at the same time. It uses argparse to read in command line options and a file that contains a list of paths of images to run VIDA on. For example the bbimage_extractor.jl file will use Julia's Distributed package to split a job amoung a number of cores. To run this you need a file that contains the paths to the list of fits images you would like to run VIDA on. Then you you can select the filter you want using the `--filter` command line option. For the other options type `-h`. If for example you wanted to fit the images with an asymmetric Gaussian filter then you would type:
```bash
julia -p ncores bbimages_extractor.jl list_of_files --filter Asymb
```

Now this script doesn't yet include the cosine ring filter for a slightly technical reason that I still need to fix. Instead if you want to use the CosineRing{N,M} filter you would type

```bash
julia -p ncores bbimages_extractor.jl list_of_files --filter N M
```
where `N` and `M` is the order of the cosine expansion in thickness and brightness. Both of these scripts have been tested in clusters with thousands of cores but if you have any problems please open an issue!
