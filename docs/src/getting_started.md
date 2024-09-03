# Getting Started

## A minimal example of extracting ring features

The basic VIDA program mirrors the following structure

```julia
using VIDA, CairoMakie
# load the image and plot it
image = load_fits("example/data/elliptical_gaussian_rot-0.00.fits")
imageviz(image)
# Build the divergence we want to fit
bh = Bhattacharyya(image)
# Create the template to use
template(θ) = SlashedGaussianRing(θ.r0, θ.σ, θ.s, θ.ξ, θ.x0, θ.y0) + θ.floor*Constant(μas2rad(100.0))

#Define our bounds
lower = (r0 = μas2rad(5.0), σ=μas2rad(1.0), 
         s=0.001, ξ=-1π, 
         x0=-μas2rad(60.0), y0 = -μas2rad(60.0), 
         floor=1e-6)
upper = (r0 = μas2rad(30.0), σ=μas2rad(15.0), 
         s=0.999, ξ=1π, 
         x0=μas2rad(60.0), y0 = μas2rad(60.0), 
         floor=100.0)

prob = VIDAProblem(bh, template, lower, upper)
# Load your optimizer and run VIDA
using OptimizationBBO
xopt, opt_temp, divmin = vida(prob, BBO_adaptive_de_rand_1_bin(); maxiters=50_000)

#plot the results
fig = triptic(image, opt_temp)
fig
```

## Idea behind VIDA

`VIDA` is based on the idea of interpreting the image as a probability distribution. Namely since any image is integrable, the space of images is in one-to-one correspondence with a probability distribution.

Therefore, our idea is very close to variational inference, hence the name *(the) Variational Image Domain Analysis*. Namely, where we view the image as a distribution and we aim to find a approximation of the distribution given some parametric family ``f_\theta(x,y)``, which for our purposes we will typically call a *template*.

The choice of template, depends on the problem of interest, namely what features we are interested in. Typically for the Event Horizon Telescope (EHT) where the images tend to be rings, we are interested in

- Radius r₀
- Width or half width σ
- Structural asymmetry τ
- Brightness asymmetry s
- Position angle ξ

`VIDA` then defines a series of templates parameterize these features.

### Templates

Currently we have 6 templates defined, although they all belong to the same family. For an example on how to see the process for defining your own template please see the [Adding a Custom Template](@ref).

### Divergences

In order to extract features we first need a cost function that penalized our parameterized distributions ``f_\theta(x,y)``. Since we are considering the image as a probability distribution, one cost function would be the distance or **divergence** between two distributions. A probability divergence is just a functional that takes in two probability distributions p,q and is minimized iff ``p\equiv q``.

Divergences are defined by the abstract type [`VIDA.AbstractDivergence`](@ref). Implementations of the this type are also expected to implement a functor that evaluates the divergence on some template.

The current recommended default template is the [`VIDA.Bhattacharyya`](@ref) divergence although all the template give similar answers.
