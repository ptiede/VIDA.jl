# Getting Started
The best way to learn how to use VIDA is to look at some of the example notebooks provided.

## Installation
`VIDA` is a registered Julia package to install
```julia
using Pkg; Pkg.add("VIDA")
```
Or go to the repl and simply type `]add VIDA`. Note that we require a Julia version >= 1.4.

Some additional dependencies that enable full functionality can be added with
```julia
Pkg.add.(["Plots","ArgParse"])
```
Plots.jl is required to use some of the plotting recipes defined in the package and ArgParse is used for some of the scripts in the example folder.

## Idea behind VIDA
`VIDA` is based on the idea of interpreting the image as a probability distribution. Namely since any image is integrable, the space of images is in one-to-one correspondence with a probability distribution.

Therefore, our idea is very close to variational inference, hence the name *(the) Variational Image Domain Analysis*. Namely, where we view the image as a distribution and we aim to find a approximation of the distribution given some parametric family ``f_\theta(x,y)``, which for our purposes we will typically call a *filter*.

The choice of filter, depends on the problem of interest, namely what features we are interested in. Typically for the Event Horizon Telescope (EHT) where the images tend to be rings, we are interested in

 - Radius r₀
 - Width or half width σ
 - Structural asymmetry τ
 - Flux asymmetry s
 - Position angle ξ

`VIDA` then defines a series of filters parameterize these features.

### Filters
Currently we have 6 filters defined, although they all belong to the same family. For an example on how to see the process for defining your own filter please see the [readme](https://github.com/ptiede/VIDA.jl/blob/master/README.md).

The filters implemented are:

 - `CosineRing{N,M}` which defines a ring with a cosine expansion in azimuthal thickness (order N) and brightness (order M).
 - `GaussianRing` which is a symmetric and circular Gaussian ring.
 - `SlashedGaussianRing` which is a circular Gaussian ring with a flux gradient across its emission.
 - `EllipticalGaussianRing` symmetric Gaussian elliptical ring, where the emission is constant across the ring, unlike with the SlashedGaussianRing.
 - `GeneralGaussianRing` A combination of the two above where the ring is allowed to be elliptical and have a intensity gradient.
 - `TIDAGaussianRing` The GeneralGaussianRing, but where the asymmetry and flux orienation are fixed relative to one another.
 - `AsymGaussian` A asymmetric Gaussian blob. This can be useful if you image has a strong non-ring component in it.
 - `Constant` Adds a constant flux floor to the image. This is very helpful for image reconstructions that tend to add small scale flux everywhere in the image.

### Divergences
In order to extract features we first need a cost function that penalized our parameterized distributions ``f_\theta(x,y)``. Since we are considering the image as a probability distribution, one cost function would be the distance or **divergence** between two distributions. A probability divergence is just a functional that takes in two probability distributions p,q and is minimized iff ``p\equiv q``. 

Divergences are defined by the abstract type `AbstractDivergence`. Implementations of the this type are also expected to implement a functor that evaluates the divergence on some filter. For an example see the implementation of the divergences in 

Currently we have two divergences implemented in `VIDA`
 - Bhattacharyya divergence (Bh)

```math
 Bh(f_\theta|I) = \int \sqrt{f_\theta(x,y)I(x,y)} dxdy.
```

 - KL divergence 
```math
 KL(f_\theta|I) = \int f_\theta(x,y)\log\left(\frac{f_\theta(x,y)}{I(x,y)}\right)dxdy. 
```
Both divergences give very similar answers, although we found the BH to be easier to maximize.



For an example of how all this all works is given in the examples folder.
