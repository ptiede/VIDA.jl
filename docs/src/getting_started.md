# Getting Started

## A minimal example of extracting ring features

We have provided a minimal example of how to run the template in examples using command line arguments.
A simpler example is

```julia
  using VIDA
  #load the image and plot it
  image = load_fits("examples/data/elliptical_gaussian_rot-0.00.fits")
  plot(image)

  #Create the template to use
  template = GaussianRing(20.0, #20 μas Gaussian ring
                            5.0,  #std dev is 5.0 μas
                            0.0,  #RA (x) location of ring center in μas
                            0.0   #DEC (y) locatin of ring center in μas
                           )
  #Plot the template
  plot(template)

  #make the measure you can choose from :KL or :Bh currently.
  bh = Bhattacharyya(image)
  
  #To call the function bh
  bh(template)

  #Now lets define our extraction problem, i.e. the template, bounds, and divergence
  lower = GaussianRing(0.1, 0.01, -60.0, -60.0)
  upper = GaussianRing(40.0, 20.0, 60.0, 60.0)
  prob = ExtractProblem(bh, template, lower, upper)
  #Now we run the extractor with the BlackBoxOptim optimizer
  opt_template, divmin = extractor(prob, BBO())

  #Now lets run 8 extractors using threads and the CMAES optimizer
  opt_template, divmin = threaded_extractor(prob, CMAES(ftol=1e-20, cov_scale=10))
  
  #plot the results
  triptic(image, opt_template)
```

## Idea behind VIDA

`VIDA` is based on the idea of interpreting the image as a probability distribution. Namely since any image is integrable, the space of images is in one-to-one correspondence with a probability distribution.

Therefore, our idea is very close to variational inference, hence the name *(the) Variational Image Domain Analysis*. Namely, where we view the image as a distribution and we aim to find a approximation of the distribution given some parametric family ``f_\theta(x,y)``, which for our purposes we will typically call a *template*.

The choice of template, depends on the problem of interest, namely what features we are interested in. Typically for the Event Horizon Telescope (EHT) where the images tend to be rings, we are interested in

- Radius r₀
- Width or half width σ
- Structural asymmetry τ
- Flux asymmetry s
- Position angle ξ

`VIDA` then defines a series of templates parameterize these features.

### Templates

Currently we have 6 templates defined, although they all belong to the same family. For an example on how to see the process for defining your own template please see the [Adding a Custom Template](@ref).

The templates implemented are:

- `CosineRing{N,M}` which defines a ring with a cosine expansion in azimuthal thickness (order N) and brightness (order M).
- `GaussianRing` which is a symmetric and circular Gaussian ring.
- `SlashedGaussianRing` which is a circular Gaussian ring with a flux gradient across its emission.
- `EllipticalGaussianRing` symmetric Gaussian elliptical ring, where the emission is constant across the ring, unlike with the SlashedGaussianRing.
- `GeneralGaussianRing` A combination of the two above where the ring is allowed to be elliptical and have a intensity gradient.
- `TIDAGaussianRing` The GeneralGaussianRing, but where the asymmetry and flux orientation are fixed relative to one another.
- `AsymGaussian` A asymmetric Gaussian blob. This can be useful if you image has a strong non-ring component in it.
- `Constant` Adds a constant flux floor to the image. This is very helpful for image reconstructions that tend to add small scale flux everywhere in the image.

### Divergences

In order to extract features we first need a cost function that penalized our parameterized distributions ``f_\theta(x,y)``. Since we are considering the image as a probability distribution, one cost function would be the distance or **divergence** between two distributions. A probability divergence is just a functional that takes in two probability distributions p,q and is minimized iff ``p\equiv q``.

Divergences are defined by the abstract type `AbstractDivergence`. Implementations of the this type are also expected to implement a functor that evaluates the divergence on some template.

Currently we have two divergences implemented in `VIDA`

- Bhattacharyya divergence (Bh)

```math
 Bh(f_\theta|I) = \int \sqrt{f_\theta(x,y)I(x,y)} dxdy.
```

- KL divergence

```math
 KL(f_\theta|I) = \int f_\theta(x,y)\log\left(\frac{f_\theta(x,y)}{I(x,y)}\right)dxdy. 
```

Both divergences give very similar answers, although we found the Bh to be easier to maximize.
