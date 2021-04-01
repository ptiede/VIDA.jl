```@meta
EditURL = "<unknown>/example/introduction.jl"
```

# Introduction to VIDA

Using VIDA is based on constructing three items:
 1. Data, i.e. an image that you want to extract features from.
 2. Cost function, i.e. pick if you want to use the KL or BH divergence
 3. Template, i.e. construct the family of distributions or templates that you will use to approximate the image.
Then all you need to do is minimize the divergence and you will have extracted you image features.

Now lets runs through how that works

## Getting started
To load VIDA we follow the typical Julia flow. Note that to include plotting functionality
you need to include Plots as well

```@example introduction
using Plots
using VIDA
using InteractiveUtils
```

### Step 1 Read in Data
`VIDA` currently only works with fits images. THe fits header is based off of what
[eht-imaging](https://github.com/achael/eht-imaging) outputs. So as long as you stick
to that standard you should be fine.

To read in an image we just use the `load_fits` function
which should work with any fits image from ehtim and clean

```@example introduction
img = load_fits(joinpath(dirname(pathof(VIDA)),"../example/data/example_image.fits"));
nothing #hide
```

To see what this img is lets print the type

```@example introduction
println(typeof(img))
```

To plot the image we can just call plot. This uses recipes and the Plots.jl framework

```@example introduction
plot(img)
```

So from the output we see that img is a EHTImage type. The information in the curly
brackets defines the parametric type information. What this means is that the image
that is constrained in the EHTImage type is a Matrix whose elements are Float64's.

Julia isn't a traditional class OOP language. Namely, methods/functions are first class
and aren't members of a class. Instead how a function behaves is dependent on the type
of argument inputs. This is known as something called multimethods or *multiple dispatch*
where at run-time the type of functions called is determined by the arguments.

In some sense OOP is just a multiple dispatch type language where the our type of
dispatch only depends on the class, i.e. the self argument in Python classes.

Now because of the lack of classes sometimes it can be difficult to figure out which
functions will act of our datatypes, e.g. the EHTImage type. Fortunately, Julia has some
convenience functions that let you see which functions/methods can act of the EHTImage
type

To see what functions can act on an EHTImage object just call

```@example introduction
methodswith(EHTImage)
```

From this list we see there are several methods that can act on EHTImage types.
To see what a certain function does you can type `?inertia` in the terminal to see the help for the inertia method.

```@example introduction
# Creating a divergence
```

In order to find the optimal template you need to first decide on your objective or
cost function. In VIDA we use probaility divergences to measure differences between the
template and image. A divergence is defined as an abstract type `AbstractDivergence`.
In VIDA a divergence is a `functor`. A functor is a type that has an anonymous function
attached to it. That means it is both a type and a function. For instance we create a
divergence by

```@example introduction
 bh = Bhattacharyya(img);
 kl = KullbackLeibler(img);
nothing #hide
```

Now to evaluate the divergence we need to pass it a template.
This can be any template your choose. The great thing about julia is that bh will use
multiple dispatch to figure out which template is being passed to the divergence.

For instance lets create a few different templates

```@example introduction
gr = GaussianRing(r0=20.0, σ=5.0, x0=0.0, y0=0.0)
ggr = GeneralGaussianRing(r0=20.0,
                          σ = 5.0,
                          τ = 0.2,
                          ξτ = 0.78,
                          s = 0.5,
                          ξs = 0.78,
                          x0=0.0,
                          y0=0.0
                        )
```

We can also plot both templates

```@example introduction
a = plot(gr, title="GaussianRing")
b = plot(ggr, title="GeneralGaussianRing")
plot(a, b, layout=(1,2), size=(600,300))
```

VIDA has a number of templates defined. These are all subtypes of the AbstractTemplate type.
To see which templates are implemented you can use the subtype method:

```@example introduction
subtypes(VIDA.AbstractTemplate)
```

Note that the AddTemplate and MulTemplate are internal templates that allow the user to easily combine two templates, for example:

```@example introduction
add = gr + 1.0*ggr
```

To evaluate the divergence between our template and image we then just evaluate the
divergence on the template

```@example introduction
@show bh(gr);
@show bh(ggr);
@show bh(add);
nothing #hide
```

Now neither template is really a great approximation to the true image. For instance
visually they look quite different, which can be checked with the `triptic` function

```@example introduction
a = triptic(img, gr)
b = triptic(img, ggr)
c = triptic(img, add)
plot(a,b,c, layout=(3,1), size=(800,800))
```

## Extracting the Optimal Template
To extract the optimal template the first thing you need to do is define your
`ExtractProblem`. This requires your divergence, initial template, and bounds.

```@example introduction
lower = GaussianRing(r0=0.1, σ=0.01, x0=-60.0, y0=-60.0);
upper = GaussianRing(r0=60.0, σ=20.0, x0=60.0, y0=60.0);
initial = GaussianRing(r0=20.0, σ=5.0, x0=0.0, y0=0.0);

prob = ExtractProblem(bh, initial, lower, upper);
nothing #hide
```

Now to run the optimizers you just need to select which optimizer to use. Currently VIDA
has three families of optimizers installed. Each one is a subtype of the VIDA.Optimizer
abstract type

```@example introduction
subtypes(VIDA.Optimizer)
```

Of the three implemented optimizers my suggestion would be to try the BBO one first.
This uses the [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) package.
BBO has a number of options that can be changed. To see them please use the julia
`?` mode, or see the documentation.

However, for this tutorial I am going to use the `CMAES` optimizer since it is faster
albeit less robust to local minima

To optimize all you need to do is run the extractor function.

```@example introduction
optfilt, divmin = extractor(prob, CMAES())
triptic(img, optfilt)
```

Well that seemed to do a terrible job. The reason is that a lot of these
images tend to have some low level flux throughout the image.
To account for this the template tends to get very big to absorb some of
this flux. To combat this you can add a constant background template to
the problem.

```@example introduction
lower = GaussianRing(r0=0.1, σ=0.01, x0=-60.0, y0=-60.0) + 1e-10*Constant();
upper = GaussianRing(r0=60.0, σ=20.0, x0=60.0, y0=60.0) + 1.0*Constant();
initial = GaussianRing(r0=20.0, σ=5.0, x0=0.0, y0=0.0) + 0.1*Constant();

prob = ExtractProblem(bh, initial, lower, upper);
optfilt, divmin = extractor(prob, CMAES())
triptic(img, optfilt)
```

We can also run multple instances of the extractor, and use Julia's Threads interface

```@example introduction
optfilt, divmin = threaded_extractor(4, prob, CMAES())
```

This will run `4` instances of the extractor function using the available threads in the
Julia session.

That's much better! Now if you wanted to capture the asymmetry in the ring you can use
other templates, for example the CosineRing template. Note that this template tends to be
a little harder to fit.

```@example introduction
lower = CosineRing{1,4}(r0=0.1,
                        σ=[0.1, -1.0], ξσ = [-π],
                        τ = 0.01, ξτ = -π,
                        s = [0.01, -1.0, -1.0, -1.0],
                        ξs = [-π,-π,-π,-π],
                        x0=-60.0, y0=-60.0
                       ) + 1e-10*Constant();
upper = CosineRing{1,4}(r0=40.0,
                        σ=[20.0, 1.0], ξσ = [π],
                        τ = 0.999, ξτ = π,
                        s = [0.999, 1.0, 1.0, 1.0],
                        ξs = [π,π,π,π],
                        x0=60.0, y0=60.0
                       ) + 1.0*Constant();
initial = CosineRing{1,4}(r0=20.0,
                        σ=[5.0, 0.1], ξσ = [0.0],
                        τ = 0.1, ξτ = 0.0,
                        s = [0.1, 0.0, 0.0, 0.0],
                        ξs = [0.0,0.0,0.0,0.0],
                        x0=0.0, y0=0.0
                       ) + 1e-2*Constant();

prob = ExtractProblem(bh, initial, lower, upper);
optfilt, divmin = extractor(prob, CMAES(verbosity=0));
triptic(img, optfilt)
```

Now looks pretty great! To see how to add a custom template see the [Adding a Custom Template](@ref) page.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

