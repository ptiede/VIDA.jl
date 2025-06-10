# # Introduction to VIDA

# Using VIDA is based on constructing three items:
#  1. Data, i.e. an image that you want to extract features from.
#  2. Cost function, i.e. pick if you want to use the KL or BH divergence
#  3. Template, i.e. construct the family of distributions or templates that you will use to approximate the image.
# Then all you need to do is minimize the divergence and you will have extracted you image features.
#
# Now lets runs through how that works

# ## Getting started
# To load VIDA we follow the typical Julia flow. Note that to include plotting functionality
# you need to include Plots as well
using VIDA

using Pkg #hide
Pkg.activate(joinpath(dirname(pathof(VIDA)), "..", "example")) #hide


using CairoMakie
using InteractiveUtils


# ### Step 1 Read in Data
# `VIDA` currently only works with fits images. The fits header is based off of what
# [eht-imaging](https://github.com/achael/eht-imaging) outputs. So as long as you stick
# to that standard you should be fine.

# To read in an image we just use the [`VIDA.load_image`](@ref) function
# which should work with any fits image from ehtim and clean

img = load_image(joinpath(dirname(pathof(VIDA)), "../example/data/example_image.fits"));

# To see what this img is lets print the type
println(typeof(img))

# To plot the image we use the `VLBISkyModels.imageviz` function.
imageviz(img)

# So from the output we see that img is a `IntensityMap` type. This follows ComradeBase and
# VLBISkyModels, see their [docs](https://ehtjulia.github.io/VLBISkyModels.jl/stable/) for
# more information.


# Julia isn't a traditional class OOP language. Namely, methods/functions are first class
# and aren't members of a class. Instead how a function behaves is dependent on the type
# of argument inputs. This is known as something called multimethods or *multiple dispatch*
# where at run-time the type of functions called is determined by the arguments.

# In some sense OOP is just a multiple dispatch type language where the our type of
# dispatch only depends on the class, i.e. the self argument in Python classes.


# Now because of the lack of classes sometimes it can be difficult to figure out which
# functions will act of our datatypes, e.g. the `IntensityMap` type. Fortunately, Julia has some
# convenience functions that let you see which functions/methods can act of the `IntensityMap`
# type
methodswith(IntensityMap)


## Creating a divergence
# In order to find the optimal template you need to first decide on your objective or
# cost function. In VIDA we use probaility divergences to measure differences between the
# template and image. A divergence is defined as an abstract type [`VIDA.AbstractDivergence`](@ref).

bh = Bhattacharyya(img);
kl = KullbackLeibler(img);

# Now to evaluate the divergence we need a template to compare the image to.
# For instance lets create a few different templates
gr = GaussianRing(μas2rad(20.0), μas2rad(5.0), 0.0, 0.0)
ggr = EllipticalSlashedGaussianRing(
    μas2rad(20.0), # r0
    μas2rad(5.0), # σ
    0.2, # τ
    0.78, # ξτ
    0.5, # s
    0.78, # ξs
    μas2rad(-10.0), # x0
    0.0 # y0
)


# We can also plot both templates
g = imagepixels(μas2rad(128.0), μas2rad(128.0), 128, 128)
fig = Figure(; size = (600, 300))
axis = axis = (aspect = 1, xticklabelsvisible = false, yticklabelsvisible = false)
image(fig[1, 1], intensitymap(gr, g), colormap = :afmhot, axis = merge(axis, (title = "GaussianRing",)))
image(fig[1, 2], intensitymap(ggr, g), colormap = :afmhot, axis = merge(axis, (title = "GeneralGaussianRing",)))
fig


# VIDA has a number of templates defined. These are all subtypes of the [`AbstractTemplate`](@ref)
# type. To see which templates are implemented you can use the subtype method:
subtypes(VLBISkyModels.AbstractImageTemplate)

# Additionally as of VIDA 0.11 we can also use any VLBISkyModels model and
# any model that satisfies the interface described [here](https://ehtjulia.github.io/VLBISkyModels.jl/dev/interface/#Model-Interface).

# Using VLBISkyModels interface we can also combine templates together
add = gr + 2.0 * shifted(gr, μas2rad(-10.0), μas2rad(10.0))

# To evaluate the divergence between our template and image we then just evaluate the
# [`VIDA.divergence`](@ref) on the template
@show divergence(kl, add);
@show divergence(kl, ggr);
@show divergence(kl, add);

# Now neither template is really a great approximation to the true image. For instance
# visually they look quite different, which can be checked with the [`VIDA.triptic`](@ref) function

figa, ax = triptic(img, gr)
figa
#-
figb, ax = triptic(img, ggr)
figb
#-
figc, ax = triptic(img, add)
figc

# ## Extracting the Optimal Template
# To extract the optimal template the first thing you need to do is define your
# a function that construct the template and parameterization you will consider
function gr_temp(θ)
    return GaussianRing(θ.r0, θ.σ, θ.x0, θ.y0)
end

# We also want to select the domain that we want to search over
lower = map(μas2rad, (r0 = 5.0, σ = 0.01, x0 = -60.0, y0 = -60.0))
upper = map(μas2rad, (r0 = 60.0, σ = 20.0, x0 = 60.0, y0 = 60.0))

# Given that we can form the [`VIDA.VIDAProblem`](@ref) type.
prob = VIDAProblem(bh, gr_temp, lower, upper);


# Now we need to optimize. VIDA uses the [Optimization.jl](https://github.com/SciML/Optimization.jl)
# meta package for optimization. That means that we can use any optimization package that
# works with optimization. For information about possible optimizers see their [docs](https://docs.sciml.ai/Optimization/stable/).

# For VIDA the classic optimizer is using the [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl).
# To use BlackBox optim we need to load the required package
using OptimizationBBO

# However, for this tutorial we will use the `BlackBoxOptim` optimizer.

# To optimize all you need to do is run the [`VIDA.vida`](@ref) function.

xopt, optfilt, divmin = vida(prob, BBO_de_rand_1_bin_radiuslimited(); maxiters = 100_000)
fig, ax = triptic(img, optfilt)
fig

# Well that seemed to do a terrible job. The reason is that a lot of these
# images tend to have some low level flux throughout the image.
# To account for this the template tends to get very big to absorb some of
# this flux. To combat this you can add a constant background template to
# the problem.
gr_temp_cont(θ) = GaussianRing(θ.r0, θ.σ, θ.x0, θ.y0) + θ.f * VLBISkyModels.Constant((μas2rad(100.0)))
lower = (r0 = μas2rad(5.0), σ = μas2rad(0.01), x0 = μas2rad(-60.0), y0 = μas2rad(-60.0), f = 1.0e-6)
upper = (r0 = μas2rad(60.0), σ = μas2rad(20.0), x0 = μas2rad(60.0), y0 = μas2rad(60.0), f = 10.0)
prob = VIDAProblem(bh, gr_temp_cont, lower, upper);
xopt, optfilt, divmin = vida(prob, BBO_de_rand_1_bin_radiuslimited(); maxiters = 50_000)
#-
fig, ax = triptic(img, optfilt)
fig

# That's much better! Now if you wanted to capture the asymmetry in the ring you can use
# other templates, for example the CosineRing template. Note that this template tends to be
# a little harder to fit.
cos_temp(θ) = EllipticalSlashedGaussianRing(θ.r0, θ.σ, θ.τ, θ.ξτ, θ.s, θ.ξs, θ.x0, θ.y0) + θ.f * θ.f * VLBISkyModels.Constant(μas2rad(100.0))
lower = (r0 = μas2rad(1.0), σ = μas2rad(0.01), τ = 0.0, ξτ = -π / 2, s = 0.001, ξs = -1π, x0 = μas2rad(-60.0), y0 = μas2rad(-60.0), f = 1.0e-6)
upper = (r0 = μas2rad(60.0), σ = μas2rad(20.0), τ = 0.5, ξτ = π / 2, s = 0.999, ξs = 1π, x0 = μas2rad(60.0), y0 = μas2rad(60.0), f = 10.0)
prob = VIDAProblem(bh, cos_temp, lower, upper);
xopt, optfilt, divmin = vida(prob, BBO_de_rand_1_bin_radiuslimited(); maxiters = 50_000);
#-
fig, ax = triptic(img, optfilt)
fig

# Now looks pretty great! To see how to add a custom template see the [Adding a Custom Template](@ref) page.
