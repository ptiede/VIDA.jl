# # Adding a Custom Template

# If you want to add your own template you just need to define a new:
# - template type
# -
# For example to add a symmetric gaussian template we can use:
using VIDA

struct SlashedExponentialRing{T} <: VIDA.AbstractImageTemplate
    αouter::T
    αinner::T #standard deviation of the Gaussian
    s::T
end

function VIDA.intensity_point(m::SlashedExponentialRing, p)
    (;X, Y) = p
    (;αinner, αouter, s) = m
    r = hypot(X, Y)
    ϕ = atan(X, -Y)

    n = (1-s*cos(ϕ))

    return n*r^αinner/(1 + r^(αouter + αinner+1))
end

# We can also add a convienence constructor
function SlashedExponentialRing(r0, αouter, αinner, s, ξs, x0, y0)
    return modify(SlashedExponentialRing(αouter, αinner, s),
            Stretch(r0, r0), Rotate(ξs), Shift(x0, y0))
end


# Then you can simply call the same optimizing functions and
# plotting functions. For example lets create a fake image and fit it

template = SlashedExponentialRing(μas2rad(20.0), 3.0, 4.0, 0.5, π/2, 0.0, 0.0)

# VIDA uses [`ComradeBase`](https://github.com/ptiede/ComradeBase.jl) and `VLBISkyModels`
# interface. We can create an image using `intensitymap`
img = intensitymap(template, μas2rad(128.0), μas2rad(128.0), 64, 64)

# We can also plot the image
using Plots
plot(img)

# Now lets see if we can get the correct parameters
bh = Bhattacharyya(img);

# To fit we need to define a fitting function. For this our template function
# needs to accept a named tuple.
temp(θ) = SlashedExponentialRing(θ.r0, θ.αout, θ.αin, θ.s, θ.ξs, θ.x0, θ.y0)
# Additionally we need to define the search region for our template extraction
upper = (r0=μas2rad(40.0), αout=10.0, αin = 10.0, s=0.999, ξs=1π, x0= μas2rad(60.0), y0= μas2rad(60.0))
lower = (r0=μas2rad(5.0), αout=1.0, αin = 0.0, s=0.001, ξs=-1π, x0= -μas2rad(60.0), y0= -μas2rad(60.0))

# We can now create our problem
prob = VIDAProblem(bh, temp, lower, upper)

# The vida method can use any optimizer that works with [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/)
# For this work we will use [`CMAEvolutionStrategy`](https://github.com/jbrea/CMAEvolutionStrategy.jl).
using OptimizationBBO
xopt, θopt, divmin = vida(prob, BBO_adaptive_de_rand_1_bin(); maxiters=50_000)

@show θopt

# Let's also plot the results
triptic(img, template)


# Now with all of this said this template actually already exists in VIDA using
# the flexible [`RingTemplate`](@ref).
rad = RadialDblPower(xopt.αin, xopt.αout)
azi = AzimuthalCosine(xopt.s, xopt.ξs)
t   = modify(RingTemplate(rad, azi), Stretch(xopt.r0), Shift(xopt.x0, xopt.y0))
#-

triptic(img, t)
