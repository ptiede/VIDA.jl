# # Adding a Custom Template

# If you want to add your own template you just need to define a new:
# - template type
# - size method
# - imagetemplate method for that type of template.
#
# For example to add a symmetric gaussian template we can use:
using VIDA

Base.@kwdef struct SymGaussian <: VIDA.AbstractImageTemplate
   σ::Float64 #standard deviation of the Gaussian
   x0::Float64 #x location of mean
   y0::Float64 #y location of mean
end

SymGaussian(p) = SymGaussian(p[1],p[2],p[3])

Base.size(::Type{SymGaussian}) = 3

function (θ::SymGaussian)(x,y)
  z2 = ((x-θ.x0)^2 + (y-θ.y0)^2)/(2.0*θ.σ^2)
  return 1.0/(2.0*π*θ.σ^2)*exp(-z2)
end

# Then you can simply call the same optimizing functions and
# plotting functions. For example lets create a fake image and fit it

template = SymGaussian(σ=20.0, x0=0.0, y0=0.0)

# Now I will use a utility function to convert a template to an EHTImage
img = VIDA.make_image(template, 64, (-60.0,60.0), (-60.0,60.0));


# Now lets see if we can get the correct parameters
bh = Bhattacharyya(img);

# Define the starting point for the optimization and the bounds
start = SymGaussian(rand(3))
upper = SymGaussian(σ=40.0, x0=60.0, y0=60.0)
lower = SymGaussian(σ=0.001, x0=-60.0, y0=-60.0)

prob = ExtractProblem(bh, start, lower, upper)

θopt, divmin = extractor(prob, CMAES())

@show θopt

# Let's also plot the results

triptic(img, template)
