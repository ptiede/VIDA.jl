using VIDA
using Plots #for plotting

#load the images
files = readdir("example/data/", join=true)
filter!(x->endswith(x, ".fits"), files)

#Load the images using Julia's . fusing syntax
images = load_ehtimfits.(files)

#Lets plots a random image
plot(images[8])

#Now constuct a divergence
bh = Bhattacharyya(images[8])
kl = KullbackLeibler(images[8])

#This accepts a filter as its argument
θ = GaussianRing(r0=20.0, σ=10.0, x0=0.0, y0=0.0) +
        0.1*Constant()

@show bh(θ)
#Plot the filter
plot(θ)
#compare the filter and the image
triptic(images[8], θ)

#Now we want to extract the optimal filter
#bbextract needs the divergence, filter, and bounds
θlower = GaussianRing(r0=0.1, σ=0.01, x0 = -60.0, y0=-60.0) + 1e-6*Constant()
θupper = GaussianRing(r0=40.0,σ=20.0, x0 = 60.0,  y0= 60.0) + 1.0*Constant()
θmin, divmin, _, _ = bbextract(bh, θ, θlower, θupper, TraceMode=:compact)

#Plot the best filter
triptic(images[8], θmin)

#What if we want a different filter?
#Easy! just change the filter function
θgen = GeneralGaussianRing(r0=20.0, σ=2.0, τ=0.1, ξτ = 0.0, s= 0.5, ξs = π/4, x0=0.0, y0=0.0) +
       0.5*Constant()
plot(θgen)

#Need new bounds
θlowergen = GeneralGaussianRing(r0 = 0.01,
                               σ = 0.1,
                               τ = 0.001, ξτ = 0.0,
                               s = 0.001, ξs = -π,
                               x0 = -60.0, y0 = -60.0) + 1e-6*Constant()
θuppergen = GeneralGaussianRing(r0 = 40.0,
                               σ = 20.0,
                               τ = 0.999, ξτ = π,
                               s = 0.999, ξs = π,
                               x0 = 60.0, y0 = 60.0) + 1.0*Constant()
θgenmin, divmin, _, _ = bbextract(bh, θgen, θlowergen, θuppergen, TraceMode=:compact)

triptic(images[8], θgenmin)

#Awesome what if we want to resolve the lumpy structure of the ring
θ16 = CosineRing{1,6}(r0=20.0,
                      σ=[2.0, -0.1], ξσ=[0.0],
                      τ = 0.2, ξτ=π/4,
                      s = [0.5,0.1,0.1,0.3,0.0,0.0], ξs = [π/4,0.0,0.0,0.0,0.0,0.0],
                      x0=0.0, y0=0.0
                      ) + 1.0*Constant()


θlower16 = CosineRing{1,6}(r0=1.0,
             σ = [0.01, -5.0],
             ξσ = [-π],
             τ = 0.01, ξτ = 0.0,
             s = [0.01,-0.99,-0.99,-0.99,-0.99,-0.99],
             ξs = [-π, -π, -π, -π,-π,-π],
             x0 = -60.0, y0 = -60.0) + 1e-6*Constant()
θupper16 = CosineRing{1,6}(r0=40.0,
            σ = [10.0,5.0],
            ξσ = [π],
            τ = 0.99, ξτ = π,
            s = [0.99,0.99,0.99,0.99,0.99,0.99],
            ξs = [π, π, π, π, π, π],
            x0 = 60.0, y0 = 60.0) + 1.0*Constant()

θ16min,divmin,_,_ = bbextract(bh, θ16, θlower16, θupper16; TraceMode=:compact, MaxFuncEvals=30000)
triptic(images[8], θ16min)

#That's cool now what if we want to fit all those files?
#Well it would be great if we could thread everything!
#and thanks to Julia this is actually pretty easy
results = []
bhs = Bhattacharyya.(images)
Threads.@threads for i in 1:length(images)
    res = bbextract(bhs[i], θ, lower, upper, MaxFuncEvals=30000)
    push!(results, res[1])
end



#To see all the types of implemented filters we can use the subtype function
@show subtypes(VIDA.AbstractFilter)
