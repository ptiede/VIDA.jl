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
lower = [0.1, 0.01, -60.0, -60.0 ,1e-6]
upper = [40.0, 20.0, 60.0, 60.0, 1.0]
θmin, divmin, _, _ = bbextract(bh, θ, lower, upper, TraceMode=:compact)

#Plot the best filter
triptic(images[8], θmin)

#What if we want a different filter?
#Easy! just change the filter function
θgen = GeneralGaussianRing(r0=20.0, σ=2.0, τ=0.1, ξτ = 0.0, s= 0.5, ξs = π/4, x0=0.0, y0=0.0) +
       0.5*Constant()
plot(θgen)

#Need new bounds
lowergen = [0.01, 0.1, 0.001, 0.0, 0.001, -π, -60.0, -60.0, 1e-6]
uppergen = [40.0, 20.1, 0.999, π, 0.999, π, 60.0, 60.0, 1]
θgenmin, divmin, _, _ = bbextract(bh, θgen, lowergen, uppergen, TraceMode=:compact)

triptic(images[8], θgenmin)

#Awesome what if we want to resolve the lumpy structure of the ring
θ26 = CosineRing{2,6}(r0=20.0,
                      σ=[2.0, -0.1], ξσ=[0.0],
                      τ = 0.2, ξτ=π/4,
                      s = [0.5,0.1,0.1,0.3,0.0,0.0], ξs = [π/4,0.0,0.0,0.0,0.0,0.0],
                      x0=0.0, y0=0.0
                      ) + 1.0*Constant()
lower26 = [1.0,
             0.01, -5.0,
             -π,
             0.01, 0.0,
             0.01,-0.99,-0.99,-0.99,-0.99,-0.99,
             -π, -π, -π, -π,-π,-π,
             -60.0, -60.0,
             1e-6]
upper26 = [40.0,
            10.0,5.0,
            π,
            0.99,π,
            0.99,0.99,0.99,0.99,0.99,0.99,
            π,π,π,π,π,π,
            60.0, 60.0,
            1]
θ24min, divmin,_,_ = bbextract(bh,θ24, lower24, upper24, TraceMode=:compact, MaxFuncEvals=30000)
triptic(images[8], θ24min)

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
