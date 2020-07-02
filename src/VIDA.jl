__precompile__()
"""
    VIDA
Is a image feature extraction tool for use with EHT images of black holes.
It assumes that the image is close to one of the filters we have implemented
and then tries to extract that feature from the image using one of the probability
divergences implemented.
"""
module VIDA

export
    #make the divergences to use for optimization
    Bhattacharyya, KullbackLeibler,
    #plotting stuff I created
    plot, plot_triptic,
    #All the filters I created
    GaussianRing,SlashedGaussianRing,EllipticalGaussianRing,
    TIDAGaussianRing,GeneralGaussianRing, Constant, AsymGaussian,
    cat,split,unpack,
    #Imaging stuff
    EHTImage, load_ehtimfits, clipimage, save_ehtimfits,
    flux, centroid, inertia,
    #Feature extraction stuff
    extract, bbextract

using BlackBoxOptim
using DataFrames
using DocStringExtensions
using FITSIO
using LaTeXStrings
using Optim
using PhysicalConstants.CODATA2018: c_0, h, k_B
using Parameters
using Plots
using Plots.PlotMeasures: mm
using Random: seed!,rand, GLOBAL_RNG, AbstractRNG
using RecipesBase
using SpecialFunctions:erf

#Load the images
include("images.jl")
#Load the visualization stuff
include("filters.jl")
#Load the divergence functions
include("divergences.jl")
#Visualization tools
include("visualizations.jl")
#filter extractor
include("extractor.jl")


end #end the module
