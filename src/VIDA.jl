#__precompile__()
"""
    VIDA
Is a image feature extraction tool for use with EHT images of black holes.
It assumes that the image is close to one of the filters we have implemented
and then tries to extract that feature from the image using one of the probability
divergences implemented.
"""
module VIDA



using BlackBoxOptim
using CMAEvolutionStrategy
using DataFrames
using DocStringExtensions
using FITSIO
using HDF5
using ImageFiltering: imfilter, Kernel.gaussian, Fill, Algorithm.FFT
using Interpolations
using LaTeXStrings
using Optim
using Parameters
using Random: seed!,rand, GLOBAL_RNG, AbstractRNG
using RecipesBase
using Requires
using SpecialFunctions:erf

export
    #make the divergences to use for optimization
    Bhattacharyya, KullbackLeibler,
    #Filters
    GaussianRing,SlashedGaussianRing,EllipticalGaussianRing,
    TIDAGaussianRing,GeneralGaussianRing, Constant, AsymGaussian,
    CosineRing,Disk,ImageFilter,LogSpiral,
    #Filter helper functions
    stack,split,unpack,
    #Image functions
    EHTImage, load_ehtimfits, load_fits, clipimage, save_fits,
    flux, centroid, inertia, rescale_image, get_radec,
    rescale, blur,
    #Movie functions
    EHTMovie, load_hdf5, save_hdf5,
    get_image, get_frames, get_times, join_frames,
    #Optimizers
    extractor, threaded_extractor, BBO, CMAES, Opt,
    #ExtractionProblem
    ExtractProblem


const C0 = 299792458
const KB = 1.38064852e-23

#Load the images
include("images.jl")
#Load the movies
include("movies.jl")
#Load the visualization stuff
include("filters.jl")
#Load the divergence functions
include("divergences.jl")
#Visualization tools
function __init__()
#    @require BlackBoxOptim="a134a8b2-14d6-55f6-9291-3336d3ab0209"
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("visualizations.jl")
end
#include("visualizations.jl")
#filter extractor
include("extractor.jl")
include("utils.jl")



end #end the module
