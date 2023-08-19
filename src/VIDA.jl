#__precompile__()
"""
    VIDA
Is a image feature extraction tool for use with EHT images of black holes.
It assumes that the image is close to one of the templates we have implemented
and then tries to extract that feature from the image using one of the probability
divergences implemented.
"""
module VIDA



import ComradeBase as CB
using DocStringExtensions
using FITSIO
using HDF5
using ImageFiltering: imfilter, Kernel.gaussian, Fill, Algorithm.FFT
using Interpolations
using LaTeXStrings
using Random: seed!,rand, GLOBAL_RNG, AbstractRNG
using RecipesBase
using Requires
using SpecialFunctions: erf
using Reexport
@reexport using ComradeBase
@reexport using VLBISkyModels

using ComradeBase: SpatialIntensityMap
import ComradeBase: intensity_point

export
    #make the divergences to use for optimization
    Bhattacharyya, KullbackLeibler, LeastSquares,Renyi,
    #Templates
    GaussianRing,SlashedGaussianRing,EllipticalGaussianRing,
    EllipticalSlashedGaussianRing, EllipticalCosineRing, Constant, AsymGaussian,
    SymCosineRing, SymCosineRingwFloor, SymCosineRingwGFloor, CosineRing,
    GaussianDisk, LogSpiral,
    #Image functions
    EHTImage, load_image, clipimage,
    inertia, blur,
    #Movie functions
    VIDAMovie, load_hdf5, save_hdf5,
    get_image, get_frames, get_times, join_frames,
    #Optimizers
    # extractor, threaded_extractor, BBO, CMAES, Opt,
    #ExtractionProblem
    ExtractProblem


const C0 = 299792458
const KB = 1.38064852e-23

#Load the images
include("images.jl")
#Load the movies
include("movies.jl")
#io stuff
include("io.jl")
#Load the visualization stuff
include("templates/templates.jl")
#Load the divergence functions
include("divergences.jl")
#Visualization tools
function __init__()
#    @require BlackBoxOptim="a134a8b2-14d6-55f6-9291-3336d3ab0209"
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("visualizations.jl")
end
#include("visualizations.jl")
#template extractor
include("extractor.jl")
# include("utils.jl")



end #end the module
