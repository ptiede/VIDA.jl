#__precompile__()
"""
    VIDA
Is a image feature extraction tool for use with EHT images of black holes.
It assumes that the image is close to one of the templates we have implemented
and then tries to extract that feature from the image using one of the probability
divergences implemented.
"""
module VIDA


using ChainRulesCore
using DocStringExtensions
using FITSIO
using HDF5
using Interpolations: extrapolate, interpolate, NoInterp, Gridded, Linear,
    AbstractExtrapolation, Flat
using NamedTupleTools
using Random: seed!, rand, GLOBAL_RNG, AbstractRNG
using SpecialFunctions: erf
using StatsBase
using StructArrays: StructArray
using Reexport
@reexport using ComradeBase
@reexport using VLBISkyModels

import ComradeBase as CB
using ComradeBase: SpatialIntensityMap
import ComradeBase: intensity_point

export
    #make the divergences to use for optimization
    Bhattacharyya, KullbackLeibler, LeastSquares, Renyi, NxCorr,
    #Templates
    GaussianRing, SlashedGaussianRing, EllipticalGaussianRing,
    EllipticalSlashedGaussianRing, EllipticalCosineRing, Constant,
    CosineRing, CosineRingwFloor, CosineRingwGFloor, CosineRing,
    GaussianDisk, LogSpiral,
    #Image functions
    load_image, clipimage,
    #Movie functions
    VIDAMovie, load_hdf5, save_hdf5,
    get_image, get_frames, get_times, join_frames

const C0 = 299792458
const KB = 1.38064852e-23

#Load the images
include("images.jl")
#Load the movies
include("movies.jl")
#io stuff
include("io.jl")
#Load the divergence functions
include("divergences.jl")
#Visualization tools


include("extractor.jl")
include("visualizations.jl")
# include("utils.jl")


end #end the module
