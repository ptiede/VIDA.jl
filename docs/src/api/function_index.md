# API

```@meta
CurrentModule = VIDA
```

## Contents

```@contents
Pages = ["function_index.md"]
```

## Index

```@index
Pages = ["function_index.md"]
```



## Filters

These are the various filter types and helper functions.

```@docs
VIDA.AbstractFilter
VIDA.AbstractImageFilter
ImageFilter
LogSpiral
Constant
Disk
AsymGaussian
GaussianRing
SlashedGaussianRing
EllipticalGaussianRing
TIDAGaussianRing
GeneralGaussianRing
CosineRing
AddFilter
MulFilter
Base.size(f::T) where {T<:AbstractFilter}
unpack
stack
split
VIDA.filter_image
```

## Images

VIDA has an image interface that reads in images using the FITS standard.

```@docs
AbstractImage
AbstractFitsImage
EHTImage
load_fits
save_fits
centroid
inertia
pixelloc
flux(::EHTImage)
rescale(::EHTImage, ::Any, ::Any, ::Any)
blur(::EHTImage, ::Any)
Base.size(img::EHTImage)
clipimage
field_of_view
```

## Movies

VIDA also has a movie interface using hdf5. Note that movies are
more than just a list of images. We also use an interpolation between frames.

```@docs
AbstractMovie
EHTMovie
load_hdf5
save_hdf5
join_frames
get_times
get_frames
get_image
flux(::EHTMovie, ::Any)
blur(::EHTMovie, ::Any)
rescale(::EHTMovie,::Any,::Any,::Any)
```

## Divergences

```@docs
AbstractDivergence
Bhattacharyya
KullbackLeibler
```


## Extractor

This defines the interface to the optimizers that can find
the optimal filter for a given image.

```@docs
Optimizer
BBO
CMAES
Opt
ExtractProblem
extractor
threaded_extractor
```

## Utilities

```@docs
make_image
```

## Visualizations

```@docs
triptic
```
