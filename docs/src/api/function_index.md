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

## Templates

These are the various template types in VIDA and VLBISkyModels.

```@docs
LogSpiral
Constant
GaussDisk
RingTemplate
RadialGaussian
RadialDblPower
RadialTruncExp
AzimuthalUniform
AzimuthalCosine
GaussianRing
CosineRing
```

In addition in VIDA we define a number of helper functions for commonly used templates
```@docs
SlashedGaussianRing
EllipticalGaussianRing
CosineRingwFloor
CosineRingwGFloor
EllipticalCosineRing
```

## Images

VIDA has an image interface that reads in images using the FITS standard.

```@docs
load_image
VIDA.blur(::ComradeBase.IntensityMap, ::Any)
clipimage
```

## Movies

VIDA also has a movie interface using hdf5. Note that movies are
more than just a list of images. We also use an interpolation between frames.

```@docs
VIDA.VIDAMovie
VIDA.load_hdf5
VIDA.save_hdf5
VIDA.join_frames
VIDA.get_times
VIDA.get_frames
VIDA.get_image
VIDA.flux(::VIDAMovie, ::Any)
VIDA.blur(::VIDAMovie, ::Any)
VLBISkyModels.regrid(::VIDAMovie)
```

## Divergences

```@docs
AbstractDivergence
Bhattacharyya
KullbackLeibler
Renyi
LeastSquares
```

## Extractor

This defines the interface to the optimizers that can find
the optimal template for a given image.

```@docs
VIDAProblem
vida
```


## Visualizations

```@docs
triptic
```
