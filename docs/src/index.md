# VIDA.jl

Variational image domain analysis for the EHT.
This package is for extracting features, such as ring from
image reconstruction of EHT data. Currently these images must be in fits format although other types may be included in the future.

## Installation

`VIDA` is a registered Julia package to install

```julia
using Pkg; Pkg.add("VIDA")
```

or go to the repl and simply type `]add VIDA`. Note that we require a Julia version >= 1.4.

Some additional dependencies that enable full functionality can be added with

```julia
Pkg.add.(["Plots","ArgParse"])
```

`Plots.jl` is required to use some of the plotting recipes defined in the package and ArgParse is used for some of the scripts in the example folder.

To extract features we use probability divergences to characterize differences between our image
and some approximation. For the divergences implemented see the [Getting Started](@ref) page. The idea is then very similar to variational inferences where we pick a parametric family of distributions which we call *templates* and then try to find the template that minimizes the divergence. For the templates that are currently implemented please see the page.

See the [API](@ref) for the complete list of documented functions and types.

## Outline

```@contents
Pages = [
  "index.md",
  "getting_started.md",
  "interface.md",
  "generated/introduction.md",
  "generated/custom_template.md",
  "api/function_index.md"
]
```
