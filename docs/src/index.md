# VIDA.jl
Variational image domain analysis for the EHT. 
This package is for extracting features, such as ring from
image reconstruction of EHT data. Currently these images must be in fits format although other types 
my be included in the future.

In order to extract features we use probability divergences to characterize differences between our image
and some approximation. For the divergences implemented see the [Getting Started](@ref) page. The idea is then very 
similar to variational inferences where we pick a parametric family of distributions which we call *filters* and
then try to find the filter that minimizes the divergence. For the filters that are currently implemented 
please see the page.

See the for the complete list of documented functions and types.


## Outline
```@contents
Pages = [
  "index.md",
  "getting_started.md",
  "function_index.md"
]
```
