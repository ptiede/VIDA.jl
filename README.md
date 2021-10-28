# VIDA

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ptiede.github.io/VIDA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ptiede.github.io/VIDA.jl/dev)
[![Build Status](https://github.com/ptiede/VIDA.jl/workflows/CI/badge.svg)](https://github.com/ptiede/VIDA.jl/actions)
[![Coverage](https://codecov.io/gh/ptiede/VIDA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ptiede/VIDA.jl)


`VIDA.jl` or the *Variational Image Domain Analysis* provides a interface to extracting features from fits images created for the EHT, using the notion of probability divergences similar to variational inference, hence the name. The currently implemented divergences are the Bhattacharyya distance/divergence as well as the Kullback-Leibler divergence. These are used to extract ring-like features from image reconstructions of black holes such as from M87. 

For more information on how to use VIDA please refer to the [documentation](https://ptiede.github.io/VIDA.jl/stable)

## Installation

`VIDA.jl` is registered so to install it just go into the REPL and type `]add VIDA`.


## Distributed scripts

In the examples folder we have a complete script that shows how to use VIDA on a cluster to extract image features from multiple images at the same time. Please see the [README](https://github.com/ptiede/VIDA.jl/tree/master/scripts) in the scripts folder for more information.
