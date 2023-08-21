# VIDA

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ptiede.github.io/VIDA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ptiede.github.io/VIDA.jl/dev)
[![Build Status](https://github.com/ptiede/VIDA.jl/workflows/CI/badge.svg)](https://github.com/ptiede/VIDA.jl/actions)
[![Coverage](https://codecov.io/gh/ptiede/VIDA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ptiede/VIDA.jl)
[![DOI](https://zenodo.org/badge/271097921.svg)](https://zenodo.org/badge/latestdoi/271097921)


`VIDA.jl` or the *Variational Image Domain Analysis* provides a interface to extracting features from fits images created for the EHT, using the notion of probability divergences similar to variational inference. This package has been used extensively in the most recent EHT publications of [Sgr A*](https://iopscience.iop.org/journal/2041-8205/page/Focus_on_First_Sgr_A_Results).

As of **version 0.11** VIDA now uses the [`ComradeBase`](https://ptiede.github.io/Comrade.jl/dev/base_api/) and [`VLBISkyModels`](https://ehtjulia.github.io/VLBISkyModels.jl/dev/interface/) interface. This implies that any model that satisifies this interface can 
now be used as a VIDA template.

For more information on how to use VIDA please refer to the [documentation](https://ptiede.github.io/VIDA.jl/stable)

## Installation

`VIDA.jl` is registered so to install it just go into the REPL and type `]add VIDA`.


## Distributed scripts

In the `example` folder we have a complete script that shows how to use VIDA on a cluster to extract image features from multiple images at the same time. Please see the [README](https://github.com/ptiede/VIDA.jl/tree/master/scripts) in the scripts folder for more information.

# Citation

If you use this work in a paper please use the following citation:

```bibtex
@ARTICLE{VIDA,
       author = {{Tiede}, Paul and {Broderick}, Avery E. and {Palumbo}, Daniel C.~M.},
        title = "{Variational Image Feature Extraction for the Event Horizon Telescope}",
      journal = {ApJ},
     keywords = {1338, 1647, 16, 162, 1855, 293},
         year = 2022,
        month = feb,
       volume = {925},
       number = {2},
          eid = {122},
        pages = {122},
          doi = {10.3847/1538-4357/ac3a6b},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022ApJ...925..122T},
}
```
