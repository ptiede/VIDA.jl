# Scripts Help

These scripts use [ArgParse.jl](https://github.com/carlobaldassi/ArgParse.jl) to read in command line options and a file that contains a list of paths of images to run VIDA on. For example the bbimage_extractor.jl file will use Julia's Distributed package to split a job amoung a number of cores. To run this you need a file that contains the paths to the list of fits images you would like to run VIDA on. Then you you can select the filter you want using the `--filter` command line option. For the other options type `-h`. If for example you wanted to fit the images with an asymmetric Gaussian filter then you would type:
```bash
julia -p ncores bbimages_extractor.jl list_of_files --filter Asymb
```

Now this script doesn't yet include the cosine ring filter for a slightly technical reason that I still need to fix. Instead if you want to use the CosineRing{N,M} filter you would type

```bash
julia -p ncores bbimages_extractor.jl list_of_files --filter N M
```
where `N` and `M` is the order of the cosine expansion in thickness and brightness. Both of these scripts have been tested in clusters with thousands of cores but if you have any problems please open an issue!