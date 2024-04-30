# Scripts

These scripts use [ArgParse.jl](https://github.com/carlobaldassi/ArgParse.jl) to read in command line options and a file that contains a list of paths of images to run VIDA on. For example the `image_extractor.jl` script will use Julia's Distributed package to split a job among a number of cores. To run this you need a file that contains the paths to the list of fits images you would like to run VIDA on. Then you you can select the template you want using the `--template` command line option. For the other options type `-h`. If for example you wanted to fit the images with 1 third order mring and a gaussian and a disk you
would type

To correctly install the script please first run `setup.jl` to install the necessary packages. Then you can run the script with the following command
```bash
julia setup.jl
```

```bash
julia -p ncores image_extractor.jl list_of_files --template stretchmring_0_3 gauss_1 disk_1
```