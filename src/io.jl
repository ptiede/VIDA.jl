
"""
$(SIGNATURES)

where `fname` is a image file. This will call a method based on the file extension.
Curerntly we have .fits, and .h5 implemented.
# Details
This reads in a fits file that is more robust to the various imaging algorithms
in the EHT, i.e. is works with clean, smili, eht-imaging.

The function returns an `IntensityMap` object that contains the relevant image and parameters
extracted from the fits file. It also ensures that we are astronomers and that the image
using sky-left coordinates.
"""
function load_image(fname)
    _, ext = splitext(fname)

    if ext == ".fits"
        return ComradeBase.load(fname, IntensityMap)
    elseif ext == ".h5" || ext == ".hdf5"
        return load_im_h5(fname)
    else
        throw("$(ext) is not a valid file extension.")
    end
end


"""
$(SIGNATURES)

where `fname` should be a hdf5 image file generated using illinois hdf5 process
# Details

The function returns an IntensityMap object that contains the relevant image and parameters
extracted from the fits file. It also ensures that we are astronomers and that the image
using sky-left coordinates.
"""
function load_im_h5(fname::String)
    img = h5open(fname, "r") do fid
        header = fid["header"]
        dsource = read(header["dsource"])
        jyscale = read(header["scale"])
        rf = read(header["freqcgs"])
        tunit = read(header["units"]["T_unit"])
        lunit = read(header["units"]["L_unit"])
        dx = read(header["camera"]["dx"])
        nx = Int(read(header["camera"]["nx"]))
        time = read(header["t"])*tunit/3600
        image = collect(fid["pol"][1,:,:]')[end:-1:1,:]


        # Now convert everything to IntensityMap
        image = image.*jyscale
        src = "SgrA"
        ra = 17.7611
        dec = -28.992

        # convert to μas
        fov = μas2rad(dx/dsource*lunit*2.06265e11)

        mjd = 53005
        ComradeBase.MinimalHeader(src, ra, dec, mjd, rf)
        g = imagepixels(fov, fov, nx, nx)

        return IntensityMap(image, g)
    end
    return img
end


@doc """
    $(SIGNATURES)

Loads an hdf5 file where `filename` should be a HDF5 file.
# Details
This reads in a hdf5 file and outputs and EHTMovie object.

# Notes
Currently this only works with movies created by *ehtim*. SMILI uses a different
format, as does Illinois, and every other group...
"""
function load_hdf5(filename; style=:ehtim)
    if style == :ehtim
        return _load_ehtimhdf5(filename)
    else
        throw("hdf5 files not from ehtim aren't implemented")
    end
end

@doc """
    $(SIGNATURES)
Saves and hdf5 file where `filename` is the write out location.
Currently style only works with ehtim, namely we save HDF5 files
that only work with ehtim.
"""
function save_hdf5(filename, mov; style=:ehtim)
    # Copy this so I don't manipulate the movie itself
    I = copy(mov.frames)
    # Now because HDF5 uses row major I need to permute the dims around
    # so that ehtim reads this in correctly.
    I = ComradeBase.baseimage(I)[end:-1:1,end:-1:1,:]
    times = mov.frames.T
    head = header(mov.frames)

    if head isa ComradeBase.NoHeader
        head = (ra = 180.0, dec = 0.0, mjd = 5000, frequency = 230e9, source="M87")
    end

    #Open and output in a safe manner
    h5open(filename, "w") do file
        #Write the Intensity
        @write file I
        #Create the header dataset and write it
        file["header"] = ""
        header = file["header"]
        #Now I write the header as attributes since this is what ehtim does
        attributes(header)["dec"] = string(head.dec)
        attributes(header)["mjd"] = string(Int(head.mjd))
        attributes(header)["pol_prim"] = "I"
        attributes(header)["polrep"] = "stokes"
        attributes(header)["psize"] = string(step(mov.frames.X))
        attributes(header)["ra"] = string(head.ra)
        attributes(header)["rf"] = string(head.frequency)
        attributes(header)["source"] = string(head.source)
        #Now write the times
        @write file times
    end
    return nothing
end

function _load_ehtimhdf5(filename)
    #Open the hdf5 file
    img = h5open(filename, "r") do fid
        header = fid["header"]
        # Read the images TODO: Make this lazy to not nuke my memory
        images = read(fid["I"])
        images = images[end:-1:1,end:-1:1,:]
        psize = parse(Float64, read(header["psize"]))
        fov = psize*size(images,1)
        times = read(fid["times"])
        source = String(read(header["source"]))
        ra = parse(Float64, read(header["ra"]))
        dec = parse(Float64, read(header["dec"]))
        mjd = parse(Float64, read(header["mjd"]))
        rf = parse(Float64, read(header["rf"]))
        psize = parse(Float64, read(header["psize"]))*3600*1e6*180.0/π
        header = ComradeBase.MinimalHeader(source, ra, dec, mjd, rf)
        g = imagepixels(fov, fov, size(images, 1), size(images, 2); header)
        gt = RectiGrid((X=g.X, Y=g.Y, T = times))
        img = IntensityMap(images, gt)
        return img
    end
    return VIDAMovie(img)
end
