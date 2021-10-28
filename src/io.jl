
"""
$(SIGNATURES)

where `fname` is a image file. This will call a method based on the file extension.
Curerntly we have .fits, and .h5 implemented.
# Details
This reads in a fits file that is more robust to the various imaging algorithms
in the EHT, i.e. is works with clean, smili, eht-imaging.

The function returns an EHTImage object that contains the relevant image and parameters
extracted from the fits file. It also ensures that we are astronomers and that the image
using sky-left coordinates.
"""

function load_image(fname)
    _, ext = splitext(fname)

    if ext == ".fits"
        return load_fits(fname)
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

The function returns an EHTImage object that contains the relevant image and parameters
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
        image = fid["pol"][1,:,:]

        # Now convert everything to EHTImage
        image = image.*jyscale
        src = "SgrA"
        ra = 17.7611
        dec = -28.992

        # convert to μas
        fovμas = dx/dsource*lunit*2.06265e11
        psize_x = fovμas/nx
        psize_y = fovμas/nx

        mjd = 53005.0
        return EHTImage(nx, nx,
                        -psize_x, psize_y,
                        src,
                        ra, dec,
                        C0/rf,
                        mjd,
                        image
                    )
    end
    return img
end


"""
$(SIGNATURES)

where `fits_name` should be a fits file generated using ehtim
# Details
This reads in a fits file that is more robust to the various imaging algorithms
in the EHT, i.e. is works with clean, smili, eht-imaging.

The function returns an EHTImage object that contains the relevant image and parameters
extracted from the fits file. It also ensures that we are astronomers and that the image
using sky-left coordinates.
"""
function load_fits(fits_name::String)
    #load the fits
    f = FITS(fits_name)
    #@assert ndims(f[1]) == 2 "load_image: First element is expected to be an ImageHDU so ndims is expected to be 2"

    #Check if there are stokes parameters (don't load them right now)
    if ndims(f[1]) == 2
        image = deepcopy(Matrix{Float64}(read(f[1])'))
    elseif ndims(f[1]) == 4
        @warn "Only stokes I will be loaded. Polarization not implemented yet."
        image = Matrix{Float64}(read(f[1])[:,:,1,1]')
    end


    header = read_header(f[1])
    #Read image dimensions
    nx = Int(header["NAXIS1"])
    ny = Int(header["NAXIS2"])
    psize_x = -abs(float(header["CDELT1"])*3600*1e6)
    psize_y = abs(float(header["CDELT2"]))*3600*1e6


    source = string(header["OBJECT"])
    ra = float(header["OBSRA"])
    dec = float(header["OBSDEC"])
    #Get frequency
    freq = 0.0
    if haskey(header, "FREQ")
        freq = parse(Float64, string(header["FREQ"]))
    elseif "CRVAL3" in keys(header)
        freq = float(header["CRVAL3"])
    end

    mjd = 0.0
    if haskey(header, "MJD")
        mjd = parse(Float64, string(header["MJD"]))
    end

    source = "NA"
    if haskey(header,"OBJECT")
        source = string(header["OBJECT"])
    end

    #Now renormalize the images if not using Jy/pixel
    bmaj = 1.0 #Nominal values
    bmin = 1.0
    if haskey(header, "BUNIT")
        if header["BUNIT"] == "JY/BEAM"
            println("Converting Jy/Beam => Jy/pixel")
            try
                bmaj = header["BMAJ"]
                bmin = header["BMIN"]
            catch
                @warn "No beam found in header using nominal values"
            end
            beamarea = (2.0*π*bmaj*bmin)/(8*log(2))
            image .= image.*(header["CDELT2"]^2/beamarea)
        end
    end

    close(f)

    return EHTImage(nx, ny,
                    psize_x, psize_y,
                    source,
                    ra, dec,
                    C0/freq,
                    mjd,
                    image)
end

@doc  """
    save_fits(image::EHTImage, fname::String)
Save the `image` as a fits object with filename `fname`
"""
function save_fits(image::EHTImage, fname::String)
    headerkeys = ["SIMPLE",
                  "BITPIX",
                  "NAXIS",
                  "NAXIS1",
                  "NAXIS2",
                  "OBJECT",
                  "CTYPE1",
                  "CTYPE2",
                  "CDELT1",
                  "CDELT2",
                  "OBSRA",
                  "OBSDEC",
                  "FREQ",
                  "CRPIX1",
                  "CRPIX2",
                  "MJD",
                  "TELESCOP",
                  "BUNIT",
                  "STOKES"]
    values = [true,
              -64,
              2,
              image.ny,
              image.nx,
              image.source,
              "RA---SIN",
              "DEC---SIN",
              image.psize_y/3600/1e6,
              image.psize_x/3600/1e6,
              image.ra,
              image.dec,
              C0/image.wavelength,
              image.nx/2+0.5,
              image.ny/2+0.5,
              image.mjd,
              "VLBI",
              "JY/PIXEL",
              "STOKES"]
    comments = ["conforms to FITS standard",
                "array data type",
                "number of array dimensions",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                ""]
    hdu = FITS(fname, "w")
    hduheader = FITSHeader(headerkeys, values, comments)
    img = copy(image.img')
    write(hdu, img, header=hduheader)
    close(hdu)
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
    I = deepcopy(reshape(mov.frames.itp.coefs, mov.nx, mov.ny, length(mov.frames.itp.knots[2])))
    # Now because HDF5 uses row major I need to permute the dims around
    # so that ehtim reads this in correctly.
    I = permutedims(I, [2,1,3])[:,end:-1:1,:]
    times = mov.frames.itp.knots[2]

    #Open and output in a safe manner
    h5open(filename, "w") do file
        #Write the Intensity
        @write file I
        #Create the header dataset and write it
        file["header"] = ""
        header = file["header"]
        #Now I write the header as attributes since this is what ehtim does
        attributes(header)["dec"] = string(mov.dec)
        attributes(header)["mjd"] = string(Int(mov.mjd))
        attributes(header)["pol_prim"] = "I"
        attributes(header)["polrep"] = "stokes"
        attributes(header)["psize"] = string(mov.psize_y/(3600.0*1e6*180.0/π))
        attributes(header)["ra"] = string(mov.ra)
        attributes(header)["rf"] = string(C0/mov.wavelength)
        attributes(header)["source"] = string(mov.source)
        #Now write the times
        @write file times
    end
    return nothing
end

function _load_ehtimhdf5(filename)
    #Open the hdf5 file
    fid = h5open(filename, "r")
    try
        header = fid["header"]
        # Read the images TODO: Make this lazy to not nuke my memory
        images = read(fid["I"])
        images = permutedims(images, [2,1,3])[end:-1:1,:,:]
        npix = Base.size(images)[1]
        times = read(fid["times"])
        source = String(read(header["source"]))
        ra = parse(Float64, read(header["ra"]))
        dec = parse(Float64, read(header["dec"]))
        mjd = parse(Float64, read(header["mjd"]))
        rf = parse(Float64, read(header["rf"]))
        psize = parse(Float64, read(header["psize"]))*3600*1e6*180.0/π
        close(fid)
        return EHTMovie(npix, npix,
                        -psize, psize,
                        source,
                        ra, dec,
                        C0/rf,
                        mjd,
                        times, images
                    )
    catch
        close(fid)
        println("Unable to open $filename")
        return nothing
    end
end
