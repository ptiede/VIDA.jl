abstract type AbstractMovie end

"""
$(TYPEDEF)

# Details
The type is to hold a EHT movie. The dimension of the movie array
is assumed to be in the form DEC,RA,Time.

`nx` is the number of pixels in the x or RA direction
`ny` is the number of pixels in the y or DEC direction
`psize_x`, `psize_y` are the pixel sizes in the x and y direction
`source` is the source we are looking at e.g. M87
`ra`,`dec` are the sources RA and DEC in J2000 coordinates using degrees
`wavelength` is the wavelength of the image.
`mjd` is the Modified Julian Date of the observation.
`frames` is the interpolation object that hold the movie frames
"""
struct EHTMovie{T<:Interpolations.GriddedInterpolation} <: AbstractMovie
    nx::Int #Number of pixel in x direction
    ny::Int #Number of pixels in y direction
    psize_x::Float64 #pixel size in μas
    psize_y::Float64 #pixel size in μas
    source::String
    ra::Float64
    dec::Float64
    wavelength::Float64 #wavelength of image in cm
    mjd::Float64 #modified julian date of observation
    frames::T
end

function EHTMovie(nx,
                  ny,
                  psize_x,
                  psize_y,
                  source,
                  ra,
                  dec,
                  wavelength,
                  mjd,
                  times,
                  images::T) where {T<:AbstractArray{Float64,3}}
    #Create the interpolation object for the movie
    #This does not need equal times
    fimages = reshape(images, nx*ny, length(times))
    sitp = interpolate((collect(1.0:(nx*ny)), times),
                        fimages,
                        (NoInterp(), Gridded(Linear()))
                      )
    return EHTMovie(nx, ny,
                    psize_x, psize_y,
                    source,
                    ra, dec,
                    wavelength,
                    mjd,
                    sitp)
end


@doc """
    $(SIGNATURES)
Joins an array of EHTImages at specified times to form an EHTMovie object.

## Inputs
 - times: An array of times that the image was created at
 - images: An array of EHTImage objects

## Outputs
EHTMovie object
"""
function join_frames(times, images::Vector{T}) where {T<:EHTImage}
    nx,ny = images[1].nx, images[1].ny
    nt = length(times)

    #Allocate the image array and fill it
    imarr = zeros(ny,nx,nt)
    for i in 1:nt
        imarr[:,:,i] = images[i].img
    end

    return EHTMovie(nx,ny,
                    images[1].psize_x, images[1].psize_y,
                    images[1].source,
                    images[1].ra, images[1].dec,
                    images[1].wavelength,
                    images[1].mjd,
                    times, imarr)

end

@doc """
    $(SIGNATURES)
Returns the times that the movie object `mov` was created at. This does not
have to be uniform in time.
"""
function get_times(mov::EHTMovie)
    return mov.frames.knots[2]
end


"""
    $(SIGNATURES)
Gets the frame of the movie object `mov` at the time t. This returns an `EHTImage`
object at the requested time. The returned object is found by linear interpolation.
"""
function get_image(mov::EHTMovie, t)
    img = reshape(mov.frames.(1:(mov.nx*mov.ny), Ref(t)), mov.ny, mov.nx)
    return EHTImage(mov.nx,
                    mov.ny,
                    mov.psize_x,
                    mov.psize_y,
                    mov.source,
                    mov.ra, mov.dec,
                    mov.wavelength,
                    mov.mjd,
                    img
                    )
end

"""
    $(SIGNATURES)
Gets all the frames of the movie object `mov`. This returns a array of `EHTImage`
objects.
"""
function get_frames(mov::EHTMovie)
    images = sizehint!(EHTImage[], length(get_times(mov)))
    for i in 1:length(get_times(mov))
        imvec = @view mov.frames.coefs[:,i]
        img2 = reshape(imvec, mov.ny, mov.nx)
        ehtimg = EHTImage(mov.nx,
                        mov.ny,
                        mov.psize_x,
                        mov.psize_y,
                        mov.source,
                        mov.ra, mov.dec,
                        mov.wavelength,
                        mov.mjd,
                        img2
                      )
        push!(images, ehtimg)
    end
    return images
end


@doc """
    $(SIGNATURES)
Blurs the `mov` with a gaussian kernel with fwhm in μas. If `fwhm` is a scalar
then the kernel is assumed to be symmetric, otherwise you
the first entry is the fwhm in the EW direction and second
the NS direction.

Returns the blurred movie.
"""
function blur(mov::EHTMovie, fwhm)
    times = get_times(mov)
    frames = get_frames(mov)
    bframes = blur.(frames, Ref(fwhm))
    return join_frames(times, frames)
end

@doc """
    rescale(img::EHTMovie, npix, xlim, ylim)
# Inputs
 - img::EHTImage : Image you want to rescale
 - npix : Number of pixels in x and y direction
 - xlim : Tuple with the limits of the image in the RA
 - ylim : Tuple with the limits of the image in DEC
"""
function rescale(mov::EHTMovie, npix, xlim, ylim)
    # Get the times and frames and apply the image method to each
    times = get_times(mov)
    frames = get_frames(mov)
    # Isn't broadcasting the best?
    rframes = rescale.(frames, Ref(npix), Ref(xlim), Ref(ylim))
    return join_frames(times, rframes)
end


"""
$(SIGNATURES)

where `filename` should be a HDF5 file.
# Details
This reads in a hdf5 file and outputs and EHTMovie object.

# Notes
Currently this only works with movies created by *ehtim*. SMILI uses a different
format, as does Illinois, and every other group.
"""
function load_hdf5(filename; style=:ehtim)
    if style == :ehtim
        return _load_ehtimhdf5(filename)
    else
        throw("hdf5 files not from ehtim aren't implemented")
    end
end

function save_hdf5(filename, mov; style=:ehtim)
    # Copy this so I don't manipulate the movie itself
    I = deepcopy(reshape(mov.frames.coefs, mov.nx, mov.ny, length(mov.frames.knots[2])))
    # Now because HDF5 uses row major I need to permute the dims around
    # so that ehtim reads this in correctly.
    I = permutedims(I, [2,1,3])[:,end:-1:1,:]
    times = mov.frames.knots[2]

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
