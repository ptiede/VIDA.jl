"""
    $(TYPEDEF)
Creates an abstract movie class
"""
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
struct EHTMovie{F,T<:Interpolations.AbstractExtrapolation{F}} <: AbstractMovie
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
                  images::T) where {F,T<:AbstractArray{F,3}}
    #Create the interpolation object for the movie
    #This does not need equal times
    fimages = reshape(images, nx*ny, length(times))
    sitp = extrapolate(interpolate((collect(1.0:(nx*ny)), times),
                        fimages,
                        (NoInterp(), Gridded(Linear()))),
                        (Interpolations.Flat(), Interpolations.Flat()))
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
    return mov.frames.itp.knots[2]
end


@doc """
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

@doc """
    $(SIGNATURES)
Gets all the frames of the movie object `mov`. This returns a array of `EHTImage`
objects.
"""
function get_frames(mov::EHTMovie)
    return get_image.(Ref(mov), get_times(mov))
end

@doc """
    $(SIGNATURES)
Returns the flux of the `mov` at the times `time` in fractional hours
"""
function flux(mov, t)
    img = get_image(mov, t)
    return flux(img)
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
    return join_frames(times, bframes)
end

@doc """
    rescale(mov::EHTMovie, npix, xlim, ylim)
# Inputs
 - mov::EHTMovie : Movie you want to rescale
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
