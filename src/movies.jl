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
"""
struct VIDAMovie{F<:AbstractVector{<:IntensityMap}, I<:Interpolations.AbstractExtrapolation} <: AbstractMovie
    frames::F
    itp::I
end


"""
    VIDAMovie(mov::IntensityMap{T,3})

Create an VIDAMovie class for easy interpolation between frames.

# Arguments

 - `mov`: An IntensityMap with axes (:X, :Y, :T) that represent the frames of a movie.
    Note that the time dimension does not have to be equi-spaced.

# Returns
An `VIDAMovie` object that behaves like `IntensityMap` but lets you interpolate between frames
with [`get_image(vida_movie, time)`](@ref).
"""
function VIDAMovie(
    mov::IntensityMap{T, 3}
    ) where {F,T<:SpatialIntensityMap{F,3}}
    #Create the interpolation object for the movie
    #This does not need equal times
    @assert propertynames(mov) == (:X, :Y, :T) "Array must have dimension X, Y, T"
    nx = size(img, :X)
    ny = size(img, :Y)
    nt = size(img, :T)
    fimages = reshape(mov, nx*ny, nt)
    sitp = extrapolate(interpolate((collect(1.0:(nx*ny)), mov.T),
                        ComradeBase.baseimage(fimages),
                        (NoInterp(), Gridded(Linear()))),
                        (Interpolations.Flat(), Interpolations.Flat()))
    return VIDAMovie(mov, sitp)
end

VIDAMovie(times, images::Vector{<:SpatialIntensityMap}) = VIDAMovie(_join_frames(times, images))

function _join_frames(times, images)
    g = axiskeys(first(images))
    mimg = stack(images)
    gt = GriddedKeys((X=g.X, Y=g.Y, T=times))
    return IntensityMap(mimg, gt)
end


@doc """
    $(SIGNATURES)
Joins an array of EHTImages at specified times to form an VIDAMovie object.

## Inputs
 - times: An array of times that the image was created at
 - images: An array of EHTImage objects

## Outputs
VIDAMovie object
"""
function join_frames(times, images::Vector{T}) where {T<:SpatialIntensityMap}
    return VIDAMovie(_join_frames(times, images))
end

@doc """
    $(SIGNATURES)
Returns the times that the movie object `mov` was created at. This does not
have to be uniform in time.
"""
function get_times(mov::VIDAMovie)
    return mov.frames.T
end


@doc """
    $(SIGNATURES)
Gets the frame of the movie object `mov` at the time t. This returns an `EHTImage`
object at the requested time. The returned object is found by linear interpolation.
"""
function get_image(mov::VIDAMovie, t)
    img = reshape(mov.frames.(1:(mov.nx*mov.ny), Ref(t)), mov.ny, mov.nx)
    return IntensityMap(img, GriddedKeys((X=mov.frames.X, Y=mov.frame.Y, T=t:t)))
end

@doc """
    $(SIGNATURES)
Gets all the frames of the movie object `mov`. This returns a array of `EHTImage`
objects.
"""
function get_frames(mov::VIDAMovie)
    return mov.frames
end

@doc """
    $(SIGNATURES)
Returns the flux of the `mov` at the times `time` in fractional hours
"""
function CB.flux(mov, t)
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
function blur(mov::VIDAMovie, fwhm)
    frames = get_frames(mov)
    bframes = map(x->convolve(x, modify(Gaussian(), Stretch(fwhm))), eachslice(frames, dims=:T))
    return bframes
end

@doc """
    regrid(mov::VIDAMovie, npix, xlim, ylim)
# Inputs
 - mov::VIDAMovie : Movie you want to regrid
 - npix : Number of pixels in x and y direction
 - xlim : Tuple with the limits of the image in the RA
 - ylim : Tuple with the limits of the image in DEC
"""
function regrid(mov::VIDAMovie, g::GriddedKeys{(:X, :Y)})
    # Get the times and frames and apply the image method to each
    frames = get_frames(mov)
    # Isn't broadcasting the best?
    rframes = map(eachslice(frames; dims=(:X, :Y))) do I
        fimg = VLBISkyModels.InterpolatedImage(I)
        img = intensitymap(fimg, g)
        return img
    end
    return rframes
end
