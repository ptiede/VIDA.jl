clipvalue(c,x) = x < c ? zero(eltype(x)) : x

@doc """
    $(SIGNATURES)
Clips the image `im` according to the value clip.
There are two modes for image clipping:
    - `:relative` which zeros the pixels whose intensity are below `clip` relative to the max.
    - `:absolute` which zeros the pixels whose intensity is below `clip` in Jy/pixel
"""
function clipimage(clip, image::SpatialIntensityMap, mode=:relative)
    cimg = copy(image)
    return clipimage!(clip, cimg, mode)
end

function clipimage!(clip, image::SpatialIntensityMap, mode=:relative)
    if mode == :absolute
        image .= clipvalue.(clip, image)
    elseif mode == :relative
        maxim = maximum(image)
        image .= clipvalue.(clip*maxim, image)
    else
        @assert false "clipimage: Mode must be one of :absolute or :relative where :absolute cuts on value and :relative on fraction of maximum intensity"
    end
    return image
end



@doc """
    regrid(img::EHTImage, nx, ny, xlim, ylim)
# Inputs
 - img::EHTImage : Image you want to regrid
 - npix : Number of pixels in x and y direction
 - xlim : Tuple with the limits of the image in the RA in μas
 - ylim : Tuple with the limits of the image in DEC in μas
"""
function regrid(img::SpatialIntensityMap, g::GriddedKeys)
    fimg = VLBISkyModels.InterpolatedImage(img)
    return intensitymap(fimg, g)
end

@doc """
    $(SIGNATURES)
Blurs the `img` with a gaussian kernel with fwhm in μas. If `fwhm` is a scalar
then the kernel is assumed to be symmetric, otherwise you
the first entry is the fwhm in the EW direction and second
the NS direction.

Returns the blurred image.
"""
function blur(img::SpatialIntensityMap, fwhm)
    σ = fwhm./(2*sqrt(2*log(2)))
    return convolve(img, modify(Gaussian(), stretch(σ)))
end
