






clipvalue(c,x) = x < c ? zero(eltype(x)) : x

@doc """
    $(SIGNATURES)
Clips the image `im` according to the value clip.
There are two modes for image clipping:
    - `:relative` which zeros the pixels whose intensity are below `clip` relative to the max.
    - `:absolute` which zeros the pixels whose intensity is below `clip` in Jy/pixel
"""
function clipimage(clip, image::IntensityMap, mode=:relative)
    cimg = copy(image)
    return clipimage!(clip, cimg, mode)
end

function clipimage!(clip, image::IntensityMap, mode=:relative)
    if mode == :absolute
        map!(x->clipvalue(clip,x),image, image)
    elseif mode == :relative
        maxim = maximum(image)
        map!(x->clipvalue(clip*maxim,x),image, image)
    else
        @assert false "clipimage: Mode must be one of :absolute or :relative where :absolute cuts on value and :relative on fraction of maximum intensity"
    end
    return image
end



@doc """
    rescale(img::EHTImage, nx, ny, xlim, ylim)
# Inputs
 - img::EHTImage : Image you want to rescale
 - npix : Number of pixels in x and y direction
 - xlim : Tuple with the limits of the image in the RA in μas
 - ylim : Tuple with the limits of the image in DEC in μas
"""
function rescale(img::IntensityMap, nx, ny, xlim, ylim)
    x_itr,y_itr = imagepixels(img)
    itp = interpolate(img.im/(-step(x_itr)*step(y_itr)), BSpline(Cubic(Line(OnGrid()))))
	sitp = scale(itp, reverse(x_itr), y_itr)
    etp = extrapolate(sitp, 0)

    #Create grid for new image
    fovy_new = (ylim[2]-ylim[1])
    psize_y = fovy_new/(npix)
    fovx_new = (xlim[2]-xlim[1])
    psize_x = fovx_new/(npix)
    x_itr_new = (fovx_new/2 - psize_x/2):-psize_x:(-fovx_new/2 + psize_x/2)
    y_itr_new = (-fovy_new/2 + psize_y/2):psize_y:(fovy_new/2 - psize_y/2)
    #Create new image
    img_new = etp(reverse(x_itr_new), y_itr_new)*psize_x*psize_y
    return EHTImage(npix, npix, -psize_x, psize_y, img.source, img.ra, img.dec,
                    img.wavelength, img.mjd, img_new)
end

@doc """
    $(SIGNATURES)
Blurs the `img` with a gaussian kernel with fwhm in μas. If `fwhm` is a scalar
then the kernel is assumed to be symmetric, otherwise you
the first entry is the fwhm in the EW direction and second
the NS direction.

Returns the blurred image.
"""
function blur(img::EHTImage, fwhm)
    # Using image template function which blurs on the pixel scale.
    # So first transform from μas to pixel number and standard deviation
    σ_px = fwhm./(2*sqrt(2*log(2)))./abs(img.psize_x)
    σ_py = fwhm./(2*sqrt(2*log(2)))./abs(img.psize_y)

    # Now I need to pick my kernel size. I am going out to 5σ for the
    # gaussian kernel. I have to add one for the convolution to play nice
    nkern = Int(floor(σ_px)*10 + 1)
    # Note I tried to use IIRGaussian but it wasn't accurate enough for us.
    bimg = imfilter(img,
                    gaussian((σ_py, σ_px),(nkern,nkern)),
                    Fill(0.0, img),
                    FFT()
                )

    #Now return the updated image
    #TODO can I find a way to do this without having to allocate a new image?
    return bimg
end
