"""
An abstact type that acts as a wrapper for image objects used in astronomy.

This is the top of the castle for images and will be rarely touched. Basically
unless you don't want to use fits images this will not be used
"""
abstract type AbstractImage end

"""
An absract image that will hold a fits image after being created or parsed in.
This will form the basis for most astronomical images that are defined.
"""
abstract type AbstractFitsImage{T<:AbstractArray} <: AbstractImage end

"""
$(TYPEDEF)

# Details
The trait is to hold a EHT image tyically in matrix form. Namely the trait will
typically be Matrix{Float64}.

`nx` is the number of pixels in the x or RA direction
`ny` is the number of pixels in the y or DEC direction
`psize_x`, `psize_y` are the pixel sizes in the x and y direction
`source` is the source we are looking at e.g. M87
`ra`,`dec` are the sources RA and DEC in J2000 coordinates using degrees
`wavelength` is the wavelength of the image.
`mjd` is the Modified Julian Date of the observation.
`img` is the actual pixeled image in Jy/pixel
"""
struct EHTImage{T} <: AbstractFitsImage{T}
    nx::Int #Number of pixel in x direction
    ny::Int #Number of pixels in y direction
    psize_x::Float64 #pixel size in μas
    psize_y::Float64 #pixel size in μas
    source::String
    #Source RA and Dec in J2000 in degrees
    ra::Float64
    dec::Float64
    wavelength::Float64 #wavelength of image in cm
    mjd::Float64 #modified julian date of observation
    #Image matrix stored in Jy/px
    img::T
end

@doc """
    get_radec(img::T) <: AbstractFitsImage
Returns two iterators (ra,dec) that give the locations
of the `img` pixels.
"""
function get_radec(img::T) where {T<: AbstractFitsImage}
    ra = range((-img.nx*img.psize_x + img.psize_x)/2.0,
               step=img.psize_x,
               length=img.nx
    )
    dec = range((-img.ny*img.psize_y + img.psize_y)/2.0,
                step=img.psize_y,
                length=img.ny
    )

    return ra,dec
end


clipvalue(c,x) = x < c ? zero(eltype(x)) : x
"""
    $(SIGNATURES)
Clips the image `im` according to the value clip, which
can either be an absolute flux in Jy/px or the intensity
relative to the maximum.
"""
function clipimage(clip, im::EHTImage, mode=:relative)
    image = im.img
    if mode == :absolute
        image = map(x->clipvalue(clip,x), image)
    elseif mode == :relative
        maxim = maximum(image)
        image = map(x->clipvalue(clip*maxim,x), image)
    else
        @assert false "clipimage: Mode must be one of :absolute or :relative where :absolute cuts on value and :relative on fraction of maximum intensity"
    end
    return EHTImage(im.nx, im.ny,
                    im.psize_x, im.psize_y,
                    im.source, im.ra, im.dec,
                    im.wavelength, im.mjd,
                    image)
end

"""
    $(SIGNATURES)
Down samples the image `im` by a `factor`

# Details
Given an image `im` with ``Nx\\times Ny`` pixels, `downsample`
converts that to an image with ``Nx/factor \\times Ny/factor`` pixels.
This can be useful when the image has much higher resolutions than is needed for feature extraction.
"""
function downsample(factor::Int, im::EHTImage)
    im_arr = im.img
    dim_arr = im_arr[1:factor:end, 1:factor:end]
    ny,nx = Base.size(dim_arr)
    psize_x = im.psize_x*factor
    psize_y = im.psize_y*factor
    return EHTImage(nx,ny,
                    psize_x,psize_y,
                    im.source, im.ra, im.dec,
                    im.wavelength, im.mjd,
                    dim_arr*factor^2)
end

"""
    $(SIGNATURES)
Given an image `im` it selects a new field of view for the image given by
`domain_x`, `domain_y`.
"""
function window_image(domain_x, domain_y, im::EHTImage)
    fovx = -im.psize_x*im.nx
    fovy = im.psize_y*im.ny
    xarr = collect((fovx+im.psize_x)/2:im.psize_x:-fovx/2)
    yarr = collect((-fovy+im.psize_y)/2:im.psize_y:fovy/2)

    imin = argmin(abs.(xarr .- domain_x[1]))
    imax = argmin(abs.(xarr .- domain_x[2]))
    jmin = argmin(abs.(yarr .- domain_y[1]))
    jmax = argmin(abs.(yarr .- domain_y[2]))
    window_image = im.img[jmin:jmax,imax:imin]
    ny,nx = Base.size(window_image)
    return EHTImage(nx,ny,
                    im.psize_x, im.psize_y,
                    im.source, im.ra, im.dec,
                    im.wavelength, im.mjd,
                    window_image)
end




"""
$(SIGNATURES)

where `fits_name` should be a fits file generated using ehtim
# Details
This reads in a fits file created using ehtim. This is because ehtim only outputs
the image and not a separate HDU for the field, so the usual fits reader doesn't work
properly.

The function returns an EHTImage object that contains the relevant image and parameters
extracted from the fits file. It also ensures that we are astronomers and that the image
using sky-left coordinates.
"""
function load_ehtimfits(fits_name::String)
    #load the fits
    f = FITS(fits_name)
    @assert ndims(f[1]) == 2 "load_image: First element is expected to be an ImageHDU so ndims is expected to be 2"
    header = read_header(f[1])
    nx = Int(header["NAXIS1"])
    ny = Int(header["NAXIS2"])
    image = deepcopy(Matrix{Float64}(read(f[1])'))
    source = string(header["OBJECT"])
    ra = float(header["OBSRA"])
    dec = float(header["OBSDEC"])
    freq = float(header["FREQ"])
    mjd = float(header["MJD"])
    psize_x = -abs(float(header["CDELT1"])*3600*1e6)
    psize_y = abs(float(header["CDELT2"]))*3600*1e6
    close(f)

    return EHTImage(nx, ny, psize_x, psize_y, source, ra, dec, C0/freq, mjd, image)
end

"""
    $(SIGNATURES)
Finds the centroid or center of light of the `img` in μas.
"""
function centroid(img::EHTImage)
    dx = abs(img.psize_x)
    dy = abs(img.psize_y)
    fovx = dx*img.nx
    fovy = dx*img.ny
    xstart = (fovx-dx)/2.0
    ystart = (-fovy+dy)/2.0
    inorm = zero(eltype(img.img))
    xcent = zero(eltype(img.img))
    ycent = zero(eltype(img.img))
    @inbounds for i in 1:img.nx
        @inbounds for j in 1:img.ny
            x = xstart - dx*(i-1)
            y = ystart + dy*(j-1)
            xcent += x*img.img[j,i]
            ycent += y*img.img[j,i]
            inorm += img.img[j,i]
        end
    end
    return xcent/inorm,ycent/inorm
end

"""
    $(SIGNATURES)
Find the image moment of inertia or **second moment**

### Notes
If `center=true` then we find the central second moment, or the second
cumulant of the image.
"""
function inertia(img::EHTImage, center=false)
    moments = Matrix{eltype(img.img)}(undef,2,2)
    dx = abs(img.psize_x)
    dy = abs(img.psize_y)
    fovx = dx*img.nx
    fovy = dx*img.ny
    xstart = (fovx-dx)/2.0
    ystart = (-fovy+dy)/2.0
    inorm = zero(eltype(img.img))
    xx = zero(eltype(img.img))
    xy = zero(eltype(img.img))
    yy = zero(eltype(img.img))
    @inbounds for i in 1:img.nx
        @inbounds for j in 1:img.ny
            x = xstart - dx*(i-1)
            y = ystart + dy*(j-1)
            xx += x*x*img.img[j,i]
            yy += y*y*img.img[j,i]
            xy += x*y*img.img[j,i]
            inorm += img.img[j,i]
        end
    end
    if center
        xcent,ycent = centroid(img)
        xx = xx - xcent*xcent
        yy = yy - ycent*ycent
        xy = xy - xcent*ycent
    end
    moments[1,1] = xx/inorm
    moments[2,2] = yy/inorm
    moments[1,2] = moments[2,1] = xy/inorm
    return moments
end

"""
    $(SIGNATURES)
Finds the image flux of an EHTImage `img`
"""
function flux(img::EHTImage)
    return sum(img.img)
end

@doc """
    rescale_image(img::EHTImage, npix, xlim, ylim)
# Inputs
 - img::EHTImage : Image you want to rescale
 - npix : Number of pixels in x and y direction
 - xlim : Tuple with the limits of the image in the RA
 - ylim : Tuple with the limits of the image in DEC
"""
function rescale_image(img::EHTImage, npix, xlim, ylim)
    fovx, fovy = field_of_view(img)
    x_itr = (fovx/2 - img.psize_x/2):-img.psize_x:(-fovx/2 + img.psize_x/2)
	y_itr = (-fovy/2 + img.psize_y/2):img.psize_y:(fovy/2 - img.psize_y/2)
	itp = interpolate(img.img, BSpline(Cubic(Line(OnGrid()))))
	etp = extrapolate(itp, 0)
	sitp = scale(etp, x_itr, y_itr)

    #Create grid for new image
    fovy_new = (ylim[2]-ylim[1])
    psize_y = fovy_new/npix
    fovx_new = (xlim[2]-xlim[1])
    psize_x = fovx_new/npix
    x_itr_new = (fovx_new/2 - psize_x/2):-psize_x:(-fovx_new/2 + psize_x/2)
    y_itr_new = (-fovy_new/2 + psize_y/2):psize_y:(fovy_new/2 - psize_y/2)
    #Create new image
    img_new = sitp(reverse(x_itr_new), y_itr_new)
    return EHTImage(npix, npix, -psize_x, psize_y, img.source, img.ra, img.dec,
                    img.wavelength, img.mjd, img_new)
end

@doc """
    field_of_view(img::EHTImage)
Finds the field of view of an EHTImage. Return a w element tuple with the
field of view in the x and y direction
"""
function field_of_view(img::EHTImage)
    return img.nx*img.psize_x, img.ny*img.psize_y
end


function save_ehtimfits(image::EHTImage, fname::String)
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
