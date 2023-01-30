"""
An abstact type that acts as a wrapper for image objects used in astronomy.

This is the top of the castle for images and will be rarely touched. Basically
unless you don't want to use fits images this will not be used
"""
abstract type AbstractImage{F<:Number} <: AbstractMatrix{F} end


"""
An absract image that will hold a fits image after being created or parsed in.
This will form the basis for most astronomical images that are defined.

The defines an interface to the FITS image methods below. In the future
I may change this to use the traits interfaces.
"""
abstract type AbstractFitsImage{F<:Number,T<:AbstractArray{F,2}} <: AbstractImage{F} end
Base.IndexStyle(::Type{<:AbstractFitsImage{T,S}}) where {T,S} = Base.IndexStyle(S)
Base.size(im::AbstractFitsImage) = size(im.img)
Base.getindex(im::AbstractFitsImage, i::Int) = getindex(im.img, i)
Base.setindex!(im::AbstractFitsImage, x, i::Int) = setindex!(im.img, x, i)



@doc """
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
struct EHTImage{F,T<:AbstractArray{F,2}} <: AbstractFitsImage{F,T}
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

function Base.similar(img::EHTImage, ::Type{T}) where{T}
    sim = similar(img.img, T)
    return EHTImage(img.nx, img.ny,
                    img.psize_x, img.psize_y,
                    img.source,
                    img.ra, img.dec,
                    img.wavelength, img.mjd,
                    sim)
end

function Base.similar(img::EHTImage, ::Type{T}, dims::Dims) where {T}
    fovx = img.psize_x*last(dims)
    fovy = img.psize_y*first(dims)
    sim = similar(img.img, T, dims)
    return EHTImage(dims[2],dims[1],
                    img.psize_x, img.psize_y,
                    img.source,
                    img.ra, img.dec,
                    img.wavelength, img.mjd,
                    sim)
end

struct FitsImageStyle <: Broadcast.AbstractArrayStyle{2} end
FitsImageStyle(::Val{2}) = FitsImageStyle()

Base.BroadcastStyle(::Type{<:AbstractFitsImage}) = FitsImageStyle()
function Base.similar(bc::Broadcast.Broadcasted{FitsImageStyle}, ::Type{ElType}) where ElType
    #Scan inputs for StokesImage
    #print(bc.args)
    img = _find_sim(bc)
    #print(Im)
    #fovxs = getproperty.(Ims, Ref(:fovx))
    #fovys = getproperty.(Ims, Ref(:fovy))
    #@assert all(i->i==first(fovxs), fovxs) "StokesImage fov must be equal to add"
    #@assert all(i->i==first(fovys), fovys) "StokesImage fov must be equal to add"
    return EHTImage(img.nx, img.ny,
                    img.psize_x, img.psize_y,
                    img.source,
                    img.ra, img.dec,
                    img.wavelength, img.mjd,
                    similar(Array{ElType}, axes(bc)))
end

#Finds the first StokesImage and uses that as the base
#TODO: If multiple StokesImages maybe I should make sure they are consistent?
_find_sim(bc::Base.Broadcast.Broadcasted) = _find_sim(bc.args)
_find_sim(args::Tuple) = _find_sim(_find_sim(args[1]), Base.tail(args))
_find_sim(x) = x
_find_sim(::Tuple{}) = nothing
_find_sim(a::EHTImage, rest) = a
_find_sim(::Any, rest) = _find_sim(rest)

#Guards to prevent someone from adding two Images with different FOV's
function Base.:+(x::AbstractFitsImage, y::AbstractFitsImage)
    @assert field_of_view(x) == field_of_view(y) "EHTImage must share same field of view"
    return x .+ y
end



@doc """
    $(SIGNATURES)
Creates an image interpolation functor from an `img` and
Interpolations.jl `interp` type.

## Notes
This pads the image with zeros for extrapolation.
"""
function image_interpolate(img::EHTImage, interp)
    fovx, fovy = field_of_view(img)
    x_itr = (fovx/2 - img.psize_x/2):-img.psize_x:(-fovx/2 + img.psize_x/2)
	y_itr = (-fovy/2 + img.psize_y/2):img.psize_y:(fovy/2 - img.psize_y/2)
    itp = interpolate(img', interp)
    etp = extrapolate(itp, 0)
    sitp = scale(etp, reverse(x_itr), y_itr)
    return sitp
end

@doc """
    $(SIGNATURES)
Returns two iterators (ra,dec) that give the locations
of the `img` pixels.
"""
function pixelloc(img::T) where {T<: AbstractFitsImage}
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

@doc """
    $(SIGNATURES)
Clips the image `im` according to the value clip.
There are two modes for image clipping:
    - `:relative` which zeros the pixels whose intensity are below `clip` relative to the max.
    - `:absolute` which zeros the pixels whose intensity is below `clip` in Jy/pixel
"""
function clipimage(clip, image::EHTImage, mode=:relative)
    cimg = copy(image)
    return clipimage!(clip, cimg, mode)
end

function clipimage!(clip, image::EHTImage, mode=:relative)
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
            xcent += x*img[j,i]
            ycent += y*img[j,i]
            inorm += img[j,i]
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
            xx += x*x*img[j,i]
            yy += y*y*img[j,i]
            xy += x*y*img[j,i]
            inorm += img[j,i]
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

@doc """
    $(SIGNATURES)
Finds the image flux of an EHTImage `img`
"""
function flux(img::EHTImage)
    return sum(img)
end

@doc """
    $(SIGNATURES)
Finds the field of view of an EHTImage. Return a w element tuple with the
field of view in the x and y direction
"""
function field_of_view(img::EHTImage)
    return img.nx*img.psize_x, img.ny*img.psize_y
end


@doc """
    rescale(img::EHTImage, npix, xlim, ylim)
# Inputs
 - img::EHTImage : Image you want to rescale
 - npix : Number of pixels in x and y direction
 - xlim : Tuple with the limits of the image in the RA in μas
 - ylim : Tuple with the limits of the image in DEC in μas
"""
function rescale(img::EHTImage, npix, xlim, ylim)
    x_itr,y_itr = pixelloc(img)
    itp = interpolate(img.img/(-step(x_itr)*step(y_itr)), BSpline(Cubic(Line(OnGrid()))))
	sitp = scale(itp, y_itr, reverse(x_itr))
    etp = extrapolate(sitp, 0)

    #Create grid for new image
    fovy_new = (ylim[2]-ylim[1])
    psize_y = fovy_new/(npix)
    fovx_new = (xlim[2]-xlim[1])
    psize_x = fovx_new/(npix)
    x_itr_new = (fovx_new/2 - psize_x/2):-psize_x:(-fovx_new/2 + psize_x/2)
    y_itr_new = (-fovy_new/2 + psize_y/2):psize_y:(fovy_new/2 - psize_y/2)
    #Create new image
    img_new = etp(y_itr_new, reverse(x_itr_new))*psize_x*psize_y
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
