@doc """

    $(SIGNATURES)
Creates an EHTImage type from the filter type. The number
of pixels in the image are given by `npix` and the field
of view in μas in the x and y direction are given by `xlim`
and `ylim`. The rest of the options are the default image
characteristics

"""
function make_image(θ::AbstractFilter,
    npix::Int, xlim, ylim;
    intensity=1.0,
    source="M87", wavelength=0.001320260,
    ra=187.7059307575226,
    dec = 12.391123223919932,
    mjd=57854.0
    )
xitr,yitr,img = filter_image(θ, npix, xlim, ylim)
scale_norm = intensity/sum(img)
X = collect(xitr)
Y = collect(yitr)
psize_x = X[2]-X[1]
psize_y = Y[2]-Y[1]
imgeht = EHTImage(npix, npix, psize_x, psize_y,
   source, ra, dec, wavelength,
   mjd, img*scale_norm)
return imgeht
end

@doc """

    $(SIGNATURES)
Creates an EHTImage type from the filter type. The number
of pixels in the image are given by `npix` and the field
of view in μas in the x and y direction are given by `xlim`
and `ylim`. To define the source we use an `source_img`.

"""
function make_image(θ::AbstractFilter,
    npix::Int, xlim, ylim,
    source_img::EHTImage;
    intensity=1.0
    )
xitr,yitr,img = filter_image(θ, npix, xlim, ylim)
scale_norm = 1.0/sum(img)
X = collect(xitr)
Y = collect(yitr)
psize_x = X[2]-X[1]
psize_y = Y[2]-Y[1]
imgeht = EHTImage(npix, npix, psize_x, psize_y,
   source_img.source, source_img.ra,
   source_img.dec, source_img.wavelength,
   source_img.mjd, img*scale_norm)
return imgeht
end
