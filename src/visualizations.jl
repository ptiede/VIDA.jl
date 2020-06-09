

function setfonts(fontsize=12)
    PyPlot.matplotlib.rc("mathtext",fontset="cm")        #computer modern font
    PyPlot.matplotlib.rc("font",family="serif",size=fontsize)  #font similar to LaTeX
    PyPlot.matplotlib.rc("text", usetex=true)
end

"""
Plots the image being read in.
"""
function PyPlot.plot(image::AbstractImage)
    setfonts()
    fig,ax = PyPlot.subplots(1)
    ax.set_aspect(:equal)
    ax.set_xlabel("RA (pixel)")
    ax.set_ylabel("DEC (pixel)")
    pim = ax.imshow(image.img, cmap=:afmhot,
                   vmin=0.0, vmax=maximum(image.img), origin="lower")
    fig.colorbar(pim, ax=ax, label="Jy/pixel")

    return fig
end



"""
$(SIGNATURES)

where `image` is templated off of EHTImage struct.

# Details
This was created to be close to the ehtim display object. It takes an
EHTImage object and plots it using PyPlot, returning the figure.

Note that is does not save the figure.

"""
function PyPlot.plot(image::T; xlim=(-50,50), ylim=(-50,50)) where {T<:EHTImage}
    setfonts()
    μas2rad = π/180.0/3600/1e6
    #Construct the image grid in μas
    psizeuas_x = -image.psize_x
    psizeuas_y = image.psize_y
    fovx = psizeuas_x*image.nx
    fovy = psizeuas_y*image.ny
    brightness_temp = image.wavelength^2/(2*k_B.val)/1e26/abs(μas2rad^2*image.psize_x*image.psize_y)
    fig,ax = PyPlot.subplots(1)
    pim = ax.imshow(image.img*brightness_temp/1e9, cmap=:afmhot,
                    extent = [fovx/2,-fovx/2, -fovy/2, fovy/2],
                    vmin=0.0, vmax=maximum(image.img)*brightness_temp/1e9,
                    origin="lower")
    ax.set_xlabel(L"$\Delta$RA  ($\mu$as)")
    ax.set_ylabel(L"$\Delta$DEC ($\mu$as)")
    fig.colorbar(pim,ax=ax, label=L"Brightness Temperature ($10^9$ K)")
    title_leg = image.source*" "*string(round(c_0.val/image.wavelength/1e9))*" GHz"
    ax.set_title(title_leg)
    return fig
end

function make_ehtimage(θ::AbstractFilter,
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


function make_ehtimage(θ::AbstractFilter,
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
    println(psize_x," ",psize_y)
    imgeht = EHTImage(npix, npix, psize_x, psize_y,
                      source_img.source, source_img.ra,
                      source_img.dec, source_img.wavelength,
                      source_img.mjd, img*scale_norm)
    return imgeht
end




function PyPlot.plot(θ::T, npix::Int=60,
                     xlim=[-40,40], ylim=[-40,40];
                     kwargs...) where {T<:AbstractFilter}
    fig,ax = PyPlot.subplots(1)
    setfonts()

    _,_,img = filter_image(θ, npix, xlim, ylim)
    fovx = (xlim[2]-xlim[1])
    fovy = (ylim[2]-ylim[1])
    pim = ax.imshow(img, cmap=:afmhot,
                    extent = [fovx/2,-fovx/2, -fovy/2, fovy/2],
                    vmin=0.0, vmax=maximum(img),
                    origin="lower"; kwargs...)
    ax.set_xlabel(L"$\Delta$RA  ($\mu$as)")
    ax.set_ylabel(L"$\Delta$DEC ($\mu$as)")

    return fig
end

"""
    $(SIGNATURES)
A plotting function used for benchmarking.
"""
function plot_filter(measure, θ::T, p1, p2, n) where {T<:AbstractFilter}
    ℓ = Matrix{Float64}(undef,n,n)
    filter_benchmark!(ℓ, measure, θ, p1, p2, n)

    println("Maximum likelihood: ", maximum(ℓ))
    fig, ax = PyPlot.subplots(1)
    pim = ax.imshow(ℓ, extent=[p1[1],p1[2],p2[1],p2[2]],
                    origin="lower", cmap=:viridis)
    ax.set_aspect(:equal,:box)
    ax.set_xlabel("r0")
    ax.set_ylabel("ξ")
    ax.axvline(θ.r0, color="white")
    ax.axhline(θ.ξ, color="white")
    fig.colorbar(pim,ax=ax)
    savefig("test.png")
    #PyPlot.savefig("test.png")
    return ℓ

end

"""
    $(SIGNATURES)
Plots a triptic plot with the `image`, i.e. an EHT image and the filter `\theta`
used for parameter estimate. The first two panels are the image and filter and
the third are RA and DEC cross-sections of both the filter and image.
"""
function plot_triptic(image::EHTImage, θ::AbstractFilter)
    setfonts()
    fig,ax = PyPlot.subplots(1,3, figsize=(7.9,3))
    fig.subplots_adjust(left=0.01,right=0.98, top=0.98,
                        bottom=0.16, wspace=0.05, hspace=0.01)

    #Construct the image grid in μas
    psizeuas_x = -image.psize_x
    psizeuas_y = image.psize_y
    fovx = psizeuas_x*image.nx
    fovy = psizeuas_y*image.ny
    dataim = image.img/sum(image.img)
    pim = ax[1].imshow(dataim, cmap=:afmhot,
                    extent = [fovx/2,-fovx/2, -fovy/2, fovy/2],
                    origin="lower", vmin=0.0, vmax=maximum(dataim))
    ax[1].text(0.15, 0.95, "Image", horizontalalignment=:center,
               verticalalignment=:center, transform=ax[1].transAxes,
               color=:white)
    ax[1].set_yticks([]);ax[1].set_xticks([])

    xitr,yitr,img = filter_image(θ, image.nx, [-fovx/2,fovx/2], [-fovy/2,fovy/2])
    img .= img*abs(psizeuas_y*psizeuas_x)
    img .= img/sum(img)
    pim = ax[2].imshow(img, cmap=:afmhot,
                    extent = [fovx/2,-fovx/2, -fovy/2, fovy/2],
                    vmin=0.0, vmax=maximum(dataim),
                    origin="lower")
    ax[2].set_yticks([]);ax[2].set_xticks([])
    ax[2].text(0.15, 0.95, "Filter", horizontalalignment=:center,
               verticalalignment=:center, transform=ax[2].transAxes,
               color=:white)
    xx = collect(xitr)
    yy = collect(yitr)
    xcol,ycol = center_of_light(xitr,yitr,img)
    imin = argmin(abs.(xx .- xcol))
    jmin = argmin(abs.(yy .- ycol))
    ax[1].axhline(xcol, color=:cornflowerblue)
    ax[1].axvline(xcol, color=:C3)
    ax[2].axhline(ycol, color=:cornflowerblue, linestyle="--")
    ax[2].axvline(ycol, color=:C3, linestyle="--")
    ax[3].plot(xx,@view(dataim[jmin,image.ny:-1:1]),
               label = "RA chord", color=:cornflowerblue)
    ax[3].plot(yy,@view(dataim[:,imin]),
               label = "DEC chord", color=:C3)
    ax[3].plot(xx,@view(img[jmin,image.ny:-1:1]), linestyle="--",
               label = "RA chord", color=:cornflowerblue)
    ax[3].plot(yy,@view(img[:,imin]), linestyle="--",
               label = "DEC chord", color=:C3)
    ax[3].set_yticks([])
    ax[3].set_xlabel(L"X,Y chords ($\mu$as)")

    #make the scale bar
    size = 40#μas
    startx = fovx/2*(1-0.1)
    starty = -fovy/2*(1-0.2)
    scale_bar1 = PyPlot.matplotlib.lines.Line2D((startx, startx-size),
                       (starty,starty), color=:white, alpha=1.0,
                       linewidth=1, markeredgecolor=:white)
    scale_bar2 = PyPlot.matplotlib.lines.Line2D((startx, startx-size),
                       (starty,starty), color=:white, alpha=1.0,
                       linewidth=1, markeredgecolor=:white)
    ax[1].add_artist(scale_bar1)
    ax[1].text(0.1, 0.12, latexstring(size)*L"\,$\mu\mathrm{as}$",
               transform=ax[1].transAxes, color="white")
    ax[2].add_artist(scale_bar2)
    ax[2].text(0.1, 0.12, latexstring(size)*L"\,$\mu\mathrm{as}$",
               transform=ax[2].transAxes, color="white")

    return fig
end
