module VIDAMakieExt

using VIDA
if isdefined(Base, :get_extension)
    using Makie
    using DimensionalData
    using ComradeBase: basedim
else
    using ..Makie
    using ..DimensionalData
    using ..ComradeBase: basedim
end

function add_scalebar!(ax, img, scale_length, color)
    fovx, fovy = fieldofview(img)
    x0 = (last(img.X))
    y0 = (first(img.Y))

    sl = (scale_length)
    barx = [x0 - fovx/32, x0 - fovx/32 - sl]
    bary = fill(y0 + fovy/32, 2)

    lines!(ax, barx, bary, color=color)
    text!(ax, (barx[1] + (barx[2]-barx[1])/2), bary[1]+fovy/64;
            text = "$(round(Int, sl)) μas",
            align=(:center, :bottom), color=color)
end


function _triptic(image::SpatialIntensityMap, θ::ComradeBase.AbstractModel;
                 figure=(size=(825, 300),), axis=(xreversed=true, aspect=DataAspect()),
                 colormap=:afmhot, scale_length=fieldofview(image).X/4)
    fig = Figure(;figure...)
    ax1 = Axis(fig[1,1]; axis...)
    ax2 = Axis(fig[1,2]; axis...)
    ax3 = Axis(fig[1,3], aspect=1, yticklabelsvisible=false)

    #Construct the image grid in μas
    g = axisdims(image)

    (;X, Y) = g
    Xitr = map(rad2μas, X)
    Yitr = map(rad2μas, Y)
    fovx, fovy = map(rad2μas, values(fieldofview(image)))

    dataim = IntensityMap(baseimage(image./flux(image)), (;X=Xitr, Y=Yitr))

    #Construct the template image
    template_img = intensitymap(θ, g)
    fimg = IntensityMap(baseimage(template_img/flux(template_img)), (;X=Xitr, Y=Yitr))


    #Get scale bar and slice data.
    size = 40#μas
    startx = fovx/2*(1-0.1)
    starty = -fovy/2*(1-0.2)
    xseg = range(startx, startx-size, length=4)
    yseg = starty*ones(4)
    #Finally the slices
    xx = collect(Xitr)
    yy = collect(Yitr)
    xcol,ycol = centroid(template_img)
    imin = argmin(abs.(xx .- xcol))
    jmin = argmin(abs.(yy .- ycol))

    color = Makie.to_colormap(colormap)[end]

    heatmap!(ax1, dataim, colormap=colormap)
    hlines!(ax1, [ycol], color=:cornflowerblue, linewidth=2, linestyle=:solid)
    vlines!(ax1, [xcol], color=:red, linewidth=2, linestyle=:solid)
    add_scalebar!(ax1, IntensityMap(parent(image), RectiGrid((X=rad2μas(image.X), Y=rad2μas(image.Y)))), rad2μas(scale_length), color)

    heatmap!(ax2, fimg, colormap=colormap)
    hlines!(ax2, [ycol], color=:cornflowerblue, linewidth=2, linestyle=:dash)
    vlines!(ax2, [xcol], color=:red, linewidth=2, linestyle=:dash)

    lines!(ax3, reverse(Xitr), dataim[jmin, :], color=:cornflowerblue)
    lines!(ax3, Yitr, dataim[:, imin], color=:red)

    lines!(ax3, reverse(Xitr), fimg[jmin, :], color=:cornflowerblue, linestyle=:dash)
    lines!(ax3, Yitr, fimg[:, imin], color=:red, linestyle=:dash)
    ax3.xlabel = "RA, DEC chords (μas)"

    return fig, (ax1, ax2, ax3)
end

VIDA.triptic(img, template; kwargs...) = _triptic(img, template; kwargs...)

end
