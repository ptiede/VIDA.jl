module VIDAPlotsExt

using VIDA
if isdefined(Base, :get_extension)
    using Plots
    using RecipesBase
    using Measures
else
    using ..Plots
    using ..RecipesBase
    using ..Measures
end

font = Plots.font





@userplot Triptic
"""
    triptic(image::SpatialIntensityMap, θ::T) where {T<:AbstractTemplate}
A plot recipe for a triptic plot with the `image`, i.e. an EHT image and the template `\theta`
used for parameter estimate. The first two panels are the image and template and
the third are RA and DEC cross-sections of both the template and image centered at the
center of light.
"""
@recipe function f(h::Triptic)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: SpatialIntensityMap) ||
        !(typeof(h.args[2]) <: ComradeBase.AbstractModel)
        error("Triptic should be given a image and template.  Got: $(typeof(h.args))")
    end
    image, θ = h.args
    #link := :both
    widen := false
    xflip := true
    framestyle := [:none,:none,:axes]
    grid := false
    size --> (825,300)
    layout := (1,3)
    left_margin --> -2mm
    right_margin --> 2mm
    bottom_margin --> 3mm
    tickfontsize --> 10


    #Construct the image grid in μas
    g = axisdims(image)
    dataim = ComradeBase.baseimage(image./flux(image))'

    #Construct the template image
    template_img = intensitymap(θ, g)
    fimg = ComradeBase.baseimage(template_img/flux(template_img))'

    (;X, Y) = g
    Xitr = map(rad2μas, X)
    Yitr = map(rad2μas, Y)
    fovx, fovy = map(rad2μas, values(fieldofview(image)))


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

    #Image plot
    @series begin
        widen := false
        xflip --> true
        seriestype := :heatmap
        xlims --> (-fovx/2, fovx/2)
        ylims --> (-fovy/2, fovy/2)
        framestyle := :none
        seriescolor --> :afmhot
        colorbar := false
        subplot := 1
        #size := (300,300)
        aspect_ratio := 1
        xflip := true
        x = collect(Xitr)
        y = collect(Yitr)
        z = dataim
        x,y,z
    end

    @series begin
        subplot := 1
        seriestype := :hline
        legend := false
        linestyle := :solid
        framestyle := :none
        seriescolor := :cornflowerblue
        y = [ycol]
        widen := false
        aspect_ratio := 1
        linewidth := 2
        y
    end

    @series begin
        subplot := 1
        seriestype := :vline
        legend := false
        linestyle := :solid
        framestyle := :none
        seriescolor := :red
        y = [xcol]
        widen := false
        aspect_ratio := 1
        linewidth := 2
        y
    end


    @series begin
        fontsize --> 10
        aspect_ratio := 1
        subplot := 1
        seriestype := :line
        seriescolor := :white
        linewidth := 1
        annotations := [(startx-size/2,starty-4, "40 μas", font(10, color=:white)),
                        (startx-size/4, -starty+4, "Image", font(10, color=:white))]
        fontcolor := :white
        xlims = (-fovx/2,fovx/2)
        ylims = (-fovy/2,fovy/2)
        x = collect(xseg)
        y = yseg
        legend := false
        framestyle := :none
        left_margin := -2mm
        right_margin := -2mm
        x,y
    end

    #Now we move onto the template Image
    @series begin
        widen := false
        xflip --> true
        seriestype := :heatmap
        xlims --> (-fovx/2, fovx/2)
        ylims --> (-fovy/2, fovy/2)
        framestyle := :none
        seriescolor --> :afmhot
        colorbar := false
        subplot := 2
        aspect_ratio := 1
        x = collect(Xitr)
        y = collect(Yitr)
        z = fimg
        x,y,z
    end

    @series begin
        subplot := 2
        seriestype := :hline
        legend := false
        linestyle := :dash
        framestyle := :none
        aspect_ratio := 1
        seriescolor := :cornflowerblue
        y = [ycol]
        widen := false
        linewidth := 2
        y
    end

    @series begin
        subplot := 2
        seriestype := :vline
        legend := false
        linestyle := :dash
        framestyle := :none
        xflip := true
        aspect_ratio := 1
        seriescolor := :red
        y = [xcol]
        widen := false
        linewidth := 2
        y
    end

    #Now do the scale bar
    @series begin
        fontsize --> 10
        subplot := 2
        seriestype := :line
        seriescolor := :white
        linewidth := 1
        aspect_ratio := 1
        annotations := [(startx-size/2,starty-4, "40 μas", font(10, color=:white)),
                        (startx-size/4-2.5, -starty+4, "Template", font(10, color=:white))]
        fontcolor := :white
        xlims = (-fovx/2,fovx/2)
        ylims = (-fovy/2,fovy/2)
        x = collect(xseg)
        y = yseg
        xflip := true
        legend := false
        framestyle := :none
        widen := false
        x,y
    end

    #Now the chord plot
    @series begin
        subplot := 3
        seriestype := :line
        linestyle := :solid
        seriescolor := :cornflowerblue
        x = collect(reverse(Xitr))
        y = dataim[jmin,end:-1:1]
        legend := false
        x,y
    end
    #Now the chord plot
    @series begin
        subplot := 3
        seriestype := :line
        linestyle := :solid
        seriescolor := :red
        x = collect(Yitr)
        y = @view dataim[:, imin]
        xguide := "RA, DEC chords (μas)"
        yticks := false
        legend := false
        x,y
    end
    #Now the chord plot
    @series begin
        subplot := 3
        seriestype := :line
        linestyle := :dash
        seriescolor := :cornflowerblue
        x = collect(Xitr)
        y = @view fimg[jmin, :]
        legend := false
        x,y
    end
    #Now the chord plot
    @series begin
        subplot := 3
        seriestype := :line
        seriescolor := :red
        linestyle := :dash
        x = collect(Yitr)
        y = @view fimg[:, imin]
        xguide := "RA, DEC chords (μas)"
        yticks := false
        legend := false
        x,y
    end

end


VIDA.triptic(args...) = triptic(args...)


end
