using RecipesBase
using Measures
font = Plots.font



"""
    plot(image::EHTImage)

where `image` is templated off of EHTImage struct.

# Details
This was created to be close to the ehtim display object. It takes an
EHTImage object and plots it according to EHT conventions.

Note that is does not save the figure.
"""
@recipe function f(image::EHTImage)

    #Define some constants
    μas2rad = π/180.0/3600/1e6
    #Construct the image grid in μas
    psizeuas_x = -image.psize_x
    psizeuas_y = image.psize_y
    fovx = psizeuas_x*image.nx
    fovy = psizeuas_y*image.ny

    brightness_temp = image.wavelength^2/(2*KB)/1e26/
                      abs(μas2rad^2*image.psize_x*image.psize_y)

    #brightness_temp = image.wavelength^2/(2*k_B.val)/1e26/
    #                  abs(μas2rad^2*image.psize_x*image.psize_y)

    tickfontsize --> 11
    guidefontsize --> 14
    size --> (500,400)
    xaxis --> "ΔRA  (μas)"
    yaxis --> "ΔDEC (μas)"
    seriescolor --> :afmhot
    aspect_ratio --> 1
    bar_width --> 0
    xlims --> (-fovx/2,fovx/2)
    ylims --> (-fovy/2,fovy/2)
    #left_margin --> -2mm
    right_margin --> 5mm
    x = collect(range(-fovx/2, fovx/2; length=image.nx))
    y = collect(range(-fovy/2, fovy/2; length=image.ny))
    z = image.img[:,end:-1:1]*brightness_temp/1e9
    seriestype := :heatmap
    fontfamily --> "sans serif"
    colorbar_title --> "Brightness Temperature (10⁹ K)"
    xflip --> true
    widen := false
    #framestyle --> :box
    title --> image.source*" "*string(round(C0/image.wavelength/1e9))*" GHz"
    linecolor-->:black
    tick_direction --> :out
    x,y,z
end

@doc """
    triptic(img::EHTImage, template)
Plots a triptic where the left panel is the `img` middle
the `template` and the right a two cross-sections of the image
and template
"""
function triptic end

@userplot Triptic
"""
    triptic(image::EHTImage, θ::T) where {T<:AbstractTemplate}
A plot recipe for a triptic plot with the `image`, i.e. an EHT image and the template `\theta`
used for parameter estimate. The first two panels are the image and template and
the third are RA and DEC cross-sections of both the template and image centered at the
center of light.
"""
@recipe function f(h::Triptic)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: EHTImage) ||
        !(typeof(h.args[2]) <: AbstractTemplate)
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
    psizeuas_x = -image.psize_x
    psizeuas_y = image.psize_y
    fovx = psizeuas_x*image.nx
    fovy = psizeuas_y*image.ny
    xitr = range(-fovx/2, fovx/2; length=image.nx)
    yitr = range(-fovy/2, fovy/2; length=image.ny)
    dataim = @view(image.img[:,end:-1:1])/sum(image.img)

    #Construct the template image
    template_img = make_image(θ, image.nx, [-fovx/2,fovx/2], [-fovy/2,fovy/2], image;
                                    intensity=1)
    fimg = template_img.img


    #Get scale bar and slice data.
    size = 40#μas
    startx = fovx/2*(1-0.1)
    starty = -fovy/2*(1-0.2)
    xseg = range(startx, startx-size, length=4)
    yseg = starty*ones(4)
    #Finally the slices
    xx = collect(xitr)
    yy = collect(yitr)
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
        x = collect(xitr)
        y = collect(yitr)
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
        x = collect(xitr)
        y = collect(yitr)
        z = @view(fimg[:,end:-1:1])
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
                        (startx-size/4, -starty+4, "Template", font(10, color=:white))]
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
        x = collect(reverse(xitr))
        y = @view dataim[jmin,image.ny:-1:1]
        legend := false
        x,y
    end
    #Now the chord plot
    @series begin
        subplot := 3
        seriestype := :line
        linestyle := :solid
        seriescolor := :red
        x = collect(yitr)
        y = @view dataim[:, image.nx-imin]
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
        x = collect(xitr)
        y = @view fimg[jmin,image.ny:-1:1]
        legend := false
        x,y
    end
    #Now the chord plot
    @series begin
        subplot := 3
        seriestype := :line
        seriescolor := :red
        linestyle := :dash
        x = collect(yitr)
        y = @view fimg[:, image.nx-imin]
        xguide := "RA, DEC chords (μas)"
        yticks := false
        legend := false
        x,y
    end

end

"""
    plot(θ::AbstractTemplate)
    plot(θ::AbstractTemplate; npix, fovx, fovy)
Plots the template `θ` using the usual EHT conventions.
The default image will use 128x128 pixels with a 120x120 field of view.
"""
@recipe function f(θ::AbstractTemplate; npix=128, fovx=120, fovy=120)

    tickfontsize --> 11
    guidefontsize --> 14
    size --> (400,400)
    xaxis --> "ΔRA  (μas)"
    yaxis --> "ΔDEC (μas)"
    seriescolor --> :afmhot
    aspect_ratio --> 1
    bar_width --> 0
    xlims --> (-fovx/2,fovx/2)
    ylims --> (-fovy/2,fovy/2)
    xitr,yitr,img = template_image(θ, npix, [-fovx/2,fovx/2],[-fovy/2,fovy/2])
    #left_margin --> -2mm
    right_margin --> 5mm
    x = collect(reverse(xitr))
    y = collect(yitr)
    z = img[:,end:-1:1]/sum(img)
    seriestype := :heatmap
    xflip --> true
    widen := false
    #framestyle --> :box
    linecolor-->:black
    tick_direction --> :out
    colorbar --> false
    x,y,z
end

export plot, plot_triptic
