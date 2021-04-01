"""
    $(SIGNATURES)
Creates an npix×npix rasterized image of the template `θ` with
limits `xlim` and `ylim`

Returns the tuple (xitr,yitr,image) where xitr,yitr are the iterators
defining the pixel locations (which are centered) and the rasterized image,
 in Jy/μas^2.

# Note
I use the pixel size definition field_of_view/npix, but the image is evaluated
at the pixel centers.

We also use the astronomer orientation and ordering.

"""
function template_image(θ::AbstractTemplate,
                      npix::Int, xlim, ylim)
    fovx = xlim[2]-xlim[1]
    fovy = ylim[2]-ylim[1]

    px = fovx/(npix)
    py = fovy/(npix)

    xitr = (fovx/2-px/2):-px:(-fovx/2)
    yitr = (-fovy/2+py/2):py:(fovy/2)
    img = Matrix{Float64}(undef,npix,npix)
    for (i,x) in enumerate(xitr)
        for (j,y) in enumerate(yitr)
            img[j,i] = θ(x,y)
        end
    end
    return (xitr,yitr,img)
end

#
#    $(SIGNATURES)
#Find the minimum square distance between an ellipse centered at (0,0) with semi-major
#axis `a` and semi-minor axis `b` and the point (`x`, `y`).
#Uses an iterative method with accuracy `ϵ` which defaults to 1e-6.
## Credit:
#Algorithm taken from
#https://github.com/0xfaded/ellipse_demo/issues/1#issuecomment-405078823
#except written in Julia and using an adaptive termination condition for accuracy
#"""
#@
@inline function ellipse_sqdist(x,y, a, b, ϵ=1e-6)
    #For simplicity we will only look at the positive quadrant
    px = abs(x)
    py = abs(y)

    #initial guess
    tx = 0.707
    ty = 0.707
    err = 1.0
    n = 0
    while err > ϵ && n < 10
        tx′, ty′ = dist_ellipse_unit(px,py,tx,ty,a,b)
        #err = hypot((tx′-tx), ty′-ty)
        err = sqrt((tx-tx′)*(tx-tx′) + (ty-ty′)*(ty-ty′))
        tx = tx′
        ty = ty′
        n+=1
    end
    return (a*tx-px)*(a*tx-px) + (b*ty-py)*(b*ty-py)
end

#Finds the closest point on the ellipse. This is an internal function
#used for ellipse_distance.
#
 @inline function dist_ellipse_unit(px, py, tx, ty, a, b)
    x′ = a*tx
    y′ = b*ty

    ex = (a*a - b*b)*tx*tx*tx/a
    ey = (b*b - a*a)*ty*ty*ty/b

    rx = x′ - ex
    ry = y′ - ey

    qx = px - ex
    qy = py - ey

    r = sqrt(rx*rx + ry*ry)
    q = sqrt(qx*qx + qy*qy)

    xx = min(1.0, max(0.0, (qx*r/q + ex)/a) )
    yy = min(1.0, max(0.0, (qy*r/q + ey)/b) )
    #t = hypot(xx,yy)
    t = sqrt(xx*xx + yy*yy)
    xx /= t
    yy /= t

    return xx,yy
end


#Rotates our points. Note that we use astronomer conventions
@inline function rotate(x,y,ξ)
    s,c = sincos(π-ξ)
    x′ = c*x - s*y
    y′ = s*x + c*y

    return (x′ ,y′)
end
