"""
    $(SIGNATURES)

Creates the specfied divergence base on the image, and divergence type.
The current options for the divergence are:
 - :KL for the Kullback-Leiber divergence
 - :Bh for the Bhattacharyya divergence, which is related to the Hellinger distance.
# Details
This returns a close that you can use to calculate the divergence given some image filter θ.
"""
function make_div(image::T, div_type, lreg=1.0) where {T<:EHTImage}
  if  div_type == :KL
    return make_kl(image, lreg)
  elseif div_type == :Bh
    return make_bh(image, lreg)
  else
    ArgumentError("div_type must be implemented.\n The current methods are: :KL, :Bh")
  end
end



"""
    $(SIGNATURES)

Closure that makes a Bhattacharyya divergence function based on the `imagefilter`
and intrinsic image, `image`.

# Details
This takes an imagefilter function and my intrinsic image and creates an
function that return the Bhattacharyya divergence of the filter and the image.

"""
function make_bh(image::T, lreg=1.0) where {T<:EHTImage}
    dx = abs(image.psize_x)
    dy = abs(image.psize_y)
    fovx = dx*image.nx
    fovy = dx*image.ny
    xstart = (fovx-dx)/2.0
    ystart = (-fovy+dy)/2.0
    intensity_norm = sum(image.img)

    @assert intensity_norm > 0 "make_bh: Something is wrong with the image. The image flux is zero"
    #Normalize the image to make it a probability dist
    norm_im = image.img/abs(intensity_norm*dx*dy)
    @fastmath function (θ::AbstractFilter)
        filter_norm = zero(eltype(norm_im))
        bsum = zero(eltype(norm_im))
        for i in 1:image.nx
            for j in 1:image.ny
                x = xstart + image.psize_x*(i-1)
                y = ystart + image.psize_y*(j-1)
                filter_value = θ(x,y)+1e-50
                @inbounds bsum += sqrt(filter_value*norm_im[j,i])
                filter_norm += filter_value
            end
        end
        return -lreg*log(bsum/sqrt(filter_norm/(dx*dy)))
    end
end

"""
    $(SIGNATURES)

Closure that makes a Kullback-Leibler divergence based on the `imagefilter`
and intrinsic image, `image`. That is it compute
            D_{KL}(filter||image)

# Details
This takes an imagefilter function and my intrinsic image and creates an
function that returns the KL-divergence of the image, where the image is the
distribution we are calculating the KL-divergence relative to.

"""
function make_kl(image::T, lreg=1.0) where {T<:EHTImage}
    dx = abs(image.psize_x)
    dy = abs(image.psize_y)
    fovx = dx*image.nx
    fovy = dx*image.ny
    xstart = (fovx-dx)/2.0
    ystart = (-fovy+dy)/2.0
    intensity_norm = sum(image.img)

    @assert intensity_norm > 0 "make_kl: Something is wrong with the image. The image flux is zero"
    #Normalize the image to make it a probability dist
    norm_im = image.img/abs(intensity_norm*dx*dy)
    function (θ::AbstractFilter)
        filter_norm = zero(eltype(norm_im))
        klsum = zero(eltype(norm_im))
        @inbounds for i in 1:image.nx
            @inbounds for j in 1:image.ny
                x = xstart + image.psize_x*(i-1)
                y = ystart + image.psize_y*(j-1)
                filter_value = θ(x,y)+1e-12
                klsum += filter_value*log(filter_value/(norm_im[j,i]+1e-12))
                filter_norm += filter_value
            end
        end
        return lreg*(klsum/filter_norm - log(filter_norm*dx*dx))
    end
end
