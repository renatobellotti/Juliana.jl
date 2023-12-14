# using CUDA
# using CUDAKernels
# using KernelAbstractions
# using Zygote


# """
# Check whether the voxel dose is above the threshold.

# Instead of a boolean, a piecewise linear response is given for a window approach:
# if dose <= threshold - R_low
#     return 0
# elif dose >= threshold + R_heigh
#     return 1
# else
#     return linear slope between the threshold and the closer margin
# """
# function dose_above_threshold(dose::T; threshold::T, R_low::T, R_high::T) where {T}
#     if dose <= threshold - R_low
#         return zero(T)
#     elseif dose >= threshold + R_high
#         return convert(T, 1.)
#     else
#         return (dose - (threshold - R_low)) / (R_high + R_low)
#     end
# end


# Zygote.@adjoint function dose_above_threshold(dose::T; threshold::T, R_low::T, R_high::T) where {T}
#     value = dose_above_threshold(dose, threshold=threshold, R_low=R_low, R_high=R_high)

#     function ddose_above_threshold(ddose::T)
#         dvalue = zero(T)
#         if ((threshold - R_low) <= dose) && (dose <= (threshold + R_high))
#             dvalue += ddose / (R_high + R_low)
#         end
#         return (dvalue,)
#     end

#     return value, ddose_above_threshold
# end


# @kernel function elementwise_above_threshold_kernel(output::AbstractArray{T, N}, dose::AbstractArray{T, N}, mask::AbstractArray{T, N}, d::T, v::T, R_low::T, R_high::T) where {T, N}
#     i = @index(Global)
#     if mask[i] == 1
#         output[i] = Juliana.dose_above_threshold(dose[i], threshold=d, R_low=R_low, R_high=R_high)
#     else
#         output[i] = zero(T)
#     end
# end


# device = CUDADevice()
# dvh_v_loss_kernel = elementwise_above_threshold_kernel(device, 32)


# function dvh_v_loss(dose::AbstractArray{T, N}; mask::AbstractArray{T, N}, d::T, v::T, R_low::T, R_high::T) where {T, N}
#     n_structure_voxels = sum(mask)

#     voxels_above_threshold = similar(dose);

#     event = dvh_v_loss_kernel(
#         voxels_above_threshold,
#         dose,
#         mask,
#         d,
#         v,
#         R_low,
#         R_high,
#         ndrange=size(dose),
#         dependencies=Event(device),
#     )
#     wait(event)
#     n_fulfilled = dot(mask, voxels_above_threshold)

#     return max(zero(T), v * n_structure_voxels - n_fulfilled)
# end


# Zygote.@adjoint function dvh_v_loss(dose::AbstractArray{T, N}; mask::AbstractArray{T, N}, d::T, v::T, R_low::T, R_high::T) where {T, N}
#     n_structure_voxels = sum(mask)

#     voxels_above_threshold = similar(dose);

#     event = dvh_v_loss_kernel(
#         voxels_above_threshold,
#         dose,
#         mask,
#         d,
#         v,
#         R_low,
#         R_high,
#         ndrange=size(dose),
#         dependencies=Event(device),
#     )
#     wait(event)
#     n_fulfilled = dot(mask, voxels_above_threshold)

#     value =  max(zero(T), v * n_structure_voxels - n_fulfilled)
    
#     # Copy to CPU for fast elementwise access.
#     voxels_above_threshold = collect(voxels_above_threshold)

#     function dvh_v_loss_adjoint(ddose::AbstractArray{T, 1}) where {T}
#         # We know the gradient exactly.
#         dval = zeros(T, length(dose))
#         if value > zero(T)
#             for i in 1:length(dose)
#                 if (voxels_above_threshold[i] > 0) && (voxels_above_threshold[i] < 1)
#                     # We have a non-vanishing gradient.
#                     dval[i] = -1 / (R_low + R_high) * ddose
#                 end
#             end
#         end
#         return (dval,)
#     end
    
#     return value, dvh_v_loss_adjoint
# end


using Zygote


"""
Check whether the voxel dose is above the threshold.

Instead of a boolean, a piecewise linear response is given for a window approach:
if dose <= threshold - R_low
    return 0
elif dose >= threshold + R_heigh
    return 1
else
    return linear slope between the threshold and the closer margin
"""
function dose_above_threshold(dose::T; threshold::T, R_low::T, R_high::T) where {T}
    if dose <= threshold - R_low
        return zero(T)
    elseif dose >= threshold + R_high
        return convert(T, 1.)
    else
        return (dose - (threshold - R_low)) / (R_high + R_low)
    end
end


Zygote.@adjoint function dose_above_threshold(dose::T; threshold::T, R_low::T, R_high::T) where {T}
    value = dose_above_threshold(dose, threshold=threshold, R_low=R_low, R_high=R_high)

    function ddose_above_threshold(ddose::T)
        dvalue = zero(T)
        if ((threshold - R_low) <= dose) && (dose <= (threshold + R_high))
            dvalue += ddose / (R_high + R_low)
        end
        return (dvalue,)
    end

    return value, ddose_above_threshold
end


# The types are not a typo (no pun intended).
# For some reason the Zygote gradient of this function does not compile with CUDA arrays...
function dvh_v_loss(dose::AbstractArray{T, N}; mask::AbstractArray{T, N}, d::T, v::T, R_low::T, R_high::T) where {T, N}
    n_structure_voxels = sum(mask)

    # collect() is needed to be able to use the ". operator".
    # The second collect() is needed because we cannot mix CUDA arrays and CPU arrays.
    voxels_above_threshold = Juliana.dose_above_threshold.(collect(dose), threshold=d, R_low=R_low, R_high=R_high)
    n_fulfilled = dot(collect(mask), voxels_above_threshold)

    return max(zero(T), v * n_structure_voxels - n_fulfilled)
end
