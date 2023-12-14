using CUDA
using KernelAbstractions
using CUDAKernels
using Interpolations
using Statistics
using Zygote

# Units: Gy always!!!


struct OptimisationConfiguration{NumberType, VolumeType, DijType}
    normalisationDose::NumberType
    normalisationStructureMask::VolumeType
    ct::VolumeType
    idealDose::VolumeType
    importance::VolumeType
    minimumDose::VolumeType
    Dij::DijType
    structures::Dict{String, VolumeType}
    prescriptions::Juliana.Prescriptions
end


function mean_dose(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}) where {T, N}
    if sum(mask) == 0
        return zero(T)
    else
        # Copy to the CPU for the element-wise access.
        structure_dose = dose[mask .== one(T)]
        return Statistics.mean(structure_dose)
    end
end


function minimum_dose(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}) where {T, N}
    # Copy to the CPU to allow for fast indexing.
    # mask = Array(mask)
    # structure_dose = Array(dose)[mask .== one(T)]
    structure_dose = dose[mask .== one(T)]
    return minimum(structure_dose)
end


function variance_dose(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}) where {T, N}
    μ = mean_dose(dose, mask)
    structure_dose = dose[mask .== one(T)]
    n_voxels_in_structure = size(structure_dose, 1)
    if n_voxels_in_structure > zero(T)
        return sum((structure_dose .- μ).^2) ./ n_voxels_in_structure
    else
        return zero(T)
    end
end


function maximum_dose(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}) where {T, N}
    return maximum(dose .* mask)
end


function my_percentile(x::AbstractArray{T, 1}, percentage::T) where {T}
    x = sort(x)
    i = max(round(percentage / convert(T, 100) * length(x)), one(T))
    i = convert(Integer, i)
    return x[i]
end

function my_percentile(x::AbstractArray{T, 1}, percentage::AbstractArray{T, 1}) where {T}
    return [my_percentile(x, p) for p in percentage]
end


function dvh_d(dose, mask, v::T) where {T}
    """
    Calculates D_<v>%, i. e. the dose to the coldest v% of the voxels in the
    given structure.
    """
    # We copy to the host because quantile needs access by element,
    # which is very slow on the GPU.
    mask_array = Array(mask)
    structure_dose = reduce(vcat, Array(dose)[mask_array .== 1])
    p = my_percentile(structure_dose, v)
    return p
end


# Need to reimplement linear interpolation because Zygote cannot deal with it.
function linear_interpolation(x::AbstractArray{T, N}, y::AbstractArray{T, N}, x_eval::T) where {T, N}
    if x_eval <= x[1]
        return y[1]
    end
    if x_eval >= x[end]
        return y[end]
    end
    i = 1
    while x_eval > x[i]
        i = i + 1
    end

    x0 = x[i-1]
    x1 = x[i]
    y0 = y[i-1]
    y1 = y[i]

    return y0 + (x_eval - x0) * (y1 - y0) / (x1 - x0 + eps(T))
end


function linear_interpolation(x::AbstractArray{T, N}, y::AbstractArray{T, N}, x_eval::AbstractArray{T}) where {T, N}
    out = zero(T, x_eval)
    for i in size(x)[1]
        if x_eval[i] <= x[1]
            out[i] = y[1]
            continue
        end
        if x_eval[i] >= x[end]
            out[i] = y[end]
            continue
        end
        k = 1
        while x_eval[i] > x[k]
            k = k + 1
        end

        x0 = x[k-1]
        x1 = x[k]
        y0 = y[k-1]
        y1 = y[k]

        out[i] = y0 + (x_eval - x0) * (y1 - y0) / (x1 - x0 + eps(T))
    end
    return out
end


function dvh_v(dose, mask, d)
    """
    Calculates an approximation to V_<d>Gy, i. e. the volume that receives
    at least <d>Gy dose.
    """
    volume_fractions = collect(LinRange(0.f0, 100.f0, 401))
    dose_values = dvh_d(dose, mask, volume_fractions)

    f = Interpolations.linear_interpolation(dose_values, volume_fractions[end:-1:1], extrapolation_bc=Interpolations.Flat())
    return f(d)
end


function normalise_dose(dose, mask, normalisation_dose)
    return dose .* (normalisation_dose / mean_dose(dose, mask))
end


function normalisation_loss(dose::AbstractArray{T, N}, normalisationStructureMask::AbstractArray{T, N}, normalisationDose::T) where {T, N}
    return (mean_dose(dose, normalisationStructureMask) - normalisationDose)^2
end


function hotspot_loss(dose::AbstractArray{T, N}, threshold::T) where {T, N}
    return sum(max.(dose .- threshold, zero(T)))
end


function hotspot_loss(dose::AbstractArray{T, N}, threshold::AbstractArray{T, N}) where {T, N}
    # Voxel-wise maximum value.
    return sum(max.(dose .- threshold, zero(T)))
end


function coldspot_loss(dose::AbstractArray{T, N}, threshold::AbstractArray{T, N}) where {T, N}
    # Voxel-wise minimum value.
    return sum(max.(threshold .- dose, zero(T)))
end


# function dvh_d_loss(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}, d::T, v::T) where {T, N}
#     # Objective for ensuring that Dv% <= d.
#     # "The coldest (100-v)% of structure voxels receive no more than dose d."
#     # That means that a fraction of v% of voxels is allowed to violate the constraint.
#     # If we clip the achieved dose values from above at d, resulting in d_clip,
#     # then a fulfilling configuration will have sum(d_clip) <= v% * d


#     # This means that a fraction v of the voxels have at most dose d.
#     # In other words, the sum of the dose in target voxels clamped from below at d
#     # should be n_target_voxels * (1 - v) or smaller. Everything above means we
#     # do not fulfill the condition.
#     structure_dose = dose .* mask
#     n_target_voxels = sum(mask)

#     excess_dose = max.(zero(T), structure_dose .- d)
#     loss = max.(0., deviation)
#     # structure_dose_clamped = clamp.(structure_dose, d, Inf)

#     # loss = max(
#     #     zero(T),
#     #     sum(structure_dose_clamped) - (1 - convert(T, v)) * n_target_voxels * d,
#     # )

"""Loss function that ensures Dmax <= threshold for the given structure."""
function Dmax_loss(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}, threshold::T) where {T, N}
    structure_dose = dose .* mask
    excess_dose = max.(zero(T), structure_dose .- threshold)
    loss = sum(excess_dose)
end


"""Loss function that ensures Dmin <= threshold for the given structure."""
    function Dmin_loss(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}, threshold::T) where {T, N}
        lacking_dose = max.(zero(T), threshold .- dose)
        lacking_dose_structure = lacking_dose .* mask
        loss = sum(lacking_dose_structure)
    end


function build_loss_parts(dose::AbstractArray{T, N}, config::OptimisationConfiguration, safety_margin::T) where {T, N}
    loss_parts = Dict{String, T}()
    # hottest_target_dose = maximum([dose for (target_name, dose) in config.prescriptions.target_doses])

    # Normalise the dose distribution to make sure that a fulfilled prescription is
    # really fulfilled even after renormalisation.
    dose = dose .* (config.normalisationDose / mean_dose(dose, config.normalisationStructureMask))

    ideal_dose_loss = dot(
        config.importance,
        (config.idealDose .- dose).^2,
    ) / sum(config.importance) / convert(T, 25)
    loss_parts["ideal_dose_loss"] = ideal_dose_loss

    maximum_loss = hotspot_loss(dose, convert(T, 1.05) .* config.normalisationDose) * convert(T, 10.)
    loss_parts["maximum_loss"] = maximum_loss

    # for (target_name, target_dose) in config.prescriptions.target_doses
    #     target_mask = config.structures[target_name]
    #     # V95_loss = dvh_v_loss(
    #     #     dose,
    #     #     mask=target_mask,
    #     #     d=convert(T, 0.95) * target_dose,
    #     #     v=convert(T, 0.95),
    #     #     # This constant regulates how much we approximate a voxel being "above" or "below" a threshold.
    #     #     # Go 15% below the threshold, but not above.
    #     #     R_low=convert(T, 0.15*target_dose),
    #     #     R_high=zero(T),
    #     # )
    #     # loss_parts["V95_loss_$target_name"] = V95_loss

        # Ensure minimum target coverage: V80% = 100%,
        # i. e. no target voxel should receive less than 80% dose.
        # minimum_loss = Dmin_loss(
        #     dose,
        #     target_mask,
        #     convert(T, 0.90) * target_dose,
        # )
        # loss_parts["minimum_loss_$target_name"] = minimum_loss / convert(T, 500.)
    # end

    loss_parts["minimum_loss"] = sum(max.(config.minimumDose .- dose, 0)) / 100

    loss_parts["normalisation_variance"] = variance_dose(
        dose,
        config.normalisationStructureMask,
    )

    # TODO: Add this for the highest dose target and cut the lower dose target.
    # mean_loss = (mean_dose(dose, target_mask) - target_dose)^2

    # Ensure homogeneity.
    # homogeneity_loss = 0.1 * variance_dose(dose, config.normalisationStructureMask)
    # loss_parts["homogeneity_loss"] = homogeneity_loss

    # OAR loss parts.
    for constraint in config.prescriptions.constraints
        if constraint.priority == soft
            continue
        end

        if constraint.kind == constraint_mean
            mean = mean_dose(dose, config.structures[constraint.structure_name])
            l = max(mean - (constraint.dose - safety_margin), zero(T))
            loss_parts["$(constraint.structure_name)_mean_loss"] = l
        elseif (constraint.kind == constraint_max) || (((constraint.kind == constraint_dvh_d) || (constraint.kind == constraint_dvh_v)) && (constraint.volume ≈ convert(T, 0.02)))
            l = Dmax_loss(
                dose,
                config.structures[constraint.structure_name],
                constraint.dose - safety_margin,
            )
            loss_parts["$(constraint.structure_name)_max_loss"] = l
        end
    end
    
    return loss_parts
end


# Version with logging.
function dose_loss!(dose::AbstractArray{T, N}, config::OptimisationConfiguration, loss_parts::Dict{String, T}, subloss_weights::Dict{String, T}) where {T, N}
    # OAR safety margin --> Make sure the dose threshold are kept by at least this margin.
    safety_margin = convert(T, 0.2)
    partial_losses = build_loss_parts(dose, config, safety_margin)

    for (key, value) in partial_losses
        loss_parts[key] = value
    end

    loss = zero(T)
    for (name, value) in loss_parts
        loss += subloss_weights[name] * value
    end

    return loss
end

# Version without logging.
function dose_loss(dose::AbstractArray{T, N}, config::OptimisationConfiguration, subloss_weights::Dict{String, T}) where {T, N}
    # OAR safety margin --> Make sure the dose threshold are kept by at least this margin.
    safety_margin = convert(T, 0.2)
    partial_losses = build_loss_parts(dose, config, safety_margin)

    loss = zero(T)
    for (name, value) in partial_losses
        loss += subloss_weights[name] * value
    end

    return loss
end


function dose_loss_gradient(dose, config::OptimisationConfiguration)
    return Zygote.gradient(dose -> dose_loss(dose, config), dose)[1]
end


function loss!(w, config::OptimisationConfiguration, loss_parts, subloss_weights::Dict{String, T}) where {T}
    return dose_loss!(config.Dij * w, config, loss_parts, subloss_weights)
end


function loss(w, config::OptimisationConfiguration, subloss_weights::Dict{String, T}) where {T}
    return dose_loss(config.Dij * w, config, subloss_weights)
end


function loss_gradient(w, config::OptimisationConfiguration, subloss_weights::Dict{String, T}) where {T}
    return Zygote.gradient(w -> loss(w, config, subloss_weights), w)[1]
end;


################################################
# Define default types.
################################################
# This is the type to be used for all numerical computations.
DefaultType = Float32
# The volume type is for all quantities that are defined on the optimisations
# points, e. g. CT, dose distribution, binary masks, WED, ...
VolumeType = CuArray{DefaultType, 1, CUDA.Mem.DeviceBuffer}
DijType = CuArray{DefaultType, 2, CUDA.Mem.DeviceBuffer}
OptimConfigType = OptimisationConfiguration{DefaultType, VolumeType, DijType}

# Precompile to speed up code loading.
# precompile(loss, (VolumeType, OptimConfigType))
# precompile(loss!, (VolumeType, OptimConfigType, Dict{String, DefaultType}))
# precompile(loss_gradient, (VolumeType, OptimConfigType))
