using CUDA

# Units: Gy always!!!


function mean_dose(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}) where {T, N}
    # @assert sum(mask) > zero(T)
    sum(dose .* mask) / sum(mask)
end


function variance_dose(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}) where {T, N}
    μ = mean_dose(dose, mask)
    structure_dose = dose[mask .== one(T)]
    n_voxels_in_structure = size(structure_dose, 1)
    @assert n_voxels_in_structure > zero(T)
    return sum((structure_dose .- μ).^2) ./ n_voxels_in_structure
end


function variance_dose_gradient(dose, mask)
    N_target = sum(mask)
    return 2 / N_target .* (dose .- mean_dose(dose, mask)) .* mask
end


function variance_dose_gradient(dose::CuArray, mask::CuArray)
    return cu(variance_dose_gradient(collect(dose), collect(mask)))
end


function hotspot_loss(dose::AbstractArray{T, N}, threshold::T) where {T, N}
    return sum(max.(dose .- threshold, zero(T)))
end


# function hotspot_loss(dose::CuArray{T, N}, threshold::T) where {T, N}
#     # collect(): Transfer to the CPU to avoid excessive memory allocations.
#     return sum(max.(collect(dose) .- threshold, zero(T)))
# end


function hotspot_loss(dose::AbstractArray{T, N}, threshold::AbstractArray{T, N}) where {T, N}
    # Voxel-wise maximum value.
    return sum(max.(dose .- threshold, zero(T)))
end


# function hotspot_loss(dose::CuArray{T, N}, threshold::AbstractArray{T, N}) where {T, N}
#     # collect(): Transfer to the CPU to avoid excessive memory allocations.
#     # Voxel-wise maximum value.
#     return sum(max.(collect(dose) .- collect(threshold), zero(T)))
# end


"""Loss function that ensures Dmax <= threshold for the given structure."""
function Dmax_loss(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}, threshold::T) where {T, N}
    structure_dose = dose .* mask
    excess_dose = max.(zero(T), structure_dose .- threshold)
    loss = sum(excess_dose)
end


function Dmax_loss_gradient(dose::AbstractArray{T, N}, mask::AbstractArray{T, N}, threshold::T) where {T, N}
    return convert.(T, dose .* mask .>= threshold)
end


function Dmax_loss_gradient(dose::CuArray{T, N}, mask::CuArray{T, N}, threshold::T) where {T, N}
    return cu(Dmax_loss_gradient(collect(dose), collect(mask), threshold))
end


function ideal_dose_loss(dose::AbstractArray{T, N}, importance, ideal_dose) where {T, N}
    return dot(
        importance,
        (ideal_dose .- dose).^2,
    ) / sum(importance) / convert(T, 25)
end


function ideal_dose_loss_gradient(dose::AbstractArray{T, N}, importance, ideal_dose) where {T, N}
    return 2 .* importance .* (dose .- ideal_dose) ./ sum(importance) ./ convert(T, 25)
end


function ideal_dose_loss_gradient(dose::CuArray{T, N}, importance::CuArray{T, N}, ideal_dose::CuArray{T, N}) where {T, N}
    return cu(ideal_dose_loss_gradient(
        collect(dose),
        collect(importance),
        collect(ideal_dose),
    ))
end


function minimum_loss(dose::AbstractArray{T, N}, minimum_dose::AbstractArray{T, N}) where {T, N}
    return sum(max.(minimum_dose .- dose, 0)) / 100
end


function minimum_loss_gradient(dose::AbstractArray{T, N}, minimum_dose::AbstractArray{T, N}) where {T, N}
    grad = similar(dose)
    for i in 1:length(dose)
        grad[i] = dose[i] < minimum_dose[i] ? convert(T, -1 / 100) : convert(T, 0)
    end
    return grad
end


function minimum_loss_gradient(dose::CuArray{T, N}, minimum_dose::CuArray{T, N}) where {T, N}
    return cu(minimum_loss_gradient(collect(dose), collect(minimum_dose)))
end


function maximum_distribution_loss(dose::AbstractArray{T, N}, maximum_dose) where {T, N}
    return sum(max.(dose .- maximum_dose, 0)) ./ 100
end


function maximum_distribution_loss_gradient(dose::AbstractArray{T, N}, maximum_dose) where {T, N}
    return convert.(Float32, dose .> maximum_dose) ./ 100
end


function maximum_distribution_loss_gradient(dose::CuArray{T, N}, maximum_dose) where {T, N}
    return cu(maximum_distribution_loss_gradient(collect(dose), collect(maximum_dose)))
end


function maximum_loss(dose::AbstractArray{T, N}, normalisation_dose) where {T, N}
    return hotspot_loss(dose, convert(T, 1.05) .* normalisation_dose) * convert(T, 10.)
end


function maximum_loss_gradient(dose::AbstractArray{T, N}, normalisation_dose) where {T, N}
    grad = similar(dose)
    for i in 1:length(dose)
        grad[i] = dose[i] > convert(T, 1.05) * normalisation_dose ? convert(T, 10) : convert(T, 0)
    end
    return grad
end


function maximum_loss_gradient(dose::CuArray{T, N}, normalisation_dose) where {T, N}
    return cu(maximum_loss_gradient(collect(dose), normalisation_dose))
end


function oar_mean_loss(dose::AbstractArray{T, N}, mask, threshold) where {T, N}
    mean = mean_dose(dose, mask)
    return max(mean - threshold, zero(T))
end


function oar_mean_loss_gradient(dose::AbstractArray{T, N}, mask, threshold) where {T, N}
    mean = mean_dose(dose, mask)
    if mean >= threshold
        return mask ./ sum(mask)
    else
        return zeros(T, length(dose))
    end
end


function oar_mean_loss_gradient(dose::CuArray{T, N}, mask, threshold) where {T, N}
    return cu(oar_mean_loss_gradient(collect(dose), collect(mask), threshold))
end


function normalisation_loss(dose, normalisation_mask, normalisation_dose)
    return (mean_dose(dose, normalisation_mask) - normalisation_dose)^2
end


function normalisation_loss_gradient(dose, normalisation_mask, normalisation_dose)
    return (2 / sum(normalisation_mask) * (mean_dose(dose, normalisation_mask) - normalisation_dose)) .* normalisation_mask
end


function build_loss_parts(dose::AbstractArray{T, N},
                          config::OptimisationConfiguration,
                          safety_margin::T) where {T, N}
    loss_parts = Dict{String, T}()

    # Normalise the dose distribution to make sure that a fulfilled prescription is
    # really fulfilled even after renormalisation.
    # dose = dose .* (config.normalisationDose / mean_dose(dose, config.normalisationStructureMask))

    # L_ideal
    loss_parts["ideal_dose_loss"] = ideal_dose_loss(dose, config.importance, config.idealDose)

    # L_hom
    loss_parts["normalisation_variance"] = variance_dose(
        dose,
        config.normalisationStructureMask,
    )

    # L_norm
    loss_parts["normalisation_loss"] = normalisation_loss(
        dose,
        config.normalisationStructureMask,
        config.normalisationDose,
    )

    # L_min
    loss_parts["minimum_loss"] = minimum_loss(dose, config.minimumDose)

    # L_max
    loss_parts["maximum_loss"] = maximum_loss(dose, config.normalisationDose)
    loss_parts["maximum_distribution_loss"] = maximum_distribution_loss(
        dose,
        config.maximumDose,
    )

    # OAR loss parts.
    for constraint in config.prescriptions.constraints
        if constraint.priority == soft
            continue
        end

        if constraint.kind == constraint_mean
            mask = config.structures[constraint.structure_name]
            threshold = constraint.dose - safety_margin
            loss_parts["$(constraint.structure_name)_mean_loss"] = oar_mean_loss(
                dose,
                mask,
                threshold,
            )
        elseif is_maximum_constraint(constraint)
            mask = config.structures[constraint.structure_name]
            threshold = constraint.dose - safety_margin
            loss_parts["$(constraint.structure_name)_max_loss"] = Dmax_loss(
                dose,
                mask,
                threshold,
            )
        end
    end
    
    return loss_parts
end


function build_loss_gradient_parts(dose::AbstractArray{T, N},
                                   config::OptimisationConfiguration,
                                   safety_margin::T) where {T, N}
    gradient_parts = Dict{String, typeof(dose)}()

    gradient_parts["ideal_dose_loss"] = ideal_dose_loss_gradient(dose, config.importance, config.idealDose)
    gradient_parts["normalisation_variance"] = variance_dose_gradient(
        dose,
        config.normalisationStructureMask,
    )
    gradient_parts["minimum_loss"] = minimum_loss_gradient(dose, config.minimumDose)
    gradient_parts["maximum_loss"] = maximum_loss_gradient(dose, config.normalisationDose)
    gradient_parts["maximum_distribution_loss"] = maximum_distribution_loss_gradient(
        dose,
        config.maximumDose,
    )
    gradient_parts["normalisation_loss"] = normalisation_loss_gradient(dose, config.normalisationStructureMask, config.normalisationDose)

    # OAR loss parts.
    for constraint in config.prescriptions.constraints
        if constraint.priority == soft
            continue
        end

        if constraint.kind == constraint_mean
            mask = config.structures[constraint.structure_name]
            threshold = constraint.dose - safety_margin
            gradient_parts["$(constraint.structure_name)_mean_loss"] = oar_mean_loss_gradient(
                dose,
                mask,
                threshold,
            )
        elseif is_maximum_constraint(constraint)
            mask = config.structures[constraint.structure_name]
            threshold = constraint.dose - safety_margin
            gradient_parts["$(constraint.structure_name)_max_loss"] = Dmax_loss_gradient(
                dose,
                mask,
                threshold,
            )
        end
    end

    return gradient_parts
end





# Version with logging.
function dose_loss!(dose::AbstractArray{T, N},
                    config::OptimisationConfiguration,
                    loss_parts::Dict{String, T},
                    subloss_weights::Dict{String, T},
                    oar_safety_margin::T) where {T, N}
    partial_losses = build_loss_parts(dose, config, oar_safety_margin)

    for (key, value) in partial_losses
        loss_parts[key] = value
    end

    loss = zero(T)
    for (name, value) in loss_parts
        if name in keys(subloss_weights)
            loss += subloss_weights[name] * value
        end
    end

    return loss
end

# Version without logging.
function dose_loss(dose::AbstractArray{T, N},
                   config::OptimisationConfiguration,
                   subloss_weights::Dict{String, T},
                   oar_safety_margin::T) where {T, N}
    partial_losses = build_loss_parts(dose, config, oar_safety_margin)

    loss = zero(T)
    for (name, value) in partial_losses
        if name in keys(subloss_weights)
            loss += subloss_weights[name] * value
        end
    end

    return loss
end


function dose_loss_gradient(dose,
                            config::OptimisationConfiguration,
                            subloss_weights::Dict{String, T},
                            oar_safety_margin::T) where {T}
    partial_gradients = build_loss_gradient_parts(dose, config, oar_safety_margin)

    gradient = zeros(T, length(dose))
    for (name, value) in partial_gradients
        if name in keys(subloss_weights)
            gradient .+= collect(subloss_weights[name]) .* value
        end
    end

    return gradient
end


function dose_loss_gradient(dose::CuArray,
                            config::OptimisationConfiguration,
                            subloss_weights::Dict{String, T},
                            oar_safety_margin::T) where {T}
    partial_gradients = build_loss_gradient_parts(dose, config, oar_safety_margin)

    gradient = CUDA.zeros(T, length(dose))
    for (name, value) in partial_gradients
        if name in keys(subloss_weights)
            gradient .+= subloss_weights[name] .* value
        end
    end

    return gradient
end


function loss!(w,
               config::OptimisationConfiguration,
               loss_parts,
               subloss_weights::Dict{String, T},
               oar_safety_margin::T) where {T}
    dose = reproducible_sparse_mv(config.Dij, w)
    return dose_loss!(dose, config, loss_parts, subloss_weights, oar_safety_margin)
end


function loss(w,
              config::OptimisationConfiguration,
              subloss_weights::Dict{String, T},
              oar_safety_margin::T) where {T}
    dose = reproducible_sparse_mv(config.Dij, w)
    return dose_loss(dose, config, subloss_weights, oar_safety_margin)
end


function loss_gradient(w,
                       config::OptimisationConfiguration,
                       subloss_weights::Dict{String, T},
                       oar_safety_margin::T) where {T}
    dose = reproducible_sparse_mv(config.Dij, w)
    grad = dose_loss_gradient(dose, config, subloss_weights, oar_safety_margin)
    return reproducible_sparse_mv(config.Dij_T, grad)
end


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
