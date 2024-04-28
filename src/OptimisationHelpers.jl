using StaticArrays


function get_optimisation_configuration(ct::Juliana.ScalarGrid,
                                        prescriptions::Juliana.Prescriptions,
                                        structures::Dict{String, Juliana.Structure},
                                        Dij,
                                        optimisation_point_indices,
                                        oar_safety_margin::Real=1.0f0,
                                        target_hotspot_percent::Real=103f0)

    # For spot placement.
    coldest_target_name, coldest_dose = Juliana.coldest_target(prescriptions)
    # For normalisation.
    normalisation_target, normalisation_dose = Juliana.hottest_target(prescriptions)

    # Only use the CT, structure masks, ideal dose etc. at the optimisation points.
    normalisation_mask = calculate_normalisation_mask(
        prescriptions,
        structures,
    )

    ideal_dose, importance = Juliana.calculate_ideal_dose_distribution(
        ct,
        prescriptions.target_doses,
        structures,
    );
    # Normalise the ideal dose distribution to make it more realistic.
    # mean_dose = sum(ideal_dose .* normalisation_mask) / sum(normalisation_mask)
    # ideal_dose .*= normalisation_dose / mean_dose

    ideal_dose_prescription_correction!(
        ideal_dose,
        prescriptions,
        structures,
        oar_safety_margin,
    )

    minimum_dose = calculate_minimum_dose_distribution(
        prescriptions,
        structures,
        ideal_dose,
        oar_safety_margin,
    )

    maximum_dose = calculate_maximum_dose_distribution(
        prescriptions,
        structures,
        ideal_dose,
        target_hotspot_percent,
    )

    opt_target_mask = convert.(Float32, Juliana.volume_at_indices(
        normalisation_mask,
        optimisation_point_indices,
    ))

    opt_structures = Dict{String, Vector{Float32}}()
    for (name, structure) in structures
        opt_structure = convert.(Float32, Juliana.volume_at_indices(
            structure.mask,
            optimisation_point_indices,
        ))
        opt_structures[name] = opt_structure
    end

    opt_ct = convert.(Float32, Juliana.volume_at_indices(
        ct.data,
        optimisation_point_indices,
    ))

    opt_ideal_dose = convert.(Float32, Juliana.volume_at_indices(
        ideal_dose,
        optimisation_point_indices,
    ))

    opt_minimum_dose = convert.(Float32, Juliana.volume_at_indices(
        minimum_dose,
        optimisation_point_indices,
    ))

    opt_maximum_dose = convert.(Float32, Juliana.volume_at_indices(
        maximum_dose,
        optimisation_point_indices,
    ))

    opt_importance = convert.(Float32, Juliana.volume_at_indices(
        importance,
        optimisation_point_indices,
    ))

    opt_Dij = Dij

    # Explicit type is mandatory to get type stable code.
    # Otherwise the VolumeType will be CuArray{Float32, 1} instead of
    # CuArray{Float32, 1, CUDA.Mem.DeviceBuffer}.
    # (I don't really understand why, though...)
    OptimisationConfiguration{typeof(normalisation_dose), typeof(opt_ct), typeof(opt_Dij)}(
        normalisation_dose,
        opt_target_mask,
        opt_ct,
        opt_ideal_dose,
        opt_importance,
        opt_minimum_dose,
        opt_maximum_dose,
        opt_Dij,
        opt_Dij',
        opt_structures,
        prescriptions,
    )
end


function is_maximum_constraint(constraint::Juliana.Constraint)
    dvh_constraints = [
        Juliana.constraint_dvh_d,
        Juliana.constraint_dvh_v,
    ]
    
    if constraint.kind == Juliana.constraint_max
        return true
    else
        return (constraint.kind in dvh_constraints) &&
            !isnothing(constraint.volume) &&
            (constraint.volume ≈ 0.02)
    end
end


function ideal_dose_prescription_correction!(ideal_dose_distribution::VolumeType,
                                             prescriptions::Juliana.Prescriptions,
                                             structures::Dict{String, Juliana.Structure},
                                             oar_safety_margin) where {VolumeType}
    for constraint in prescriptions.constraints
        if is_maximum_constraint(constraint)
            threshold_for_structure = constraint.dose - oar_safety_margin
            mask = structures[constraint.structure_name].mask
            ideal_dose_distribution[mask .== 1] .= min.(
                ideal_dose_distribution[mask .== 1],
                threshold_for_structure,
            )
        end
    end
end


function calculate_normalisation_mask(prescriptions, structures; OAR_expansion=0.3)
    hottest_target_name, normalisation_dose = Juliana.hottest_target(prescriptions)

    normalisation_mask = copy(structures[hottest_target_name].mask)
    normalisation_mask_full = copy(normalisation_mask)

    for constraint in prescriptions.constraints
        oar = structures[constraint.structure_name]
        oar_mask = (oar.distanceFromStructure .<= OAR_expansion)
        normalisation_mask = normalisation_mask .&& (.!oar_mask)
    end

    # We do not cut away the OARs if that would remove > 50% of the structure.
    if sum(normalisation_mask) < 0.5 * sum(normalisation_mask_full)
        normalisation_mask = normalisation_mask_full
    end
    
    return normalisation_mask
end


function calculate_minimum_dose_distribution(prescriptions, structures, ideal_dose_distribution, oar_safety_margin, ignore_soft_constraints::Bool=true)
    hottest_target_name, _ = Juliana.hottest_target(prescriptions)
    target_mask = copy(structures[hottest_target_name].mask)
    minimum_dose = zeros(Float32, size(ideal_dose_distribution)...)
    for (name, dose) in prescriptions.target_doses
        target_mask .= max.(target_mask, structures[name].mask)
        minimum_dose .= max.(minimum_dose, structures[name].mask .* dose)
    end

    for constraint in prescriptions.constraints
        if ignore_soft_constraints && (constraint.priority == soft)
            continue
        end

        if !is_maximum_constraint(constraint)
            continue
        end

        oar = structures[constraint.structure_name]
        oar_mask = oar.mask .&& (target_mask .== 1)
        minimum_dose[oar_mask] .= min.(
            constraint.dose .- oar_safety_margin,
            minimum_dose[oar_mask],
        )
    end
    return clamp.(minimum_dose, 0, Inf)
end


function calculate_maximum_dose_distribution(prescriptions, structures, ideal_dose_distribution, hotspot_percent::Real)
    hottest_target_name, hottest_target_dose = Juliana.hottest_target(prescriptions)
    target_mask = copy(structures[hottest_target_name].mask)
    for (name, _) in prescriptions.target_doses
        target_mask .= max.(target_mask, structures[name].mask)
    end

    maximum_dose = copy(ideal_dose_distribution)
    maximum_dose[target_mask] .= ideal_dose_distribution[target_mask] .* (hotspot_percent / 100)

    return maximum_dose
end


function normalise_dose(dose, mask, normalisation_dose)
    return dose .* (normalisation_dose / mean_dose(dose, mask))
end


function get_optimisation_grid(optimisation_points, grid)
    p_min = reshape(minimum(optimisation_points, dims=2), (:,))
    p_max = reshape(maximum(optimisation_points, dims=2), (:,))

    shape = convert.(Int64, ceil.((p_max .- p_min) ./ grid.spacing)) .+ 1

    return Juliana.Grid(
        SVector{3}(grid.spacing),
        SVector{3}(p_min),
        SVector{3}(shape),
    )
end


"""
    get_optimisation_points_from_prescription(
        grid::Juliana.Grid,
        prescriptions::Juliana.Prescriptions,
        structures::Dict{String, Juliana.Structure};
        [interest_distance::Real=2],
        [checkerboard_skip_n_inslice::Integer=2]
        [checkerboard_skip_n_slices::Integer=2]
    ) -> BitArray{3}, Matrix{T:<Real}, Matrix{S:<Integer}

Select which voxels to use as optimisation points.

# Returns
- optimisation_mask<br />
    A binary mask indicating whether the voxel is an optimisation point.
- optimisation_points<br />
    XYZ-coordinates of the optimisation points in the provided grid's
    coordinate system.
- optimisation_point_indices<br />
    3D indices of the optimisation points in the provided grid.
"""
function get_optimisation_points_from_prescription(grid::Juliana.Grid,
                                                   prescriptions::Juliana.Prescriptions,
                                                   structures::Dict{String, Juliana.Structure};
                                                   interest_distance::Real=2,
                                                   checkerboard_skip_n_inslice::Integer=5,
                                                   checkerboard_skip_n_slices::Integer=1,
                                                   margin_skip_n_slices::Integer=1,
                                                   small_oar_voxel_threshold::Integer=2000)
    # Select points that are within a distance of interest_distance [cm]
    # from any target, plus part of any OAR with a hard mean dose constraint
    # to avoid oversparing mean-constraint OARs.
    distance_from_roi = Array{Float32, 3}(undef, grid.size[1], grid.size[2], grid.size[3])
    fill!(distance_from_roi, Inf)

    for (name, value) in prescriptions.target_doses
        distance_from_roi .= min.(
            distance_from_roi,
            structures[name].distanceFromStructure,
        )
    end

    mean_constraint_oar_mask = distance_from_roi .< -Inf
    fill!(mean_constraint_oar_mask, false)
    for constraint in prescriptions.constraints
        if (constraint.priority != soft) && (constraint.kind == constraint_mean)
            mean_constraint_oar_mask .= mean_constraint_oar_mask .|| structures[constraint.structure_name].mask
        end
    end

    optimisation_roi_mask = (distance_from_roi .<= interest_distance) .|| mean_constraint_oar_mask
    optimisation_roi_mask .= optimisation_roi_mask .&& Juliana.build_checker_board_mask(grid, checkerboard_skip_n_inslice, checkerboard_skip_n_slices)

    # Make sure to select all points that are both part of an OAR with a
    # non-soft maximum dose constraint and a target.
    target_mask = distance_from_roi .<= minimum(grid.spacing[1:2])

    oar_mask = similar(optimisation_roi_mask)
    fill!(oar_mask, false)
    for constraint in prescriptions.constraints
        if (constraint.priority != soft) && Juliana.is_maximum_constraint(constraint)
            oar_mask .= oar_mask .|| structures[constraint.structure_name].mask
        end
    end

    overlap_mask = target_mask .&& oar_mask

    # Make sure that the margin voxels of the targets are included,
    # unless they intersect with an OAR that has a non-soft mean dose constraint.
    margin_mask = similar(overlap_mask)
    fill!(margin_mask, false)
    for (name, value) in prescriptions.target_doses
        margin_mask .= margin_mask .|| Juliana.get_margin_mask(structures[name])
    end

    n_points_per_slice = vec(sum(margin_mask, dims=(1, 2)))
    start = findfirst(n_points_per_slice .> 0)
    stop = findlast(n_points_per_slice .> 0)
    for iz in start:stop-1
        if (mod(iz-start, margin_skip_n_slices+1) > 0)
            margin_mask[:, :, iz] .= 0
        end
    end

    margin_mask = margin_mask .&& (.!mean_constraint_oar_mask)

    # Make sure to include the binary masks of OARs that have a non-soft
    # mean dose constraint and only a few voxels.
    small_oars = similar(overlap_mask)
    fill!(small_oars, false)
    for constraint in prescriptions.constraints
        s = structures[constraint.structure_name]
        if (constraint.priority != soft) && (sum(s.mask) <= small_oar_voxel_threshold)
            small_oars .= small_oars .|| structures[constraint.structure_name].mask
        end
    end

    # Build the full mask.
    optimisation_roi_mask = optimisation_roi_mask .|| overlap_mask .|| margin_mask .|| small_oars

    # Convert the optimisation point mask to indices and positions.
    raw_indices = [Tuple(p) for p in findall(optimisation_roi_mask)]

    optimisation_point_indices = Matrix{Int32}(undef, 3, length(raw_indices));
    optimisation_points = Matrix{Float32}(undef, 3, length(raw_indices))
    for (i, index) in enumerate(raw_indices)
        optimisation_point_indices[:, i] .= index
        optimisation_points[:, i] .= Juliana.index_to_xyz(index, grid)
    end

    return optimisation_roi_mask, optimisation_points, optimisation_point_indices
end


function get_default_subloss_weights(prescriptions, structures; oar_weights::Float32=1f0)
    subloss_weights = Dict{String, Float32}(
        "ideal_dose_loss" => 1.f0,
        "maximum_loss" => 1.f0,
        "minimum_loss" => 1.f0,
        "normalisation_variance" => 1.f0,
        # "normalisation_loss" => 10.f0,
        # "normalisation_loss" => 50.f0,
    )

    for constraint in prescriptions.constraints
        if constraint.priority == Juliana.soft
            continue
        end

        # Skip OARs for which there is no optimisation point.
        if (constraint.structure_name ∉ keys(structures)) || sum(structures[constraint.structure_name])== 0
            continue
        end

        if constraint.kind == Juliana.constraint_mean
            subloss_weights["$(constraint.structure_name)_mean_loss"] = oar_weights
        elseif Juliana.is_maximum_constraint(constraint)
            subloss_weights["$(constraint.structure_name)_max_loss"] = oar_weights
        end
    end

    return subloss_weights
end
