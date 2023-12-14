using CUDA
using PyCall


struct PatientInfo
    patient
    plan
    ct_path::String

    PatientInfo(data_dir::String, plan_file::String, patient_ID::String) = begin
        data = pyimport("pyftpp.data")
        plan = pyimport("pyftpp.treatment_plan")

        patient = data.FionaPatient(data_dir, patient_ID, 0)
        plan = plan.TreatmentPlan(patient, plan_file)
        ct_path = "$data_dir/CTs/$(patient_ID)_0.dat"

        new(patient, plan, ct_path)
    end
end


struct PatientData
    ct::ScalarGrid
    structures::Dict{String, Structure}
    prescriptions::Prescriptions
end


function volume_at_indices(volume::AbstractArray{T, 3}, indices::AbstractArray) where {T}
    return [volume[index...] for index in eachrow(indices)]
end


function load_patient_data(data_dir::String, patient_ID::String; series=0)
    ct_path = "$(data_dir)/CTs/$(patient_ID)_$series.dat"
    ct = Juliana.load_ct_dat_file(ct_path);

    structure_names = Juliana.structure_names(data_dir, patient_ID)
    structures_to_load = [name for name in structure_names if name != "AUTO_FULL_BODY"]
    structures = Juliana.load_structures(
        data_dir,
        patient_ID,
        ct.grid,
        which=structures_to_load,
    );

    prescriptions = Juliana.load_prescriptions(data_dir, patient_ID, structures);

    return ct_path, PatientData(
        ct,
        structures,
        prescriptions,
    )
end


function get_optimisation_configuration(ct::ScalarGrid,
                                        prescriptions::Prescriptions,
                                        structures::Dict{String, Structure},
                                        Dij,
                                        optimisation_point_indices;
                                        series=0,
                                        optim_grid_spacing=0.35)

    # For spot placement.
    coldest_target_name, coldest_dose = coldest_target(prescriptions)
    # For normalisation.
    normalisation_target, normalisation_dose = hottest_target(prescriptions)

    optimisation_grid = Juliana.FionaStandalone.calculate_optimisation_grid(
        ct.grid,
        prescriptions,
        structures,
        optim_grid_spacing,
    )

    # Only use the CT, structure masks, ideal dose etc. at the optimisation points.
    normalisation_mask = calculate_normalisation_mask(
        prescriptions,
        structures,
    )

    ideal_dose, importance = calculate_ideal_dose_distribution(
        ct,
        prescriptions.target_doses,
        structures,
    );
    # Normalise the ideal dose distribution to make it more realistic.
    mean_dose = sum(ideal_dose .* normalisation_mask) / sum(normalisation_mask)
    ideal_dose .*= normalisation_dose / mean_dose

    ideal_dose_prescription_correction!(
        ideal_dose,
        prescriptions,
        structures,
    )

    minimum_dose = calculate_minimum_dose_distribution(
        prescriptions,
        structures,
        ideal_dose,
    )

    opt_target_mask = cu(convert.(Float32, volume_at_indices(
        normalisation_mask,
        optimisation_point_indices,
    )))

    opt_structures = Dict{String, CuArray{Float32, 1}}()
    for (name, structure) in structures
        opt_structure = convert.(Float32, volume_at_indices(
            structure.mask,
            optimisation_point_indices,
        ))
        opt_structures[name] = cu(opt_structure)
    end

    opt_ct = cu(convert.(Float32, volume_at_indices(
        ct.data,
        optimisation_point_indices,
    )))

    opt_ideal_dose = cu(convert.(Float32, volume_at_indices(
        ideal_dose,
        optimisation_point_indices,
    )))

    opt_minimum_dose = cu(convert.(Float32, volume_at_indices(
        minimum_dose,
        optimisation_point_indices,
    )))

    opt_importance = cu(convert.(Float32, volume_at_indices(
        importance,
        optimisation_point_indices,
    )))

    OptimisationConfiguration(
        normalisation_dose,
        opt_target_mask,
        opt_ct,
        opt_ideal_dose,
        opt_importance,
        opt_minimum_dose,
        cu(collect(Dij)),
        opt_structures,
        prescriptions,
    )
end


# function getOptimisationConfiguration(type, info::PatientInfo, angles, nozzle_extraction, spot_spacing, optimisation_grid_resolution, dose_grid_resolution, output_dir, debugging::Bool)
#     patient = info.patient
#     plan = info.plan
#     ct_path = info.ct_path

#     spot_placement_target = patient.prescriptions.coldest_target()[1]
#     normalisation_target, normalisation_dose = patient.prescriptions.hottest_target()
#     normalisation_dose = convert(type, normalisation_dose)

#     dij_config = DijCalculationConfig(
#         patient,
#         spot_placement_target,
#         normalisation_target,
#         normalisation_dose,
#         angles,
#         nozzle_extraction,
#         spot_spacing,
#         dose_grid_resolution,
#         optimisation_grid_resolution,
#     )
#     optimisation = calculate_Dij_Fiona(dij_config, ct_path, output_dir, debugging)

#     # Dij = cu(sparse(optimisation.Dij.toarray()))
#     Dij = cu(optimisation.Dij.toarray())
#     target_mask = get_structure(optimisation, optimisation.normalisation_structure_name)
#     target_dose = convert(type, optimisation.normalisation_dose)

#     MaskType = typeof(target_mask)
#     structures = Dict{String, MaskType}()
#     for name in patient.structures.structure_names
#         if name == "AUTO_FULL_BODY"
#             continue
#         end
#         structures[name] = get_structure(optimisation, name)
#     end

#     optimisation.initialise_optimisation()
#     w = cu(reshape(optimisation.w, :))

#     targets = get_targets(
#         patient,
#         optimisation,
#         normalisation_target,
#         normalisation_dose,
#     )
#     ideal_dose, importance = calculate_flattened_ideal_dose_distribution(info, optimisation);
#     ideal_dose = cu(ideal_dose)

#     # TODO: Only load data once, preferably in Julia.
#     # Correct the ideal dose distribution for the OAR Dmax prescriptions.
#     data_dir = info.patient._data_dir
#     patient_ID = info.patient.patient_ID
#     ct = Juliana.load_ct_dat_file("$(data_dir)/CTs/$(patient_ID)_0.dat");

#     structure_names = Juliana.structure_names(data_dir, patient_ID)
#     structures_to_load = [name for name in structure_names if name != "AUTO_FULL_BODY"]
#     structures_3d = Juliana.load_structures(data_dir, patient_ID, ct.grid, which=structures_to_load);

#     prescriptions = Juliana.load_prescriptions(data_dir, patient_ID, structures_3d);

#     ideal_dose_prescription_correction!(ideal_dose, prescriptions, structures)

#     return OptimisationConfiguration(
#         target_dose,
#         target_mask,
#         targets,
#         cu(Array{type}(reshape(patient.ct.data, :))),
#         cu(Array{type}(reshape(ideal_dose, :))),
#         cu(Array{type}(reshape(importance, :))),
#         Dij,
#         structures,
#         prescriptions,
#     ), w, optimisation
# end


# function get_structure(optimisation, name)
#     mask = optimisation.structure_mask(name)
#     return cu(reshape(mask, :))
# end;


# function calculate_Dij_Fiona(config::DijCalculationConfig, ct_path::String, output_dir::String, debugging)
#     optim = pyimport("pyftpp.optimise.optimisation")
#     optimisation = optim.TensorflowSpotOptimisation(
#         config.patient,
#         config.spot_placement_target,
#         config.normalisation_target,
#         config.normalisation_dose,
#         ct_path,
#         output_dir
#     )
#     optimisation.build_Dij(
#         config.angles,
#         config.nozzle_extraction,
#         config.spot_spacing,
#         dose_grid_resolution=config.dose_grid_resolution,
#         optimisation_grid_resolution=config.optimisation_grid_resolution,
#         # TODO: Change in production!
#         debugging=debugging,
#     )
#     optimisation.load_Dij()
    
#     return optimisation
# end


# function get_targets(patient, optimisation, normalisation_target, normalisation_dose)
#     n_targets = length(patient.prescriptions.target_doses)

#     @assert n_targets > 0
#     hottest_mask = get_structure(optimisation, normalisation_target)
#     hottest_mask_bool = collect(.!iszero.(hottest_mask))
#     targets = [(normalisation_target, hottest_mask, normalisation_dose)]

#     if n_targets > 1
#         for (name, dose) in patient.prescriptions.target_doses
#             if name == normalisation_target
#                 continue
#             end
#             mask = get_structure(optimisation, name)
#             mask_bool = collect(.!iszero.(mask))
#             mask = mask_bool .&& .!hottest_mask_bool
#             mask = cu(convert(AbstractArray{Float32}, mask))
#             dose = convert(Float32, dose)
#             push!(targets, (name, mask, dose))
#         end
#     end
    
#     return targets
# end


function is_maximum_constraint(constraint::Constraint)
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
                                             prescriptions::Prescriptions,
                                             structures::Dict{String, Structure}) where {VolumeType}
    for constraint in prescriptions.constraints
        if is_maximum_constraint(constraint)
            threshold_for_structure = constraint.dose
            mask = structures[constraint.structure_name].mask
            ideal_dose_distribution[mask .== 1] .= min.(
                ideal_dose_distribution[mask .== 1],
                threshold_for_structure,
            )
        end
    end
end


function calculate_normalisation_mask(prescriptions, structures; OAR_expansion=0.3)
    hottest_target_name, normalisation_dose = hottest_target(prescriptions)

    normalisation_mask = structures[hottest_target_name].mask
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


function calculate_minimum_dose_distribution(prescriptions, structures, ideal_dose_distribution; OAR_expansion=0.3, safety_margin=0.5)
    hottest_target_name, _ = hottest_target(prescriptions)
    target_mask = structures[hottest_target_name].mask
    for (name, _) in prescriptions.target_doses
        target_mask .= max.(target_mask, structures[name].mask)
    end

    minimum_dose = (ideal_dose_distribution .- safety_margin) .* target_mask

    for constraint in prescriptions.constraints
        oar = structures[constraint.structure_name]
        oar_mask = (oar.distanceFromStructure .<= OAR_expansion) .&& (target_mask .== 1)
        minimum_dose[oar_mask] .= min.(
            constraint.dose .- safety_margin,
            minimum_dose[oar_mask],
        )
    end
    return clamp.(minimum_dose, 0, Inf)
end
