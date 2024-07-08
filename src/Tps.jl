using SparseArrays


function load_juliana_tps(fiona_standalone_bin_path="$(@__DIR__)/../data", nozzle_extraction=0)
    huToSpPathFile = "$(fiona_standalone_bin_path)/huToSp.json"
    convert_to_sp = hu_to_sp_factory(huToSpPathFile)

    depth_dose_curves, sigma_mcs_curves, phase_space_no_preabsorber, phase_space_with_preabsorber = load_machine_parameters(fiona_standalone_bin_path, nozzle_extraction)

    return JulianaTps(
        convert_to_sp,
        depth_dose_curves,
        sigma_mcs_curves,
        phase_space_no_preabsorber,
        phase_space_with_preabsorber,
    )
end


function to_gpu(tps::JulianaTps)
    d_depth_dose_curves = cu(tps.d_depth_dose_curves)
    d_sigma_mcs_curves = cu(tps.d_sigma_mcs_curves)
    d_phase_space_no_preabsorber = cu(tps.d_phase_space_no_preabsorber);
    d_phase_space_with_preabsorber = cu(tps.d_phase_space_with_preabsorber);

    return JulianaTps(
        tps.convert_to_sp,
        d_depth_dose_curves,
        d_sigma_mcs_curves,
        d_phase_space_no_preabsorber,
        d_phase_space_with_preabsorber,
    )
end


function place_spots(tps::JulianaTps,
                     ct::ScalarGrid,
                     field_definitions,
                     # Structures that defines the iso-center of the fields.
                     field_structures,
                     preabsorber::String="AUTO")
    @assert preabsorber in ("IN", "OUT", "AUTO")
    # Dummy values: We don't optimise.
    target_dose = 1f0

    densities = Juliana.ScalarGrid(
        tps.convert_to_sp.(ct.data),
        ct.grid,
    )

    fields = Juliana.place_spots_on_grid(
        densities,
        field_definitions,
        field_structures,
        tps.d_depth_dose_curves,
    )
    plan = Juliana.TreatmentPlan(fields, target_dose)
    plan = Juliana.initialise_spot_weights(plan, tps)

    return plan
end


# Calculate dose using Juliana.
function calculate_Dij_matrix(tps::JulianaTps, ct, field::Juliana.FieldDefinition, optimisation_points)
    # Convert HUs to relative stopping power.
    densities = convert.(Float32, tps.convert_to_sp.(ct.data))

    # Convert gantry & couch angles to direction vectors.
    gantry_angle = convert(Float32, field.gantryAngle)
    couch_angle = convert(Float32, field.couchAngle)
    nozzle_extraction = convert(Float32, field.nozzleExtraction)
    directions = angles_to_direction([gantry_angle], [couch_angle])

    # Calculate WED.
    d_wed = cu(Vector{Float32}(undef, size(optimisation_points, 2)))
    d_densities = cu(densities)
    d_grid = cu(ct.grid)
    d_optim_points = cu(optimisation_points)
    d_directions = cu(directions)

    n_points = size(optimisation_points, 2)

    event = calculate_wed_simple(
        d_wed,
        d_densities,
        d_grid,
        d_optim_points,
        d_directions,
        ndrange=(n_points, 1),
    )
    wait(event)

    # Copy spots to the GPU.
    df = DataFrame(field.spots);
    d_spots_t = cu([s.t for s in field.spots])
    d_spots_u = cu([s.u for s in field.spots])
    d_spots_energy = cu([s.energykeV for s in field.spots])
    d_spots_absorbers = cu([s.numberOfAbsorbers for s in field.spots])
    d_field_center = cu(convert.(Float32, [
        field.fieldCenter["x"],
        field.fieldCenter["y"],
        field.fieldCenter["z"],
    ]))
    d_optimisation_points = cu(collect(optimisation_points'))
    d_tps = Juliana.to_gpu(tps)

    # Allocate memory.
    N = size(optimisation_points, 2)
    M = length(field.spots)
    Dij_juliana = cu(zeros(N, M))

    # Calculate the Dij matrix.
    runDijKernel(
        Dij_juliana,
        d_wed,
        d_field_center,
        gantry_angle,
        couch_angle,
        d_spots_t,
        d_spots_u,
        d_spots_energy,
        d_spots_absorbers,
        cu(optimisation_points'),
        d_tps.d_depth_dose_curves,
        d_tps.d_sigma_mcs_curves,
        d_tps.d_phase_space_no_preabsorber,
        d_tps.d_phase_space_with_preabsorber,
    )

    return sparse(collect(Dij_juliana))
end


function calculate_Dij_matrix(tps::JulianaTps, ct, plan, optimisation_points)
    return hcat([calculate_Dij_matrix(tps, ct, f, optimisation_points) for f in plan.fields]...)
end
