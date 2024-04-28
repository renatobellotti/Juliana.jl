using CUDA
using DataFrames
using SparseArrays


function calculate_Dij(field,
                       wed,
                       optimisation_points,
                       depth_dose_curves,
                       sigma_mcs_curves,
                       phase_space_no_preabsorber,
                       phase_space_with_preabsorber)
    # Move WED to the GPU.
    d_wed = cu(wed)
    
    # Move spot parameters to the GPU.
    df = DataFrame(field.spots);
    d_spots_t = cu(df.t)
    d_spots_u = cu(df.u)
    d_spots_energy = cu(df.energykeV)
    d_spots_absorbers = cu(df.numberOfAbsorbers)
    d_field_center = cu(convert.(Float32, [
        field.fieldCenter["x"],
        field.fieldCenter["y"],
        field.fieldCenter["z"],
    ]))

    # Allocate memory on the GPU.
    N = size(optimisation_points, 1)
    M = length(d_spots_u)
    d_Dij_juliana = cu(ones(N, M)*NaN32);
    d_optimisation_points = cu(optimisation_points)
    
    d_depth_dose_curves = cu(depth_dose_curves)
    d_sigma_mcs_curves = cu(sigma_mcs_curves)
    d_phase_space_no_preabsorber = cu(phase_space_no_preabsorber)
    d_phase_space_with_preabsorber = cu(phase_space_with_preabsorber);
    
    # Run the kernel.
    Juliana.runDijKernel(
        d_Dij_juliana,
        d_wed,
        d_field_center,
        convert(Float32, field.gantryAngle * π / 180f0),
        convert(Float32, field.couchAngle * π / 180f0),
        d_spots_t,
        d_spots_u,
        d_spots_energy,
        d_spots_absorbers,
        d_optimisation_points,
        d_depth_dose_curves,
        d_sigma_mcs_curves,
        d_phase_space_no_preabsorber,
        d_phase_space_with_preabsorber,
    );
    return d_Dij_juliana
end


function calculate_Dij(fields, wed, optimisation_points, fiona_standalone_bin_path)
    Dij_all_fields = []
    for (i, field) in enumerate(fields)
        field_wed = wed[:, i]

        depth_dose_curves, sigma_mcs_curves, phase_space_no_preabsorber, phase_space_with_preabsorber = Juliana.load_machine_parameters(
            fiona_standalone_bin_path,
            field.nozzleExtraction,
        )

        d_Dij = calculate_Dij(
            field,
            field_wed,
            optimisation_points,
            depth_dose_curves,
            sigma_mcs_curves,
            phase_space_no_preabsorber,
            phase_space_with_preabsorber,
        );
        push!(Dij_all_fields, collect(d_Dij))
    end

    return sparse(hcat(Dij_all_fields...))
end


# needed for the lookup tables
function calculate_Dij(lookup_table_directory,
                       ct,
                       beam_arrangement,
                       optimisation_points,
                       fields)
    # Calculate relative stopping powers.
    huToSpPathFile = "$(lookup_table_directory)/huToSp.json"
    convert_to_sp = Juliana.hu_to_sp_factory(huToSpPathFile)
    @time densities = Juliana.ScalarGrid(
        convert.(Float32, convert_to_sp.(ct.data)),
        ct.grid,
    );
    
    # Calculate the WED
    @time wed = Juliana.calculate_wed(
        densities,
        beam_arrangement.gantry_angles,
        beam_arrangement.couch_angles,
        optimisation_points',
    );

    # Calculate the influence matrix.
    return @time Juliana.calculate_Dij(
        fields,
        wed,
        optimisation_points,
        lookup_table_directory,
        
    )
end
