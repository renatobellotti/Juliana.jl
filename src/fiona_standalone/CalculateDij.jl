import ..Grid
import ..Prescriptions
import ..Structure


function calculate_optimisation_grid(ct_grid::Grid,
                                     prescriptions::Prescriptions,
                                     structures::Dict{String, Structure},
                                     optim_grid_spacing::Real)
    # This function is a direct translation of the pyftpp code,
    # which in turn is a direct translation of the Fiona standalone code.
    EXTENSION = 3.
    
    # Get bounding box of all relevant structures
    relevant_structure_names = String[]

    for (name, dose) in prescriptions.target_doses
        push!(relevant_structure_names, name)
    end

    for constraint in prescriptions.constraints
        push!(relevant_structure_names, constraint.structure_name)
    end

    relevant_structure_names

    p_min = [Inf, Inf, Inf]
    p_max = [-Inf, -Inf, -Inf]

    for name in relevant_structure_names
        points = structures[name].points
        p_min .= min.(p_min, reshape(maximum(points, dims=1), :))
        p_max .= max.(p_max, reshape(maximum(points, dims=1), :))
    end

    p_min .-= EXTENSION
    p_max .+= EXTENSION

    # Round the origin in order to make them lie on CT grid nodes.
    @. p_min = floor((p_min - ct_grid.origin) / ct_grid.spacing) * ct_grid.spacing + ct_grid.origin
    # Needed in Fiona because of the way in which voxels are assigned to
    # structures? I don't understand it fully...
    p_min[3] += 0.1
    
    optim_grid_size = Base.convert.(Int64, ceil.((p_max .- p_min) ./ optim_grid_spacing))
    spacing = [optim_grid_spacing, optim_grid_spacing, optim_grid_spacing]
    
    return Grid(
        Base.convert.(Float32, spacing),
        Base.convert.(Float32, p_min),
        optim_grid_size,
    )
end


function calculate_Dij(working_dir::String,
                       ct_path::String,
                       target_dose::Real,
                       target_structure::Structure, # Used for the spot placement.
                       fiona_standalone_bin_path::String,
                       fiona_jar_path::String,
                       optimisation_grid::Grid,
                       gantry_angles::Vector,
                       couch_angles::Vector,
                       nozzle_extractions::Vector;
                       log_wed::Bool=false,
                       debugging::Bool = false,
                       optimization_points = nothing)
    if !debugging
        main_config = MainConfig(
            ct_path,
            working_dir,
            target_dose,
            fiona_standalone_bin_path,
        );

        optimisation_config = OptimizationSettings(
            target_dose,
            0.9 * target_dose,
            StructureDefinition(target_structure, 0),
            Array{StructureConstraints, 1}(undef, 0),
            convert(OptimizationGrid, optimisation_grid),
        );

        spot_placement_config = SpotPlacementConfig(
            gantry_angles,
            couch_angles,
            nozzle_extractions,
            "AUTO",
            target_structure,
        );

        output = run_spot_placement(
            fiona_jar_path,
            working_dir,
            false,
            false,
            main_config,
            optimisation_config,
            spot_placement_config,
        )
        plan = read_plan_file("$(working_dir)/result_plan.json");

        if !isnothing(optimization_points)
            optimizationPointsPath = "$(working_dir)/custom_optimization_points.txt"
            write_optimisation_points(
                optimizationPointsPath,
                optimization_points,
                optimisation_grid,
            )
        else
            optimizationPointsPath = nothing
        end

        output = run_optimization(
            fiona_jar_path,
            working_dir,
            true,
            log_wed,
            main_config,
            optimisation_config,
            spot_placement_config,
            plan,
            optimizationPointsPath=optimizationPointsPath,
        )
    end

    optimisation_points = load_optimisation_points("$working_dir/optimization_points.json")
    Dij = load_Dij("$working_dir/dij_matrix.dat", size(optimisation_points, 1))

    return Dij, optimisation_points
end
