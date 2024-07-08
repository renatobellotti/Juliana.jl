using StaticArrays


function position_t_u_energy_absorber(spot::Juliana.Spot)
    return (spot.t, spot.u, spot.energykeV, spot.numberOfAbsorbers)
end


function build_stu_grid(field, placement_structure; placement_margin=0.5, lateral_spacing=0.4, depth_spacing=0.2)
    # Figure out which points are relevant for spot placement.
    # These are called placement_points.
    spot_placement_mask = placement_structure.distanceFromStructure .<= placement_margin
    placement_points, placement_indices = Juliana.mask_to_points_and_indices(
        placement_structure.grid,
        spot_placement_mask,
    )
    
    # Convert the placement point to STU coordinates to determine the extent of
    # the spot placement grid.
    iso_center = (
        field.fieldCenter["x"],
        field.fieldCenter["y"],
        field.fieldCenter["z"],
    )

    placement_points_stu = similar(placement_points)
    for (j, p) in enumerate(eachcol(placement_points))
        placement_points_stu[:, j] .= Juliana.XyzToStu(
            p[1],
            p[2],
            p[3],
            iso_center,
            field.gantryAngle,
            field.couchAngle,
        )
    end

    # Determine the STU coordinate extent.
    s_min = minimum(placement_points_stu[1, :])
    t_min = minimum(placement_points_stu[2, :])
    u_min = minimum(placement_points_stu[3, :])

    s_max = maximum(placement_points_stu[1, :])
    t_max = maximum(placement_points_stu[2, :])
    u_max = maximum(placement_points_stu[3, :])

    stu_min = convert.(Float32, (s_min, t_min, u_min))
    stu_max = convert.(Float32, (s_max, t_max, u_max))
    stu_spacing = convert.(Float32, (depth_spacing, lateral_spacing, lateral_spacing))

    # Round to grid spacing to ensure that the margin is actually kept.
    stu_min = floor.(stu_min ./ stu_spacing) .* stu_spacing
    stu_max =  ceil.(stu_max ./ stu_spacing) .* stu_spacing

    max_index = convert.(Integer, round.(stu_max ./ stu_spacing))
    min_index = convert.(Integer, round.(stu_min ./ stu_spacing))
    stu_size = max_index .- min_index

    return Juliana.Grid(
        SVector{3}(collect(stu_spacing)),
        SVector{3}(collect(stu_min)),
        SVector{3}(collect(stu_size)),
    )
end


"""
The WED cube is an array of CT shape. Its values must have been calculated at least for all points on the STU grid.
"""
function place_spots_on_grid(field::Juliana.FieldDefinition,
                             placement_structure::Juliana.Structure,
                             wed_cube::AbstractArray{T, 3},
                             depth_dose_curves::Juliana.LookUpTable;
                             placement_margin::Real=0.5,
                             lateral_spacing::Real=0.4,
                             depth_spacing::Real=0.2) where {T<:Real}
    stu_grid = build_stu_grid(
        field,
        placement_structure;
        placement_margin=placement_margin,
        lateral_spacing=lateral_spacing,
        depth_spacing=depth_spacing,
    )
    s_min, t_min, u_min = stu_grid.origin
    s_max, t_max, u_max = stu_grid.origin .+ stu_grid.spacing .* stu_grid.size

    spot_data_type = typeof(t_min)

    # Actually place the spots.
    ranges = [Juliana.get_range(depth_dose_curves, E) for E in depth_dose_curves.energies]
    iso_center = getindex.(Ref(field.fieldCenter), ("x", "y", "z"))

    current_ID = 0
    spots = Vector{Juliana.Spot{spot_data_type}}(undef, 0)
    for t in t_min:lateral_spacing:t_max
        for u in u_min:lateral_spacing:u_max
            for s in s_min:depth_spacing:s_max
                p = Juliana.StuToXyz(
                    s, t, u,
                    iso_center,
                    field.gantryAngle,
                    field.couchAngle,
                )
                ip = Juliana.xyz_to_index(p, placement_structure.grid)
                if !Juliana.index_is_valid(ip[1], ip[2], ip[3], size(wed_cube))
                    # Spot would be placed outside of the CT grid.
                    continue
                end
                depth = wed_cube[ip[1], ip[2], ip[3]]

                # Check whether we need a preabsorber.
                preabsorber = 0
                current_ranges = ranges
                if depth <= minimum(ranges)
                    preabsorber = 1
                    current_ranges = ranges .- Juliana.PREABSORBER_WED
                end

                # Get the closest energy.
                energy_index = argmin(abs.(current_ranges .- depth))
                E = depth_dose_curves.energies[energy_index]

                if placement_structure.distanceFromStructure[ip...] <= placement_margin
                    spot = Juliana.Spot{spot_data_type}(
                        current_ID,
                        convert(T, t),
                        convert(T, u),
                        convert(T, 1), # unit weight initialisation
                        convert(T, E),
                        preabsorber,
                    )
                    push!(spots, spot)
                    current_ID += 1
                end
            end
        end
    end

    # Remove duplicate spots and reindex them.
    spots = unique(position_t_u_energy_absorber, spots)

    ID = 0
    for (i, spot) in enumerate(spots)
        spots[i] = Juliana.Spot(
            ID,
            spot.t,
            spot.u,
            spot.weight,
            spot.energykeV,
            spot.numberOfAbsorbers,
        )
        ID += 1
    end

    return spots
end


function place_spots_on_grid(densities::Juliana.ScalarGrid,
                             field::Juliana.FieldDefinition,
                             placement_structure::Juliana.Structure,
                             depth_dose_curves::Juliana.LookUpTable;
                             placement_margin::Real=0.5,
                             lateral_spacing::Real=0.4,
                             depth_spacing::Real=0.2)
    # Figure out which points are relevant for spot placement.
    # These are called placement_points.
    spot_placement_mask = placement_structure.distanceFromStructure .<= placement_margin
    placement_points, placement_indices = Juliana.mask_to_points_and_indices(
        placement_structure.grid,
        spot_placement_mask,
    )

    wed = calculate_wed(
        densities,
        field.gantryAngle,
        field.couchAngle,
        placement_points,
    )
    wed_cube = Juliana.flat_vector_to_cube(
        placement_structure.grid,
        placement_indices,
        wed,
    )

    return place_spots_on_grid(
        field,
        placement_structure,
        wed_cube,
        depth_dose_curves;
        placement_margin=placement_margin,
        lateral_spacing=lateral_spacing,
        depth_spacing=depth_spacing,
    )
end


function place_spots_on_grid(densities::Juliana.ScalarGrid,
                             fields::Vector{Juliana.FieldDefinition{T}},
                             placement_structures::Vector{S},
                             depth_dose_curves::Juliana.LookUpTable;
                             placement_margin::Real=0.5,
                             lateral_spacing::Real=0.4,
                             depth_spacing::Real=0.2) where {T<:Real, S <: Juliana.Structure}
    grids = [s.grid for s in placement_structures]
    @assert all([grid == grids[1] for grid in grids])
    grid = grids[1]

    N = length(fields)
    all_spots = Vector{Vector{Juliana.Spot{T}}}(undef, N)
    for j in 1:N
        # Figure out which points are relevant for spot placement.
        # These are called placement_points.
        spot_placement_mask = placement_structures[j].distanceFromStructure .<= placement_margin
        placement_points, placement_indices = Juliana.mask_to_points_and_indices(
            grid,
            spot_placement_mask,
        )

        wed = calculate_wed(
            densities,
            fields[j].gantryAngle,
            fields[j].couchAngle,
            placement_points,
        )

        wed_cube = Juliana.flat_vector_to_cube(
            grid,
            placement_indices,
            wed[:, 1],
        )
        all_spots[j] = place_spots_on_grid(
            fields[j],
            placement_structures[j],
            wed_cube,
            depth_dose_curves;
            placement_margin=placement_margin,
            lateral_spacing=lateral_spacing,
            depth_spacing=depth_spacing,
        )
    end
    all_spots = relabel_spots(all_spots)

    fields_new = Vector{Juliana.FieldDefinition{Float32}}(undef, length(fields))
    for (i, (field_old, spots)) in enumerate(zip(fields, all_spots))
        fields_new[i] = Juliana.FieldDefinition(
            field_old.label,
            field_old.id,
            field_old.gantryAngle,
            field_old.couchAngle,
            field_old.nozzleExtraction,
            field_old.fieldCenter,
            spots,
        )
    end

    return fields_new
end


function place_spots_on_grid(densities::Juliana.ScalarGrid,
                             plan::Juliana.TreatmentPlan{T},
                             placement_structures::Vector{Juliana.Structure},
                             depth_dose_curves::Juliana.LookUpTable;
                             placement_margin::Real=0.5,
                             lateral_spacing::Real=0.4,
                             depth_spacing::Real=0.2) where {T<:Real}
    all_spots = place_spots_on_grid(
        densities,
        plan.fields,
        placement_structures,
        depth_dose_curves;
        placement_margin=placement_margin,
        lateral_spacing=lateral_spacing,
        depth_spacing=depth_spacing,
    )
    fields_new = Vector{Juliana.FieldDefinition{T}}(undef, length(plan.fields))
    for (i, field_old) in enumerate(plan.fields)
        fields_new[i] = Juliana.FieldDefinition(
            field_old.label,
            field_old.id,
            field_old.gantryAngle,
            field_old.couchAngle,
            field_old.nozzleExtraction,
            field_old.fieldCenter,
            all_spots[i],
        )
    end
    return Juliana.TreatmentPlan(fields_new, plan.prescribedDose)
end


function hexagonal_tu_placement(stu_grid::Juliana.Grid)
    start = stu_grid.origin
    stop = start .+ stu_grid.size .* stu_grid.spacing
    
    N, M = stu_grid.size[2], stu_grid.size[3]
    positions = Matrix{Float32}(undef, (N+1) * M, 2)

    i = 1
    for it in 0:N
        for iu in 0:M-1
            if mod(it, 2) == 1
                iu += 0.5
            end
            positions[i, :] .= (start[2] + it * stu_grid.spacing[2], start[3] + iu * stu_grid.spacing[3])
            i += 1
        end
    end
    
    return positions
end


"""
The WED cube is an array of CT shape. Its values must have been calculated at least for all points on the STU grid.
"""
function place_spots_on_hexagonal_grid(field::Juliana.FieldDefinition,
                                       placement_structure::Juliana.Structure,
                                       wed_cube::AbstractArray{T, 3},
                                       depth_dose_curves::Juliana.LookUpTable;
                                       first_spot_ID::Integer=0,
                                       placement_margin::Real=0.5,
                                       lateral_spacing::Real=0.4,
                                       depth_spacing::Real=0.4) where {T<:Real}
    stu_grid = Juliana.build_stu_grid(
        field,
        placement_structure;
        placement_margin=placement_margin,
        lateral_spacing=lateral_spacing,
        depth_spacing=depth_spacing,
    )
    s_min = stu_grid.origin[1]
    s_max = s_min + stu_grid.spacing[1] * stu_grid.size[1]

    spot_data_type = Float32

    # Actually place the spots.
    ranges = [Juliana.get_range(depth_dose_curves, E) for E in depth_dose_curves.energies]
    ranges_with_preabsorber = ranges .- Juliana.PREABSORBER_WED
    iso_center = (
        field.fieldCenter["x"],
        field.fieldCenter["y"],
        field.fieldCenter["z"],
    )

    tu_positions = hexagonal_tu_placement(stu_grid)

    current_ID = first_spot_ID
    spots = Vector{Juliana.Spot{spot_data_type}}(undef, 0)
    for (t, u) in eachrow(tu_positions)
        for s in s_min:depth_spacing:s_max
            p = Juliana.StuToXyz(
                s, t, u,
                iso_center,
                field.gantryAngle,
                field.couchAngle,
            )
            ip = Juliana.xyz_to_index(p, placement_structure.grid)
            if !Juliana.index_is_valid(ip[1], ip[2], ip[3], size(wed_cube))
                # Spot would be placed outside of the CT grid.
                continue
            end
            depth = wed_cube[ip[1], ip[2], ip[3]]

            # Check whether we need a preabsorber.
            preabsorber = 0
            current_ranges = ranges
            if depth <= minimum(ranges)
                preabsorber = 1
                current_ranges = ranges_with_preabsorber
            end

            # Get the closest energy.
            energy_index = argmin(abs.(current_ranges .- depth))
            E = depth_dose_curves.energies[energy_index]

            if placement_structure.distanceFromStructure[ip[1], ip[2], ip[3]] <= placement_margin
                spot = Juliana.Spot{spot_data_type}(
                    current_ID,
                    convert(spot_data_type, t),
                    convert(spot_data_type, u),
                    convert(spot_data_type, 1), # unit weight initialisation
                    convert(spot_data_type, E),
                    preabsorber,
                )
                push!(spots, spot)
                current_ID += 1
            end
        end
    end

    # Remove duplicate spots and reindex them.
    spots = unique(Juliana.position_t_u_energy_absorber, spots)

    ID = first_spot_ID
    for (i, spot) in enumerate(spots)
        spots[i] = Juliana.Spot(
            ID,
            spot.t,
            spot.u,
            spot.weight,
            spot.energykeV,
            spot.numberOfAbsorbers,
        )
        ID += 1
    end

    return spots
end


function place_spots_on_hexagonal_grid(densities::Juliana.ScalarGrid,
                                       field::Juliana.FieldDefinition,
                                       placement_structure::Juliana.Structure,
                                       depth_dose_curves::Juliana.LookUpTable;
                                       first_spot_ID::Integer=0,
                                       placement_margin::Real=0.5,
                                       lateral_spacing::Real=0.4,
                                       depth_spacing::Real=0.4)
    placement_mask = placement_structure.distanceFromStructure .<= placement_margin
    placement_points, placement_indices = Juliana.mask_to_points_and_indices(densities.grid, placement_mask)
    
    wed = Juliana.calculate_wed(
        densities,
        [field.gantryAngle],
        [field.couchAngle],
        placement_points,
    )

    wed_cube = Juliana.flat_vector_to_cube(densities.grid, placement_indices, wed)
    
    return place_spots_on_hexagonal_grid(
        field,
        placement_structure,
        wed_cube,
        depth_dose_curves;
        placement_margin=placement_margin,
        lateral_spacing=lateral_spacing,
        depth_spacing=depth_spacing,
        first_spot_ID=first_spot_ID,
    )
end


function wed_for_energy_and_preabsorber(E::Real,
                                        use_preabsorber::Bool,
                                        depth_dose_curve::Juliana.LookUpTable;
                                        preabsorber_WED=Juliana.PREABSORBER_WED)
    wed = Juliana.get_range(depth_dose_curve, E)
    if use_preabsorber
        wed -= preabsorber_WED
    end
    return wed
end


function relabel_spots(spots_per_field::Vector{Vector{Juliana.Spot{Float32}}})
    spots_new = [Vector{Juliana.Spot{Float32}}(undef, length(spots)) for spots in spots_per_field]
    spot_ID = 0
    for (i, spots) in enumerate(spots_per_field)
        for (j, spot) in enumerate(spots)
            spots_new[i][j] = Juliana.Spot(
                spot_ID,
                spot.t,
                spot.u,
                spot.weight,
                spot.energykeV,
                spot.numberOfAbsorbers,
            )
            spot_ID += 1
        end
    end
    return spots_new
end


######################################
# Initialisation of spot weights.
######################################
using LineSearches
using Optim


function ideal_dose_1D(z; zmin, zmax, σ=0.8f0)
    if z < zmin
        return exp(-(zmin - z)^2 / (2*σ^2))
    elseif z < zmax
        return 1
    else
        return exp(-(zmax - z)^2 / (2*σ^2))
    end
end


"""
    evaluate_curve(wed, E, preabsorber, ddc, preabsorber_wed=Juliana.PREABSORBER_WED)

Evaluate the depth-dose curve for the given energy and preabsorber using the lookup table ddc
at depth wed. Interpolate linearly.
"""
function evaluate_curve(wed, E, preabsorber, ddc, preabsorber_wed=Juliana.PREABSORBER_WED)
    i = findfirst(==(E), ddc.energies)
    curve = ddc.table[i, 1:ddc.firstNanIndex[i]-1]
    wed_nodes = collect(ddc.x0[i]:ddc.dx[i]:ddc.x0[i]+(length(curve)-1)*ddc.dx[i]);

    if preabsorber == 1
        wed_nodes = wed_nodes .- Juliana.PREABSORBER_WED
        i_positive = findfirst(>=(0), wed_nodes)
        wed_nodes = wed_nodes[i_positive:end]
        curve = curve[i_positive:end]
    end

    # Interpolate linearly.
    if wed > maximum(wed_nodes)
        j = -1
    else
        j = findfirst(>(wed), wed_nodes)
        if isnothing(j)
            j = -1
        else
            j -= 1
        end
    end

    if j == -1
        return curve[end]
    elseif j == 0
        return curve[1]
    else
        return curve[j] + (curve[j+1] - curve[j]) * (wed_nodes[j] - wed) / (wed_nodes[j+1] - wed_nodes[j])
    end
end

"""
    calculate_depth_dose_profile(depth::Real, ddc_functions, ddc_weights)

Evaluate the 1D depth-dose profile at the given depth for the given Bragg peak functions
and Bragg peak weights.
"""
function calculate_depth_dose_profile(depth::Real, ddc_functions, ddc_weights)
    result = 0
    for (ddc, w) in zip(ddc_functions, ddc_weights)
        result += ddc(depth) * w
    end
    return result
end


function optimise_depth_dose_profile(depths,
                                     ddc_functions,
                                     early_stopping_delta,
                                     early_stopping_patience,
                                     min_range,
                                     max_range)
    early_stopping_delta = 0.0001
    early_stopping_patience = 10
    width = 0.1

    # TODO: Change from spot ranges to target extent.
    depths = collect(min_range:0.01:max_range+1);

    # Only optimise the plateau region.
    ideal = [ideal_dose_1D(z, zmin=min_range, zmax=max_range) for z in depths]
    i_min = findfirst(ideal .> 1 - width/2)
    i_max = findlast(ideal .> 1 - width/2)

    depths = depths[i_min:i_max]
    ideal = ideal[i_min:i_max]

    Dij = Matrix{Float32}(undef, length(depths), length(ddc_functions));
    for (j, f) in enumerate(ddc_functions)
        Dij[:, j] .= f.(depths)
    end

    ###########################
    # OPTIMISE
    ###########################
    # Define the loss function and its gradient.
    function l(weights)
        dose = Dij * weights
        difference = dose .- ideal

        return sum(difference.^2)
    end

    function l_grad!(gradient, weights)
        dose = Dij * weights
        gradient .= 2 * Dij' * (dose .- ideal)
    end

    # Prepare optimisation.
    initial_iteration = 1
    maximum_iterations = 100

    early_stopping = Juliana.build_early_stopping(
        early_stopping_delta,
        early_stopping_patience,
    )

    # Optimise depth dose profile.
    w = ones(length(ddc_functions))
    w .= w ./ maximum(Dij * w)

    alg = LBFGS(linesearch=LineSearches.HagerZhang())
    options = Optim.Options()
    d = Optim.promote_objtype(alg, w, :finite, true, l, l_grad!)
    state = Optim.initial_state(alg, options, d, w)

    # Optimise.
    history = Vector{Dict{String, Float32}}()
    gradients = Array{Vector{Float32}, 1}()
    loss = Inf
    for i in initial_iteration:maximum_iterations
        # Check for convergence.
        if early_stopping(loss, i)
            break
        end

        # Iterate.
        Optim.update_state!(d, state, alg)
        Optim.update_g!(d, state, alg)
        Optim.update_h!(d, state, alg)
        clamp!(state.x, 0, Inf)

        loss = l(state.x)
    end
    w_opt = state.x
    clamp!(w_opt, 0, Inf);
end


function initialise_spot_weights(spots, tps::Juliana.JulianaTps)
    # Build a list of (energy, preabsorber) pairs.
    configs = sort(unique([(s.energykeV, s.numberOfAbsorbers) for s in spots]));
    ranges_in_plan = sort([Juliana.wed_for_energy_and_preabsorber(E, preab == 1, tps.d_depth_dose_curves) for (E, preab) in configs])
    min_range = minimum(ranges_in_plan)
    max_range = maximum(ranges_in_plan)
    is_in_target = depth -> (min_range <= depth) && (depth <= max_range)

    ddc_functions = [wed -> evaluate_curve(wed, E, preab, tps.d_depth_dose_curves, Juliana.PREABSORBER_WED) for (E, preab) in configs]

    # Optimise depth dose profile.
    early_stopping_delta = 0.01
    early_stopping_patience = 10

    # TODO: Change from spot ranges to target extent.
    depths = collect(min_range:0.01:max_range+1);

    weights = optimise_depth_dose_profile(
        depths,
        ddc_functions,
        early_stopping_delta,
        early_stopping_patience,
        min_range,
        max_range,
    )

    # Build a new spot list with the given weights.
    spots_new = similar(spots)
    for (i, spot) in enumerate(spots)
        w = weights[findfirst(==((spot.energykeV, spot.numberOfAbsorbers)), configs)]
        spots_new[i] = Juliana.Spot(
            spot.id,
            spot.t,
            spot.u,
            convert(Float32, w),
            spot.energykeV,
            spot.numberOfAbsorbers,
        )
    end

    return spots_new
end


function initialise_spot_weights(field::Juliana.FieldDefinition{T}, tps::Juliana.JulianaTps) where {T <: Real}
    spots = initialise_spot_weights(field.spots, tps)
    return Juliana.FieldDefinition{T}(
        field.label,
        field.id,
        field.gantryAngle,
        field.couchAngle,
        field.nozzleExtraction,
        field.fieldCenter,
        spots,
    )
end


function initialise_spot_weights(plan::Juliana.TreatmentPlan{T}, tps::Juliana.JulianaTps) where {T <: Real}
    fields = [initialise_spot_weights(f, tps) for f in plan.fields]
    return Juliana.TreatmentPlan{T}(fields, plan.prescribedDose)
end
