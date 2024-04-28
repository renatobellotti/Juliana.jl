using CUDA


function get_in_field_mask(d_direction, d_target_points, d_grid)
    d_in_beam = cu(zeros(Float32, collect(d_grid.size)...))
    event = Juliana.cast_rays(
        d_in_beam,
        d_target_points,
        d_direction,
        d_grid,
        ndrange=size(d_target_points, 2),
    )
    wait(event)
    return d_in_beam
end


function wed_variance_scores(gantry_angles, couch_angles, d_target_points, d_densities, d_grid)
    # Convert the beam angles to direction vectors.
    angles = vec(collect(Iterators.product(gantry_angles, couch_angles)))
    N_angles =length(angles)

    directions = Matrix{Float32}(undef, 3, N_angles)
    for (i, (gantry, couch)) in enumerate(angles)
        directions[:, i] .= Juliana.angles_to_direction(gantry, couch)
    end
    d_directions = cu(directions);

    # Calculate the WED for each target point and each direction.
    N_points = size(d_target_points, 2)
    d_wed = cu(zeros(Float32, N_points, N_angles))

    event = Juliana.calculate_wed_simple(
        d_wed,
        d_densities,
        cu(d_grid),
        d_target_points,
        d_directions,
        ndrange=(N_points, N_angles),
    )
    wait(event)
    wed = collect(d_wed)

    # Calculate the maximum WED difference for each angle.
    wed_range = vec(maximum(wed, dims=1) .- minimum(wed, dims=1));
    wed_range = reshape(
        wed_range,
        (length(gantry_angles), length(couch_angles)),
    )

    return wed_range
end


function calculate_score_map(distances, d_target_points, d_grid, gantry_angles, couch_angles; which="sum")
    scores = Array{Float32, 2}(
        undef,
        length(gantry_angles),
        length(couch_angles),
    )
    for (j, couch_angle) in enumerate(couch_angles)
        for (i, gantry_angle) in enumerate(gantry_angles)
            direction = Juliana.angles_to_direction(
                gantry_angle,
                couch_angle,
            )
            in_beam_mask = Juliana.get_in_field_mask(
                cu(direction),
                d_target_points,
                d_grid,
            )
            masked_distance = collect(in_beam_mask) .* distances
            if which == "min"
                masked_distance = replace(masked_distance, 0 => Inf)
                scores[i, j] = minimum(masked_distance)
            elseif which == "sum"
                scores[i, j] = sum(masked_distance)
            end
        end
    end
    return scores
end


"""
    find_n_largest(scores, gantry_angles, couch_angles, n)

The the n beam angles with the highest score.

Axis 1 in scores is gantry angles, axis 2 is couch angles.
"""
function find_n_largest(scores, gantry_angles, couch_angles, n)
    scores = deepcopy(scores)
    
    best_angles = Array{Tuple{Float32, Float32}, 1}()
    for i in 1:n
        i_gantry, i_couch = Tuple(argmax(scores))
        push!(best_angles, (gantry_angles[i_gantry], couch_angles[i_couch]))
        scores[i_gantry, i_couch] = -Inf
    end
    return best_angles
end


function is_shooting_from_feet(gantry_angle, couch_angle)
    gantry_forbidden = (45 <= gantry_angle) && (gantry_angle <= 135)
    couch_forbidden = (-135 <= couch_angle) && (couch_angle <= -45)
    return gantry_forbidden && couch_forbidden
end


function is_shooting_through_nasal_cavity(gantry_angle, couch_angle)
    gantry_forbidden = (-30 <= gantry_angle) && (gantry_angle <= 60)
    couch_forbidden = (-150 <= couch_angle) && (couch_angle <= -30)
    return gantry_forbidden || couch_forbidden
end
    

function is_valid_angle(gantry_angle, couch_angle)
    return !is_shooting_from_feet(gantry_angle, couch_angle) &&
           !is_shooting_through_nasal_cavity(gantry_angle, couch_angle)
end


function angle_difference(a, b)
    couch_diff = abs(a[2] - b[2])
    if couch_diff >= 180
        couch_diff = 180 - couch_diff
    end
    
    difference = [
        abs(a[1] - b[1]),
        couch_diff,
    ]
    return difference
end


"""
    select_candidates(candidates, n; MINIMUM_SEPARATION=30)

Select n entries from candidates such that the MINIMUM_SEPARATION between
them is not exceeded.
"""
function select_candidates(candidates, n; MINIMUM_SEPARATION=30)
    @assert n <= length(candidates)
    angles = Array{Float32, 2}(undef, 2, 1)
    for candidate in candidates
        if is_valid_angle(candidate...)
            angles[:, 1] .= collect(candidate)
            break
        end
    end

    for candidate in candidates
        candidate = collect(candidate)
        distances = [angle_difference(angle, candidate) for angle in eachcol(angles)]
        if is_valid_angle(candidate...) && (minimum([maximum(d) for d in distances]) > MINIMUM_SEPARATION)
            angles = hcat(angles, candidate)
        end
        if size(angles, 2) == n
            break
        end
    end
    return angles
end


function select_beam_arrangement(ct,
                                 structures,
                                 prescriptions;
                                 n_beams=4,
                                 gantry_angles=-30:10:180,
                                 couch_angles=-180:10:180)
    target_mask = Juliana.calculate_whole_target_mask(
        prescriptions,
        structures,
    )
    target_points = Juliana.mask_to_points(
        ct.grid,
        target_mask,
    )

    # Calculate a tensor that contains each voxel's distance to the
    # closest OAR.
    distances = Array{Float32}(undef, size(ct.data))
    fill!(distances, Inf)
    for constraint in prescriptions.constraints
        oar = structures[constraint.structure_name]
        distances = min.(distances, oar.distanceFromStructure)
    end
    # The slices without contours still have Inf, which will confuse the
    # optimiser. We set it to zero.
    distances = replace(distances, Inf => 0)

    d_target_points = cu(target_points)
    d_grid = cu(ct.grid)

    scores = Juliana.calculate_score_map(
        distances,
        d_target_points,
        d_grid, 
        gantry_angles,
        couch_angles,
    )

    candidates = Juliana.find_n_largest(
        scores,
        gantry_angles,
        couch_angles,
        length(scores),
    )

    return Juliana.select_candidates(candidates, 4)
end
