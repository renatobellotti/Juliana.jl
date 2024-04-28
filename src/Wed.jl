using CUDAKernels
using KernelAbstractions


function calculate_wed(px, py, pz,
                       dx, dy, dz,
                       grid,
                       densities)
    # Implementation of Siddon's algorithm.
    START_STEP = 100.f0

    start_pos_x = px + START_STEP * dx
    start_pos_y = py + START_STEP * dy
    start_pos_z = pz + START_STEP * dz

    wed = 0.f0


    # Setup the ray marching initial intersections
    # and marching directions.
    #############################################################
    # - 0.5 spacing: The CT grid values are given at the center
    # of the voxel. The beginning of the CT is half a voxel
    # to the bottom lower left (in CT coordinates).
    Δx = grid.spacing[1]
    x = grid.origin[1] - grid.spacing[1] / 2
    if start_pos_x > px
        Δx *= -1
        x += grid.size[1] * grid.spacing[1]
    end

    Δy = grid.spacing[2]
    y = grid.origin[2] - grid.spacing[2] / 2
    if start_pos_y > py
        Δy *= -1
        y += grid.size[2] * grid.spacing[2]
    end

    Δz = grid.spacing[3]
    z = grid.origin[3] - grid.spacing[3] / 2
    if start_pos_z > pz
        Δz *= -1
        z += grid.size[3] * grid.spacing[3]
    end

    # Initialise the marching: arch until you encounter a point that lies within the CT.
    tx = 100.f0
    if !(abs(dx) ≈ 0)
        tx = (start_pos_x - x) / (dx * START_STEP)
    end

    ty = 100.f0
    if !(abs(dy) ≈ 0)
        ty = (start_pos_y - y) / (dy * START_STEP)
    end

    tz = 100.f0
    if !(abs(dz) ≈ 0)
        tz = (start_pos_z - z) / (dz * START_STEP)
    end

    t = -1000.f0
    finished = false
    while !finished
        # Calculate the current ray-grid intersections
        # in terms of the ray parameter t for each direction.
        tx = 100.f0
        if !(abs(dx) ≈ 0)
            tx = (start_pos_x - x) / (dx * START_STEP)
        end

        ty = 100.f0
        if !(abs(dy) ≈ 0)
            ty = (start_pos_y - y) / (dy * START_STEP)

        end

        tz = 100.f0
        if !(abs(dz) ≈ 0)
            tz = (start_pos_z - z) / (dz * START_STEP)
        end

        # Get the next intersection point.
        if (tx <= ty) && (tx <= tz)
            t_next = tx
            x += Δx
        elseif (ty <= tx) && (ty <= tz)
            t_next = ty
            y += Δy
        elseif (tz <= tx) && (ty <= ty)
            t_next = tz
            z += Δz
        else
            return NaN
        end

        if t == 0
            t = t_next
            continue
        end

        # Check convergence.
        if t_next > 1.f0
            t_next = 1.f0
            finished = true
        end

        # Update the WED.
        p_eval_x = start_pos_x - (t + t_next) / 2 * START_STEP * dx
        p_eval_y = start_pos_y - (t + t_next) / 2 * START_STEP * dy
        p_eval_z = start_pos_z - (t + t_next) / 2 * START_STEP * dz

        ix = convert(Int32, round((p_eval_x - grid.origin[1]) / grid.spacing[1]) + 1)
        iy = convert(Int32, round((p_eval_y - grid.origin[2]) / grid.spacing[2]) + 1)
        iz = convert(Int32, round((p_eval_z - grid.origin[3]) / grid.spacing[3]) + 1)

        if Juliana.index_is_valid(ix, iy, iz, grid.size)
            # Only update once we've reached the CT entry region.
            d = (t_next - t) * START_STEP
            ρ = densities[ix, iy, iz]
            wed += d * ρ
        end
        t = t_next
    end
    return wed
end


# CPU implementation: Called when points is not a GPU matrix.
function calculate_wed(points::Matrix, gantry_angles, couch_angles, grid, densities::Array)
    directions = Juliana.angles_to_direction(gantry_angles, couch_angles)

    wed = Matrix{Float32}(undef, size(points, 2), size(directions, 2))
    for j in 1:size(directions, 2)
        @Threads.threads for i in 1:size(points, 2)
            x = points[1, i]
            y = points[2, i]
            z = points[3, i]

            dx = directions[1, j]
            dy = directions[2, j]
            dz = directions[3, j]

            point_wed = Juliana.calculate_wed(
                x, y, z,
                dx, dy, dz,
                grid,
                densities,
            )

            wed[i, j] = point_wed
        end
    end
    return wed
end


@kernel function wed_kernel(d_wed, @Const(d_densities), @Const(d_grid), @Const(d_points), @Const(d_directions))
    index = @index(Global, Cartesian)
    i = index[1]
    j = index[2]

    x = d_points[1, i]
    y = d_points[2, i]
    z = d_points[3, i]

    direction_x = d_directions[1, j]
    direction_y = d_directions[2, j]
    direction_z = d_directions[3, j]

    wed = calculate_wed(
        x, y, z,
        direction_x, direction_y, direction_z,
        d_grid,
        d_densities,
    )

    d_wed[i, j] = wed
end


calculate_wed_simple = wed_kernel(CUDADevice(), 32);


function calculate_wed(densities::Juliana.ScalarGrid,
                       gantry_angle::Real,
                       couch_angle::Real,
                       grid::Grid)
    points = Juliana.grid_to_points(grid)
    indices = Juliana.xyz_to_index(points, grid);

    wed = calculate_wed(
        densities,
        [gantry_angle],
        [couch_angle],
        points,
    )
    wed_cube = Juliana.flat_vector_to_cube(grid, indices, wed)

    return Juliana.ScalarGrid(wed_cube, grid)
end


function calculate_wed(densities::Juliana.ScalarGrid,
                       gantry_angle::Real,
                       couch_angle::Real,
                       points::AbstractMatrix{U}) where {U<:Real}
    return calculate_wed(
        densities,
        [gantry_angle],
        [couch_angle],
        points,
    )
end


function calculate_wed(densities::Juliana.ScalarGrid,
                       gantry_angles::Vector{U},
                       couch_angles::Vector{U},
                       points::AbstractMatrix{U}) where {U<:Real}
    d_densities = CuArray{eltype(densities.data), 3}(densities.data)

    M = length(gantry_angles)
    directions = Matrix{eltype(gantry_angles)}(undef, 3, M)
    for j in 1:M
        directions[:, j] .= vec(Juliana.angles_to_direction(
            gantry_angles[j],
            couch_angles[j],
        ))
    end
    d_directions = CuArray{eltype(directions), 2}(directions)
    
    d_points = CuArray{U, 2}(points)
    N = size(d_points, 2)

    d_grid = cu(densities.grid)
    d_wed = CuArray{Float64, 2}(undef, N, M)
    fill!(d_wed, Inf)

    event = Juliana.calculate_wed_simple(
        d_wed,
        d_densities,
        d_grid,
        d_points,
        d_directions,
        ndrange=(N, M),
    )
    wait(event)
    wed = collect(d_wed)
    return wed
end
