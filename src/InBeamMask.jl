using CUDA
using CUDAKernels
using KernelAbstractions


function fill_in_path_mask!(px, py, pz,
                            dx, dy, dz,
                            grid,
                            in_path_mask)
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
            in_path_mask[ix, iy, iz] = true
        end
        t = t_next
    end
    return wed
end


@kernel function calculate_in_path_mask_kernel(d_points, d_direction, d_grid, d_in_path_mask)
    j = @index(Global, Cartesian)[1]

    fill_in_path_mask!(
        d_points[1, j],
        d_points[2, j],
        d_points[3, j],
        d_direction[1],
        d_direction[2],
        d_direction[3],
        d_grid,
        d_in_path_mask,
    )
end


calculate_in_path_mask_kernel_gpu = calculate_in_path_mask_kernel(CUDADevice(), 32);


function calculate_in_beam_mask(structure_mask, gantry_angle, couch_angle, grid)
    points, indices = Juliana.mask_to_points_and_indices(grid, structure_mask)
    
    in_path_mask = Array{Bool, 3}(undef, grid.size...)
    fill!(in_path_mask, false)
    
    direction = Juliana.angles_to_direction(gantry_angle, couch_angle)

    d_points = cu(points)
    d_direction = cu(direction)
    d_grid = cu(grid)
    d_in_path_mask = cu(in_path_mask)

    event = calculate_in_path_mask_kernel_gpu(
        d_points,
        d_direction,
        d_grid,
        d_in_path_mask,
        ndrange=size(points, 2),
    )
    wait(event)

    return collect(d_in_path_mask)
end
