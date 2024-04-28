using SparseArrays


function shift_ct(ct, shift)
    HU_AIR = -1000

    # Calculate shifted grid at whose positions we need the HU values.
    shifted_grid = Juliana.Grid(
        ct.grid.spacing,
        ct.grid.origin .+ SVector{3}(shift),
        ct.grid.size,
    )

    # Add guard cells filled with air because
    # the interpolation cannot deal with boundary cells.
    n_guard_cells = convert(Int32, ceil(maximum(abs.(shift ./ ct.grid.spacing)))+1)

    guarded_size = SVector{3}(
        ct.grid.size[1]+2*n_guard_cells,
        ct.grid.size[2]+2*n_guard_cells,
        ct.grid.size[3]+2*n_guard_cells,
    )

    guarded_grid = Juliana.Grid(
        ct.grid.spacing,
        ct.grid.origin .- n_guard_cells .* ct.grid.spacing,
        guarded_size,
    )

    ct_guarded_data = Array{Int16, 3}(undef, guarded_size[1], guarded_size[2], guarded_size[3])
    fill!(ct_guarded_data, HU_AIR)
    ct_guarded_data[n_guard_cells+1:end-n_guard_cells, n_guard_cells+1:end-n_guard_cells, 1+n_guard_cells:end-n_guard_cells] .= ct.data;
    ct_guarded = Juliana.ScalarGrid(
        ct_guarded_data,
        guarded_grid,
    );

    # Points at which we need to calculate the HU values.
    points, indices = Juliana.grid_to_points_and_indices(shifted_grid);

    # Interpolate the HU values.
    hu_interpolated = convert.(
        Int16,
        round.(Juliana.interpolate_linearly(ct_guarded, points)),
    )
    hu_interpolated_cube = Juliana.flat_vector_to_cube(
        shifted_grid,
        indices,
        hu_interpolated,
    )
    # The shifted CT should be in the SAME coordinates as the original CT, i. e. use the same grid.
    return Juliana.ScalarGrid(
        convert.(eltype(ct.data), hu_interpolated_cube),
        ct.grid,
    )
end


function scale_ct(ct, factor, minimum_HU=-1000, maximum_HU=3071)
    @warn "Rescaled CTs are clamped to [$(minimum_HU), $(maximum_HU)]"
    return Juliana.ScalarGrid(
        convert.(
            eltype(ct.data),
            clamp.(
                round.(ct.data .* factor),
                minimum_HU,
                maximum_HU,
            ),
        ),
        ct.grid
    )
end


function generate_scenarios(ct::Juliana.ScalarGrid, shift::Float32=0.3f0, ct_difference_scale::Float32=0.03f0)
    scenarios = Vector{typeof(ct)}(undef, 9)
    scenarios[1] = Juliana.ScalarGrid(copy(ct.data), ct.grid)
    scenarios[2] = Juliana.shift_ct(ct, [shift, 0.f0, 0.f0])
    scenarios[3] = Juliana.shift_ct(ct, [-shift, 0.f0, 0.f0])
    scenarios[4] = Juliana.shift_ct(ct, [0.f0, shift, 0.f0])
    scenarios[5] = Juliana.shift_ct(ct, [0.f0, -shift, 0.f0])
    scenarios[6] = Juliana.shift_ct(ct, [0.f0, 0.f0, shift])
    scenarios[7] = Juliana.shift_ct(ct, [0.f0, 0.f0, -shift])
    scenarios[8] = Juliana.scale_ct(ct, (1.f0 - ct_difference_scale))
    scenarios[9] = Juliana.scale_ct(ct, (1.f0 + ct_difference_scale))
    return scenarios
end


function build_robust_Dij(Dijs, prescriptions, structures, optimisation_point_indices)
    # Define regions of interest.
    oar_mask = Juliana.build_constraint_oar_mask(
        prescriptions,
        structures,
    )
    target_mask = Juliana.calculate_whole_target_mask(
        prescriptions,
        structures,
    )

    oar_only_mask = oar_mask .&& (.!(target_mask))
    target_only_mask = target_mask .&& (.!(oar_mask))
    overlap_mask = oar_mask .&& target_mask

    # Transpose the matrix to speedup memory accesses.
    # Each column of Dij_T corresponds to an optimisation point.
    Dijs_T = [sparse(Dij') for Dij in Dijs]

    # Determine for each row to which region the associated
    # region belongs.
    is_optim_point_oar = [oar_only_mask[ix, iy, iz] for (ix, iy, iz) in eachcol(optimisation_point_indices)]
    is_optim_point_target = [target_only_mask[ix, iy, iz] for (ix, iy, iz) in eachcol(optimisation_point_indices)]
    is_optim_point_overlap = [overlap_mask[ix, iy, iz] for (ix, iy, iz) in eachcol(optimisation_point_indices)]

    # Allocate the new matrix using the same sparsity structure as the first scenario's Dij matrix.
    Dij_T = Dijs_T[1]
    Dij_T_combined = similar(Dij_T)
    fill!(Dij_T_combined, 0)

    n_spots = size(Dij_T_combined, 1)
    n_points = size(Dij_T_combined, 2)

    # Reduce the matrix row by row.
    columns = Vector{SparseVector}(undef, n_points)
    for j in 1:n_points
        if is_optim_point_target[j]
            # Target point.
            val = min.([Dij[:, j] for Dij in Dijs_T]...)
        elseif is_optim_point_overlap[j]
            # Overlap point.
            val = mean([Dij[:, j] for Dij in Dijs_T])
        else
            # OARs or normal tissue point.
            val = max.([Dij[:, j] for Dij in Dijs_T]...)
        end
        val = sparsevec(val)

        columns[j] = val
    end

    # Combine the columns to a new sparse matrix.
    n_elements = sum([length(col.nzval) for col in columns])

    col_indices = Vector{Int64}(undef, n_elements)
    row_indices = Vector{Int64}(undef, n_elements)
    values = Vector{Float64}(undef, n_elements)

    i_entry = 1
    for column_index in 1:length(columns)
        n_entries = length(columns[column_index].nzind)
        if n_entries > 0
            col_indices[i_entry:i_entry+n_entries-1] .= column_index
            row_indices[i_entry:i_entry+n_entries-1] = columns[column_index].nzind
            values[i_entry:i_entry+n_entries-1] = columns[column_index].nzval
            i_entry += n_entries
        end
    end

    Dij_full = sparse(row_indices, col_indices, values);

    # Consistency check.
    for (j, column) in enumerate(columns)
        @assert Dij_full[:, j] == column
    end

    # Transpose back to have shape (n_optimisation_points, n_spots).
    Dij_full = sparse(Dij_full')

    return Dij_full
end
