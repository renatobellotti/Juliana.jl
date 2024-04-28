"""
    interpolate_linearly(nodes::Juliana.ScalarGrid, p::Tuple{T, T, T})

Calculate the value at position p using trilinear interpolation and the given values on a regular grid.
"""
function interpolate_linearly(nodes::ScalarGrid, p::Tuple{T, T, T}) where {T<:Real}
    # Implements the following algorithm:
    # https://en.wikipedia.org/wiki/Trilinear_interpolation
    # Origin of the "interpolation neighbourhood cube cell".
    ix = convert(Int64, floor((p[1] - nodes.grid.origin[1]) / nodes.grid.spacing[1])) + 1
    iy = convert(Int64, floor((p[2] - nodes.grid.origin[2]) / nodes.grid.spacing[2])) + 1
    iz = convert(Int64, floor((p[3] - nodes.grid.origin[3]) / nodes.grid.spacing[3])) + 1

    # We need at least one boundary cell on all sides.
    if (ix < 2) || (ix >= nodes.grid.size[1]) || (iy < 2) || (iy >= nodes.grid.size[2]) || (iz < 2) || (iz >= nodes.grid.size[3])
        return  convert(T, NaN)
    end

    x0 = nodes.grid.origin[1] + (ix - 1) * nodes.grid.spacing[1]
    y0 = nodes.grid.origin[2] + (iy - 1) * nodes.grid.spacing[2]
    z0 = nodes.grid.origin[3] + (iz - 1) * nodes.grid.spacing[3]

    # Relative position of the evaluation point,
    # i. e. interpolation weights.
    xd = (p[1] - x0) / nodes.grid.spacing[1]
    yd = (p[2] - y0) / nodes.grid.spacing[2]
    zd = (p[3] - z0) / nodes.grid.spacing[3]

    # Values at the closest interpolation points.
    c000 = nodes.data[ix  , iy  , iz  ]
    c001 = nodes.data[ix  , iy  , iz+1]
    c010 = nodes.data[ix  , iy+1, iz  ]
    c011 = nodes.data[ix  , iy+1, iz+1]
    c100 = nodes.data[ix+1, iy  , iz  ]
    c101 = nodes.data[ix+1, iy  , iz+1]
    c110 = nodes.data[ix+1, iy+1, iz  ]
    c111 = nodes.data[ix+1, iy+1, iz+1]

    # Interpolate along x.
    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    # Interpolate along y.
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd
    
    # Interpolate along z.
    c = c0 * (1 - zd) + c1 * zd
    
    return c
end


"""
    interpolate_linearly(ct::ScalarGrid, points)

Calculate the value at the given points for the given regular ScalarGrid.
This function uses multi-threading to speed up the calculations.
"""
function interpolate_linearly(nodes::ScalarGrid, points)
    results = Vector{eltype(points)}(undef, size(points, 2))
    @Threads.threads for i in 1:size(points, 2)
        @inbounds p = (points[1, i], points[2, i], points[3, i])
        @inbounds results[i] = Juliana.interpolate_linearly(nodes, p)
    end
    return results
end


function interpolate_linearly(nodes::Juliana.ScalarGrid, new_grid::Juliana.Grid)
    points, indices = Juliana.grid_to_points_and_indices(new_grid)
    values = convert.(Float32, Juliana.interpolate_linearly(nodes, points))
    values[isnan.(values)] .= 0
    return Juliana.ScalarGrid(
        Juliana.flat_vector_to_cube(
            new_grid,
            indices,
            values,
        ),
        new_grid,
    )
end
