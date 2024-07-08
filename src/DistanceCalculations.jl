using Base.Threads
using LinearAlgebra
using ProgressBars
using StaticArrays


"""
    split_by_z(points)

Split an array of points by their z coordinates.

The shape of points is (N, 3). The column order is x, y, z.
"""
function split_by_z(points::Matrix{T}) where {T}
    contours = Vector{Matrix{T}}(undef, 0)
    start_of_contour_index = 1
    for i in 1:size(points, 1)
        if (points[i, 3] ≈ points[start_of_contour_index, 3])
            continue
        else
            push!(contours, points[start_of_contour_index:i-1, :])
            start_of_contour_index = i
        end
    end
    push!(contours, points[start_of_contour_index:end, :])
    return contours
end


"""
Return a vector of contour groups that share the common z coordinate.
"""
function group_by_slice(contours::Vector{Matrix{T}}) where {T}
    groups = Vector{Vector{Matrix{T}}}(undef, 0)

    current_group = Vector{Matrix{T}}(undef, 0)
    current_z = contours[1][1, 3]

    for contour in contours
        if contour[1, 3] != current_z
            # Store the current group.
            push!(groups, current_group)
            # Start a new group.
            current_group = Vector{Matrix{T}}(undef, 1)
            current_group[1] = contour
            current_z = contour[1, 3]
        else
            # Extend the current group.
            push!(current_group, contour)
        end
    end
    # Store the final group.
    push!(groups, current_group)
    
    return groups
end


"""
    split_by_contour_index(points)

Split an array of points by their contour index.

The shape of points is (N, 4). The column order is x, y, z, contour_index.
"""
function split_by_contour_index(points::AbstractArray{T, 2}) where {T}
    contours = Vector{Matrix{T}}(undef, 0)
    start_of_contour_index = 1
    for i in 1:size(points, 1)
        if (points[i, 4] == points[start_of_contour_index, 4])
            continue
        else
            push!(contours, points[start_of_contour_index:i-1, :])
            start_of_contour_index = i
        end
    end
    push!(contours, points[start_of_contour_index:end, :])
    return contours
end


function shortest_distance_to_line(x, a, b)
    d_ab = b .- a
    d_xa = x .- a
    
    # TODO: Why squared???
    t = dot(d_xa, d_ab) / norm(d_ab)^2
    
    if t < 0
        return norm(x .- a)
    elseif t < 1
        x_proj = a .+ t .* d_ab
        return norm(x .- x_proj)
    else
        return norm(x .- b)
    end
end


"""
    shortest_distance_to_polygon(x, contour)

Calculate the shortest in-slice distance of point x to the contour.

Shape of contour: (d, N), where d is 2 or 3.
If d = 3, dimension 3 will be ignored.
The in-slice distance is the distance in the xy plane.
Both inside and outside points will have positive distance values!
"""
function shortest_distance_to_polygon(x, contour::AbstractArray{T, 2}) where {T}
    N = size(contour, 2)
    minimum_distance = Inf
    for i in 1:N
        p0 = (contour[1, i], contour[2, i])
        p1 = (contour[1, mod1(i+1, N)], contour[2, mod1(i+1, N)])

        d = shortest_distance_to_line(
            x,
            p0,
            p1,
        )
        minimum_distance = min(minimum_distance, d)
    end

    return minimum_distance
end


"""
    contour_line_intersections(line_value, contour, y_spacing)

Calculate the intersection points between the contour and the curve y=line_value.

Size of contour: (3, N)
"""
function contour_line_intersections(line_value::T, contour::AbstractArray{T, 2}, y_spacing::T) where {T}
    N = size(contour, 2)
    intersection_points = Vector{Tuple{T, T}}()
    for i in 1:N
        p0 = (contour[1, i], contour[2, i])
        p1 = (contour[1, mod1(i+1, N)], contour[2, mod1(i+1, N)])

        # Intersection conting rules:
        # https://web.archive.org/web/20130126163405/http://geomalgorithms.com/a03-_inclusion.html
        # 1.) an upward edge includes its starting endpoint, and excludes its final endpoint;
        # 2.) a downward edge excludes its starting endpoint, and includes its final endpoint;
        # 3.) horizontal edges are excluded
        # 4.) the edge-ray intersection point must be strictly right of the point P.
        intersects_p0 = abs(p0[2] - line_value) < 1e-6
        intersects_p1 = abs(p1[2] - line_value) < 1e-6

        is_upward = p1[2] > p0[2]
        is_downward = p1[2] < p0[2]
        is_horizontal = (!is_upward) && (!is_downward)
        
        if is_horizontal
            continue
        end

        if intersects_p0
            if is_upward
                push!(intersection_points, p0)
            end
            continue
        end

        if intersects_p1
            if is_downward
                push!(intersection_points, p1)
            end
            continue
        end

        t = (line_value - p0[2]) / (p1[2] - p0[2])
        if (0. <= t) && (t <= 1.)
            push!(intersection_points, p0 .+ t .* (p1 .- p0))
        end
    end

    return intersection_points
end


"""
    calculate_distance_mask_for_contour(shape, spacing, contour, x_indices, y_indices)

Iterate over all voxels and calculate the in-slice minimum distance to the contour.
The even-odd rule determines the sign.

Size of contour: (3, N)
"""
function calculate_distance_mask_for_contour(shape, spacing, contour, x_indices, y_indices)
    distance = zeros(Float32, (length(x_indices), length(y_indices)))

    for (nj, j) in enumerate(y_indices)
        # -1: Julia indices start at 1
        y = (j-1)*spacing[2]
        intersection_points = contour_line_intersections(y, contour, spacing[2])

        for (ni, i) in enumerate(x_indices)
            # -1: Julia indices start at 1
            x = (i-1)*spacing[1]

            is_inside = false
            for p in intersection_points
                if p[1] <= x
                    is_inside = !is_inside
                end
            end

            p = (x, y)
            d = shortest_distance_to_polygon(p, contour)
            if is_inside
                d = -d
            end
            distance[ni, nj] = d
        end
    end
    return distance
end


"""
    calculate_distance_mask(shape, spacing, contours)

Calculate the in-slice distance from the contours on the given grid (assuming origin zero).

shape: Tuple{Integer, Integer, Integer}
spacing: Array{Real} of size (3, )
contours: Vector{Matrix{Real, 2}} whose elements are 3D point coordinates of shape (N, 3)
"""
function calculate_distance_mask(shape, spacing::AbstractArray{T, 1}, contours) where {T}
    mask = Array{T, 3}(undef, shape)
    fill!(mask, Inf)

    l = Base.ReentrantLock()
    Threads.@threads for contour in contours
        z = contour[1, 3]
        # +1: Julia indices start at one
        iz = convert(Int64, round(z / spacing[3]) + 1)
        slice_mask = calculate_distance_mask_for_contour(
            shape[1:2],
            spacing,
            contour',
            collect(1:shape[1]),
            collect(1:shape[2]),
        )
        # Need to lock because multiple contours might share the same
        # z-coordinate.
        lock(l)
        try
            mask[:, :, iz] .= min.(mask[:, :, iz], slice_mask)
        finally
            unlock(l)
        end
    end
    return mask
end


"""
    calculate_distance_mask(shape, spacing, contours, index_ranges)

Calculate the in-slice distance from the contours on the given grid (assuming origin zero).

shape: Tuple{Integer, Integer, Integer}
spacing: Array{Real} of size (3, )
points: Vector{Matrix{Real}} whose elements are point coordinates of shape (N, 4)
    (x, y, z, contour_index)
index_ranges: Array{Real} of size (2, 2) where index_range(i, j) is the lower (1)
                and upper (2) range of indices of axis j for which to calculate the distances
"""
function calculate_distance_mask(shape, spacing::AbstractArray{T, 1}, points, index_ranges) where {T}
    mask = Array{T, 3}(undef, shape)
    fill!(mask, Inf)

    x_indices = collect(index_ranges[1, 1]:index_ranges[2, 1])
    y_indices = collect(index_ranges[1, 2]:index_ranges[2, 2])

    contours = split_by_contour_index(points)
    # contour_groups = group_by_slice(contours)
    
    l = Base.ReentrantLock()
    Threads.@threads for contour in contours
        z = contour[1, 3]
        # +1: Julia indices start at one
        iz = convert(Int64, round(z / spacing[3]) + 1)
        slice_mask = calculate_distance_mask_for_contour(
            shape[1:2],
            spacing,
            contour',
            x_indices,
            y_indices,
        )
        # Need to lock because multiple contours might share the same
        # z-coordinate.
        lock(l)
        try
            mask[x_indices, y_indices, iz] .= min.(
                mask[x_indices, y_indices, iz],
                slice_mask,
            )
        finally
            unlock(l)
        end
    end
    return mask
end


# The code below is to calculate the full 3D distance (not only in-slice)
# from a structure, but only outside of the structure.
# This requires that a binary mask has been obtained in advance.
struct Rectangle{T <: Real}
    p_start::Tuple{T, T, T}
    extent::Tuple{T, T, T}
end


"""
    build_rectangles(structure)

Build a vector of rectangles that cover the binary mask of the `structure`.

This is is an exact (up to numerical errors) representation of the binary
mask, assuming piecewise constant interpolation (i. e. voxels are squares).
"""
function build_rectangles(structure::Juliana.Structure)
    grid = structure.grid
    @assert grid.origin == zeros(3)
    spacing = (grid.spacing[1], grid.spacing[2], grid.spacing[3])

    # + 1: Julia indices start at one.
    ix_min, iy_min, iz_min = convert.(
        Int32,
        floor.(vec(minimum(structure.points[:, 1:3], dims=1)) ./ spacing),
    ) .+ 1
    ix_max, iy_max, iz_max = convert.(
        Int32,
        ceil.(vec(maximum(structure.points[:, 1:3], dims=1)) ./ spacing),
    ) .+ 1;

    rectangles = Vector{Rectangle{Float32}}(undef, 0)

    for iz in iz_min:iz_max
        mask = structure.mask[:, :, iz];
        n_rectangles_in_slice = 0

        i_start = (-1f0, -1f0, -1f0)
        is_inside = false
        for iy in iy_min:iy_max
            for ix in ix_min:ix_max+1
                if mask[ix, iy]
                    if is_inside
                        continue
                    else
                        i_start = convert.(Float32, (ix, iy, iz))
                        is_inside = true
                    end
                else
                    if is_inside
                        # Close the rectangle and add it to the list.
                        i_stop = (ix-1, iy-1, iz)
                        r = Rectangle(
                            # -1: Julia indices start at 1.
                            # - 0.5: The voxel position indicates its center,
                            # so the surface is shifted by half a voxel spacing.
                            (i_start .- 1 .- 0.5f0) .* spacing,
                            # + 1: Left and right faces of the rectangle.
                            (i_stop .- i_start .+ 1) .* spacing,
                        )
                        push!(rectangles, r)
                        is_inside = false
                        n_rectangles_in_slice += 1
                    end
                end
            end
            @assert !is_inside
        end
    end
    
    return rectangles
end


"""
    distance_point_face(p, x_min, x_max, y_min, y_max, z)

Calculate the unsigned shortest distance between the rectangle spanned by
[`x_min`, `x_max`] × [`y_min`, `y_max`] at "height" `z` and point `p`.
"""
function distance_point_face(p::Tuple{T, T, T}, x_min::T, x_max::T, y_min::T, y_max::T, z::T) where {T <: Real}
    in_x_range = x_min <= p[1] <= x_max
    in_y_range = y_min <= p[2] <= y_max
    
    Δx = min(abs(p[1] - x_min), abs(p[1] - x_max))
    Δy = min(abs(p[2] - y_min), abs(p[2] - y_max))
    Δz = abs(p[3] - z)

    if in_x_range && in_y_range
        return Δz
    elseif in_x_range
        return √((Δy)^2 + (Δz)^2)
    elseif in_y_range
        return √((Δx)^2 + (Δz)^2)
    else
        return √((Δx)^2 + (Δy)^2 + (Δz)^2)
    end
    return 0
end


"""
    distance_point_rectangle(point, rectangle)

Calculate the shortest distance between the `point` and the `rectangle`.
"""
function distance_point_rectangle(p::Tuple{T, T, T}, r::Rectangle{T}) where {T <: Real}
    p_min = r.p_start
    p_max = r.p_start .+ r.extent

    # xy plane, z fix.
    distance_to_inferior_face = distance_point_face(
        (p[1], p[2], p[3]),
        p_min[1], p_max[1],
        p_min[2], p_max[2],
        p_min[3],
    )
    distance_to_superior_face = distance_point_face(
        (p[1], p[2], p[3]),
        p_min[1], p_max[1],
        p_min[2], p_max[2],
        p_max[3],
    )
    # xz plane, y fix.
    distance_to_anterior_face = distance_point_face(
        (p[1], p[3], p[2]),
        p_min[1], p_max[1],
        p_min[3], p_max[3],
        p_min[2],
    )
    distance_to_posterior_face = distance_point_face(
        (p[1], p[3], p[2]),
        p_min[1], p_max[1],
        p_min[3], p_max[3],
        p_max[2], 
    )
    # yz plane, x fix.
    distance_to_right_face = distance_point_face(
        (p[2], p[3], p[1]),
        p_min[2], p_max[2],
        p_min[3], p_max[3],
        p_min[1], 
    )
    distance_to_left_face = distance_point_face(
        (p[2], p[3], p[1]),
        p_min[2], p_max[2],
        p_min[3], p_max[3],
        p_max[1], 
    )
    
    return min(
        distance_to_inferior_face,
        distance_to_superior_face,
        distance_to_anterior_face,
        distance_to_posterior_face,
        distance_to_right_face,
        distance_to_left_face,
    )
end


"""
    distance_point_rectangles(p, rectangles)

Calculate the shortest unsigned distance between point `p` and any of the `rectangles`.
"""
function distance_point_rectangles(p::Tuple{T, T, T}, rectangles::Vector{Rectangle{T}}) where {T <: Real}
    distances = Vector{T}(undef, length(rectangles))
    for (i, r) in enumerate(rectangles)
        distances[i] = distance_point_rectangle(p, r)
    end
    return minimum(distances)
end


# expand by 2.2cm because the optimisation uses a 2cm margin by default
# to place optimisation points and we want to be on the safe side
"""
    calculate_outside_distance_3D(structure[, expansion=2.2f0])

Calculate the 3D distance of any voxel within distance `expansion` (in L-∞ norm)
to the structure. Points outside this distance will have distance `Inf`.

# Why is this function needed?

If the binary masks have been calculated based on the in-slice
distance. The masks are correct, but the distance cube is wrong and leads to
discontinuities in three dimensions.
"""
function calculate_outside_distance_3D(structure, expansion=2.2f0)
    distances_new = zeros(eltype(structure.distanceFromStructure), size(structure.distanceFromStructure)...)
    fill!(distances_new, Inf)

    @assert all(structure.grid.origin .== (0f0, 0f0, 0f0))
    spacing = (
        structure.grid.spacing[1],
        structure.grid.spacing[2],
        structure.grid.spacing[3],
    )
    Δ = convert.(Int32, ceil.(expansion ./ spacing))
    # + 1: Julia indices start at one.
    ix_min, iy_min, iz_min = convert.(
        Int32,
        floor.(vec(minimum(structure.points[:, 1:3], dims=1)) ./ spacing),
    ) .+ 1 .- Δ
    ix_min = clamp(ix_min, 1, structure.grid.size[1])
    iy_min = clamp(iy_min, 1, structure.grid.size[2])
    iz_min = clamp(iz_min, 1, structure.grid.size[3])

    ix_max, iy_max, iz_max = convert.(
        Int32,
        ceil.(vec(maximum(structure.points[:, 1:3], dims=1)) ./ spacing),
    ) .+ 1 .+ Δ
    ix_max = clamp(ix_max, 1, structure.grid.size[1])
    iy_max = clamp(iy_max, 1, structure.grid.size[2])
    iz_max = clamp(iz_max, 1, structure.grid.size[3])

    rectangles = build_rectangles(structure)

    Threads.@threads for iz in ProgressBar(iz_min:iz_max)
        for iy in iy_min:iy_max
            for ix in ix_min:ix_max
                if structure.mask[ix, iy, iz]
                    continue
                end
                distances_new[ix, iy, iz] = distance_point_rectangles(
                    (
                        (ix-1) * spacing[1],
                        (iy-1) * spacing[2],
                        (iz-1) * spacing[3],
                    ),
                    rectangles,
                )
            end
        end
    end
    
    return distances_new
end


"""
    outside_distance_from_structure_3D_consistent(structure[, expansion=2.2f0])

Make sure that the distance from the given `structure` (outside of the contours)
is consistent in 3D, i. e. do not use the in-slice distance.
"""
function outside_distance_from_structure_3D_consistent(structure::Juliana.Structure, expansion=2.2f0)
    new_distance = calculate_outside_distance_3D(structure, expansion)
    new_distance[structure.mask] .= 0

    return Juliana.Structure(
        structure.name,
        structure.points,
        structure.mask,
        new_distance,
        structure.grid,
        structure.is_target,
    )
end
