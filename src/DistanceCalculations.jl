using LinearAlgebra


function split_by_z(points)
    contours = []
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


function shortest_distance_to_line(x::AbstractArray{T, 1}, a::AbstractArray{T, 1}, b::AbstractArray{T, 1}) where {T}
    d_ab = b - a
    d_xa = x - a
    
    # TODO: Why squared???
    t = dot(d_xa, d_ab) / norm(d_ab)^convert(T, 2)
    
    if t < zero(T)
        return norm(x - a)
    elseif t < one(T)
        x_proj = a + t * d_ab
        return norm(x - x_proj)
    else
        return norm(x - b)
    end
end


function shortest_distance_to_polygon(x::AbstractArray{T, 1}, contour::AbstractArray{T, 2}) where {T}
    N = size(contour, 2)
    minimum_distance = Inf
    for i in 1:N
        p0 = contour[1:2, i]
        p1 = contour[1:2, mod1(i+1, N)]

        d = shortest_distance_to_line(x[1:2], p0, p1)
        minimum_distance = min(minimum_distance, d)
    end

    return minimum_distance
end


"""
Size of contour: (3, N)
"""
function contour_line_intersections(line_value::T, contour::AbstractArray{T, 2}, y_spacing::T) where {T}
    N = size(contour, 2)
    intersection_points = Array{Array{T, 1}}(undef, 0)
    for i in 1:N
        p0 = contour[1:2, i]
        p1 = contour[1:2, mod1(i+1, N)]

        # Avoid division by zero for parallel lines.
        if (abs(p0[2] - p1[2]) <= 1e-6)
            # Axis parallel line.
            if abs(p0[2] - line_value) <= 0.5*y_spacing
                # The line is close to the scan line.
                push!(intersection_points, p0)
                push!(intersection_points, p1)
            end
            continue
        end

        t = (line_value - p0[2]) / (p1[2] - p0[2])

        if (0. <= t) && (t <= 1.)
            push!(intersection_points, p0 + t * (p1 - p0))
        end
    end

    return intersection_points
end


"""
function calculate_distance_mask_for_contour(shape, spacing, contour, x_indices, y_indices)

    Iterate over all voxels and calculate the minimum distance to the contour.
    The even-odd rule determines the sign.

Size of contour: (3, N)
"""
function calculate_distance_mask_for_contour(shape, spacing, contour, x_indices, y_indices)
    distance = zeros((length(x_indices), length(y_indices)))

    for (nj, j) in enumerate(y_indices)
        # -1: Julia indices start at 1
        y = (j-1)*spacing[2]
        intersection_points = contour_line_intersections(y, contour, spacing[2])

        for (ni, i) in enumerate(x_indices)
            # -1: Julia indices start at 1
            x = (i-1)*spacing[1]

            crossed_intersections = [p for p in intersection_points if p[1] <= x]
            is_inside = length(crossed_intersections) % 2 == 1

            p = [x, y]
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
Calculate the distance from the contours on the given grid (assuming origin zero).

shape: Tuple{Integer, Integer, Integer}
spacing: Array{Real} of size (3, )
contours: Vector{Matrix{Real, 2}} whose elements are 3D point coordinates of shape (N, 3)
"""
function calculate_distance_mask(shape, spacing, contours)
    mask = Array{Float32, 3}(undef, shape)
    fill!(mask, Inf)
    for contour in contours
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
        mask[:, :, iz] .= slice_mask
    end
    return mask
end
