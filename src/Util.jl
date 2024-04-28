using CUDA
using DataFrames
using JSON
using LinearAlgebra
using StaticArrays


function mask_to_points(grid::Grid, mask::AbstractArray{T, 3}) where {T}
    points = Array{Float32, 2}(undef, (3, convert(Int64, sum(mask))));

    i = 1
    for iz in 1:grid.size[3]
        for iy in 1:grid.size[2]
            for ix in 1:grid.size[1]
                if mask[ix, iy, iz] > 0
                    points[:, i] .= grid.origin .+ [ix-1, iy-1, iz-1] .* grid.spacing
                    i += 1
                end
            end
        end
    end

    @assert sum(isnan.(points)) == 0
    return points
end


function mask_to_points_and_indices(grid::Grid, mask::AbstractArray{T, 3}) where {T}
    N = convert(Int64, sum(mask))
    points = Array{Float32, 2}(undef, (3, N))
    indices = Array{Int64, 2}(undef, (3, N))

    i = 1
    for iz in 1:grid.size[3]
        for iy in 1:grid.size[2]
            for ix in 1:grid.size[1]
                if mask[ix, iy, iz] > 0
                    points[:, i] .= grid.origin .+ (ix-1, iy-1, iz-1) .* grid.spacing
                    indices[:, i] .= (ix, iy, iz)
                    i += 1
                end
            end
        end
    end

    @assert sum(isnan.(points)) == 0
    return points, indices
end


function flat_vector_to_cube(grid, indices, values)
    cube = zeros(grid.size[1], grid.size[2], grid.size[3])
    for (value, index) in zip(values, eachcol(indices))
        cube[index[1], index[2], index[3]] = value
    end
    return cube
end


"""
    angles_to_direction(gantry, couch)

Convert angles to a direction vector pointing from zero towards the nozzle
(unit length).

The angles are expected to be in degree.
The implementation is taken from Fiona, SiddonXyzKernel (in gpu/dosecalc).
Look for xStart, yStart, zStart.
"""
function angles_to_direction(gantry::Number, couch::Number)
    # These factors are assuming HFS orientation.
    # Taken from the Fiona code.
    factors = [1.f0, -1.f0, -1.f0]

    direction = factors .* [
        sind(gantry) * cosd(couch),
        cosd(gantry),
        -sind(gantry) * sind(couch),
    ]
    direction ./= norm(direction)
    return direction
end


function angles_to_direction(gantry::AbstractVector{T}, couch::AbstractVector{T}) where {T<:Number}
    @assert size(gantry) == size(couch)

    directions = Matrix{Float32}(undef, 3, length(gantry))
    for (i, (gantry_angle, couch_angle)) in enumerate(zip(gantry, couch))
        directions[:, i] .= angles_to_direction(gantry_angle, couch_angle)
    end

    return directions
end


function calculate_whole_target_mask(prescriptions, structures)
    shape = size(structures[prescriptions.target_doses[1][1]].mask)
    target_mask = zeros(Bool, shape)

    for (name, _) in prescriptions.target_doses
        target_mask .= target_mask .|| structures[name].mask
    end

    return target_mask
end


function calculate_whole_target_distance(prescriptions, structures)
    return min.([structures[name].distanceFromStructure for (name, dose) in prescriptions.target_doses]...)
end


function build_constraint_oar_mask(prescriptions, structures)
    mask = copy(structures[collect(keys(structures))[1]].mask)
    fill!(mask, false)
    for constr in prescriptions.constraints
        mask .= mask .|| structures[constr.structure_name].mask
    end
    return mask
end


function index_is_valid(ix, iy, iz, shape)
    return (ix >= 1)        && (iy >= 1)        && (iz >= 1) && 
           (ix <= shape[1]) && (iy <= shape[2]) && (iz <= shape[3])
end


function build_checker_board_mask(grid)
    checker_board = zeros(Float32, Tuple(grid.size))
    
    first_in_z_is_white = true
    for iz in 1:grid.size[3]
        first_in_z_is_white = !first_in_z_is_white
        first_in_y_is_white = first_in_z_is_white
        for iy in 1:grid.size[2]
            is_white = first_in_y_is_white
            for ix in 1:grid.size[1]
                checker_board[ix, iy, iz] = is_white
                is_white = !is_white
            end
            first_in_y_is_white = !first_in_y_is_white
        end
    end
    return checker_board
end


function build_checker_board_mask_shifted(grid, n)
    checker_board = zeros(Float32, Tuple(grid.size))
    
    first_in_z_index = 0
    for iz in 1:grid.size[3]
        first_in_z_index = (first_in_z_index + 1) % (n+1)
        first_in_y_index = first_in_z_index
        for iy in 1:grid.size[2]
            is_white = first_in_y_index
            for ix in 1:grid.size[1]
                checker_board[ix, iy, iz] = is_white == n
                is_white = (is_white + 1) % (n+1)
            end
            first_in_y_index = (first_in_y_index + 1) % (n+1)
        end
    end
    return checker_board
end


function build_checker_board_mask(grid, n_inslice, n_slices)
    checker_board = zeros(Float32, Tuple(grid.size))
    for iz in 1:grid.size[3]
        for iy in 1:grid.size[2]
            for ix in 1:grid.size[1]
                checker_board[ix, iy, iz] = (mod1(ix, n_inslice) == 1) && (mod1(iy, n_inslice) == 1) && (mod1(iz, n_slices) == 1)
            end
        end
    end
    return checker_board
end


function shift_to_zero(ct::Juliana.ScalarGrid,
                       structures::Dict{String, Juliana.Structure},
                       dose_distributions;
                       calculate_distances::Bool=true)
    origin = ct.grid.origin
    new_ct = Juliana.ScalarGrid(
        ct.data,
        Juliana.Grid(ct.grid.spacing, SVector{3}(zeros(Float32, 3)), ct.grid.size),
    )

    new_structures = Dict{String, Juliana.Structure}()
    for (name, structure) in structures
        new_points = copy(structure.points)
        new_points[:, 1:3] .= (new_points[:, 1:3]' .- origin)'

        distance = copy(structure.distanceFromStructure)
        mask = copy(structure.mask)
        if calculate_distances
            distance .= calculate_distance_mask(
                Tuple(new_ct.grid.size),
                new_ct.grid.spacing,
                Juliana.split_by_contour_index(new_points),
            )
            mask .= Juliana.to_binary_mask(distance)
        end

        new_structures[name] = Juliana.Structure(
            structure.name,
            new_points,
            mask,
            distance,
            new_ct.grid,
            structure.is_target,
        )
    end

    new_dose_distributions = [Juliana.ScalarGrid(dose.data, new_ct.grid) for dose in dose_distributions]

    return new_ct, new_structures, new_dose_distributions
end


"""
    xyz_to_index(point_xyz, grid)

Convert XYZ coordinates to 3D grid indices.
"""
function xyz_to_index(point_xyz, grid)
    # +1: Indices in Julia start at 1.
    # Ex.: The origin is zero, spacing 1.
    #      The point (0.4, 0.4, 0.4) should have index (1, 1, 1).
    return convert.(Int32, round.((point_xyz .- grid.origin) ./ grid.spacing) .+ 1)
end


"""
    xyz_to_index(x, y, z, grid)

Convert XYZ coordinates to 3D grid indices.
"""
function xyz_to_index(x, y, z, grid)
    # +1: Indices in Julia start at 1.
    # Ex.: The origin is zero, spacing 1.
    #      The point (0.4, 0.4, 0.4) should have index (1, 1, 1).
    return (
        convert(Int32, round((x - grid.origin[1]) / grid.spacing[1]) + 1),
        convert(Int32, round((y - grid.origin[2]) / grid.spacing[2]) + 1),
        convert(Int32, round((z - grid.origin[3]) / grid.spacing[3]) + 1),
    )
end


"""
    index_to_xyz(index, grid::Juliana.Grid)

Convert 3D grid indices to XYZ coordinates.
"""
function index_to_xyz(index, grid::Juliana.Grid)
    # -1: Indices in Julia start at 1.
    return grid.origin .+ grid.spacing .* (index .- 1)
end


"""
    volume_at_indices(volume::AbstractArray{T, 3}, indices::AbstractArray) where {T}

Create a flat array containing the values of volume at the given indices.
"""
function volume_at_indices(volume::AbstractArray{T, 3}, indices::AbstractArray) where {T}
    values = zeros(T, size(indices, 2))
    for (i, index) in enumerate(eachcol(indices))
        values[i] = volume[index[1], index[2], index[3]]
    end
    return values
end


function count_inslice_neighbours(mask)
    cnt = Array{Int32, 3}(undef, size(mask, 1), size(mask, 2), size(mask, 3))
    fill!(cnt, 0)
    for iz in 1:size(mask, 3)
        for iy in 2:size(mask, 2)-1
            for ix in 2:size(mask, 1)-1
                cnt[ix, iy, iz] = mask[ix-1, iy, iz] + mask[ix+1, iy, iz] + mask[ix, iy-1, iz] + mask[ix, iy+1, iz]
            end
        end
    end
    return cnt
end


function get_margin_mask(structure)
    return (0 .< count_inslice_neighbours(structure.mask) .< 4) .&& structure.mask
end


function get_axis_aligned_bounding_box_grid(points,
                                            grid::Juliana.Grid,
                                            extension::Real=0)
    @assert all(grid.origin .== 0)

    pmin = minimum(points[:, 1:3], dims=1) .- extension
    pmax = maximum(points[:, 1:3], dims=1) .+ extension

    # Flatten the arrays to be vectors.
    pmin = vec(pmin)
    pmax = vec(pmax)

    # Align the extremum points to the grid.
    pmin = floor.(pmin ./ grid.spacing) .* grid.spacing
    pmax = ceil.(pmax ./ grid.spacing) .* grid.spacing

    shape = convert.(Integer, round.((pmax .- pmin) ./ grid.spacing))

    return Juliana.Grid(
        SVector{3}(grid.spacing),
        SVector{3}(pmin),
        SVector{3}(shape),
    )
end


# expand by 2.2cm because the optimisation uses a 2cm margin by default
# to place optimisation points and we want to be on the safe side
function calculate_distance_in_neighbourhood(structure, grid; expansion=2.2f0)
    points = structure.points

    # Sometimes (For the whole body and patient surface contours for G3 data),
    # some contour points lie half a voxel spacing outside of the grid.
    # We clamp to avoid these corner cases.
    points[:, 1:3] .= clamp.(points[:, 1:3], 0, Inf) 

    bounding_grid = Juliana.get_axis_aligned_bounding_box_grid(
        points,
        grid,
        expansion,
    )

    pmin = bounding_grid.origin
    pmax = bounding_grid.origin .+ bounding_grid.spacing .* bounding_grid.size

    start = convert.(Integer, floor.(pmin ./ bounding_grid.spacing))
    stop = convert.(Integer, ceil.(pmax ./ bounding_grid.spacing));

    index_range = [
        start[1] start[2]
         stop[1]  stop[2]
    ]

    return Juliana.calculate_distance_mask(
        Tuple(grid.size),
        grid.spacing,
        points,
        index_range,
    )
end


function add_distance_and_mask_in_neighbourhood(structure::Juliana.Structure,
                                                inside_distance::Real;
                                                expansion=2.2f0)
    distance_mask = calculate_distance_in_neighbourhood(
        structure,
        structure.grid,
        expansion=expansion,
    )
    binary_mask = Juliana.to_binary_mask(
        distance_mask,
        max_distance=inside_distance,
    )
    return Juliana.Structure(
        structure.name,
        structure.points,
        binary_mask,
        distance_mask,
        structure.grid,
        structure.is_target,
    )
end


function count_8_neighbours(img)
    new = Matrix{Integer}(undef, size(img))
    fill!(new, 0)
    for j in 2:size(img, 2)-1
        for i in 2:size(img, 1)-1
            new[i, j] = sum(img[i-1:i+1, j-1:j+1])
        end
    end
    return new
end

function erode(mask)
    new = convert.(Bool, similar(mask))
    for iz in 1:size(mask, 3)
        n_neighbours = count_8_neighbours(mask[:, :, iz])
        new[:, :, iz] .= n_neighbours .>= 8
    end
    return new
end

function dilate(mask)
    new = convert.(Bool, similar(mask))
    for iz in 1:size(mask, 3)
        n_neighbours = count_8_neighbours(mask[:, :, iz])
        new[:, :, iz] .= n_neighbours .> 1
    end
    return new
end


function split_contours_by_long_lines(structure; SPLIT_DISTANCE=3)
    # SPLIT_DISTANCE is in cm.
    new_points = copy(structure.points)

    contours = Juliana.split_by_contour_index(new_points)
    contour_lengths = [size(contour, 1) for contour in contours]

    for (contour_index, contour) in enumerate(contours)
        N = size(contour, 1)
        contour_distances = Vector{Float32}(undef, N-1)
        for i in 2:N
            a = contour[i  , 1:2]
            b = contour[i-1, 1:2]
            dist = √(sum((a .- b).^2))
            contour_distances[i-1] = dist
        end
    
        for (i, dist) in enumerate(contour_distances)
            if dist >= SPLIT_DISTANCE
                update_start_index = sum(contour_lengths[1:contour_index-1]) + i + 1
                new_points[update_start_index:end, 4] .+= 1
            end
        end
    end

    s_split = Juliana.Structure(
        structure.name,
        new_points,
        zeros(structure.grid.size...),
        zeros(structure.grid.size...),
        structure.grid,
        structure.is_target,
    )

    distance_mask = Juliana.calculate_distance_in_neighbourhood(
            s_split,
            structure.grid;
            expansion=2.2f0,
    )

    return Juliana.Structure(
        structure.name,
        new_points,
        Juliana.to_binary_mask(
            distance_mask,
            max_distance=0.5*minimum(structure.grid.spacing),
        ),
        distance_mask,
        structure.grid,
        structure.is_target,
    )
end


function hu_to_sp_factory(path_to_conversion_file)
    open(path_to_conversion_file) do file
        hu_to_sp_dict = JSON.parse(file)

        @assert hu_to_sp_dict["dz"] == 1
        @assert hu_to_sp_dict["z0"] == -1000
        
        dz = convert(Float32, hu_to_sp_dict["dz"])
        z0 = convert(Float32, hu_to_sp_dict["z0"])
        densities = convert.(Float32, hu_to_sp_dict["densities"]);

        function convert_to_sp(value)
            sp_index = convert(
                Int64,
                round((value - z0) / dz),
            ) + 1
            return densities[sp_index]
        end
        
        return convert_to_sp
    end
end


function get_range(depth_dose_curves::Juliana.LookUpTable, energy::Real)
    i = argmin(abs.(energy .- depth_dose_curves.energies))
    table = depth_dose_curves.table[i, :]
    table = table[.!isnan.(table)]
    return argmax(table) * depth_dose_curves.dx[i] - depth_dose_curves.x0[i]
end


function local_minima_indices(vector::AbstractVector{T}, window_size::Integer) where {T}
    @assert window_size % 2 == 1
    halfwindow = ÷(window_size, 2)
    minima = Vector{Integer}(undef, 0)
    for i in halfwindow+1:length(vector)-halfwindow
        if all((vector[i-halfwindow:i-1] .> vector[i]) .&& (vector[i+1:i+halfwindow] .> vector[i]))
            push!(minima, i)
        end
    end
    return minima
end


function get_fwhm(depth_dose_curves::Juliana.LookUpTable, energy::Real)
    i = argmin(abs.(energy .- depth_dose_curves.energies))

    # Get half maximum value.
    dose_values = copy(depth_dose_curves.table[i, :])
    dose_values[isnan.(dose_values)] .= -Inf;
    max = maximum(dose_values)
    half_max = 0.5 * max

    # Get the two values closest to the half maximum.
    difference_to_halfmax = abs.(dose_values .- half_max)
    local_minima = local_minima_indices(difference_to_halfmax, 5)
    @assert length(local_minima) == 2
    i_start, i_stop = local_minima

    # Take the difference to obtain the FWHM (full width at half maximum).
    x_start = depth_dose_curves.x0[i] + i_start * depth_dose_curves.dx[i]
    x_stop = depth_dose_curves.x0[i] + i_stop * depth_dose_curves.dx[i]

    return x_stop - x_start
end


function conformity_index(dose_distribution::AbstractArray,
                          target::Juliana.Structure,
                          target_dose::Real)
    return sum(dose_distribution .>= target_dose) / sum(target.mask)
end


function conformity_index(dose_distribution::Juliana.ScalarGrid,
                          target::Juliana.Structure,
                          target_dose::Real)
    return conformity_index(dose_distribution, target, target_dose)
end


"""
    reproducible_sparse_mv(A, x)

Implements a bitwise reproducible sparse matrix-vector (SPV) multiplication on
an NVIDIA GPU. The default SPV is NOT reproducible bitwise and therefore
unsuitable for a medical context.

WARNING: Matrix A must be in CSR format!!!
"""
function reproducible_sparse_mv(A::CUDA.CUSPARSE.CuSparseMatrixCSR{T, Int32}, x::AbstractVector{T}) where {T<:Number}
    Y = cu(zeros(eltype(A), size(A, 1)))
    CUDA.CUSPARSE.mv!('N', 1, A, x, 0, Y, 'O', CUDA.CUSPARSE.CUSPARSE_SPMV_CSR_ALG2)
    return Y
end


function reproducible_sparse_mv(A::AbstractMatrix{T}, x::AbstractVector{T}) where {T<:Number}
    return A * x
end


function grid_to_points(grid::Juliana.Grid)
    N = prod(grid.size)
    points = Matrix{eltype(grid.spacing)}(undef, 3, N)
    i = 1
    for ix in 1:grid.size[1]
        for iy in 1:grid.size[2]
            for iz in 1:grid.size[3]
                points[:, i] .= grid.origin .+ grid.spacing .* (ix-1, iy-1, iz-1)
                i += 1
            end
        end
    end
    return points
end


function grid_to_points_and_indices(grid::Juliana.Grid)
    N = prod(grid.size)
    points = Matrix{eltype(grid.spacing)}(undef, 3, N)
    indices = Matrix{eltype(grid.size)}(undef, 3, N)
    i = 1
    for ix in 1:grid.size[1]
        for iy in 1:grid.size[2]
            for iz in 1:grid.size[3]
                points[:, i] .= grid.origin .+ grid.spacing .* (ix-1, iy-1, iz-1)
                indices[:, i] .= (ix, iy, iz)
                i += 1
            end
        end
    end
    return points, indices
end


function expand_structure(s::Juliana.Structure, expansion)
    points = Matrix{Missing}(missing, size(s.points, 1), size(s.points, 2))
    return Juliana.Structure(
        s.name,
        points,
        s.distanceFromStructure .<= expansion,
        s.distanceFromStructure,
        s.grid,
        s.is_target,
    )
end


function spot_weights(plan)
    spots = vcat([DataFrame(field.spots) for field in plan.fields]...);
    # Ensure that the spots are ordered correctly.
    @assert spots.id == collect(0:size(spots, 1)-1)
    return convert.(Float32, spots.weight)
end


function get_fiona_dose_calc_grid(ct_grid, resolution=0.35f0)
    @assert ct_grid.origin == zeros(3)
    p_max = ct_grid.spacing .* ct_grid.size

    spacing = SVector{3}(resolution, resolution, resolution)

    return Grid(
        spacing,
        # +0.01: To be consistent with the way in which the dose calculation grid is
        # defined in FIonA.
        SVector{3}(0f0, 0f0, 0.01f0),
        # floor(): ensure the dose calc grid lies within the CT grid.
        SVector{3}(Base.convert.(Int32, floor.(p_max ./ spacing))),
    )
end


function update_spot_weights(plan::TreatmentPlan{T}, w::Vector{T}) where {T<:Real}
    n_fields = length(plan.fields)
    n_spots_per_field = [length(f.spots) for f in plan.fields]

    @assert sum(n_spots_per_field) == length(w)

    # Make sure that the spots are ordered by field.
    last_id = -1
    for field in plan.fields
        spot_IDs = [spot.id for spot in field.spots];
        @assert spot_IDs[1] == (last_id + 1)
        @assert sort(spot_IDs) == spot_IDs
        last_id = spot_IDs[end]
    end

    # Create a new TreatmentPlan object for the current spot weights.
    new_fields = Array{FieldDefinition{T}, 1}(undef, n_fields)
    id = 0
    for (i_field, field) in enumerate(plan.fields)
        new_spots = Array{Spot{T}, 1}(undef, length(field.spots))
        for (i_spot, spot) in enumerate(field.spots)
            new_spots[i_spot] = Spot{T}(
                spot.id,
                spot.t,
                spot.u,
                w[id+1],
                spot.energykeV,
                spot.numberOfAbsorbers,
            )
            id += 1
        end
        new_fields[i_field] = FieldDefinition{T}(
            field.label,
            field.id,
            field.gantryAngle,
            field.couchAngle,
            field.nozzleExtraction,
            field.fieldCenter,
            new_spots,
        )
    end
    new_plan = TreatmentPlan(new_fields, plan.prescribedDose)
    return new_plan
end
