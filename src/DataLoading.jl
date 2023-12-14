using CSV
using DataFrames
using JSON
using Logging
using NPZ


CACHE_DIR::String = get(ENV, "JULIANA_CACHE_DIR", "$(homedir())/juliana_cache")
mkpath(CACHE_DIR)


struct Grid{S, T}
    spacing::Array{S, 1}
    origin::Array{S, 1}
    size::Array{T, 1}
end


struct ScalarGrid{A<:AbstractArray, S, T}
    data::A
    grid::Grid{S, T}
end


struct Structure{A<:AbstractArray, B<:AbstractArray, C<:AbstractArray, S, T}
    name::String
    points::A
    mask::B
    # This is actually only the in-slice distance. We do not consider z!
    distanceFromStructure::C
    grid::Grid{S, T}
    is_target::Bool
end


@enum ConstraintType constraint_dvh_d constraint_dvh_v constraint_mean constraint_max constraint_unknown_type
@enum ConstraintDirection upper_constraint lower_constraint
@enum ConstraintPriority hard soft unknown_priority


struct Constraint{T}
    structure_name::String
    kind::ConstraintType
    dose::T
    volume::Union{Nothing, T}
    priority::ConstraintPriority
    direction::ConstraintDirection
end


struct Prescriptions{T}
    target_doses::Vector{Tuple{String, T}}
    constraints::Array{Constraint, 1}
end


# CTs and dose distributions.
function load_dat_file(T, path::String, read_orientation::Bool)
    open(path) do file
        if read_orientation
            buffer = Array{UInt8, 1}(undef, 3)
            readbytes!(file, buffer)
            orientation = String(buffer)
        end

        buffer = Array{UInt8, 1}(undef, 3*4)
        readbytes!(file, buffer)
        origin = collect(reinterpret(Float32, buffer))

        buffer = Array{UInt8, 1}(undef, 3*4)
        readbytes!(file, buffer)
        spacing = collect(reinterpret(Float32, buffer))

        buffer = Array{UInt8, 1}(undef, 3*4)
        readbytes!(file, buffer)
        data_size = tuple(collect(reinterpret(Int32, buffer))...)
        n = data_size[1] * data_size[2] * data_size[3]

        buffer = Array{UInt8, 1}(undef, n*sizeof(T))
        readbytes!(file, buffer)
        data = collect(reinterpret(T, buffer))
        data = reshape(data, data_size)

        grid = Grid(spacing, origin, collect(data_size))
        return ScalarGrid(data, grid)
    end
end


function load_ct_dat_file(path::String)
    return load_dat_file(Int16, path, true)
end


function load_dose_dat_file(path::String)
    return load_dat_file(Float32, path, false)
end


# Structures.
function to_binary_mask(distance_mask; max_distance=0.f0)
    return (distance_mask .<= max_distance)
end

function to_binary_mask(contours, grid; max_distance=0.f0)
    distance_mask = calculate_distance_mask(
        Tuple(grid.size),
        grid.spacing,
        contours,
    )
    return to_binary_mask(distance_mask; max_distance=max_distance)
end


function load_npy_structure(name::String, path::String, grid::Grid{S, T}, is_target::Bool) where {S, T}
    points = Array{Float32, 2}(npzread(path))
    contours = split_by_z(points);
    distance_mask = calculate_distance_mask(
        Tuple(grid.size),
        grid.spacing,
        contours,
    )
    binary_mask = to_binary_mask(distance_mask)
    return Structure(name, points, binary_mask, distance_mask, grid, is_target)
end


function structure_names(data_dir::String, patient_ID::String; series=0)
    structure_dir = "$(data_dir)/structures/$(patient_ID)/$(series)"
    file_names = readdir(structure_dir)
    npy_file_names = [f for f in file_names if endswith(f, ".npy")]
    return [f[1:end-4] for f in npy_file_names]
end


function load_structures(data_dir, patient_ID, mask_grid; series=0, which::Union{Nothing, Array{String, 1}}=nothing, target_names::Union{Nothing, Array{String}}=nothing, cache_dir=CACHE_DIR)
    if isnothing(which)
        names = structure_names(data_dir, patient_ID, series=series)
    else
        names = which
    end

    structure_dir = "$(data_dir)/structures/$(patient_ID)/$(series)"
    structures = Dict{String, Structure}()

    for (i, name) in enumerate(names)
        if is_cached(cache_dir, patient_ID, name, series=series)
            @debug "Reading structure $name from cache"
            structures[name] = read_from_cache(
                cache_dir,
                patient_ID,
                name,
                series=series,
            )
        else
            @debug "Reading structure $name from a .npy file"
            path = "$structure_dir/$name.npy"
            if isnothing(target_names)
                is_target = false
            else
                is_target = name in target_names
            end
            structure = load_npy_structure(name, path, mask_grid, is_target)
            structures[structure.name] = structure
            write_to_cache(cache_dir, patient_ID, structure, series=series)
        end
    end
    return structures
end


# Prescriptions.
function load_prescriptions(data_dir, patient_ID, structures)
    target_doses = load_target_doses(data_dir, patient_ID)
    target_dose = maximum(values(target_doses))
    target_doses = [(name, dose) for (name, dose) in target_doses]

    constraints_path = "$(data_dir)/prescriptions/$(patient_ID)_constraints.csv"
    df = DataFrame(CSV.File(constraints_path));

    constraints = Array{Constraint, 1}(undef, 0)
    for row in eachrow(df)
        if !(row.structure_name in keys(structures))
            @warn "Could not find a structure named $(row.structure_name)"
            continue
        end

        structure = structures[row.structure_name]
        spacing = structure.grid.spacing
        voxel_volume = spacing[1] * spacing[2] * spacing[3]
        structure_volume = sum(structure.mask) * voxel_volume

        soft_hard = row.soft_hard
        if ismissing(soft_hard)
            soft_hard = ""
        end

        constraint = build_constraint(
            row.structure_name,
            soft_hard,
            row.constraint_quantity,
            row.constraint_threshold,
            target_dose,
            structure_volume,
        )
        push!(constraints, constraint)
    end
    return Prescriptions(target_doses, constraints)
end


function load_target_doses(data_dir, patient_ID)
    target_doses_path = "$(data_dir)/dose_per_target.csv"
    target_dose_table = DataFrame(CSV.File(target_doses_path))
    patient_targets = target_dose_table[
        in.(target_dose_table.patient_ID, Ref([patient_ID])),
        ["target_name", "total_dose_gy"]
    ]
    target_doses = Dict{String, Float32}(eachrow(patient_targets))
    return target_doses
end


function is_percentage(query::AbstractString)
    return !isnothing(match(r"^\d+\.?\d*%$", query))
end

function parse_dose(dose::AbstractString, target_dose_gy::T) where {T}
    dose = lowercase(dose)
    dose = replace(dose, " " => "")
    if endswith(dose, "gy")
        dose = dose[1:end-2]
        dose = parse(T, dose)
    elseif is_percentage(dose)
        dose = replace(dose, "%" => "")
        dose = parse(T, dose) * target_dose_gy / convert(T, 100) 
    else
        error("Could not parse dose string '$dose'")
    end
    return dose
end

function parse_volume(volume::AbstractString, structure_volume::T) where {T}
    number_only_pattern = r"\d+\Z"
    if is_percentage(volume)
        return parse(T, replace(volume, "%" => "")) / convert(T, 100)
    elseif !isnothing(match(number_only_pattern, volume))
        return parse(T, volume) / convert(T, 100)
    elseif endswith(lowercase(volume), "cc")
        return parse(T, replace(volume, "cc" => "")) / structure_volume
    else
        error("Malformatted volume string '$volume'")
    end
end

function build_constraint(structure_name::AbstractString,
                          soft_hard::AbstractString,
                          constraint_quantity::AbstractString,
                          constraint_threshold::AbstractString,
                          target_dose::T,
                          structure_volume::T) where {T}
    dvh_D_pattern = r"^D\d+(\.\d+)?(cc)?(%)?\Z"
    dvh_V_pattern = r"^V\d+(\.\d+)?(gy|%)?\Z"

    if !ismissing(soft_hard) && (soft_hard == "hard")
        priority = hard
    elseif !ismissing(soft_hard) && (soft_hard == "soft")
        priority = soft
    else
        priority = unknown_priority
    end

    if constraint_quantity == "D_mean"
        kind = constraint_mean
        dose = parse_dose(constraint_threshold, target_dose)
        volume = nothing
    elseif constraint_quantity == "D_max"
        kind = constraint_max
        dose = parse_dose(constraint_threshold, target_dose)
        volume = nothing
    elseif !isnothing(match(dvh_D_pattern, constraint_quantity))
        # D<volume> < threshold
        kind = constraint_dvh_d
        # Remove the "D" at the beginning.
        dose = parse_dose(constraint_threshold, target_dose)
        volume = parse_volume(constraint_quantity[2:end], structure_volume)
    elseif !isnothing(match(dvh_V_pattern, constraint_quantity))
        # V<dose> < threshold
        kind = consqtraint_dvh_v
        # Remove the "V" at the beginning.
        dose = parse_dose(constraint_quantity[2:end], target_dose)
        volume = parse_volume(threshold, structure_volume)
    else
        error("Could not parse constraint $constraint")
    end

    constraint = Constraint(
        String(structure_name),
        kind,
        dose,
        volume,
        priority,
        upper_constraint,
    )
    
    return constraint
end


function hottest_target(prescriptions)
    hottest_name, hottest_dose = prescriptions.target_doses[1]
    for (name, dose) in prescriptions.target_doses[2:end]
        if dose > hottest_dose
            hottest_name = name
            hottest_dose = dose
        end
    end
    return hottest_name, hottest_dose
end


function coldest_target(prescriptions)
    coldest_name, coldest_dose = prescriptions.target_doses[1]
    for (name, dose) in prescriptions.target_doses[2:end]
        if dose < coldest_dose
            coldest_name = name
            coldest_dose = dose
        end
    end
    return coldest_name, coldest_dose
end


# Cache system to speedup structure loading.
function write_to_cache(cache_dir::String,
                        patient_ID::String,
                        structure::Structure; series=0)
    @debug "Writing structure $(structure.name), $(patient_ID) to cache $(cache_dir)..."

    patient_cache_dir = "$cache_dir/$(patient_ID)"
    mkpath(patient_cache_dir)
    basename = "$patient_cache_dir/$(series)_$(structure.name)"

    npzwrite("$(basename)_points.npy", structure.points)
    npzwrite("$(basename)_binary_mask.npy", structure.mask)
    npzwrite("$(basename)_distanceFromStructure.npy", structure.distanceFromStructure)
    
    open("$(basename)_grid.json", "w") do file
        write(file, JSON.json(structure.grid, 4))
    end
    
    open("$(basename)_isTarget.txt", "w") do file
        write(file, convert(String, structure.is_target ? "true" : "false"))
    end
    return nothing
end

function is_cached(cache_dir::String,
                   patient_ID::String,
                   structure_name::String; series=0)
    patient_cache_dir = "$cache_dir/$(patient_ID)"
    basename = "$patient_cache_dir/$(series)_$(structure_name)"

    return isfile("$(basename)_points.npy") &&
        isfile("$(basename)_binary_mask.npy") &&
        isfile("$(basename)_distanceFromStructure.npy") &&
        isfile("$(basename)_grid.json") &&
        isfile("$(basename)_isTarget.txt")
end


function read_from_cache(cache_dir::String,
                         patient_ID::String,
                         structure_name::String; series=0)
    patient_cache_dir = "$cache_dir/$(patient_ID)"
    basename = "$patient_cache_dir/$(series)_$(structure_name)"
    
    points = npzread("$(basename)_points.npy")
    mask = npzread("$(basename)_binary_mask.npy")
    distances = npzread("$(basename)_distanceFromStructure.npy")
    
    mask = BitArray{3}(mask)
    
    grid_dict = JSON.parsefile("$(basename)_grid.json")
    grid = Grid(
        convert(Vector{Float32}, grid_dict["spacing"]),
        convert(Vector{Float32}, grid_dict["origin"]),
        convert(Vector{Int32}, grid_dict["size"]),
    )
    is_target = false
    open("$(basename)_isTarget.txt", "r") do file
        text = readline(file)
        if text == "true"
            is_target = true
        end
    end
    return Structure(
        structure_name,
        points,
        mask,
        distances,
        grid,
        is_target,
    )
end


function xyz_to_index(point_xyz, grid) where {T}
    # +1: Indices in Julia start at 1.
    # Ex.: The origin is zero, spacing 1.
    #      The point (0.4, 0.4, 0.4) should have index (1, 1, 1).
    return convert.(Int32, round.((point_xyz .- grid.origin) ./ grid.spacing) .+ 1)
end

# function xyz_to_index(points, grid)
#     return collect(hcat([xyz_to_index(point, grid) for point in eachrow(points)]...)')
# end

function index_to_xyz(index, grid) where {T}
    # -1: Indices in Julia start at 1.
    return grid.origin .+ grid.spacing .* (index .- 1)
end

# function index_to_xyz(indices, grid)
#     return collect(hcat([index_to_xyz(index, grid) for index in eachrow(indices)]...)')
# end
