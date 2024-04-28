using CSV
using DataFrames
using JSON
using Logging
using NPZ
using StaticArrays


CACHE_DIR::String = get(ENV, "JULIANA_CACHE_DIR", "$(homedir())/juliana_cache")
mkpath(CACHE_DIR)


function parse_plan_file(path)
    plan = JSON.parsefile(path)
    fields = plan["fields"]
    return BeamArrangement(
        convert.(Float32, [f["gantry_angle"] for f in fields]),
        convert.(Float32, [f["couch_angle"] for f in fields]),
        convert.(Float32, [f["nozzle_extraction"] for f in fields]),
    )
end


"""
    load_patient_data(data_dir::String,
                      patient_ID::String;
                      series::Integer=0) -> Juliana.PatientData
"""
function load_patient_data(data_dir::String, patient_ID::String; series::Integer=0)
    ct_path = "$(data_dir)/CTs/$(patient_ID)_$series.dat"
    ct = Juliana.load_ct_dat_file(ct_path);

    structure_names = Juliana.structure_names(data_dir, patient_ID)
    structures_to_load = [name for name in structure_names if !(name in ["AUTO_FULL_BODY", "BODY"])]
    structures = Juliana.load_structures(
        data_dir,
        patient_ID,
        ct.grid,
        which=structures_to_load,
    );
    for (name, structure) in structures
        structures[name] = Structure(
            structure.name,
            structure.points[:, 1:3],
            structure.mask,
            structure.distanceFromStructure,
            structure.grid,
            structure.is_target,
        )
    end

    prescriptions = Juliana.load_prescriptions(data_dir, patient_ID, structures);

    return ct_path, PatientData(
        ct,
        structures,
        prescriptions,
    )
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

        grid = Grid(
            SVector{3}(spacing),
            SVector{3}(origin),
            SVector{3}(collect(data_size)),
        )
        return ScalarGrid(data, grid)
    end
end


function write_ct_dat_file(path::String, ct::Juliana.ScalarGrid)
    @assert eltype(typeof(ct.data)) == Int16
    open(path, "w") do file
        write(file, reinterpret(UInt8, b"\x48\x46\x53")) #HFS
        write(file, reinterpret(UInt8, convert.(Float32, collect(ct.grid.origin))))
        write(file, reinterpret(UInt8, collect(ct.grid.spacing)))
        write(file, reinterpret(UInt8, collect(convert.(Int32, size(ct.data)))))
        write(file, reinterpret(UInt8, vec(ct.data)));
    end
end


function load_ct_dat_file(path::String)
    return load_dat_file(Int16, path, true)
end


function write_dose_dat_file(path::String, dose::Juliana.ScalarGrid)
    open(path, "w") do file
        write(file, reinterpret(UInt8, convert.(Float32, dose.grid.origin)))
        write(file, reinterpret(UInt8, dose.grid.spacing))
        write(file, reinterpret(UInt8, collect(convert.(Int32, size(dose.data)))))
        write(file, reinterpret(UInt8, vec(dose.data)));
    end
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
    if size(points, 2) == 3
        new = Array{Float32, 2}(undef, size(points, 1), 4)
        new[:, 1:3] .= points
        contours = split_by_z(points);
        i = 1
        index = 1
        for contour in contours
            new[i:i+size(contour, 1)'-1, 4] .= index
            i += size(contour, 1)
            index += 1
        end
        points = new
    end

    contours = split_by_contour_index(points)
    distance_mask = calculate_distance_mask(
        Tuple(grid.size),
        grid.spacing,
        contours,
    )
    binary_mask = to_binary_mask(distance_mask, max_distance=0.5*minimum(grid.spacing))
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
function parse_oar_constraints_file(path, target_dose, structures; T=Float32)
    df = DataFrame(CSV.File(path));

    constraints = Array{Constraint{T}, 1}(undef, 0)
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

    return constraints
end


function load_prescriptions(data_dir, patient_ID, structures)
    target_doses = load_target_doses(data_dir, patient_ID)
    target_dose = maximum(values(target_doses))
    target_doses = [(name, dose) for (name, dose) in target_doses]

    T = typeof(target_dose)

    constraints_path = "$(data_dir)/prescriptions/$(patient_ID)_constraints.csv"
    constraints = parse_oar_constraints_file(
        constraints_path,
        target_dose,
        structures;
        T=T,
    )
    
    return Prescriptions{T}(target_doses, constraints)
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
    elseif dose == "unknown"
        dose = convert(T, NaN)
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
        kind = constraint_dvh_v
        # Remove the "V" at the beginning.
        dose = parse_dose(constraint_quantity[2:end], target_dose)
        volume = parse_volume(constraint_threshold, structure_volume)
    else
        kind = constraint_unknown_type
        dose = convert(T, NaN)
        volume = convert(T, NaN)
    end

    constraint = Constraint{T}(
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
        SVector{3}(convert(Vector{Float32}, grid_dict["spacing"])),
        SVector{3}(convert(Vector{Float32}, grid_dict["origin"])),
        SVector{3}(convert(Vector{Int32}, grid_dict["size"])),
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


"""
    cached_structures(cache_dir::String, patient_ID::String; series=0)

Return a set trhat contains the names of all structures that are in the cache.
This function does NOT check whether the cache entry is complete!
"""
function cached_structures(cache_dir::String, patient_ID::String; series=0)
    if !ispath("$(cache_dir)/$(patient_ID)")
        return Set{String}()
    end
    filenames = readdir("$(cache_dir)/$(patient_ID)")
    filenames = [f for f in filenames if startswith(f, string(series))]
    # Remove the "_binary_mask.npy" and so on, then concatenate.
    structure_names = Set{String}()
    to_delete = (
        "_binary_mask.npy",
        "_distanceFromStructure.npy",
        "_grid.json",
        "_isTarget.txt",
        "_points.npy",
    )
    for filename in filenames
        name = replace(filename, "$(series)_" => "", count=1)
        for s in to_delete
            name = replace(name, s => "")
        end
        push!(structure_names, name)
    end
    return structure_names
end


function write_plan_config(path::String, config::TreatmentPlan)
    open(path, "w") do file
        write(file, JSON.json(config, 4))
    end
end


function read_plan_file(path::String)
    open(path, "r") do file
        plan_dict = JSON.parse(file)
        return TreatmentPlan(plan_dict)
    end
end
