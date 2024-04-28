# This file contains code for parsing patient data from PSIPlan, the legacy
# in-house treatment planning system used at PSI.
# The functions defined here are only useful if you need to access plan data.
# CTs, structure sets and dose distributions can be converted to DICOM files
# using the DicomConverter in-house tool.


function read_int(file)
    buffer = Vector{UInt8}(undef, 4)
    readbytes!(file, buffer)
    return only(reinterpret(Int32, buffer))
end

function read_int_array(file, n)
    buffer = Vector{UInt8}(undef, n * 4)
    readbytes!(file, buffer)
    return reinterpret(Int32, buffer)
end


function load_spot(file, data_scaling, field_center)
    buffer = Vector{UInt8}(undef, 3*4)

    readbytes!(file, buffer)
    spot_coordinates = reinterpret(Int32, buffer)
    spot_coordinates = spot_coordinates * (data_scaling * 100) # cm

    s, t, u = spot_coordinates
    
    buffer = Vector{UInt8}(undef, 4)
    readbytes!(file, buffer)
    min_dist_to_target_surface = only(reinterpret(Int32, buffer))
    
    readbytes!(file, buffer)
    n_absorbers = only(reinterpret(Int32, buffer))
    
    readbytes!(file, buffer)
    energy_lut_index = only(reinterpret(Int32, buffer))
    
    readbytes!(file, buffer)
    t_index = only(reinterpret(Int32, buffer))
    
    readbytes!(file, buffer)
    u_index = only(reinterpret(Int32, buffer))
    
    readbytes!(file, buffer)
    absorber_wed = only(reinterpret(Int32, buffer))
    absorber_wed *= data_scaling * 100
    
    readbytes!(file, buffer)
    residual_water_equivalent_range_after_preabsorber = only(reinterpret(Int32, buffer))
    residual_water_equivalent_range_after_preabsorber *= data_scaling * 100

    readbytes!(file, buffer)
    energy_MeV = only(reinterpret(Int32, buffer))
    energy_MeV *= data_scaling
    
    return s, t, u, energy_MeV, n_absorbers, absorber_wed
end


"""
    load_spli2_file(path, weights)

Load the information about a field from a .SPLI2 file. The function requires
that the spot weights have been loaded from the corresponding .BTOP file before.
"""
function load_spli2_file(path, weights; spot_ID_start=0)
    open(path) do file
        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        n_bytes = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        stu_grid_size = reinterpret(Int32, buffer)

        readbytes!(file, buffer)
        stu_grid_origin = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        stu_grid_spacing = reinterpret(Int32, buffer)

        # The data type is "not presently used"...
        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        data_type = reinterpret(Int32, buffer)

        readbytes!(file, buffer)
        data_size = reinterpret(Int32, buffer)

        readbytes!(file, buffer)
        header_units = only(reinterpret(Int32, buffer))
        scaling = (10.)^(-header_units)

        stu_grid_spacing = stu_grid_spacing .* (scaling * 100) # cm
        stu_grid_origin = stu_grid_origin .* (scaling * 100)

        readbytes!(file, buffer)
        data_units = only(reinterpret(Int32, buffer))
        data_scaling = (10.)^(-data_units)

        readbytes!(file, buffer)
        gantry_angle = only(reinterpret(Int32, buffer))
        gantry_angle *= data_scaling / 100 # degrees

        readbytes!(file, buffer)
        couch_angle = only(reinterpret(Int32, buffer))
        couch_angle *= data_scaling / 100 # degrees

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        field_center = reinterpret(Int32, buffer)
        field_center = field_center * data_scaling

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        n_spots = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        pool = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        version_num = only(reinterpret(Int32, buffer))

        readbytes!(file, buffer)
        max_range = only(reinterpret(Int32, buffer))

        readbytes!(file, buffer)
        min_range = only(reinterpret(Int32, buffer))

        readbytes!(file, buffer)
        nozzle_extraction = only(reinterpret(Int32, buffer))
        nozzle_extraction *= data_scaling / 100 # cm

        buffer = Array{UInt8, 1}(undef, 156)
        readbytes!(file, buffer)
        blank = String(buffer)

        buffer = Array{UInt8, 1}(undef, 256)
        readbytes!(file, buffer)
        comment = String(buffer)

        spots = Vector{Juliana.FionaStandalone.Spot}(undef, n_spots)
        for (i, w) in enumerate(weights)
            s, t, u, energy_MeV, n_absorbers, absorber_wed = load_spot(
                file,
                data_scaling,
                field_center,
            )
            spots[i] = Juliana.FionaStandalone.Spot(
                spot_ID_start + i-1,
                t,
                u,
                w,
                energy_MeV * 1000,
                n_absorbers,
            )
        end

        return gantry_angle, couch_angle, nozzle_extraction, field_center, spots
    end
end


"""
    load_btop_file(path)

Load the spot weights from a .BTOP file.
"""
function load_btop_file(path)
    open(path) do file
        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        n_bytes = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        stu_grid_size = reinterpret(Int32, buffer)

        readbytes!(file, buffer)
        stu_grid_origin = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        stu_grid_spacing = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        data_type = reinterpret(Int32, buffer)

        readbytes!(file, buffer)
        data_size = reinterpret(Int32, buffer)

        readbytes!(file, buffer)
        header_units = only(reinterpret(Int32, buffer))
        scaling = (10.)^(-header_units)

        stu_grid_spacing = stu_grid_spacing .* (scaling * 100) # cm
        stu_grid_origin = stu_grid_origin .* (scaling * 100)

        readbytes!(file, buffer)
        data_units = only(reinterpret(Int32, buffer))
        data_scaling = (10.)^(-data_units)

        readbytes!(file, buffer)
        gantry_angle = only(reinterpret(Int32, buffer))
        gantry_angle *= data_scaling / 100 # degrees

        readbytes!(file, buffer)
        couch_angle = only(reinterpret(Int32, buffer))
        couch_angle *= data_scaling / 100 # degrees

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        field_center = reinterpret(Int32, buffer)
        field_center = field_center * (data_scaling * 100) # cm

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        n_spots = only(reinterpret(Int32, buffer))

        buffer = Array{UInt8, 1}(undef, 176)
        readbytes!(file, buffer)
        blank = String(buffer)

        buffer = Array{UInt8, 1}(undef, 256)
        readbytes!(file, buffer)
        comment = String(buffer)

        # Parse actual spot weights.
        buffer = Array{UInt8, 1}(undef, 4*n_spots)
        readbytes!(file, buffer)
        spot_weights = reinterpret(Int32, buffer)
        spot_weights *= data_scaling
        
        return spot_weights
    end
end


"""
Example of a treatment definition file::

    TREATMENT_NUMBER: 0
    NUMBER_OF_SERIES: 2;
    CT_SET: 0;
    T_SET: 3;
    PLAN: 0;
    START_DOSE: 0.00000;
    END_DOSE: 54.0000;
    END_OF_SERIES:
    CT_SET: 0;
    T_SET: 4;
    PLAN: 0;
    START_DOSE: 54.0000;
    END_DOSE: 74.0000;
    END_OF_SERIES:
"""
function load_treatment_definition_file(treatment_file)
    content = String(read(open(treatment_file)))
    lines = split(content, "\n")
    
    line = lines[1]
    line = split(line, " ")
    @assert length(line) == 2
    @assert line[1] == "TREATMENT_NUMBER:"
    treatment_ID = parse(Int32, line[2])
    
    line = lines[2]
    line = split(line, " ")
    @assert length(line) == 2
    @assert line[1] == "NUMBER_OF_SERIES:"
    # end-1: Remove ";" at the end.
    n_series = parse(Int32, replace(line[2][1:end-1]))
    
    i = 3
    identifiers = Vector{Dict{String, Real}}(undef, n_series)
    for series in 1:n_series
        identifiers[series] = Dict{String, Real}()
        while true
            line = lines[i]
            if line == "END_OF_SERIES:"
                i += 1
                break
            end
            parts = split(line, " ")
            @assert length(parts) == 2
            
            # Remove the ":" at the end.
            name = parts[1][1:end-1]

            # Remove the ";" at the end.
            value_string = parts[2][1:end-1]
            value = tryparse(Int32, value_string)
            if isnothing(value)
                value = parse(Float32, value_string)
            end

            identifiers[series][name] = value
            i += 1
        end
    end
    return identifiers
end


"""
Example PLAN_INFO file:

      1.63636
       3
      1.00000     0.984000     0.338753     0.554324 F0
      1.00000     0.986000     0.338066     0.553199 F1
      1.00000     0.986000     0.338066     0.553199 F2
"""
function load_plan_def_file(plan_def_file)
    content = String(read(open(plan_def_file)))
    lines = split(content, "\n")[3:end]
    
    field_names = [split(line)[end] for line in lines if length(line) > 0]
    return field_names
end


function load_field(patient_data_directory, patient_ID, ct, t, f, field_ID; spot_ID_start=0)
    spot_weight_path = "$(patient_data_directory)/$(patient_ID)_CT$(ct)_T$(t)_$(f).BTOP"
    spot_list_path = "$(patient_data_directory)/$(patient_ID)_CT$(ct)_T$(t)_$(f).SPLI2"

    weights = Juliana.load_btop_file(spot_weight_path);
    gantry_angle, couch_angle, nozzle_extraction, field_center, spots = Juliana.load_spli2_file(
        spot_list_path,
        weights;
        spot_ID_start=spot_ID_start,
    )

    field = Juliana.FionaStandalone.FieldDefinition(
        f,
        field_ID,
        gantry_angle,
        couch_angle,
        nozzle_extraction,
        Dict(
            "x" => field_center[1],
            "y" => field_center[2],
            "z" => field_center[3],
        ),
        spots,
    )

    return field
end


function load_plan(patient_data_directory, patient_ID, ct, t, p, prescribed_dose)
    plan_def_file = "$(patient_data_directory)/$(patient_ID)_CT$(ct)_T$(t)_P$(p).PLAN_DEF"
    field_names = Juliana.load_plan_def_file(plan_def_file)
    
    fields = Vector{Juliana.FionaStandalone.FieldDefinition}(undef, length(field_names))
    spot_ID_start = 0
    for (i, f) in enumerate(sort(field_names))
        fields[i] = load_field(
            patient_data_directory,
            patient_ID,
            ct,
            t,
            f,
            i;
            spot_ID_start=spot_ID_start,
        )
        spot_ID_start += length(fields[i].spots)
    end
    return Juliana.FionaStandalone.TreatmentPlan(
        fields,
        prescribed_dose,
    )
end


"""
series_ID must start at 0!!!
"""
function load_treatment_series(patient_data_directory, patient_ID, treatment_ID, series_ID)
    treatment_file = "$(patient_data_directory)/$(patient_ID)_TR$(treatment_ID).TREATMENT_DEF"
    
    series_infos = Juliana.load_treatment_definition_file(treatment_file)
    # +1: Julia indices start at one.
    infos = series_infos[series_ID+1]
    
    ct = infos["CT_SET"]
    t = infos["T_SET"]
    p = infos["PLAN"]
    prescribed_dose = infos["END_DOSE"] - infos["START_DOSE"]
    
    return load_plan(
        patient_data_directory,
        patient_ID,
        ct,
        t,
        p,
        prescribed_dose,
    )
end


function load_treatment(patient_data_directory, patient_ID, treatment_ID)
    treatment_file = "$(patient_data_directory)/$(patient_ID)_TR$(treatment_ID).TREATMENT_DEF"
    series_infos = Juliana.load_treatment_definition_file(treatment_file)
    
    return [Juliana.load_treatment_series(patient_data_directory, patient_ID, treatment_ID, series_ID) for series_ID in 0:length(series_infos)-1]
end


function load_ctct_file(filename)
    open(filename) do file
        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        n_bytes = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        ct_size = reinterpret(Int32, buffer)

        # Not sure what the difference is to nx, ny, nz...
        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        n_voxels = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        origin_index = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        origin = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 3*4)
        readbytes!(file, buffer)
        spacing = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        data_type = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        n_bits_per_voxel = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        header_units = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        data_units = only(reinterpret(Int32, buffer))

        # ?
        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        n_files = only(reinterpret(Int32, buffer))

        # ?
        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        file_num = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        couch = only(reinterpret(Int32, buffer))

        # Unused values.
        buffer = Vector{UInt8}(undef, 4*4)
        readbytes!(file, buffer)

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        patient_weight = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        extra_weight = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        system_weight = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        ct_type = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        gantry_type = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        study_ID = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        series_ID = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 4)
        readbytes!(file, buffer)
        couch_offset_y = only(reinterpret(Int32, buffer))

        buffer = Vector{UInt8}(undef, 116)
        readbytes!(file, buffer)
        blank = reinterpret(Int32, buffer)

        buffer = Vector{UInt8}(undef, 256)
        readbytes!(file, buffer)
        comment = reinterpret(Int32, buffer)

        if n_bits_per_voxel != 16
            @warn "Ignoring that n_bits_per_voxel != 16"
        end

        buffer = Vector{UInt8}(undef, prod(ct_size)*2)
        readbytes!(file, buffer)
        data = reinterpret(Int16, buffer)
        data = reshape(data, collect(ct_size)...)
        data = data[1:end, end:-1:1, end:-1:1]
        
        buffer = Vector{UInt8}(undef, prod(ct_size)*2)
        readbytes!(file, buffer)

        # Apply the scaling factors and build the ScalarGrid.
        header_scaling = (10f0)^(-header_units + 2)
        data_scaling = (10f0)^(-data_units)
        
        origin = origin .* header_scaling   # cm
        spacing = spacing .* header_scaling # cm
        data = convert.(Int16, (data .- 1000))
        
        return Juliana.ScalarGrid(
            data,
            Juliana.Grid(
                collect(spacing),
                collect(origin),
                collect(ct_size),
            ),
        )
    end
end


function load_dat2_file(path)
    open(path) do file
        # Read the header fields.
        n_bytes                  = read_int(file)
        number_of_curves         = read_int(file)
        max_points_per_curve     = read_int(file)
        depth_separation         = read_int(file)
        water_equ_range_peak_sep = read_int(file)
        data_type                = read_int(file) # "Not presently used."
        n_bits_per_block         = read_int(file)
        header_units             = read_int(file)
        data_units               = read_int(file)
        pool                     = read_int(file) # 1: no preabsorber, 2: preabsorber, 3: auto
        pool_version             = read_int(file)
        max_range                = read_int(file)
        min_range                = read_int(file)
        nozzle_extension         = read_int(file)

        # Blank.
        buffer = Vector{UInt8}(undef, 200)
        readbytes!(file, buffer)

        # Comment.
        buffer = Vector{UInt8}(undef, 256)
        readbytes!(file, buffer)

        # Read the actual payload.
        number_of_absorbers = read_int_array(file, number_of_curves)
        absorber_wer        = read_int_array(file, number_of_curves)
        ranges              = read_int_array(file, number_of_curves)
        energies            = read_int_array(file, number_of_curves) # Nominal energies.
        depths              = read_int_array(file, max_points_per_curve)

        buffer = Vector{UInt8}(undef, number_of_curves * max_points_per_curve * 4)
        readbytes!(file, buffer)
        dose_values = reinterpret(Int32, buffer)
        dose_values = collect(reshape(dose_values, (max_points_per_curve, number_of_curves))')

        buffer = Vector{UInt8}(undef, number_of_curves * max_points_per_curve * 4)
        readbytes!(file, buffer)
        sigma_values = reinterpret(Int32, buffer)
        sigma_values = collect(reshape(dose_values, (max_points_per_curve, number_of_curves))');

        # Convert to physical units.
        header_scaling = (10f0)^(-header_units)

        max_range        = max_range        * header_scaling # cm
        min_range        = min_range        * header_scaling # cm
        nozzle_extension = nozzle_extension * header_scaling # cm
        absorber_wer     = absorber_wer     * header_scaling # cm
        ranges           = ranges          .* header_scaling # cm
        energies         = energies        .* header_scaling # MeV
        depths           = depths          .* header_scaling # cm
        dose_values      = dose_values     .* header_scaling # nGy
        sigma_values     = sigma_values    .* header_scaling # cm

        # Build the structures to be used in JulianA.
        x0 = Vector{Float32}(undef, number_of_curves)
        fill!(x0, depths[1])
        dx = depths[2:end] .- depths[1:end-1]

        indices = sortperm(energies)

        depth_dose_curves = Juliana.LookUpTable(
            energies[indices] * 1000,
            x0[indices],
            dx[indices],
            dose_values[indices, :],
        )
        sigma_values = Juliana.LookUpTable(
            energies[indices] * 1000,
            x0[indices],
            dx[indices],
            sigma_values[indices, :],
        )

        return depth_dose_curves, sigma_values, nozzle_extension
    end
end
