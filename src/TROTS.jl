using Downloads
using ZipFile
using MAT
using StatsBase
using Printf



"""
    transform_angles_to_G2_range(gantry_angle, couch_angle)

Input ranges:
gantry_angle ∈ [0, 360]
couch_angle ∈ [0, 360]

Output ranges:
gantry_angle ∈ [-30, 180]
couch_angle ∈ [-180, 180]
"""
function transform_angles_to_G2_range(gantry_angle, couch_angle)
    @assert 0 <= gantry_angle && gantry_angle <= 360
    @assert 0 <= couch_angle && couch_angle <= 360
    if !(0 <= gantry_angle && gantry_angle <= 180)
        gantry_angle = 360 - gantry_angle
        couch_angle = couch_angle + 180
    end
    if couch_angle > 180
        couch_angle = couch_angle - 360
    end
    
    return gantry_angle, couch_angle
end



function load_TROTS_data(patient_ID="Protons1", tmp_dir=nothing)
    # Load only these structures to speedup the loading.
    all_structures_to_load = Dict(
        "Protons1" => [
            "CTV High",
            "CTV Intermediate 10 mm",
            "CTV Low Shrunk 10 mm",
            "CTV Combined Ring 0-10 mm",
            "Brainstem",
            "Oesophagus",
            "Oral Cavity",
            "Parotid (L)",
            "Parotid (R)",
            "SMG (L)",
            "SMG (R)",
            "Spinal Cord",
        ],
    )

    if !(patient_ID in keys(all_structures_to_load))
        @error ("Beam arrangement for patient ID $(patient_ID) not ready; try one of these: $(keys(all_structures_to_load))")
    end

    structures_to_load = all_structures_to_load[patient_ID]

    # Download the data.
    URL = "https://files.sebastiaanbreedveld.nl/TROTS/TROTS_in_DICOM.zip"
    PREFIX = "DICOMs/$(patient_ID)/"

    if isnothing(tmp_dir)
        tmp_dir = mktempdir()
    end

    file = Downloads.download(URL, IOBuffer());
    reader = ZipFile.Reader(file);

    # Save DICOM files to a temporary directory.
    ct_files = []
    for f in reader.files
        if startswith(f.name, PREFIX) && f.name != PREFIX
            base_name = replace(f.name, "DICOMs/$(patient_ID)/" => "")

            # Adjust file names to convention.
            if startswith(base_name, "CTSlice")
                i = parse(Int64, replace(base_name, "CTSlice" => "", ".dcm" => ""))
                base_name = "CT.$(i).dcm"
            elseif base_name == "structs.dcm"
                base_name = "RS.dcm"
            else
                continue
            end

            open("$(tmp_dir)/$(base_name)", "w") do file
                write(file, read(f))
            end

            if startswith(base_name, "CT")
                push!(ct_files, base_name)
            end
        end
    end
    sort!(
        ct_files,
        by=n -> parse(Int64, replace(n, "CT." => "", ".dcm" => "")),
    )
    ct_paths = ["$(tmp_dir)/$f" for f in ct_files]

    # Remove the plan.csv file in order to not confuse load_dicom_directory().
    # rm("$(tmp_dir)/planning.csv")

    # Load the temporary DICOM files.
    # ct, structures = load_dicom_directory(tmp_dir; structure_names=structures_to_load)
    ct = Juliana.read_dicom_ct(ct_paths);
    structures = Juliana.read_dicom_structureset("$(tmp_dir)/RS.dcm", ct.grid)

    # Shift structures to zero because JulianA expects it for most functions.
    structures = Dict{String, Juliana.Structure}(
        name => structures[name] for name in structures_to_load
    )
    ct, structures = Juliana.shift_to_zero(ct, structures, []; calculate_distances=false)

    # Shift the CT to have HU >= -1000 instead of -1024.
    ct = Juliana.ScalarGrid(
        ct.data .+ convert(Int16, 24),
        ct.grid,
    )

    # Round z-spacing to 2 digits to avoid issues with floating point precision for DICOM exports.
    # Does not influence z-coordinate by more than 0.01cm.
    @assert ct.grid.origin == zeros(3)
    z_spacing_new = round(ct.grid.spacing[3], digits=2)

    for name in keys(structures)
        s = structures[name]
        z_new = round.(s.points[:, 3] ./ ct.grid.spacing[3]) .* z_spacing_new
        @assert maximum(abs.(s.points[:, 3] .- z_new)) < 0.01
        s.points[:, 3] = z_new
    end

    grid_new = Juliana.Grid(
        SVector{3}(ct.grid.spacing[1], ct.grid.spacing[2], z_spacing_new),
        ct.grid.origin,
        ct.grid.size,
    )
    ct = Juliana.ScalarGrid(ct.data, grid_new);

    # Calculate the distance masks only after the rounding.
    ct, structures = Juliana.shift_to_zero(ct, structures, []; calculate_distances=true)

    # Assign prescriptions.
    presc = nothing
    if patient_ID == "Protons1"
        @warn("WARNING: The parotid constraints aim to spare both, but the guideline demands only one.")
        presc = Juliana.Prescriptions(
            # Target values are taken from the original papers:
            # https://doi.org/10.1088/0031-9155/58/19/6969
            # https://doi.org/10.1016/j.ijrobp.2015.01.031
            [
                ("CTV High", 66.f0),
                ("CTV Intermediate 10 mm", 66.f0),
                ("CTV Low Shrunk 10 mm", 54.f0),
            ],
            # OAR constraints are based on international guidelines.
            # https://doi.org/10.1016/j.ijrobp.2019.06.2540
            [
                Juliana.Constraint("Brainstem", Juliana.constraint_max, 54.f0, 0.f0, Juliana.hard, Juliana.upper_constraint),
                Juliana.Constraint("Spinal Cord", Juliana.constraint_max, 50.f0, 0.f0, Juliana.hard, Juliana.upper_constraint),
                Juliana.Constraint("Parotid (L)", Juliana.constraint_mean, 26.f0, 0.f0, Juliana.hard, Juliana.upper_constraint),
                Juliana.Constraint("Parotid (R)", Juliana.constraint_mean, 26.f0, 0.f0, Juliana.hard, Juliana.upper_constraint),
                Juliana.Constraint("Oesophagus", Juliana.constraint_mean, 45.f0, 0.f0, Juliana.hard, Juliana.upper_constraint),
                Juliana.Constraint("Oral Cavity", Juliana.constraint_mean, 40.f0, 0.f0, Juliana.hard, Juliana.upper_constraint),
                Juliana.Constraint("SMG (L)", Juliana.constraint_mean, 35.f0, 0.f0, Juliana.hard, Juliana.upper_constraint),
                Juliana.Constraint("SMG (R)", Juliana.constraint_mean, 35.f0, 0.f0, Juliana.hard, Juliana.upper_constraint),
            ],
        );
    else
        @warn("Prescription for patient ID $(patient_ID) not ready; try one of these: \"Protons1\"")
    end
    
    # Download .mat file for the given patient.
    # TODO: Find a way to speed up the download...
    # URL = "https://files.sebastiaanbreedveld.nl/TROTS/Protons.zip"
    # file = Downloads.download(URL, IOBuffer())
    # reader = ZipFile.Reader(file)

    # i = parse(Int64, replace(patient_ID, "Protons" => ""))
    # n = @sprintf "%02d" i
    # filename = "Protons_$(n).mat"

    # matlab_data = nothing
    # for f in reader.files
    #     if f.name == filename
    #         path = "$(tmp_dir)/$(f.name)"
    #         open(path, "w") do file
    #             write(file, read(f))
    #         end
    #         matlab_data = matread(path)
    #         break
    #     end
    # end
    
    # # Build field definitions and field targets.
    # beams = matlab_data["patient"]["Beams"]["BeamConfig"]

    # gantry_angles = beams["Gantry"]
    # @assert size(gantry_angles, 1) == 1
    # gantry_angles = vec(gantry_angles)
    # couch_angles = beams["Couch"]
    # @assert size(couch_angles, 1) == 1
    # couch_angles = vec(couch_angles)

    # Hard-code beam arrangement because the .mat files for all patients
    # are 2.2GB, which takes too long to download.
    all_gantry_angles = Dict(
        "Protons1" => (60.f0, 180.f0, 60.f0),
    )
    all_couch_angles = Dict(
        "Protons1" => (0.f0, 0.f0, 180.f0),
    )
    all_placement_targets = Dict(
        "Protons1" => structures["CTV Combined Ring 0-10 mm"]
    )

    gantry_angles = all_gantry_angles[patient_ID]
    couch_angles = all_couch_angles[patient_ID]
    placement_target = all_placement_targets[patient_ID]
    iso_center = mean(placement_target.points, dims=1)[1:3]

    fields = Vector{Juliana.FieldDefinition{Float32}}(undef, length(gantry_angles))
    for (i, (g, c)) in enumerate(zip(gantry_angles, couch_angles))
        fields[i] = Juliana.FieldDefinition(
            "F_$(i)",
            i,
            convert(Float32, g),
            convert(Float32, c),
            0.f0,
            iso_center,
        )
    end
    field_centers = [placement_target for _ in 1:length(fields)]

    return ct, structures, presc, fields, field_centers
end
