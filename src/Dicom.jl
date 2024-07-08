using DICOM
using Logging
using SHA
using StaticArrays


RELEVANT_DIGITS::Int16 = 4


function get_study_instance_uid(patient_ID)
    # Taken from pydicom:
    # https://github.com/pydicom/pydicom/blob/master/pydicom/uid.py#L460-L469
    # https://github.com/pydicom/pydicom/blob/bbd76ff3e607de3851baafc89883f314a65a52e4/pydicom/uid.py#L240
    prefix = "1.2.826.0.1.3680043.8.498."
    max_uid_len = 64
    avail = max_uid_len - length(prefix)

    postfix = join(SHA.sha512(patient_ID))[1:avail]
    return "$(prefix)$(postfix)"
end


"""
    load_dicom_directory(dicom_dir; structure_names::AbstractArray{String}=[])

structure_names: Name of the structures to be loaded.
If empty: Load all structures.

Binary mask: All voxels within less than 0.5 * minimum(ct spacing).
"""
function load_dicom_directory(dicom_dir; structure_names::AbstractArray{String}=Vector{String}(undef, 0))
    filenames = readdir(dicom_dir)

    # Load the CT.
    ct_files = [f for f in filenames if startswith(f, "CT.")]
    ct_files = sort(ct_files, by=f -> parse(Int64, String(split(f, ".")[end-1])))
    ct_paths = ["$(dicom_dir)/$(f)" for f in ct_files];
    ct = Juliana.read_dicom_ct(ct_paths)

    # Load the structures.
    structure_files = [f for f in filenames if startswith(f, "RS.")]
    @assert length(structure_files) == 1
    structureset_path = "$(dicom_dir)/$(structure_files[1])"
    structures = Juliana.read_dicom_structureset(structureset_path, ct.grid; calculate_distance=false);

    if isempty(structure_names)
        structure_names = collect(keys(structures))
    end

    selected_structures = Dict{String, Juliana.Structure}()
    for name in structure_names
        selected_structures[name] = structures[name]
    end

    ct, shifted_structures, _ = Juliana.shift_to_zero(ct, selected_structures, [], calculate_distances=false);

    for name in structure_names
        structure = shifted_structures[name]
        @debug "Calculate distance mask for $name with $(size(structure.points, 1)) points..."
        shifted_structures[name] = Juliana.add_distance_and_mask_in_neighbourhood(
            structure,
            0.5*minimum(ct.grid.spacing),
        )
    end

    # TODO: Load dose distributions
    return ct, shifted_structures
end


function dicom_export_to_directory(ct,
                                   structures::Dict{String, Juliana.Structure},
                                   output_dir,
                                   study_instance_UID,
                                   new_patient_ID,
                                   doses::Dict{String, Juliana.ScalarGrid}=Dict{String, Juliana.ScalarGrid}())
    frame_of_reference_UID = "$(study_instance_UID).0"
    ct_series_instance_UID = "$(study_instance_UID).1"
    structureset_series_instance_UID = "$(study_instance_UID).2"
    patient_name = "$(new_patient_ID)^$(new_patient_ID)"

    # Export the CT.
    dicom_cts = Juliana.ct_to_dicom(
        ct,
        study_instance_UID,
        frame_of_reference_UID,
        ct_series_instance_UID,
        new_patient_ID,
        patient_name,
    );
    for (i, dicom_ct_slice) in enumerate(dicom_cts)
        DICOM.dcm_write("$(output_dir)/CT.$(i-1).dcm", dicom_ct_slice)
    end

    # Export the structures.
    dicom_structures = Juliana.structures_to_dicom(
        collect(values(structures)),
        study_instance_UID,
        frame_of_reference_UID,
        ct_series_instance_UID,
        structureset_series_instance_UID,
        new_patient_ID,
        patient_name,
        ct,
        dicom_cts,
    );
    DICOM.dcm_write("$(output_dir)/RS.dcm", dicom_structures);

    # Export the dose distributions.
    for (i, (label, dose)) in enumerate(doses)
        dose_series_instance_UID = "$(study_instance_UID).$(3+(i-1))"
        dicom_dose = Juliana.dose_to_dicom(
            dose,
            patient_name,
            new_patient_ID,
            study_instance_UID,
            dose_series_instance_UID,
            frame_of_reference_UID,
        )
        DICOM.dcm_write("$(output_dir)/RD_plan_$(label).dcm", dicom_dose)
    end
end


##############
# CT
##############
function read_dicom_ct(paths)
    intercept = -1000
    slope = 1
    
    row_spacing = Vector{Float32}(undef, length(paths))
    col_spacing = Vector{Float32}(undef, length(paths))
    img_position = Array{Float32, 2}(undef, length(paths), 3)
    
    slices = []

    for (i, path) in enumerate(paths)
        slice = DICOM.dcm_parse(path);
        @assert slice.ImageOrientationPatient == [1, 0, 0, 0, 1, 0]
        if slice.RescaleIntercept != -1000
            @warn "slice.RescaleIntercept != -1000"
        end
        if slice.RescaleSlope != 1.
            @warn "slice.RescaleSlope != 1."
        end
        if i == 1
            intercept = slice.RescaleIntercept
            slope = slice.RescaleSlope
        else
            @assert intercept == slice.RescaleIntercept
            @assert slope == slice.RescaleSlope
        end
        row_spacing[i] = slice.PixelSpacing[1]
        col_spacing[i] = slice.PixelSpacing[2]
        img_position[i, :] = slice.ImagePositionPatient
        
        # Transpose: DICOM is row-major, Julia is column-major.
        # https://dicom.innolitics.com/ciods/ct-image/image-pixel/7fe00010
        push!(slices, slice.PixelData')
        
        if i > 1
            @assert all(row_spacing[1] == row_spacing[i])
            @assert all(col_spacing[1] == col_spacing[i])
            @assert all(img_position[1, 1:2] .== img_position[i, 1:2])
        end
    end
    
    slice_separations = Vector{Float32}(undef, length(paths)-1)
    for i in 1:length(paths)-1
        slice_separations[i] = img_position[i+1, 3] - img_position[i, 3]
    end
    if !all(slice_separations[1] .== slice_separations)
        @warn "Not all CT slices have the exactly same spacing!"
    end
    # Less than one 100th of a mm difference.
    @assert maximum(abs.(slice_separations[1] .- slice_separations)) < 1e-3
    
    origin = img_position[1, :] ./ 10.f0
    values = convert.(Int16, slope .* cat(slices..., dims=3) .+ intercept)
    spacing = [row_spacing[1], col_spacing[1], slice_separations[1]] ./ 10.f0
    
    return ScalarGrid(
        values,
        Grid(
            SVector{3}(spacing),
            SVector{3}(origin),
            SVector{3}(size(values)),
        ),
    )
end


function ct_slice_to_dicom(array, spacing, origin, orientation, slice_index; decrease_precision=true)
    dicom_ct = DICOM.DICOMData(
        Dict{Tuple{UInt16, UInt16}, Any}(),
        :little,
        false,
        Dict{Tuple{UInt16,UInt16},String}(),
    )
    dicom_ct.Modality = "CT"
    dicom_ct.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2" # CT Image Storage
    dicom_ct.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
    dicom_ct.SeriesDescription = "CT Images"
    dicom_ct.StudyID = "0"
    dicom_ct.SeriesNumber = 0
    dicom_ct.AcquisitionNumber = slice_index
    dicom_ct.InstanceNumber = slice_index

    dicom_ct.RescaleIntercept = -1000
    dicom_ct.RescaleSlope = 1.
    # Transpose: DICOM is row-major, Julia is column-major.
    # https://dicom.innolitics.com/ciods/ct-image/image-pixel/7fe00010
    dicom_ct.PixelData = convert.(
        UInt16,
        (array' .- dicom_ct.RescaleIntercept) ./ dicom_ct.RescaleSlope,
    )

    # Explanation of PixelRepresentation:
    # https://dicom.innolitics.com/ciods/ct-image/image-pixel/00280103
    dicom_ct.PixelRepresentation = 0
    dicom_ct.HighBit = 15
    dicom_ct.Rows = size(array, 2)
    dicom_ct.Columns = size(array, 1)
    dicom_ct.NumberOfFrames = 1
    dicom_ct.SamplesPerPixel = 1
    dicom_ct.PixelSpacing = [spacing[1] * 10, spacing[2] * 10]
    if decrease_precision
        dicom_ct.PixelSpacing = round.(
            dicom_ct.PixelSpacing,
            digits=RELEVANT_DIGITS,
        )
    end

    dicom_ct.SliceThickness = spacing[3] * 10
    dicom_ct.SliceLocation = (origin[3] + spacing[3] * slice_index) * 10.f0
    if decrease_precision
        dicom_ct.SliceThickness = round(dicom_ct.SliceThickness, digits=RELEVANT_DIGITS)
        dicom_ct.SliceLocation = round(dicom_ct.SliceLocation, digits=RELEVANT_DIGITS)
    end

    if decrease_precision
        dicom_ct.ImagePositionPatient = [
            round(origin[1] * 10.f0, digits=RELEVANT_DIGITS),
            round(origin[2] * 10.f0, digits=RELEVANT_DIGITS),
            round((origin[3] + spacing[3] * slice_index) * 10.f0, digits=2),
        ]
    else
        dicom_ct.ImagePositionPatient = [
            origin[1] * 10.f0,
            origin[2] * 10.f0,
            (origin[3] + spacing[3] * slice_index) * 10.f0,
        ]
    end

    # Assumption about the patient orientation!!!
    dicom_ct.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    dicom_ct.PatientPosition = orientation
    @assert orientation == "HFS" # other orientations are not tested!!!
    dicom_ct.PhotometricInterpretation = "MONOCHROME2"
    dicom_ct.BitsAllocated = 16
    dicom_ct.BitsStored = 16
    dicom_ct.PixelRepresentation = 0 # Unsigned integer.

    # Use the default "DICOM Implicit VR Little Endian Transfer Syntax"
    # https://dicom.nema.org/dicom/2013/output/chtml/part05/chapter_10.html#sect_10.1
    dicom_ct.meta[tag"TransferSyntaxUID"] = "1.2.840.10008.1.2"
    # CT Image Storage
    dicom_ct.meta[tag"MediaStorageSOPClassUID"] = "1.2.840.10008.5.1.4.1.1.2"
    
    return dicom_ct
end


function ct_to_dicom(ct::ScalarGrid,
                     study_instance_UID::String,
                     frame_of_reference_UID::String,
                     ct_series_instance_UID::String,
                     patient_ID::String,
                     patient_name::String;
                     decrease_precision=true)
    slices = []
    for iz in 1:size(ct.data, 3)
        dicom_slice = ct_slice_to_dicom(
            ct.data[:, :, iz],
            ct.grid.spacing,
            ct.grid.origin,
            "HFS",
            iz-1,
        )

        slice_ID = "$(ct_series_instance_UID).$(iz-1)"
        dicom_slice.meta[tag"MediaStorageSOPInstanceUID"] = slice_ID
        dicom_slice.PatientName = patient_name
        dicom_slice.PatientID = patient_ID
        dicom_slice.SOPInstanceUID = slice_ID
        dicom_slice.StudyInstanceUID = study_instance_UID
        dicom_slice.SeriesInstanceUID = ct_series_instance_UID
        dicom_slice.FrameOfReferenceUID = frame_of_reference_UID

        push!(slices, dicom_slice)
    end

    return slices
end


#############
# Dose
#############
function read_dicom_dose(path::String)
    dicom_data = DICOM.dcm_parse(path)
    
    data = dicom_data.PixelData * dicom_data.DoseGridScaling
    nx, ny, nz = size(data)
    data_new = Array{eltype(data), 3}(undef, ny, nx, nz)
    for iz in 1:nz
        data_new[:, :, iz] = data[:, :, iz]'
    end
    data = data_new

    
    offsets = dicom_data.GridFrameOffsetVector
    differences = Vector{Float32}(undef, length(offsets) - 1)
    for i in 2:length(offsets)
        differences[i-1] = offsets[i] - offsets[i-1]
    end
    differences .= round.(differences, digits=4)
    @assert all(differences.== differences[1])
    @assert isa(dicom_data.SliceThickness, Vector{Any}) || dicom_data.SliceThickness ≈ differences[1]
    @assert dicom_data.DoseUnits == "GY"
    @assert dicom_data.DoseType == "PHYSICAL"

    spacing = [
        dicom_data.PixelSpacing[1],
        dicom_data.PixelSpacing[2],
        differences[1],
    ] / 10. # cm
    origin = dicom_data.ImagePositionPatient / 10. # cm
    
    return ScalarGrid(
        data * 1.1, # Convert from Gy to Gy RBE.
        Grid(
            SVector{3}(spacing),
            SVector{3}(origin),
            SVector{3}(collect(size(data))),
        ),
    )
end


function dose_to_dicom(dose::ScalarGrid,
                       patient_name::String,
                       patient_ID::String,
                       study_instance_UID::String,
                       series_instance_UID::String,
                       frame_of_reference_UID::String)
    dicom_dose = DICOM.DICOMData(
        Dict{Tuple{UInt16, UInt16}, Any}(),
        :little,
        false,
        Dict{Tuple{UInt16,UInt16},String}(),
    )

    spacing = collect(dose.grid.spacing)
    origin = collect(dose.grid.origin)

    # Physical dose.
    ###################################################
    # Assumption: This dose object contains dose RBE.
    ###################################################
    data = dose.data ./ 1.1

    # DICOM is row-major:
    # https://dicom.innolitics.com/ciods/rt-dose/image-pixel/7fe00010
    nx, ny, nz = size(data)
    pixel_data = Array{eltype(data), 3}(undef, ny, nx, nz)
    for iz in 1:nz
        pixel_data[:, :, iz] = data[:, :, iz]'
    end

    # Make sure we can really fit the entire range into our ints.
    dose_scaling = (maximum(pixel_data) + 1) / typemax(UInt16)

    # Monochrome2 is only 16bit.
    # https://pydicom.github.io/pydicom/dev/guides/encoding/rle_lossless.html
    array = convert.(UInt16, round.(pixel_data / dose_scaling))

    dicom_dose.Modality = "RTDOSE"
    dicom_dose.SOPClassUID = "1.2.840.10008.5.1.4.1.1.481.2" # RT Dose Storage
    dicom_dose.PixelData = array
    dicom_dose.DoseGridScaling = dose_scaling
    # Explanation of PixelRepresentation:
    # https://dicom.innolitics.com/ciods/ct-image/image-pixel/00280103
    dicom_dose.PixelRepresentation = 0
    dicom_dose.HighBit = 15
    dicom_dose.Rows = size(data, 2)
    dicom_dose.Columns = size(data, 1)
    dicom_dose.NumberOfFrames = size(data, 3)
    dicom_dose.PixelSpacing = [
        spacing[1],
        spacing[2],
    ] .* 10. # mm
    dicom_dose.SliceThickness = spacing[3] * 10 # mm
    dicom_dose.GridFrameOffsetVector = [round(i*spacing[3], digits=RELEVANT_DIGITS) for i in 1:size(data, 3)] .* 10. # mm

    dicom_dose.ImagePositionPatient = origin .* 10. # mm
    # Assumption about the patient position!!!
    dicom_dose.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    dicom_dose.PhotometricInterpretation = "MONOCHROME2"
    dicom_dose.DoseUnits = "GY"
    dicom_dose.DoseType = "PHYSICAL"
    dicom_dose.DoseSummationType = "PLAN"
    dicom_dose.BitsAllocated = 16
    dicom_dose.BitsStored = 16
    dicom_dose.SamplesPerPixel = 1
    # Use the default "DICOM Implicit VR Little Endian Transfer Syntax"
    # https://dicom.nema.org/dicom/2013/output/chtml/part05/chapter_10.html#sect_10.1
    dicom_dose.meta[tag"TransferSyntaxUID"] = "1.2.840.10008.1.2"
    dicom_dose.meta[tag"MediaStorageSOPClassUID"] = "1.2.840.10008.5.1.4.1.1.481.2" # RT Dose Storage
    # dicom_dose.meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()

    dicom_dose.PatientName = patient_name
    dicom_dose.PatientID = patient_ID
    dicom_dose.StudyInstanceUID = study_instance_UID
    dicom_dose.SeriesInstanceUID = series_instance_UID
    dicom_dose.SOPInstanceUID = series_instance_UID
    dicom_dose.FrameOfReferenceUID = frame_of_reference_UID

    return dicom_dose
end


#############
# Structures
#############
function read_dicom_structure_points(path)
    structures = DICOM.dcm_parse(path);

    structure_points = Dict{String, Array{Float32, 2}}()
    is_target = Dict{String, Bool}()
    for (i, (roi, seq, obs)) in enumerate(zip(structures.StructureSetROISequence, structures.ROIContourSequence, structures.RTROIObservationsSequence))
        @debug "Loading structure $(roi.ROIName)..."

        if isnothing(seq.ContourSequence) || (length(seq.ContourSequence) == 0)
            @warn "No contours for $(roi.ROIName). Skipping..."
            continue
        end

        all_points = Vector{Array{Float32, 2}}(undef, 0)
        for (contour_index, contour) in enumerate(seq.ContourSequence)
            n_points = contour.NumberOfContourPoints
            # mm -> cm
            points = reshape(contour.ContourData, (3, n_points))' / 10.

            new_points = Array{Float32, 2}(undef, size(points, 1), 4)
            new_points[:, 1:3] .= points
            new_points[:, 4] .= contour_index
            push!(all_points, new_points)
        end
        all_points = vcat(all_points...)
        structure_points[roi.ROIName] = all_points
        # See:
        # https://dicom.innolitics.com/ciods/rt-structure-set/rt-roi-observations/30060080/300600a4
        is_target[roi.ROIName] = obs.RTROIInterpretedType in ["GTV", "CTV", "PTV"]
    end
    return structure_points, is_target
end


function read_dicom_structureset(path, mask_grid; calculate_distance=false)
    structures = Dict{String, Structure}()
    structure_points, structure_is_target = read_dicom_structure_points(path)
    for (name, points) in structure_points
        contours = split_by_contour_index(points);
        if calculate_distance
            distance_mask = calculate_distance_mask(
                Tuple(mask_grid.size),
                mask_grid.spacing,
                contours,
            )
            binary_mask = to_binary_mask(distance_mask)
        else
            distance_mask = Array{Float32, 3}(
                undef,
                mask_grid.size[1],
                mask_grid.size[2],
                mask_grid.size[3],
            )
            binary_mask = Array{Bool, 3}(
                undef,
                mask_grid.size[1],
                mask_grid.size[2],
                mask_grid.size[3],
            )
        end

        structures[name] = Structure(
            name,
            points,
            binary_mask,
            distance_mask,
            mask_grid,
            structure_is_target[name],
        )
    end
    return structures
end


function structures_to_dicom(structures::AbstractVector{Juliana.Structure},
                             study_instance_UID::String,
                             frame_of_reference_UID::String,
                             base_ct_series_instance_UID::String,
                             structureset_series_instance_UID::String,
                             patient_ID::String,
                             patient_name::String,
                             ct,
                             ct_datasets;
                             drop_precision::Bool=true)
    ds = DICOM.DICOMData(
        Dict{Tuple{UInt16, UInt16}, Any}(),
        :little,
        false,
        Dict{Tuple{UInt16,UInt16},String}(),
    )
    ds.Modality = "RTSTRUCT"
    ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.481.3" # RT Structure Set Storage
    ds.StructureSetLabel = "StructureSet"
    ds.StructureSetName = "StructureSet"
    ds.ApprovalStatus = "UNAPPROVED"
    ds.AccessionNumber = ""
    ds.StudyID = ""
    ds.SeriesNumber = 0

    referenced_frame_of_reference = DICOM.DICOMData(
        Dict{Tuple{UInt16, UInt16}, Any}(),
        :little,
        false,
        Dict{Tuple{UInt16,UInt16},String}(),
    )
    referenced_frame_of_reference.FrameOfReferenceUID = frame_of_reference_UID
    referenced_frame_of_reference.RTReferencedStudySequence = build_referenced_study_sequence(
        study_instance_UID,
        base_ct_series_instance_UID,
        ct_datasets,
    )
    ds.ReferencedFrameOfReferenceSequence = [referenced_frame_of_reference]

    roi_seq = []
    for (i, structure) in enumerate(structures)
        roi_ds = DICOM.DICOMData(
            Dict{Tuple{UInt16, UInt16}, Any}(),
            :little,
            false,
            Dict{Tuple{UInt16,UInt16},String}(),
        )
        # -1: Julia indices start at 1.
        roi_ds.ROINumber = i - 1
        roi_ds.ReferencedFrameOfReferenceUID = frame_of_reference_UID
        roi_ds.ROIName = structure.name
        roi_ds.ROIGenerationAlgorithm = "MANUAL"
        push!(roi_seq, roi_ds)
    end
    ds.StructureSetROISequence = roi_seq

    roi_contour_sequence = []
    for (i, structure) in enumerate(structures)
        contours = build_contour_sequence(
            structure,
            ct,
            ct_datasets,
            decrease_precision=drop_precision,
        )
        roi_contour = DICOM.DICOMData(
            Dict{Tuple{UInt16, UInt16}, Any}(),
            :little,
            false,
            Dict{Tuple{UInt16,UInt16},String}(),
        )
        # -1: Julia indices start at 1.
        roi_contour.ReferencedROINumber = i - 1
        roi_contour.ContourSequence = contours
        push!(roi_contour_sequence, roi_contour)
    end
    ds.ROIContourSequence = roi_contour_sequence

    observations = []
    for (i, structure) in enumerate(structures)
        observation = DICOM.DICOMData(
        Dict{Tuple{UInt16, UInt16}, Any}(),
            :little,
            false,
            Dict{Tuple{UInt16,UInt16},String}(),
        )
        # -1: Julia indices start at 1.
        observation.ObservationNumber = i - 1
        observation.ReferencedROINumber = i - 1
        observation.ROIObservationLabel = structure.name
        observation.ROIInterpreter = ""

        n = lowercase(structure.name)
        is_target = occursin("ctv", n) || occursin("gtv", n) || occursin("ptv", n)

        observation.RTROIInterpretedType = is_target ? "PTV" : "ORGAN"

        push!(observations, observation)
    end
    ds.RTROIObservationsSequence = observations

    # Use the default "DICOM Implicit VR Little Endian Transfer Syntax"
    # https://dicom.nema.org/dicom/2013/output/chtml/part05/chapter_10.html#sect_10.1
    ds.meta[tag"TransferSyntaxUID"] = "1.2.840.10008.1.2"
    # RT Structure Set Storage
    ds.meta[tag"MediaStorageSOPClassUID"] = "1.2.840.10008.5.1.4.1.1.481.3"
    ds.meta[tag"MediaStorageSOPInstanceUID"] = structureset_series_instance_UID

    ds.PatientName = patient_name
    ds.PatientID = patient_ID
    ds.StudyInstanceUID = study_instance_UID
    ds.SeriesInstanceUID = structureset_series_instance_UID
    ds.SOPInstanceUID = structureset_series_instance_UID
    ds.FrameOfReferenceUID = frame_of_reference_UID

    return ds
end


function build_contour_sequence(structure::Juliana.Structure,
                                ct::Juliana.ScalarGrid,
                                ct_datasets;
                                decrease_precision::Bool = true)
    if size(structure.points, 2) == 4
        slices = Juliana.split_by_contour_index(structure.points)
    elseif size(structure.points, 2) == 3
        slices = Juliana.split_by_z(structure.points)
    else
        error("Structure point matrix must have 3 or 4 columns...")
    end

    contours = []
    for points in slices
        points = points[:, 1:3]
        z = points[1, 3]
        @assert maximum(abs.(points[:, 3] .- z)) <= 1e-6
        # +1: Julia indices start at 1.
        slice_ind = convert(
            Int64,
            round((z - ct.grid.origin[3] - ct.grid.spacing[3] / 2) / ct.grid.spacing[3]) + 1,
        )

        ####################################################
        # TODO: Make sure all the incides start at 1 etc.
        ####################################################

        # Avoid explicitly closing the contour: DICOM takes care of this because
        # of the CLOSED_PLANAER geometric type (see below).
        if points[1, :] ≈ points[end, :]
            points = points[1:end-1, :]
        end

        ct_dataset = ct_datasets[slice_ind]

        contour = DICOM.DICOMData(
        Dict{Tuple{UInt16, UInt16}, Any}(),
            :little,
            false,
            Dict{Tuple{UInt16,UInt16},String}(),
        )

        contour_image = DICOM.DICOMData(
        Dict{Tuple{UInt16, UInt16}, Any}(),
            :little,
            false,
            Dict{Tuple{UInt16,UInt16},String}(),
        )
        contour_image.ReferencedSOPClassUID = ct_dataset.SOPClassUID
        contour_image.ReferencedSOPInstanceUID = ct_dataset.SOPInstanceUID

        contour.ContourImageSequence = [contour_image]
        contour.ContourGeometricType = "CLOSED_PLANAR"
        contour.NumberOfContourPoints = size(points, 1)
        # We use cm because the Fiona standalone does, but DICOM expects mm.
        points = points * 10
        # DICOM export to Fiona makes problems if the numbers are not exactly
        # the slice spacing.
        if decrease_precision
            points = round.(points, digits=Juliana.RELEVANT_DIGITS)
        end
        contour.ContourData = vec(points')
        push!(contours, contour)
    end

    return contours
end


function build_referenced_study_sequence(study_instance_UID::String,
                                         base_ct_series_instance_UID::String,
                                         ct_datasets)
    ref_study = DICOM.DICOMData(
        Dict{Tuple{UInt16, UInt16}, Any}(),
        :little,
        false,
        Dict{Tuple{UInt16,UInt16},String}(),
    )
    ref_study.ReferencedSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.3" # RT Structure Set Storage
    ref_study.ReferencedSOPInstanceUID = study_instance_UID
    ref_study.RTReferencedSeriesSequence = build_referenced_series_sequence(
        base_ct_series_instance_UID,
        ct_datasets,
    )
    return [ref_study]
end


function build_referenced_series_sequence(base_ct_series_instance_UID::String, ct_datasets)
    ref_series = DICOM.DICOMData(
        Dict{Tuple{UInt16, UInt16}, Any}(),
        :little,
        false,
        Dict{Tuple{UInt16,UInt16},String}(),
    )
    ref_series.SeriesInstanceUID = base_ct_series_instance_UID

    contour_image_sequence = []
    for ct_ds in ct_datasets
        ds = DICOM.DICOMData(
            Dict{Tuple{UInt16, UInt16}, Any}(),
            :little,
            false,
            Dict{Tuple{UInt16,UInt16},String}(),
        )
        ds.ReferencedSOPClassUID = ct_ds.SOPClassUID
        ds.ReferencedSOPInstanceUID = ct_ds.SOPInstanceUID
        push!(contour_image_sequence, ds)
    end
    ref_series.ContourImageSequence = contour_image_sequence
    return [ref_series]
end


##############
# RT Ion Plan
##############
function parse_dicom_beam(beam, spot_index_start)
    T = Float32
    gantry_angle = convert(T, beam.IonControlPointSequence[1].GantryAngle)
    couch_angle = convert(T, beam.IonControlPointSequence[1].PatientSupportAngle)
    nozzle_extraction = convert(T, beam.IonControlPointSequence[1].SnoutPosition)
    field_center = convert.(T, beam.IonControlPointSequence[1].IsocenterPosition)
    tune_ID = beam.IonControlPointSequence[1].ScanSpotTuneID

    n_spots = 0
    for point in beam.IonControlPointSequence
        n_spots += point.NumberOfScanSpotPositions
    end
    spots = Vector{Juliana.FionaStandalone.Spot{T}}(undef, n_spots)
    spot_index = 0
    for point in beam.IonControlPointSequence
        for i in 0:point.NumberOfScanSpotPositions-1
            E = point.NominalBeamEnergy * 1000 # keV
            u = point.ScanSpotPositionMap[i*2+1] / 10 # cm
            t = point.ScanSpotPositionMap[i*2+2] / 10 # cm
            w = point.ScanSpotMetersetWeights[i+1]

            @assert isnothing(point.GantryAngle) || (point.GantryAngle == gantry_angle)
            @assert isnothing(point.PatientSupportAngle) || point.PatientSupportAngle == couch_angle
            @assert isnothing(point.SnoutPosition) || point.SnoutPosition == nozzle_extraction
            @assert isnothing(point.IsocenterPosition) || convert.(T, point.IsocenterPosition) == field_center
            @assert isnothing(point.ScanSpotTuneID) || point.ScanSpotTuneID == tune_ID

            spots[spot_index+1] = Juliana.FionaStandalone.Spot(
                spot_index_start + spot_index,
                convert(T, t),
                convert(T, u),
                convert(T, w),
                convert(T, E),
                beam.NumberOfRangeShifters, # Number of absorbers??? TODO
            )
            spot_index += 1
        end
    end

    return Juliana.FionaStandalone.FieldDefinition(
        String(beam.BeamName),
        beam.BeamNumber,
        gantry_angle,
        couch_angle,
        nozzle_extraction / 10, # cm
        Dict(
            "x" => field_center[1] / 10, # cm
            "y" => field_center[2] / 10, # cm
            "z" => field_center[3] / 10, # cm
        ),
        spots,
    ), tune_ID
end


function load_dicom_treatment_plan(path, target_dose)
    dcm = dcm_parse(path)
    
    fields = Vector{Juliana.FionaStandalone.FieldDefinition{Float32}}(undef, 0)
    tune_IDs = []
    n_spots_before = 0
    for beam in dcm.IonBeamSequence
        field, tune_ID = parse_dicom_beam(beam, n_spots_before)
        n_spots_before += length(field.spots)
        push!(fields, field)
        push!(tune_IDs, tune_ID)
    end
    
    plan = Juliana.FionaStandalone.TreatmentPlan(fields, target_dose)
    
    return plan, tune_IDs
end
