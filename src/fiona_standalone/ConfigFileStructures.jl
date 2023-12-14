import ..Structure
import ..split_by_z
import ..Grid


struct MainConfig
    # Paths to data and config files.
    ct::String
    optimizationSettings::String
    spotConfiguration::String
    doseResults::String
    # Path to which the plan JSON file will be written. (output!)
    planResults::String
    prescribedDose::Real
    # Path to a plan JSON (input). When using the RUN command,
    # this is overwritten during thefield definition step. However, useful if
    # only a dose calculation is needed.
    plan::String
    preabsorberWed::Real
    cutoffFactor::Real
    doseResolution::Real
    huToSp::String
    beamdata::String
    beamlinePhasespace::String
    nozzlePhasespace::String
    absorberPhasespace::String
end

function MainConfig(ct::String,
                    optimizationSettings::String,
                    spotConfiguration::String,
                    doseResults::String,
                    planResults::String,
                    prescribedDose::Real,
                    binPath::String,
                    plan::String)
    MainConfig(
        ct,
        optimizationSettings,
        spotConfiguration,
        doseResults,
        planResults,
        prescribedDose,
        plan,
        4.173,                   # preabsorberWed
        5.0,                     # cutoffFactor
        0.35,                    # doseResolution
        "$binPath/huToSp.json",  # huToSp
        "$binPath/beamdata.xml", # beamdata
        "$binPath/beamline.xml", # beamlinePhasespace
        "$binPath/nozzle.xml",   # nozzlePhasespace
        "$binPath/nozzle.xml",   # absorberPhasespace
    )
end

function MainConfig(ct_path::String, output_dir::String, prescribed_dose::Real, bin_path::String)
    spot_config_path  = "$output_dir/spot_config.json"
    optim_config_path = "$output_dir/optim_config.json"
    result_dose_path  = "$output_dir/result_dose.dat"
    result_plan_path  = "$output_dir/result_plan.json"
    properties_path   = "$output_dir/config.json"

    MainConfig(
        ct_path,
        optim_config_path,
        spot_config_path,
        result_dose_path,
        result_plan_path,
        prescribed_dose,
        bin_path,
        result_plan_path,
    )
end


FionaContour = Dict{String, Vector{Dict{String, Real}}}

function convert(FionaContour, contour::Array{T, 2}) where {T}
    N_points = size(contour, 1)
    point_list = Array{Dict{String, T}, 1}(undef, N_points)
    for i in 1:N_points
        point_list[i] = Dict{String, T}(
            "x" => contour[i, 1],
            "y" => contour[i, 2],
            "z" => contour[i, 3],
        )
    end
    return Dict{String, Array{Dict{String, T}, 1}}(
        "points" => point_list,
    )
end


struct StructureDefinition
    label::String
    id::Integer
    contours::Array{FionaContour, 1}
    zMin::Real
    zMax::Real
end

function StructureDefinition(structure::Structure, id::Integer)
    contours = split_by_z(structure.points)
    StructureDefinition(
        structure.name,
        id,
        [convert(FionaContour, contour) for contour in contours],
        minimum(structure.points[:, 3]),
        maximum(structure.points[:, 3]),
    )
end


struct Constraint
    id::Integer
    importance::Real
    dose::Real
    volume::Real
    # Type can be either 'DOSE_VOLUME' or 'MEAN'.
    type::String
    enforced::Bool
    # limitType is either "UPPER" or "LOWER".
    limitType::String
end


struct StructureConstraints
    structure::StructureDefinition
    constraints::Array{Constraint, 1}
end


struct OptimizationGrid
    size::Tuple{Integer, Integer, Integer}
    spacing::Tuple{Real, Real, Real}
    origin::Tuple{Real, Real, Real}
end

function convert(OptimizationGrid, grid::Grid)
    OptimizationGrid(
        Tuple(grid.size),
        Tuple(grid.spacing),
        Tuple(grid.origin),
    )
end


struct OptimizationSettings
    targetDose::Real
    targetDoseAtEdge::Real
    target::StructureDefinition
    structureConstraints::Array{StructureConstraints, 1}
    optimizationGrid::OptimizationGrid
    numberOfIterations::Integer
    sigmaFalloff::Real
    importanceExponent::Real
    falloffDistance::Real
    inTargetDistance::Real
    minWeight::Real
    hotSpotEnabled::Bool
    hotSpotDose::Real
    hotSpotImportance::Real
    rrEnabled::Bool
    rrScaleValue::Real
    rrScenarioImportance::Real
end

# Default arguments via outer constructor.
function OptimizationSettings(targetDose::Real,
                              targetDoseAtEdge::Real,
                              target::StructureDefinition,
                              structureConstraints::Array{StructureConstraints, 1},
                              optimizationGrid::OptimizationGrid)
    OptimizationSettings(
        targetDose,
        targetDoseAtEdge,
        target,
        structureConstraints,
        optimizationGrid,
        1,              # numberOfIterations
        0.8,            # sigmaFalloff
        10.0,           # importanceExponent
        -0.1672349,     # falloffDistance
        3.0,            # inTargetDistance
        70.0,           # minWeight
        false,          # hotSpotEnabled
        107.0,          # hotSpotDose
        3.0,            # hotSpotImportance
        false,          # rrEnabled
        3.0,            # rrScaleValue
        1.0,            # rrScenarioImportance
    )
end


struct SpotPlacementFieldDefinition
    targetStructureId::Integer
    gantryAngle::Real
    couchAngle::Real
    nozzleExtraction::Real
    preabsorberSetting::String
    # Spot placement margin and spacing.
    margin::Real
    xSpacing::Real
    ySpacing::Real
end

function SpotPlacementFieldDefinition(targetStructureId::Integer,
                         gantryAngle::Real,
                         couchAngle::Real,
                         nozzleExtraction::Real,
                         preabsorberSetting::String)
    SpotPlacementFieldDefinition(
        targetStructureId,
        gantryAngle,
        couchAngle,
        nozzleExtraction,
        preabsorberSetting,
        0.5,
        0.4,
        0.4,
    )
end


struct SpotPlacementConfig
    fields::Array{SpotPlacementFieldDefinition, 1}
    structures::Array{StructureDefinition}
end

function SpotPlacementConfig(gantry_angles::Vector{T},
                             couch_angles::Vector{T},
                             nozzle_extractions::Vector{T},
                             preabsorber_setting::String,
                             spot_placement_structure::Structure) where {T}
    field_definitions = Vector{SpotPlacementFieldDefinition}(undef, length(gantry_angles))
    for (i, (gantry_angle, couch_angle, nozzle_extraction)) in enumerate(zip(gantry_angles, couch_angles, nozzle_extractions))
        field_definitions[i] = SpotPlacementFieldDefinition(
            0,
            gantry_angle,
            couch_angle,
            nozzle_extraction,
            preabsorber_setting,
        )
    end
    return SpotPlacementConfig(
        field_definitions,
        [StructureDefinition(spot_placement_structure, 0)],
    )
end


struct Spot
    id::Integer
    t::Real
    u::Real
    weight::Real
    energykeV::Real
    # 0 or 1
    numberOfAbsorbers::Integer
end


struct FieldDefinition
    label::String
    id::Integer
    gantryAngle::Real
    couchAngle::Real
    nozzleExtraction::Real
    fieldCenter::Dict{String, Real}
    spots::Vector{Spot}
end


struct TreatmentPlan
    fields::Vector{FieldDefinition}
    prescribedDose::Real
end


# Constructors to build the structures from a JSON representation.
function Spot(spot_dict::Dict)
    Spot(
        spot_dict["id"],
        spot_dict["t"],
        spot_dict["u"],
        spot_dict["weight"],
        spot_dict["energykeV"],
        spot_dict["numberOfAbsorbers"],
    )
end


function FieldDefinition(field_dict::Dict)
    spots = [Spot(spot_dict) for spot_dict in field_dict["spots"]];
    FieldDefinition(
        field_dict["label"],
        field_dict["id"],
        field_dict["gantryAngle"],
        field_dict["couchAngle"],
        field_dict["nozzleExtraction"],
        field_dict["fieldCenter"],
        spots,
    )
end


function TreatmentPlan(plan_dict::Dict)
    fields = [FieldDefinition(field_dict) for field_dict in plan_dict["fields"]]
    TreatmentPlan(fields, plan_dict["prescribedDose"])
end
