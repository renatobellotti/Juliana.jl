using Adapt
using StaticArrays


struct Grid{S<:Real, T<:Integer}
    spacing::SVector{3, S}
    origin::SVector{3, S}
    size::SVector{3, T}
end


struct GpuGrid{T, U}
    spacing::T
    origin::T
    size::U
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
    constraints::Array{Constraint{T}, 1}
end


# TPS = Treatment Planning System
abstract type AbstractTps end


# function Adapt.adapt_structure(to, grid::Grid{Vector{U}, Vector{V}}) where {U<:Real, V<:Integer}
#     d_spacing = Adapt.adapt_structure(to, grid.spacing)
#     d_origin  = Adapt.adapt_structure(to, grid.origin)
#     d_size    = Adapt.adapt_structure(to, grid.size)
    
#     return Grid{CuArray{U, 1, CUDA.Mem.DeviceBuffer}, CuArray{V, 1, CUDA.Mem.DeviceBuffer}}(
#         d_spacing,
#         d_origin,
#         d_size,
#     )
# end


function Adapt.adapt_structure(to, grid::Grid)
    d_spacing = Adapt.adapt_structure(to, collect(grid.spacing))
    d_origin  = Adapt.adapt_structure(to, collect(grid.origin))
    d_size    = Adapt.adapt_structure(to, collect(grid.size))
    
    return GpuGrid(
        d_spacing,
        d_origin,
        d_size,
    )
end


function Adapt.adapt_structure(to, grid::Juliana.GpuGrid)
    return Juliana.GpuGrid(
        Adapt.adapt_structure(to, grid.spacing),
        Adapt.adapt_structure(to, grid.origin),
        Adapt.adapt_structure(to, grid.size),
    )
end


struct PatientData{T}
    ct::ScalarGrid
    structures::Dict{String, Structure}
    prescriptions::Prescriptions{T}
end


struct BeamArrangement{T<:Real, V<:AbstractVector{T}}
    gantry_angles::V
    couch_angles::V
    nozzle_extractions::V
end


struct OptimisationConfiguration{NumberType, VolumeType, DijType}
    normalisationDose::NumberType
    normalisationStructureMask::VolumeType
    ct::VolumeType
    idealDose::VolumeType
    importance::VolumeType
    minimumDose::VolumeType
    maximumDose::VolumeType
    Dij::DijType
    Dij_T::DijType # Needed because we cannot efficiently transpose GPU arrays.
    structures::Dict{String, VolumeType}
    prescriptions::Prescriptions{NumberType}
end


function to_gpu(config::Juliana.OptimisationConfiguration)
    opt_target_mask = cu(config.normalisationStructureMask)
    opt_ct = cu(config.ct)
    opt_ideal_dose = cu(config.idealDose)
    opt_importance = cu(config.importance)
    opt_minimum_dose = cu(config.minimumDose)
    opt_maximum_dose = cu(config.maximumDose)
    opt_Dij = CUDA.CUSPARSE.CuSparseMatrixCSR(config.Dij)
    opt_Dij_T = CUDA.CUSPARSE.CuSparseMatrixCSR(config.Dij')
    opt_structures = Dict{String, typeof(opt_ideal_dose)}()
    for (name, structure) in config.structures
        opt_structures[name] = cu(structure)
    end
    prescriptions = config.prescriptions
    return Juliana.OptimisationConfiguration{typeof(config.normalisationDose), typeof(opt_ct), typeof(opt_Dij)}(
        config.normalisationDose,
        opt_target_mask,
        opt_ct,
        opt_ideal_dose,
        opt_importance,
        opt_minimum_dose,
        opt_maximum_dose,
        opt_Dij,
        opt_Dij_T,
        opt_structures,
        config.prescriptions,
    )
end


struct TargetCoverage{T<:Real}
    target_name::String
    V95::T       #  %
    Dmax::T      # Gy
    Dmean::T     # Gy
    Dmin::T      # Gy
    HI_D5_D95::T # no units
end

function TargetCoverage(dose::Juliana.ScalarGrid,
                        structure::Juliana.Structure,
                        prescribed_dose::Real)
    dose_data = convert.(Float32, dose.data)
    V95 = Juliana.dvh_v(
        dose_data,
        structure.mask,
        0.95f0 * prescribed_dose,
    )
    return TargetCoverage(
        structure.name,
        V95,
        maximum(dose_data .* structure.mask),
        sum(dose_data .* structure.mask) / sum(structure.mask),
        minimum(dose_data[structure.mask]),
        homogeneity(dose_data, structure.mask),
    )
end


function calculate_constraint_fulfillment(prescriptions, structures, dose_distribution)
    T = Float32
    constraint_results = Vector{Tuple{Juliana.Constraint{T}, T}}(
        undef,
        length(prescriptions.constraints),
    )
    for (i, constr) in enumerate(prescriptions.constraints)
        constraint_results[i] = constr, Juliana.evaluate_constraint(
            constr,
            structures[constr.structure_name].mask,
            dose_distribution,
        )
    end
    return constraint_results
end


struct Report
    label
    dose_distribution
    constraint_fulfillments
    target_coverages
    Dmax
end

function Report(label::String,
                dose_distribution::Juliana.ScalarGrid,
                prescriptions::Juliana.Prescriptions,
                structures)
    # Calculate how much the OAR constraints are fulfilled.
    constraint_results = calculate_constraint_fulfillment(
        prescriptions,
        structures,
        dose_distribution,
    )

    # Calculate target coverage for each target.
    target_coverages = Vector{TargetCoverage}(
        undef,
        length(prescriptions.target_doses),
    )
    for (i, (target_name, target_dose)) in enumerate(prescriptions.target_doses)
        target_coverages[i] = TargetCoverage(
            dose_distribution,
            structures[target_name],
            target_dose,
        )
    end

    # Build the report.
    return Report(
        label,
        dose_distribution,
        constraint_results,
        target_coverages,
        maximum(dose_distribution.data),
    )
end


struct RobustReport
    reports
end


function RobustReport(doses::Dict{String, Juliana.ScalarGrid}, prescriptions, structures)
    reports = Vector{Juliana.Report}(undef, length(doses))
    for (i, (label, dose)) in enumerate(doses)
        reports[i] = Juliana.Report(
            label,
            dose,
            prescriptions,
            structures,
        )
    end
    return RobustReport(reports)
end


function evaluate_constraint(constraint::Juliana.Constraint,
                             structure_mask,
                             dose::Juliana.ScalarGrid)
    if Juliana.is_maximum_constraint(constraint)
        return maximum(dose.data .* structure_mask)
    elseif constraint.kind == Juliana.constraint_mean
        return sum(dose.data .* structure_mask) / sum(structure_mask)
    elseif constraint.kind == Juliana.constraint_dvh_d
        return Juliana.dvh_d(
            convert.(Float32, dose.data),
            structure_mask,
            constraint.volume,
        )
    else
        @error "Unknown constraint type $(constraint.kind); $(constraint)"
    end
end


# Equality operators for some Juliana structs.
function Base.:(==)(a::Grid{S, T}, b::Grid{S, T}) where {S<:AbstractArray, T<:AbstractArray}
    return (a.origin == b.origin) && (a.spacing == b.spacing) && (a.size == b.size)
end

function Base.:(==)(a::ScalarGrid{A, S, T}, b::ScalarGrid{A, S, T}) where {A<:AbstractArray, S, T}
    return (a.data == b.data) && (a.grid == b.grid)
end


function Base.:(==)(a::Constraint{T}, b::Constraint{T}) where {T}
    return (a.structure_name == b.structure_name) &&
        (a.kind == b.kind) &&
        (a.dose == b.dose) &&
        (a.volume == b.volume) &&
        (a.priority == b.priority) &&
        (a.direction == b.direction)
end


function Base.:(==)(a::Prescriptions{T}, b::Prescriptions{T}) where {T}
    return (a.target_doses == b.target_doses) && (a.constraints == b.constraints)
end


function Base.:(==)(a::Structure{A, B, C, S, T}, b::Structure{A, B, C, S, T}) where {A<:AbstractArray, B<:AbstractArray, C<:AbstractArray, S, T}
    return (a.name == b.name) &&
        (a.points == b.points) &&
        (a.mask == b.mask) &&
        (a.distanceFromStructure == b.distanceFromStructure) &&
        (a.grid == b.grid) &&
        (a.is_target == b.is_target)
end


# Copy operators.
function Base.copy(s::Juliana.Structure)
    return Juliana.Structure(
        s.name,
        copy(s.points),
        copy(s.mask),
        copy(s.distanceFromStructure),
        s.grid,
        s.is_target,
    )
end


# Machine parameter data types.
struct LookUpTable{T<:Real, VectorType<:AbstractArray{T}, IntVectorType<:AbstractArray, Table<:AbstractMatrix{T}}
    energies::VectorType
    x0::VectorType
    dx::VectorType
    table::Table
    firstNanIndex::IntVectorType
end

struct PhaseSpace{T}
    energies::T

    a0t::T
    a1t::T
    a2t::T

    a0u::T
    a1u::T
    a2u::T

    is_symmetric::Bool
end

"""
    function PhaseSpace(energies, a0, a1, a2)

Constructor for symmetric phase spaces, e. g. for the nozzle and the preabsorber.
"""
function PhaseSpace(energies::T, a0::T, a1::T, a2::T) where {T}
    return PhaseSpace(
        energies,
        a0,
        a1,
        a2,
        a0,
        a1,
        a2,
        true,
    )
end


function Adapt.adapt_structure(to, ddc::LookUpTable)
    return LookUpTable(
        Adapt.adapt_structure(to, ddc.energies),
        Adapt.adapt_structure(to, ddc.x0),
        Adapt.adapt_structure(to, ddc.dx),
        Adapt.adapt_structure(to, ddc.table),
        Adapt.adapt_structure(to, ddc.firstNanIndex),
    )
end


function Adapt.adapt_structure(to, phase_space::PhaseSpace)
    return PhaseSpace(
        Adapt.adapt_structure(to, phase_space.energies),
        Adapt.adapt_structure(to, phase_space.a0t),
        Adapt.adapt_structure(to, phase_space.a1t),
        Adapt.adapt_structure(to, phase_space.a2t),
        Adapt.adapt_structure(to, phase_space.a0u),
        Adapt.adapt_structure(to, phase_space.a1u),
        Adapt.adapt_structure(to, phase_space.a2u),
        phase_space.is_symmetric,
    )
end


struct JulianaTps
    convert_to_sp
    d_depth_dose_curves
    d_sigma_mcs_curves
    d_phase_space_no_preabsorber
    d_phase_space_with_preabsorber
end


struct Spot{T<:Real}
    id::Integer
    t::T
    u::T
    weight::T
    energykeV::T
    # 0 or 1
    numberOfAbsorbers::Int32
end


struct FieldDefinition{T<:Real}
    label::String
    id::Integer
    gantryAngle::T
    couchAngle::T
    nozzleExtraction::T
    fieldCenter::Dict{String, T}
    spots::Vector{Spot{T}}
end


struct TreatmentPlan{T<:Real}
    fields::Vector{FieldDefinition{T}}
    prescribedDose::T
end


# Constructors to build the structures from a JSON representation.
function Spot(spot_dict::Dict)
    T = Float32
    Spot{T}(
        Base.convert(Int32, spot_dict["id"]),
        Base.convert(T, spot_dict["t"]),
        Base.convert(T, spot_dict["u"]),
        Base.convert(T, spot_dict["weight"]),
        Base.convert(T, spot_dict["energykeV"]),
        Base.convert(Int32, spot_dict["numberOfAbsorbers"]),
    )
end


function FieldDefinition(field_dict::Dict)
    T = Float32
    spots = [Spot(spot_dict) for spot_dict in field_dict["spots"]];
    FieldDefinition{T}(
        field_dict["label"],
        field_dict["id"],
        convert(T, field_dict["gantryAngle"]),
        convert(T, field_dict["couchAngle"]),
        convert(T, field_dict["nozzleExtraction"]),
        Base.convert(Dict{String, T}, field_dict["fieldCenter"]),
        spots,
    )
end


function FieldDefinition(label::String,
                         id::Integer,
                         gantryAngle::T,
                         couchAngle::T,
                         nozzleExtraction::T,
                         fieldCenter,
                         spots::Vector{Spot{T}}) where {T<:Real}
    return FieldDefinition(
        label,
        id,
        gantryAngle,
        couchAngle,
        nozzleExtraction,
        Dict{String, T}(
            "x" => fieldCenter[1],
            "y" => fieldCenter[2],
            "z" => fieldCenter[3],
        ),
        spots,
    )
end


function FieldDefinition(label::String,
                         id::Integer,
                         gantryAngle::T,
                         couchAngle::T,
                         nozzleExtraction::T,
                         fieldCenter) where {T<:Real}
    return FieldDefinition(
        label,
        id,
        gantryAngle,
        couchAngle,
        nozzleExtraction,
        fieldCenter,
        Vector{Spot{T}}(undef, 0),
    )
end


function TreatmentPlan(plan_dict::Dict)
    T = Float32
    fields = [FieldDefinition(field_dict) for field_dict in plan_dict["fields"]]
    TreatmentPlan(fields, convert(T, plan_dict["prescribedDose"]))
end


import Base.==

function ==(f::FieldDefinition, f2::FieldDefinition)
    return (f.gantryAngle == f2.gantryAngle) && (f.couchAngle == f2.couchAngle) && (f.nozzleExtraction == f2.nozzleExtraction) && (f.id == f2.id) && (f.label == f2.label) && f.fieldCenter == f2.fieldCenter && (f.spots == f2.spots)
end


function ==(plan::TreatmentPlan, plan2::TreatmentPlan)
    return (plan.prescribedDose == plan2.prescribedDose) && (plan.fields == plan2.fields)
end
