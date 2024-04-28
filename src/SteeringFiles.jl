function write_field_block(file, field, coordinate_type, field_dose; default_tune_path=raw"$TCS_HOME/MachineData/ClinicalDefault.tuneTds")
    n_spots = length(field.spots)
    BEV_lift = 1 # [cm]
    BEV_pan = 0  # [°]
    
    table_coordinate = round.((
        field.fieldCenter["x"],
        field.fieldCenter["y"],
        field.fieldCenter["z"],
    ), digits=2);

    write(file, "FIST[Version=0,Number=0]:")
    write(file, "NAME[STR=\"$(field.label)\"]")
    write(file, ",")
    write(file, "DESC[STR=\"\"]")
    write(file, ",")
    write(file, "FISP[INT=$(n_spots)]")
    write(file, ",")
    write(file, "FIDO[FLT=$(field_dose)]")
    write(file, ",")
    write(file, "NOEX[FLT=$(field.nozzleExtraction)]")
    write(file, ",")
    write(file, "TUFN[STR=\"$(default_tune_path)\"]")
    write(file, ",")
    write(file, "COOR[STR=\"$(coordinate_type)\"]")
    write(file, ",")
    write(file, "STPO[FLTARR=$(field.gantryAngle),$(field.couchAngle),$(table_coordinate[1]),$(table_coordinate[2]),$(table_coordinate[3])]")
    write(file, ",")
    write(file, "BELI[FLT=$(BEV_lift)]")
    write(file, ",")
    write(file, "BEPA[FLT=$(BEV_pan)]")
    write(file, ",")
    # From documentation:
    # S_Nozzle: Nozzle extension in TPS in [cm]
    # This is only used for steering file verification, it is not used for any
    # calculation in SFGen
    write(file, "SNOZ[FLT=$(field.nozzleExtraction)]")

    write(file, "\n")

    # Write the spot blocks.
    for (i, spot) in enumerate(field.spots)
        spot_size = "medium"

        spot_sizes = ["small", "medium", "large"]
        @assert spot_size in spot_sizes

        spot_size_ID = findfirst(==(spot_size), spot_sizes)

        # -1: Julia indices start at 1, but SFGen expects to start at 0.
        write(file, "DOSP[Version=0,Number=$(i-1)]:")
        write(file, "DODA[FLTARR=$(spot.energykeV / 1000),$(Juliana.PREABSORBER_WED * spot.numberOfAbsorbers),$(spot_size_ID),$(spot.weight),$(spot.t),$(spot.u)]")
        write(file, "\n")
    end
end


"""
Specification:
https://codebeamer.psi.ch/svn/PROSCAN_TCS/trunk/Tools/SteeringFileGenerator/Generator/Tools/SFGen/doc/TPS-SFGen-API-Specification.pdf
"""
function write_sfgen_text_input(plan, path, fraction_dose, patient_ID;
                                patient_name::String="",
                                treatment_name="research",
                                # Weight of patient alone in kg.
                                # This is a dummy value right now.
                                patient_weight::Real=10,
                                # Total weight of patient with equipment on table in kg.
                                # This is a dummy value right now.
                                system_weight::Real=10,
                                # If true: Spot dose is in MU.
                                # If false: Spot dose is in number of protons.
                                spot_dose_in_MU::Bool=true,
                                couch_type::String="PatientCouchShort",
                                nozzle_cover_type::String="Undefined",
                                coordinate_type::String="PAT",
                                plan_type::String="Conventional")
    if length(patient_name) == 0
        patient_name = "$(patient_ID)^$(patient_ID)"
    end

    @assert couch_type in [
        "PatientCouchLong",
        "PatientCouchShort",
        "PhantomCouch",
        "DosimetryCouch",
        "MeasuringCouchLong",
        "MeasuringCouchShort",
    ]

    @assert nozzle_cover_type in [
        "NoCover",
        "ShortCover",
        "LongCover",
        "Undefined",
    ]

    @assert coordinate_type in [
        "BEAM",
        "GEO",
        "MACH",
        "PAT",
    ]

    @assert plan_type in [
        "Conventional",
        "Conventional_4D",
        "DAPT_Daily",
        "DAPT_Template",
        "DAPT_Fallback",
    ]

    open(path, "w") do file
        field = plan.fields[1];
        N_spots_fraction = sum([length(f.spots) for f in plan.fields]);
        fraction_number = "1";

        table_coordinate = round.((
            field.fieldCenter["x"],
            field.fieldCenter["y"],
            field.fieldCenter["z"],
        ), digits=2);

        # Write fraction start block.
        # Number: Block ID (starts at zero)
        # Version: Block format version
        write(file, "FRST[Version=0,Number=0]:")
        write(file, "PAID[STR=\"$(patient_ID)\"]")
        write(file, ",")
        write(file, "PANA[STR=\"$(patient_name)\"]")
        write(file, ",")
        write(file, "PLVE[STR=\"$(treatment_name)\"]")
        write(file, ",")
        write(file, "AREA[STR=\"G2\"]")
        write(file, ",")
        write(file, "FRSP[INT=$(N_spots_fraction)]");
        write(file, ",")
        write(file, "FRDO[FLT=$(fraction_dose)]");
        write(file, ",")
        write(file, "SDMU[INT=$(convert(Int32, spot_dose_in_MU))]");
        write(file, ",")
        write(file, "PAWE[FLT=$(patient_weight)]");
        write(file, ",")
        write(file, "SYWE[FLT=$(system_weight)]");
        write(file, ",")
        write(file, "COTY[STR=\"$(couch_type)\"]");
        write(file, ",")
        write(file, "NCTY[STR=\"$(nozzle_cover_type)\"]");
        write(file, ",")
        write(file, "COOR[STR=\"$(coordinate_type)\"]");
        write(file, ",")
        write(file, "FRNB[STR=\"$(fraction_number)\"]");
        write(file, ",")
        write(file, "BEV0[FLTARR=0.0,$(field.couchAngle),$(table_coordinate[1]),$(table_coordinate[2]),$(table_coordinate[3])]");
        write(file, ",")
        write(file, "BEV1[FLTARR=90.0,$(field.couchAngle),$(table_coordinate[1]),$(table_coordinate[2]),$(table_coordinate[3])]");
        write(file, "\n")

        # Write a field block.
        for field in plan.fields
            field_dose = fraction_dose / length(plan.fields)
            write_field_block(file, field, coordinate_type, field_dose)
        end

        # Write fraction end block.
        # Apparently it does not matter if the file length and CRC checksum are wrong,
        # these values are adjusted when converting from TXT to TPS.
        # There is a typo in the documentation: The attribute is actually called "FILD",
        # not "FILE".
        write(file, "FREN[Version=0,Number=0]:FILD[INT=1],CRCD[INT=1747304150]");
    end;
    
    return nothing
end
