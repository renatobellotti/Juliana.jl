using DataFrames


function update_spot_weights(plan::TreatmentPlan, w::Vector{T}) where {T}
    n_fields = length(plan.fields)
    n_spots_per_field = Array{Int64, 1}(undef, n_fields)

    # Make sure that the spots are ordered by field.
    last_id = -1
    for field in plan.fields
        spot_df = DataFrame(field.spots);
        @assert spot_df[1, "id"] == (last_id + 1)
        @assert sort(spot_df[!, "id"]) == spot_df[!, "id"]
        last_id = spot_df[end, "id"]
    end

    # Create a new TreatmentPlan object for the current spot weights.
    new_fields = Array{FieldDefinition, 1}(undef, n_fields)
    id = 0
    for (i_field, field) in enumerate(plan.fields)
        new_spots = Array{Spot, 1}(undef, length(field.spots))
        for (i_spot, spot) in enumerate(field.spots)
            new_spots[i_spot] = Spot(
                spot.id,
                spot.t,
                spot.u,
                w[id+1],
                spot.energykeV,
                spot.numberOfAbsorbers,
            )
            id += 1
        end
        new_fields[i_field] = FieldDefinition(
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
