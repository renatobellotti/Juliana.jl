"""
Calculate the Fiona loss function without OAR sparing.

The magic numbers in the parameters come from Gabriel Meier. They are probably
inspired by the physics of proton therapy.

The original Fiona implementation is in the class
ch.psi.ftpp.gpu.optimization.model.OptimizationTargetPoints in the Fiona repo.
"""
function calculate_ideal_dose_distribution(ct::ScalarGrid,
                                           target_doses::Vector{Tuple{String, T}},
                                           structures::Dict{String, Juliana.Structure};
                                           σ=0.8f0,
                                           falloff_distance=-0.1672349,
                                           target_dose_at_edge_factor=0.9f0,
                                           importance_exponent=10f0) where {T<:Real}
    ideal_dose = zeros(T, size(ct.data))
    importance = zeros(T, size(ct.data))
    target_distr = nothing

    for (target_name, target_dose) in target_doses
        target_dose_at_edge = target_dose * target_dose_at_edge_factor

        # We want to ensure target coverage, so we spill a bit outside of the
        # target. Unit: cm.
        distances = copy(structures[target_name].distanceFromStructure)
        distances .-= falloff_distance

        # Calculate the ideal distribution of this target.
        clamp!(distances, 0, Inf)
        target_distr = target_dose .* exp.(-distances.^2 ./ (2*σ^2))
        target_importance = min.(target_distr ./ target_dose_at_edge, 1f0).^importance_exponent

        # Combine the ideal dose distributions for each target.
        ideal_dose .= max.(ideal_dose, target_distr)
        importance .= max.(importance, target_importance)
    end

    # No importance to air voxels.
    # Here we differ from Fiona: Fiona divides the importance of those voxels by
    # a factor of 10, we set the importance to zero.
    # importance[info.patient.ct.data .< -950] .= 0f0
    return ideal_dose, importance
end