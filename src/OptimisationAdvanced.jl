using DataFrames
using NPZ


function is_fulfilled(D_in, ideal_dose_distribution, D_presc, d_margin, inside_trans_mask, outside_mask, target_dose, σ=0.8f0)
    A = exp(- d_margin^2 / (2 * σ^2))
    
    N_in = sum(inside_trans_mask)
    N_out = sum(outside_mask)
    lhs = N_in * Juliana.mean_dose(ideal_dose_distribution .* (D_in / target_dose), convert.(Float32, inside_trans_mask)) + D_in * A * N_out
    return lhs <= (N_in + N_out) * D_presc
end


function find_Din(ideal_dose_distribution, D_presc, d_margin, inside_trans_mask, outside_mask, target_dose, σ=0.8f0)
    for D_in in target_dose:-0.5f0:0
        if is_fulfilled(D_in, ideal_dose_distribution, D_presc, d_margin, inside_trans_mask, outside_mask, target_dose, σ)
            return D_in
        end
    end
    @error "Could not find a reasonable D_in!!!"
end


function calculate_prescribed_outside_dose(ideal_dose_distribution, D_presc, d_margin, inside_trans_mask, outside_mask, target_dose, σ=0.8f0)
    D_in = find_Din(
        ideal_dose_distribution,
        D_presc,
        d_margin,
        inside_trans_mask,
        outside_mask,
        target_dose,
        σ,
    )
    A = exp(- d_margin^2 / (2 * σ^2))
    return A * D_in
end


@doc raw"""
    add_helper_constraints!(ct,
                            prescriptions,
                            structures,
                            OUTSIDE_FIRST_OPTIM_DISTANCE,
                            AUTO_STRUCTURE_PREFIX)

Add a helper structure for the OAR region outside of the target (plus OUTSIDE_FIRST_OPTIM_DISTANCE)
to help the optimiser fulfill mean dose constraints without unnecessarily compromising
target coverage.


**Dose model for mean dose constraint structures**

There are three regions:

- Inside
- transition (i. e. outside the target but inside the margin around it)
- outside (i. e. outside the target plus margin).

In the inside and transition region, we use the ideal dose distribution normalised to $D_{in}$.

Outside, we aim to impose mean dose of $x$.

$$\begin{align}
 D_{mean} &= \frac{1}{N} \sum_i d_i\\
          &= \frac{1}{N} \left( N_{in, trans} \mu_{in, trans} + N_{out} \mu_{out}\right)\\
          &= \frac{1}{N} \left( N_{in, trans} \mu_{in, trans} + N_{out} x(D_{in})\right)\\
          &\leq D_{presc}
\end{align}$$

In the transition region, we impose a Gaussian falloff:

$$\begin{align}
    f(d) = D_{in} exp\left( -\frac{d^2}{2\sigma^2} \right)
\end{align}$$

We impose boundary conditions:

$$\begin{align}
    f(d_{margin}) = D_{in} exp\left( -\frac{d_{margin}^2}{2\sigma^2} \right) = x
\end{align}$$

Renaming constants, we obtain $x = A D_{in}$, where $A := exp\left( -\frac{d_{margin}^2}{2\sigma^2} \right)$.

Rearranging terms in the first inequality, we obtain the final bound to be fulfilled:

$$N_{in, trans} \mu_{in, trans} + D_{in} A N_{out} \leq N D_{presc}$$

Higher D_{in} means better target coverage, therefore we would like to have $D_{in}$ as large as possible.

$$max_{D_{in} \in [0, 100\%]} D_{in} \quad \mathrm{s.~t.~the~condition~is~fulfilled}$$
"""
function add_helper_constraints!(ct,
                                 prescriptions,
                                 structures,
                                 OUTSIDE_FIRST_OPTIM_DISTANCE,
                                 AUTO_STRUCTURE_PREFIX)
    ideal_dose_distribution, _ = Juliana.calculate_ideal_dose_distribution(
        ct,
        prescriptions.target_doses,
        structures,
    )

    target_distances = Juliana.calculate_whole_target_distance(prescriptions, structures)
    target_mask_expanded = target_distances .<= OUTSIDE_FIRST_OPTIM_DISTANCE
    target_dose = Juliana.coldest_target(prescriptions)[2]

    for constr in prescriptions.constraints
        if (constr.priority != Juliana.soft) && (constr.kind == Juliana.constraint_mean)
            # Split the OAR into in inside and an outside/transition region.
            structure = structures[constr.structure_name]
            inside_trans_mask = structure.mask .&& convert.(Bool, target_mask_expanded)
            outside_mask = structure.mask .&& .!(convert.(Bool, target_mask_expanded))
            if sum(outside_mask) == 0
                continue
            end

            # Create an empty distance mask for the optimisation structure.
            nx, ny, nz = structure.grid.size
            dist = Array{Float32, 3}(undef, nx, ny, nz)
            fill!(dist, NaN)

            # Build the optimisation structure.
            outside_part = Juliana.Structure(
                "$(AUTO_STRUCTURE_PREFIX)$(structure.name)",
                Array{Float32, 3}(undef, 0, 0, 0),
                outside_mask,
                dist,
                structure.grid,
                structure.is_target,
            )

            # Build the constraint for the optimisation structure.
            if sum(outside_mask) != sum(structure.mask)
                println("Add helper constraint for outside part of: $(structure.name)")
                structures[outside_part.name] = outside_part
                dose = calculate_prescribed_outside_dose(
                    ideal_dose_distribution,
                    constr.dose,
                    OUTSIDE_FIRST_OPTIM_DISTANCE,
                    inside_trans_mask,
                    outside_mask,
                    target_dose,
                )
                println("Original constraint: $(constr.dose) New: $(dose)")
                outside_constr = Juliana.Constraint(
                    outside_part.name,
                    constr.kind,
                    dose,
                    constr.volume,
                    constr.priority,
                    constr.direction,
                )
                push!(prescriptions.constraints, outside_constr)
                println(outside_constr)
            end
        end
    end
end


function optimise_head_and_neck(config,
                                output_dir,
                                w_init,
                                oar_safety_margin,
                                early_stopping_delta,
                                early_stopping_patience,
                                weight_log_period,
                                maximum_iterations,
                                AUTO_STRUCTURE_PREFIX)
    weight_log_dir = "$(output_dir)/weights"
    mkpath(weight_log_dir)
    
    w = copy(w_init)

    ##############################################################
    # Phase 1: Optimise target coverage using prior knowledge
    ##############################################################
    # Optimise for Gaussian dose falloff.
    subloss_weights = Dict{String, Float32}(
        "ideal_dose_loss" => 1.f0,
        "maximum_loss" => 0.f0,
        "maximum_distribution_loss" => 0.f0,
        "minimum_loss" => 0.f0,
        "normalisation_variance" => 0.f0,
    )

    for constraint in config.prescriptions.constraints
        if constraint.priority == Juliana.soft
            continue
        end

        # Skip OARs for which there is no optimisation point.
        if sum(config.structures[constraint.structure_name])== 0
            continue
        end

        if constraint.kind == Juliana.constraint_mean
            subloss_weights["$(constraint.structure_name)_mean_loss"] = 0f0
        elseif Juliana.is_maximum_constraint(constraint)
            subloss_weights["$(constraint.structure_name)_max_loss"] = 0f0
        end
    end

    w_opt, history, gradients = Juliana.optimise(
        w,
        config,
        subloss_weights,
        oar_safety_margin=oar_safety_margin,
        early_stopping_delta=early_stopping_delta,
        early_stopping_patience=early_stopping_patience,
        weight_log_directory=weight_log_dir,
        weight_log_period=weight_log_period,
        maximum_iterations=maximum_iterations,
    )
    
    n_iterations_so_far = size(history, 1)

    open("$(output_dir)/n_iterations_phase1.txt", "w") do file
        write(file, string(n_iterations_so_far))
    end
    
    w_opt_phase1 = collect(copy(w_opt))
    npzwrite("$output_dir/w_opt_phase1.npy", w_opt_phase1)
    npzwrite("$output_dir/weights/weights_iteration_$(n_iterations_so_far).npy", w_opt_phase1)
    
    ###################################################################################
    # Phase 2: Optimise also Dmax constraints and Dmean parts outside of the target.
    ###################################################################################
    subloss_weights["maximum_loss"] = 1.f0
    subloss_weights["maximum_distribution_loss"] = 1.f0
    subloss_weights["minimum_loss"] = 1.f0
    subloss_weights["normalisation_variance"] = 1.f0

    for constraint in config.prescriptions.constraints
        if constraint.priority == Juliana.soft
            continue
        end

        # Skip OARs for which there is no optimisation point.
        if sum(config.structures[constraint.structure_name])== 0
            continue
        end

        if (constraint.kind == Juliana.constraint_mean) .&& startswith(constraint.structure_name, AUTO_STRUCTURE_PREFIX)
            subloss_weights["$(constraint.structure_name)_mean_loss"] = 30f0
        elseif Juliana.is_maximum_constraint(constraint)
            subloss_weights["$(constraint.structure_name)_max_loss"] = 30f0
        end
    end

    w_opt, history2, gradients2 = Juliana.optimise(
        w_opt,
        config,
        subloss_weights,
        oar_safety_margin=oar_safety_margin,
        early_stopping_delta=early_stopping_delta,
        early_stopping_patience=early_stopping_patience,
        weight_log_directory=weight_log_dir,
        weight_log_period=weight_log_period,
        maximum_iterations=maximum_iterations,
        initial_iteration=n_iterations_so_far,
    )

    n_iterations_so_far += size(history2, 1)

    open("$(output_dir)/n_iterations_phase2.txt", "w") do file
        write(file, string(n_iterations_so_far))
    end
    
    w_opt_phase2 = collect(copy(w_opt))
    npzwrite("$output_dir/w_opt_phase2.npy", w_opt_phase2)
    npzwrite("$output_dir/weights/weights_iteration_$(n_iterations_so_far).npy", w_opt_phase2)

    ##########################################################
    # Phase 3: Optimise for all Dmax and Dmean constraints.
    ##########################################################
    for constraint in config.prescriptions.constraints
        if constraint.priority == Juliana.soft
            continue
        end

        # Skip OARs for which there is no optimisation point.
        if sum(config.structures[constraint.structure_name])== 0
            continue
        end

        if (constraint.kind == Juliana.constraint_mean)
            subloss_weights["$(constraint.structure_name)_mean_loss"] = 30f0
        end
    end

    w_opt, history3, gradients3 = Juliana.optimise(
        cu(w_opt_phase2),
        config,
        subloss_weights,
        oar_safety_margin=oar_safety_margin,
        early_stopping_delta=early_stopping_delta,
        early_stopping_patience=early_stopping_patience,
        weight_log_directory=weight_log_dir,
        weight_log_period=weight_log_period,
        maximum_iterations=maximum_iterations,
        initial_iteration=n_iterations_so_far,
    )

    n_iterations_so_far += size(history3, 1)

    open("$(output_dir)/n_iterations_phase3.txt", "w") do file
        write(file, string(n_iterations_so_far))
    end

    w_opt_phase3 = collect(copy(w_opt))
    npzwrite("$output_dir/w_opt_phase3.npy", w_opt_phase3)
    npzwrite("$output_dir/weights/weights_iteration_$(n_iterations_so_far).npy", w_opt_phase3)

    #########################
    # Postprocess results.
    #########################
    # Ensure range of spot weights and copy to the CPU.
    w_opt = clamp!(w_opt, 0., typemax(Float32))
    

    gradients = collect(hcat(gradients...)')

    # Code taken from: https://stackoverflow.com/a/54170025
    loss_df = vcat(DataFrame.(history)...)
    loss_df2 = vcat(DataFrame.(history2)...)
    loss_df3 = vcat(DataFrame.(history3)...)
    loss_df = vcat(loss_df, loss_df2, loss_df3)
    
    return w_opt, gradients, loss_df
end


function expand_Dmax_constraint_structures!(structures, prescriptions, oar_expansion)
    # Expand structures with Dmax constraints to avoid single-voxel Dmax violations.
    to_expand = Set()
    for constr in prescriptions.constraints
        if (constr.priority != Juliana.soft) && Juliana.is_maximum_constraint(constr)
            push!(to_expand, constr.structure_name)
        end
    end
    for name in to_expand
        structures[name] = Juliana.expand_structure(structures[name], oar_expansion)
    end

    # println("The following structures have been expanded by $(oar_expansion)cm:")
    # println(to_expand)
end


function optimise_head_and_neck(ct,
                                structures,
                                prescriptions,
                                plan,
                                tps,
                                output_dir;
                                robust=false,
                                oar_safety_margin=1.0f0, # Gy
                                oar_expansion=0.1, # mm
                                target_hotspot_percent=103f0,
                                early_stopping_delta=0.25f0,
                                early_stopping_patience=25,
                                weight_log_period=10,
                                maximum_iterations=3000)
    AUTO_STRUCTURE_PREFIX = "AUTO_OPTIM_"
    OUTSIDE_FIRST_OPTIM_DISTANCE = 0.2f0

    # Preprocess the structures and prescriptions.
    prescriptions = deepcopy(prescriptions)
    structures = deepcopy(structures)

    add_helper_constraints!(
        ct,
        prescriptions,
        structures,
        OUTSIDE_FIRST_OPTIM_DISTANCE,
        AUTO_STRUCTURE_PREFIX,
    )
    expand_Dmax_constraint_structures!(
        structures,
        prescriptions,
        oar_expansion,
    )

    # Determine optimisation points.
    optimisation_mask, optimisation_points, optimisation_point_indices = Juliana.get_optimisation_points_from_prescription(
        ct.grid,
        prescriptions,
        structures,
        checkerboard_skip_n_inslice=5,
        checkerboard_skip_n_slices=1,
        margin_skip_n_slices=1,
    )
    npzwrite("$(output_dir)/optimisation_mask.npy", optimisation_mask)

    println("Optimising for ", sum(optimisation_mask))

    # Calculate the Dij matrix.
    start = time()
    if robust
        @time ct_scenarios = Juliana.generate_scenarios(ct);

        @info "Calculate Dij for all $(length(ct_scenarios)) scenarios..."
        Dij_scenarios = [
            Juliana.calculate_Dij_matrix(
                tps,
                ct,
                plan,
                optimisation_points,
            ) for ct in ct_scenarios
        ]

        @info "Build the robust Dij matrix from the scenario Dij matrices."
        Dij = Juliana.build_robust_Dij(
            Dij_scenarios,
            prescriptions,
            structures,
            optimisation_point_indices,
        )

        @info "Releasing memory for the robust scenarios."
        ct_scenarios = nothing
        Dij_scenarios = nothing
        GC.gc()
    else
        Dij = Juliana.calculate_Dij_matrix(
            tps,
            ct,
            plan,
            optimisation_points,
        )
    end
    stop = time()

    open("$(output_dir)/time_Dij_calculation.txt", "w") do file
        write(file, "$(stop - start)s")
    end

    # Prepare the optimisation.
    config = Juliana.get_optimisation_configuration(
        ct,
        prescriptions,
        structures,
        Dij,
        optimisation_point_indices,
        oar_safety_margin,
        target_hotspot_percent,
    )
    config = Juliana.to_gpu(config)

    # Initialise the spot weights.
    spots = vcat([DataFrame(field.spots) for field in plan.fields]...)
    # Ensure that the spots are ordered correctly.
    @assert spots.id == collect(0:size(spots, 1)-1)
    initial_weights = convert.(Float32, spots.weight)

    w = cu(copy(initial_weights))
    tmp = Juliana.reproducible_sparse_mv(config.Dij, w)
    mean_dose = sum(collect(tmp) .* collect(config.normalisationStructureMask)) / sum(config.normalisationStructureMask)
    w *= config.normalisationDose / mean_dose;

    npzwrite("$(output_dir)/initial_weights.npy", initial_weights)

    # Optimise.
    w_opt, gradients, loss_df = optimise_head_and_neck(
        config,
        output_dir,
        w,
        oar_safety_margin,
        early_stopping_delta,
        early_stopping_patience,
        weight_log_period,
        maximum_iterations,
        AUTO_STRUCTURE_PREFIX,
    )

    return collect(w_opt), gradients, loss_df
end
