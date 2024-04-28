using LineSearches
using Optim
using NPZ


"""
    build_early_stopping(delta::T, patience)

Stops the iteration if the function value has not decreased by more than delta
in the last patience iterations.
"""
function build_early_stopping(delta::T, patience) where {T}
    previous_best = typemax(T)
    previous_best_iteration = 1

    function early_stopping(value, iteration)
        if value <= (previous_best - delta)
            previous_best = value
            previous_best_iteration = iteration
        end
        
        return (iteration - previous_best_iteration) > patience
    end

    return early_stopping
end


function optimise(w::AbstractArray{T, N},
                  config::Juliana.OptimisationConfiguration,
                  subloss_weights::Dict{String, T};
                  oar_safety_margin::T=1.f0,
                  minimum_spot_weight::T=13f0,
                  early_stopping_delta=0.5f0,
                  early_stopping_patience=25,
                  weight_log_directory::String="",
                  weight_log_period=0,
                  maximum_iterations=3000,
                  initial_iteration=1) where {T, N}
    function my_loss(w)
        clamp!(w, minimum_spot_weight, typemax(T))
        loss = Juliana.loss!(
            w,
            config,
            Dict{String, Float32}(),
            subloss_weights,
            oar_safety_margin,
        )

        return loss
    end

    function my_loss_gradient!(gradient, w)
        clamp!(w, minimum_spot_weight, typemax(T))

        grad = Juliana.loss_gradient(w, config, subloss_weights, oar_safety_margin)
        gradient[:] = grad[:]
    end

    # We don't call the Optim.optimize function directly
    # because we want to log the sublosses.
    # See the following issue:
    # https://github.com/JuliaNLSolvers/Optim.jl/issues/1024
    early_stopping = build_early_stopping(
        early_stopping_delta,
        early_stopping_patience,
    )

    alg = LBFGS(linesearch=LineSearches.HagerZhang())
    options = Optim.Options()
    d = Optim.promote_objtype(alg, w, :finite, true, my_loss, my_loss_gradient!)
    state = Optim.initial_state(alg, options, d, w)

    history = Vector{Dict{String, T}}()
    gradients = Array{Vector{T}, 1}()
    loss = Inf
    for i in initial_iteration:maximum_iterations
        # Check for convergence.
        if early_stopping(loss, i)
            break
        end

        # Iterate.
        Optim.update_state!(d, state, alg)
        Optim.update_g!(d, state, alg)
        Optim.update_h!(d, state, alg)

        # Log sublosses.
        loss_parts = Dict{String, T}()
        w = state.x
        clamp!(w, minimum_spot_weight, typemax(T))
        loss = Juliana.loss!(w, config, loss_parts, subloss_weights, oar_safety_margin)
        println(loss)
        loss_parts["total_loss"] = loss
        push!(history, copy(loss_parts))
        grad = Juliana.loss_gradient(w, config, subloss_weights, oar_safety_margin)
        push!(gradients, grad)

        # Dump weights to disk.
        if (weight_log_period > 0) && (mod(i, weight_log_period) == 1)
            npzwrite("$(weight_log_directory)/weights_iteration_$(i).npy", collect(w))
        end
    end
    w_opt = state.x
    clamp!(w_opt, minimum_spot_weight, typemax(T))

    return w_opt, history, gradients
end


function optimise_masked(w::AbstractArray{T, N},
                         config::Juliana.OptimisationConfiguration,
                         subloss_weights::Dict{String, T},
                         weight_mask::AbstractArray{T, N};
                         oar_safety_margin::T=1.f0,
                         minimum_spot_weight::T=13f0,
                         early_stopping_delta=0.5f0,
                         early_stopping_patience=25,
                         weight_log_directory::String="",
                         weight_log_period=0,
                         maximum_iterations=3000,
                         initial_iteration=1) where {T, N}
    function my_loss(w)
        clamp!(w, minimum_spot_weight, typemax(T))
        loss = Juliana.loss!(
            w,
            config,
            Dict{String, Float32}(),
            subloss_weights,
            oar_safety_margin,
        )

        return loss
    end

    function my_loss_gradient!(gradient, w)
        clamp!(w, minimum_spot_weight, typemax(T))

        grad = Juliana.loss_gradient(w, config, subloss_weights, oar_safety_margin)
        gradient[:] = grad[:]
    end

    # We don't call the Optim.optimize function directly
    # because we want to log the sublosses.
    # See the following issue:
    # https://github.com/JuliaNLSolvers/Optim.jl/issues/1024
    early_stopping = build_early_stopping(
        early_stopping_delta,
        early_stopping_patience,
    )

    alg = LBFGS(linesearch=LineSearches.HagerZhang())
    options = Optim.Options()
    d = Optim.promote_objtype(alg, w, :finite, true, my_loss, my_loss_gradient!)
    state = Optim.initial_state(alg, options, d, w)

    history = Vector{Dict{String, T}}()
    gradients = Array{Vector{T}, 1}()
    loss = Inf
    for i in initial_iteration:maximum_iterations
        # Check for convergence.
        if early_stopping(loss, i)
            break
        end

        # Iterate.
        Optim.update_state!(d, state, alg)
        Optim.update_g!(d, state, alg)
        Optim.update_h!(d, state, alg)

        # Log sublosses.
        loss_parts = Dict{String, T}()
        w = state.x
        clamp!(w, minimum_spot_weight, typemax(T))
        w .*= weight_mask
        loss = Juliana.loss!(w, config, loss_parts, subloss_weights, oar_safety_margin)
        loss_parts["total_loss"] = loss
        push!(history, copy(loss_parts))
        grad = Juliana.loss_gradient(w, config, subloss_weights, oar_safety_margin)
        push!(gradients, grad)

        # Dump weights to disk.
        if (weight_log_period > 0) && (mod(i, weight_log_period) == 1)
            npzwrite("$(weight_log_directory)/weights_iteration_$(i).npy", collect(w))
        end
    end
    w_opt = state.x
    clamp!(w_opt, minimum_spot_weight, typemax(T))

    return w_opt, history, gradients
end


"""
    optimise_juliana(tps, ct, structures, prescriptions, plan)

Optimise the spot weights of the spots in the given plan using
the JulianA autoplanning algorithm (https://doi.org/10.48550/arXiv.2305.10211).

This function is meant to be a wrapper for researchers who do not care about the
optimisation algorithm at all and just want an easy to use blackbox function.
"""
function optimise_juliana(tps, ct, structures, prescriptions, plan)
    # Select optimisation points.
    optimisation_mask, optimisation_points, optimisation_point_indices = Juliana.get_optimisation_points_from_prescription(
        ct.grid,
        prescriptions,
        structures,
    )
    # Calculate influence matrix.
    Dij = Juliana.calculate_Dij_matrix(
        tps,
        ct,
        plan,
        optimisation_points,
    )
    # Prepare optimisation.
    initial_weights = Juliana.spot_weights(plan)
    config = Juliana.get_optimisation_configuration(
        ct,
        prescriptions,
        structures,
        Dij,
        optimisation_point_indices,
    )
    w = cu(initial_weights)
    config = Juliana.to_gpu(config)
    subloss_weights = Juliana.get_default_subloss_weights(prescriptions, config.structures)
    # Optimise.
    w_opt, history, gradients = Juliana.optimise(
        w,
        config,
        subloss_weights;
        early_stopping_patience=1,
    )
    # Build the optimised treatment plan.
    w_opt_cpu = collect(w_opt);
    return Juliana.update_spot_weights(plan, w_opt_cpu)
end