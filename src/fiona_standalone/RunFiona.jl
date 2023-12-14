function run_standalone(fiona_jar_path::String,
                        output_dir::String,
                        log_dij::Bool,
                        log_wed::Bool,
                        path_to_properties::String,
                        mode::String;
                        optimizationPointsPath=nothing)
    if log_dij
        log_dij_command = "-dijLogDir $output_dir"
    else
        log_dij_command = ""
    end

    if log_wed
        log_wed_command = "-wedLogDir $output_dir"
    else
        log_wed_command = ""
    end

    if isnothing(optimizationPointsPath)
        optim_point_command = ""
    else
        optim_point_command = "-optimizationPointsFile $optimizationPointsPath"
    end

    cmd = "java -jar $(fiona_jar_path) -c $(mode) $(log_dij_command) $(log_wed_command) $(optim_point_command) $(path_to_properties)"
    @info cmd
    parts = String.(split(cmd, " "))
    parts = [p for p in parts if !isempty(p)]
    return read(Cmd(parts), String)
end


function write_config_files(output_dir::String,
                            main_config::MainConfig,
                            optimisation_config::OptimizationSettings,
                            spot_placement_config::SpotPlacementConfig)
    write_properties("$output_dir/config.properties", main_config)
    write_optimisation_config("$output_dir/optim_config.json", optimisation_config)
    write_spot_placement_config("$output_dir/spot_config.json", spot_placement_config)
end

function run_spot_placement(fiona_jar_path::String,
                            output_dir::String,
                            log_dij::Bool,
                            log_wed::Bool,
                            main_config::MainConfig,
                            optimisation_config::OptimizationSettings,
                            spot_placement_config::SpotPlacementConfig)
    write_config_files(
        output_dir,
        main_config,
        optimisation_config,
        spot_placement_config,
    )
    open("$output_dir/result_plan.json", "w") do file
        write(file, "{}")
    end
    return run_standalone(
        fiona_jar_path,
        output_dir,
        log_dij,
        log_wed,
        "$output_dir/config.properties",
        "GENERATE_PLAN",
    )
end


function run_optimization(fiona_jar_path::String,
                          output_dir::String,
                          log_dij::Bool,
                          log_wed::Bool,
                          plan::TreatmentPlan;
                          optimizationPointsPath=nothing)
    write_plan_config("$output_dir/result_plan.json", plan)

    return run_standalone(
        fiona_jar_path,
        output_dir,
        log_dij,
        log_wed,
        "$output_dir/config.properties",
        "OPTIMIZE",
        optimizationPointsPath=optimizationPointsPath,
    )
end

function run_optimization(fiona_jar_path::String,
                          output_dir::String,
                          log_dij::Bool,
                          log_wed::Bool,
                          main_config::MainConfig,
                          optimisation_config::OptimizationSettings,
                          spot_placement_config::SpotPlacementConfig,
                          plan::TreatmentPlan;
                          optimizationPointsPath=nothing)
    write_config_files(
        output_dir,
        main_config,
        optimisation_config,
        spot_placement_config,
    )
    write_plan_config("$output_dir/result_plan.json", plan)
    run_optimization(
        fiona_jar_path,
        output_dir,
        log_dij,
        log_wed,
        plan,
        optimizationPointsPath=optimizationPointsPath,
    )
end

function run_dose_calculation(fiona_jar_path::String,
          output_dir::String,
          log_dij::Bool,
          log_wed::Bool,
          main_config::MainConfig,
          optimisation_config::OptimizationSettings,
          spot_placement_config::SpotPlacementConfig)
    write_config_files(
        output_dir,
        main_config,
        optimisation_config,
        spot_placement_config,
    )
    run_dose_calculation(
        fiona_jar_path,
        output_dir,
        log_dij,
        log_wed,
    )
end

function run_dose_calculation(fiona_jar_path::String,
                              output_dir::String,
                              log_dij::Bool,
                              log_wed::Bool)
    return run_standalone(
        fiona_jar_path,
        output_dir,
        log_dij,
        log_wed,
        "$output_dir/config.properties",
        "CALCULATE_DOSE",
    )
end
