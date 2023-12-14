using JSON
import ..xyz_to_index


function write_properties(path::String, config::MainConfig)
    open(path, "w") do file
        for property in fieldnames(MainConfig)
            write(file, "$property=$(getfield(config, property))\n")
        end
    end
end

function write_optimisation_config(path::String, config::OptimizationSettings)
    open(path, "w") do file
        write(file, JSON.json(config, 4))
    end
end

function write_spot_placement_config(path::String, config::SpotPlacementConfig)
    open(path, "w") do file
        write(file, JSON.json(config, 4))
    end
end

function write_plan_config(path::String, config::TreatmentPlan)
    open(path, "w") do file
        write(file, JSON.json(config, 4))
    end
end


function read_plan_file(path::String)
    open(path, "r") do file
        plan_dict = JSON.parse(file)
        return TreatmentPlan(plan_dict)
    end
end


function write_optimisation_points(path, points, grid)
    open(path, "w") do file
        for point in eachrow(points)
            index = xyz_to_index(point, grid)
            # Julia indices start at one, but Fiona standalone starts at zero.
            ix = index[1] - 1
            iy = index[2] - 1
            iz = index[3] - 1
            # We assume column-major layout for Fiona Standalone.
            index_flat = ix + iy * grid.size[1] + iz * grid.size[1] * grid.size[2]
            write(file, "$(point[1])")
            write(file, " ")
            write(file, "$(point[2])")
            write(file, " ")
            write(file, "$(point[3])")
            write(file, " ")
            write(file, "$(index_flat)")
            write(file, "\n")
        end
    end
end