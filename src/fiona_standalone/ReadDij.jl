using SparseArrays


function load_optimisation_points(path::String)
    open(path, "r") do file
        content = JSON.parse(file)
        x_coords = content["xCoordinates"]
        y_coords = content["yCoordinates"]
        z_coords = content["zCoordinates"]

        return hcat(x_coords, y_coords, z_coords)
    end
end


function point_to_grid_index(point::AbstractArray{T, 1}, grid::Grid) where {T}
    return Base.convert.(Int64, round.((point .- grid.origin) ./ grid.spacing))
end


function point_to_grid_index(points::AbstractArray{T, 2}, grid::Grid) where {T}
    indices = Array{Int64, 2}(undef, size(points))
    for (i, point) in enumerate(eachrow(points))
        indices[i, :] .= point_to_grid_index(point, grid)
    end
    return indices
end


function load_Dij(path::String, n_optimisation_points)
    open(path, "r") do file
        buffer = Array{UInt8, 1}(undef, 4)
        readbytes!(file, buffer)
        n_entries = collect(reinterpret(Int32, buffer))[1]

        buffer = Array{UInt8, 1}(undef, 4)
        readbytes!(file, buffer)
        n_spots = collect(reinterpret(Int32, buffer))[1]

        buffer = Array{UInt8, 1}(undef, n_entries * 4)
        readbytes!(file, buffer)
        entries = collect(reinterpret(Float32, buffer))

        buffer = Array{UInt8, 1}(undef, n_entries * 4)
        readbytes!(file, buffer)
        index = collect(reinterpret(Int32, buffer))

        buffer = Array{UInt8, 1}(undef, (n_spots+1) * 4)
        readbytes!(file, buffer)
        spot_stop = collect(reinterpret(Int32, buffer))

        # Sometimes the format is corrupted...
        # Perhaps they implicitly assume a past-the-end spot pointer?
        # TODO: Find out what is going on.
        # For now simply ignore excess entries.
        if (spot_stop[end] - 1) != length(index)
        new_length = spot_stop[end]
        index = index[1:new_length]
        entries = entries[1:new_length]
        end

        Dij = SparseMatrixCSC(
            n_optimisation_points,
            n_spots,
            # Java indexing (Fiona) starts at 0, Julia starts at 1.
            spot_stop.+1,
            index.+1,
            entries,
        )
        return Dij
    end
end
