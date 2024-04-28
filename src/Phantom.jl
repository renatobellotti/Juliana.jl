using StaticArrays


function build_water_slab()
    ct_data = zeros(Int16, 100, 100, 20)
    fill!(ct_data, -1000)
    ct_data[40:80, 20:80, 5:15] .= 0

    ct = Juliana.ScalarGrid(
        ct_data,
        Juliana.Grid(
            SVector{3}([0.1f0, 0.1f0, 0.1f0]),
            SVector{3}([0.f0, 0.f0, 0.f0]),
            SVector{3}(collect(size(ct_data))),
        ),
    )

    points = convert.(Float32, vcat([[
        3.9 1.9 (i-1)*0.1 i-5 i-5
        7.9 1.9 (i-1)*0.1 i-5 i-5
        7.9 7.9 (i-1)*0.1 i-5 i-5
        3.9 7.9 (i-1)*0.1 i-5 i-5
    ] for i in 5:15]...))

    distance = Juliana.calculate_distance_mask(
        Tuple(ct.grid.size),
        ct.grid.spacing,
        Juliana.split_by_z(points),
    )

    whole_body = Juliana.Structure(
        "whole_body",
        points, 
        Juliana.to_binary_mask(distance),
        distance,
        ct.grid,
        true,
    )
    
    return ct, whole_body
end


function build_alternating_bone_water()
    ct_data = zeros(Int16, 100, 100, 20)
    fill!(ct_data, -1000)

    for iy in 20:80
        if iy % 2 == 0
            ct_data[40:80, iy, 5:15] .= 0
        else
            # http://www.radclass.mudr.org/content/hounsfield-units-scale-hu-ct-numbers
            ct_data[40:80, iy, 5:15] .= 700
        end
    end

    ct = Juliana.ScalarGrid(
        ct_data,
        Juliana.Grid(
            SVector{3}([0.1f0, 0.1f0, 0.1f0]),
            SVector{3}([0.f0, 0.f0, 0.f0]),
            SVector{3}(collect(size(ct_data))),
        ),
    )

    points = convert.(Float32, vcat([[
        4. 2 i*0.1 i-5
        8. 2 i*0.1 i-5
        8. 8 i*0.1 i-5
        4. 8 i*0.1 i-5
    ] for i in 5:15]...))

    distance = Juliana.calculate_distance_mask(
        Tuple(ct.grid.size),
        ct.grid.spacing,
        Juliana.split_by_z(points),
    )

    whole_body = Juliana.Structure(
        "whole_body",
        points, 
        Juliana.to_binary_mask(distance),
        distance,
        ct.grid,
        true,
    )
    
    return ct, whole_body
end


function build_bone_gradient()
    ct_data = zeros(Int16, 100, 100, 20)
    fill!(ct_data, -1000)

    for iy in 20:80
        ct_data[40:80, iy, 5:15] .= convert(
            Int16,
            round(3000. * exp(-(iy - 20.) / 60.)),
        )
    end

    ct = Juliana.ScalarGrid(
        ct_data,
        Juliana.Grid(
            SVector{3}([0.1f0, 0.1f0, 0.1f0]),
            SVector{3}([0.f0, 0.f0, 0.f0]),
            SVector{3}(collect(size(ct_data))),
        ),
    )

    points = convert.(Float32, vcat([[
        4. 2 i*0.1 i-5
        8. 2 i*0.1 i-5
        8. 8 i*0.1 i-5
        4. 8 i*0.1 i-5
    ] for i in 5:15]...))

    distance = Juliana.calculate_distance_mask(
        Tuple(ct.grid.size),
        ct.grid.spacing,
        Juliana.split_by_z(points),
    )

    whole_body = Juliana.Structure(
        "whole_body",
        points, 
        Juliana.to_binary_mask(distance),
        distance,
        ct.grid,
        true,
    )
    
    return ct, whole_body
end


function build_bone_gradient_long()
    ct_data = zeros(Int16, 100, 100, 20)
    fill!(ct_data, -1000)

    for iy in 20:80
        ct_data[40:80, iy, 5:15] .= convert(
            Int16,
            round(3000. * exp(-(iy - 20.) / 60.)),
        )
    end

    ct = Juliana.ScalarGrid(
        ct_data,
        Juliana.Grid(
            SVector{3}([0.2f0, 0.2f0, 0.2f0]),
            SVector{3}([0.f0, 0.f0, 0.f0]),
            SVector{3}(collect(size(ct_data))),
        ),
    )

    points = convert.(Float32, vcat([[
        8.  4 i*0.2 i-5
       16.  4 i*0.2 i-5
       16. 16 i*0.2 i-5
        8. 16 i*0.2 i-5
    ] for i in 5:15]...))

    distance = Juliana.calculate_distance_mask(
        Tuple(ct.grid.size),
        ct.grid.spacing,
        Juliana.split_by_z(points),
    )

    whole_body = Juliana.Structure(
        "whole_body",
        points, 
        Juliana.to_binary_mask(distance),
        distance,
        ct.grid,
        true,
    )
    
    return ct, whole_body
end


function build_water_slab_with_cavity()
    ct_data = zeros(Int16, 200, 200, 60)
    fill!(ct_data, -1000)
    ct_data[80:120, 60:120, 30:35] .= 0
    ct_data[90:110, 80:100, 30:35] .= -1000

    ct = Juliana.ScalarGrid(
        ct_data,
        Juliana.Grid(
            SVector{3}{[0.1f0, 0.1f0, 0.1f0]},
            SVector{3}([0.f0, 0.f0, 0.f0]),
            SVector{3}(collect(size(ct_data))),
        ),
    )

    points = convert.(Float32, vcat([[
         8.  6 i*0.1 i-5
        12.  6 i*0.1 i-5
        12. 12 i*0.1 i-5
         8. 12 i*0.1 i-5
    ] for i in 30:35]...))

    distance = Juliana.calculate_distance_mask(
        Tuple(ct.grid.size),
        ct.grid.spacing,
        Juliana.split_by_z(points),
    )

    whole_body = Juliana.Structure(
        "whole_body",
        points, 
        Juliana.to_binary_mask(distance),
        distance,
        ct.grid,
        true,
    )
    
    return ct, whole_body
end
