using Juliana
using Test
using Random
Random.seed!(530479)

@testset "DistanceCalculation.jl" begin
    @testset "shortest_distance_to_line()" begin
        # Distance calculation mask.
        rhombus = convert.(Float64, [
            4 0 0;
            0 2 0;
            4 4 0;
            8 2 0;
            4 0 0;
        ])

        p0 = [4., 0.];
        p1 = [0., 2.];
        x = [2., 1.];

        @test isapprox(
            Juliana.shortest_distance_to_line(x, p0, p1),
            0,
            atol=1e-12,
        )

        x = [1., 0.5]
        @test isapprox(
            Juliana.shortest_distance_to_line(x, p0, p1),
            2 * sin(atan(1, 2)),
            atol=1e-12,
        )

        x = [3., 2.]
        @test isapprox(
            Juliana.shortest_distance_to_line(x, p0, p1),
            3 * sin(atan(1.5, 3.)),
            atol=1e-12,
        )
    end

    @testset "shortest_distance_to_polygon()" begin
        # Distance calculation mask.
        rectangle = convert.(Float64, [
            0 0 0;
            2 0 0;
            2 2 0;
            0 2 0;
        ])

        x = [0., 1., 0.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            0,
            atol=1e-12,
        )

        x = [1., 0., 0.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            0,
            atol=1e-12,
        )

        x = [1., 1., 0.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            1.,
            atol=1e-12,
        )

        x = [-1., 0., 0.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            1.,
            atol=1e-12,
        )

        x = [-1., 1., 0.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            1.,
            atol=1e-12,
        )

        x = [-1., -1., 0.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            √(2),
            atol=1e-12,
        )

        x = [-1., -1., 1.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            √(2),
            atol=1e-12,
        )

        x = [-1., -1., 3.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            √(2),
            atol=1e-12,
        )

        x = [3., 1., 0.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            1,
            atol=1e-12,
        )

        x = [3., 3., 0.];
        @test isapprox(
            Juliana.shortest_distance_to_polygon(x, rectangle'),
            √(2),
            atol=1e-12,
        )
    end
end


@testset "CoordinateTransformations.jl" begin
    p = (0.f0, 1.f0, 1.f0)
    iso = (1f0, 2f0, 3f0)

    gantry_angle = convert(Float32, rand() * 180)
    couch_angle = convert(Float32, rand() * 360 - 180)

    s, t, u = Juliana.XyzToStu(p[1], p[2], p[3], iso, gantry_angle, couch_angle)
    p_new = Juliana.StuToXyz(s, t, u, iso, gantry_angle, couch_angle)

    println(all(abs.(p .- p_new) .< 1e-6))
end