using Juliana
using Test

@testset "Juliana.jl" begin
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
