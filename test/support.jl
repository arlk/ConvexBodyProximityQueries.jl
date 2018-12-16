using GeometryTypes: Point, Vec
using GeometryTypes: HyperSphere, HyperRectangle, HyperCube
using ConvexBodyProximityQueries: support

@testset "support functions" begin
    @testset "array of vertices" begin
        vertices = SMatrix{2, 3}([3.0 5.0 7.0; 13.0 15.0 11.0])
        dir = SVector{2}([1.0, 2.0])
        @test all(support(vertices, dir) .≈ [5.0, 15.0])
        dir = SVector{2}([-1.0, 2.0])
        @test all(support(vertices, dir) .≈ [5.0, 15.0])
        dir = SVector{2}([-1.0, -2.0])
        @test all(support(vertices, dir) .≈ [3.0, 13.0])
        dir = SVector{2}([1.0, -2.0])
        @test all(support(vertices, dir) .≈ [7.0, 11.0])

        vertices = SMatrix{3, 4}([3.0 5.0 7.0 5.0; 13.0 15.0 11.0 13.0; 1.0 2.0 3.0 4.0])
        dir = SVector{3}([1.0, 1.0, 2.0])
        @test all(support(vertices, dir) .≈ [5.0, 13.0, 4.0])
        dir = SVector{3}([1.0, 1.0, -2.0])
        @test all(support(vertices, dir) .≈ [5.0, 15.0, 2.0])
        dir = SVector{3}([1.0, -1.0, -2.0])
        @test all(support(vertices, dir) .≈ [7.0, 11.0, 3.0])
        dir = SVector{3}([-1.0, -1.0, -2.0])
        @test all(support(vertices, dir) .≈ [3.0, 13.0, 1.0])
    end # array of vertices
end # support functions
