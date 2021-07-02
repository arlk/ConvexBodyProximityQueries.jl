using ConvexBodyProximityQueries: support

@testset "support functions" begin
    @testset "array of vertices" begin
        vertices = SMatrix{2, 3}([3.0 5.0 7.0; 13.0 15.0 11.0])
        dir = SVector{2}([1.0, 2.0])
        @test isapprox(support(vertices, dir), [5.0, 15.0], rtol=1e-3)
        dir = SVector{2}([-1.0, 2.0])
        @test isapprox(support(vertices, dir), [5.0, 15.0], rtol=1e-3)
        dir = SVector{2}([-1.0, -2.0])
        @test isapprox(support(vertices, dir), [3.0, 13.0], rtol=1e-3)
        dir = SVector{2}([1.0, -2.0])
        @test isapprox(support(vertices, dir), [7.0, 11.0], rtol=1e-3)

        vertices = SMatrix{3, 4}([3.0 5.0 7.0 5.0; 13.0 15.0 11.0 13.0; 1.0 2.0 3.0 4.0])
        dir = SVector{3}([1.0, 1.0, 2.0])
        @test isapprox(support(vertices, dir), [5.0, 13.0, 4.0], rtol=1e-3)
        dir = SVector{3}([1.0, 1.0, -2.0])
        @test isapprox(support(vertices, dir), [5.0, 15.0, 2.0], rtol=1e-3)
        dir = SVector{3}([1.0, -1.0, -2.0])
        @test isapprox(support(vertices, dir), [7.0, 11.0, 3.0], rtol=1e-3)
        dir = SVector{3}([-1.0, -1.0, -2.0])
        @test isapprox(support(vertices, dir), [3.0, 13.0, 1.0], rtol=1e-3)
    end # array of vertices
end # support functions
