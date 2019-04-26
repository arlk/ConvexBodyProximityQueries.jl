@testset "Obstacles" begin
    @test typeof(randpoly(rand(2), n=3, scale=rand())) <: ConvexPolygon{2, 3}
    @test typeof(randpoly(rand(2), n=8, scale=π)) <: ConvexPolygon{2, 8}
    @test typeof(randpoly(rand(2), n=15, scale=rand())) <: ConvexPolygon{2, 15}
    @test typeof(@point rand(2)) <: ConvexPolygon{2, 1}
    @test typeof(@point rand(3)) <: ConvexPolygon{3, 1}
    @test typeof(@line rand(2), rand(2)) <: ConvexPolygon{2, 2}
    @test typeof(@line rand(3), rand(3)) <: ConvexPolygon{3, 2}
    @test typeof(@rect rand(2), rand(2)) <: ConvexPolygon{2, 4}
    @test typeof(@rect rand(3), rand(3)) <: ConvexPolygon{3, 8}
    @test typeof(@square rand(2), π) <: ConvexPolygon{2, 4}
    @test typeof(@square rand(3), rand()) <: ConvexPolygon{3, 8}
end
