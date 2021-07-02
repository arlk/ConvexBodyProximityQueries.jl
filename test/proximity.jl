using GeometryBasics: Point, Vec
using GeometryBasics: HyperSphere, HyperRectangle

import ConvexBodyProximityQueries.support

function ConvexBodyProximityQueries.support(sphere::HyperSphere{N, T}, dir::AbstractVector) where {N, T}
    SVector{N}(sphere.center + sphere.r*normalize(dir, 2))
end

@generated function ConvexBodyProximityQueries.support(rect::HyperRectangle{N, T}, dir::AbstractVector) where {N, T}
    exprs = Array{Expr}(undef, (N,))
    for i = 1:N
        exprs[i] = :(rect.widths[$i]*(dir[$i] ≥ 0.0 ? 1.0 : -1.0)/2.0 + rect.origin[$i])
    end

    return quote
        Base.@_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return SVector{N, T}(elements)
    end
end

@testset "gjk (2D)" begin
    vec2 = SVector{2}([1.0, 0.0])
    @testset "array of vertices" begin
        rect2 = SMatrix{3,4}([0.0 1.0 1.0 0.0; 0.0 0.0 2.0 2.0; 1.0 1.0 1.0 1.0])
        transform2(θ,x,y) = SMatrix{2,3}([cos(θ) -sin(θ) x; sin(θ) cos(θ) y])

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(0,3,0) * rect2
        ret = closest_points(polyA, polyB, vec2)
        @test collision_detection(polyA, polyB, vec2) == false
        @test isapprox(ret[1], [1.0, 0.0], rtol=1e-3)
        @test isapprox(ret[2], [3.0, 0.0], rtol=1e-3)

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(0,1,0) * rect2
        ret = closest_points(polyA, polyB, vec2)
        @test collision_detection(polyA, polyB, vec2) == true

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(0,0.5,0) * rect2
        ret = closest_points(polyA, polyB, vec2)
        @test collision_detection(polyA, polyB, vec2) == true

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(π/4,3,0) * rect2
        ret = closest_points(polyA, polyB, vec2)
        @test collision_detection(polyA, polyB, vec2) == false
        @test isapprox(ret[1], [1.0, 1.4142135623730951], rtol=1e-3)
        @test isapprox(ret[2], [1.585786437626905, 1.4142135623730954], rtol=1e-3)

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(π/4,1,0) * rect2
        ret = closest_points(polyA, polyB, vec2)
        @test collision_detection(polyA, polyB, vec2) == true

        polyA = transform2(0,0,0) * rect2
        polyB = transform2(π/4,0.5,0) * rect2
        ret = closest_points(polyA, polyB, vec2)
        @test collision_detection(polyA, polyB, vec2) == true

        polyA = transform2(-π/4,0,0) * rect2
        polyB = transform2(π/4,3,0) * rect2
        ret = closest_points(polyA, polyB, vec2)
        @test collision_detection(polyA, polyB, vec2) == false
        @test isapprox(ret[1], [1.914213562373095, 0.9142135623730953], rtol=1e-3)
        @test isapprox(ret[2], [2.0, 1.0], rtol=1e-3)

        polyA = transform2(-π/4,0,0) * rect2
        polyB = transform2(π/4,1,0) * rect2
        ret = closest_points(polyA, polyB, vec2)
        @test collision_detection(polyA, polyB, vec2) == true

        polyA = transform2(-π/4,0,0) * rect2
        polyB = transform2(π/4,0.5,0) * rect2
        ret = closest_points(polyA, polyB, vec2)
        @test collision_detection(polyA, polyB, vec2) == true
    end # array of vertices

    @testset "sphere" begin
        sphereA = HyperSphere(Point(0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(3.0, 0.0), 1.0)
        ret = closest_points(sphereA, sphereB, vec2)
        @test collision_detection(sphereA, sphereB, vec2) == false
        @test isapprox(ret[1], [1.0, 0.0], rtol=1e-3)
        @test isapprox(ret[2], [2.0, 0.0], rtol=1e-3)

        sphereA = HyperSphere(Point(0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(2.0, 0.0), 1.0)
        ret = closest_points(sphereA, sphereB, vec2)
        @test collision_detection(sphereA, sphereB, vec2) == true

        sphereA = HyperSphere(Point(0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(1.0, 0.0), 1.0)
        ret = closest_points(sphereA, sphereB, vec2)
        @test collision_detection(sphereA, sphereB, vec2) == true
    end # sphere

    @testset "rectangle" begin
        rectangleA = HyperRectangle(Vec(0.0, 0.0), Vec(1.0, 2.0))
        rectangleB = HyperRectangle(Vec(3.0, 0.0), Vec(1.0, 2.0))
        ret = closest_points(rectangleA, rectangleB, vec2)
        @test collision_detection(rectangleA, rectangleB, vec2) == false
        @test isapprox(ret[1], [0.5, 1.0], rtol=1e-3)
        @test isapprox(ret[2], [2.5, 1.0], rtol=1e-3)

        rectangleA = HyperRectangle(Vec(0.0, 0.0), Vec(1.0, 2.0))
        rectangleB = HyperRectangle(Vec(1.0, 0.0), Vec(1.0, 2.0))
        ret = closest_points(rectangleA, rectangleB, vec2)
        @test collision_detection(rectangleA, rectangleB, vec2) == true

        rectangleA = HyperRectangle(Vec(0.0, 0.0), Vec(1.0, 2.0))
        rectangleB = HyperRectangle(Vec(0.5, 0.0), Vec(1.0, 2.0))
        ret = closest_points(rectangleA, rectangleB, vec2)
        @test collision_detection(rectangleA, rectangleB, vec2) == true
    end # rectangle
end # gjk

@testset "gjk (3D)" begin
    vec3 = SVector{3}([1.0, 0.0, 0.0])
    @testset "array of vertices" begin
        # TODO (add several tetrahedron configurations)
    end # array of vertices

    @testset "sphere" begin
        sphereA = HyperSphere(Point(0.0, 0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(3.0, 3.0, 0.0), 1.0)
        ret = closest_points(sphereA, sphereB, vec3)
        @test collision_detection(sphereA, sphereB, vec3) == false
        @test isapprox(ret[1], [0.7071067841796789, 0.7071067781934159, 0.0], rtol=1e-3)
        @test isapprox(ret[2], [2.2928932158203210, 2.2928932218065840, 0.0], rtol=1e-3)

        sphereA = HyperSphere(Point(0.0, 0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(2.0, 0.0, 0.0), 1.0)
        ret = closest_points(sphereA, sphereB, vec3)
        @test collision_detection(sphereA, sphereB, vec3) == true

        sphereA = HyperSphere(Point(0.0, 0.0, 0.0), 1.0)
        sphereB = HyperSphere(Point(1.0, 0.0, 0.0), 1.0)
        ret = closest_points(sphereA, sphereB, vec3)
        @test collision_detection(sphereA, sphereB, vec3) == true
    end # sphere

    @testset "rectangle" begin
        rectangleA = HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(1.0, 2.0, 3.0))
        rectangleB = HyperRectangle(Vec(3.0, 3.0, 0.0), Vec(1.0, 2.0, 3.0))
        ret = closest_points(rectangleA, rectangleB, vec3)
        @test collision_detection(rectangleA, rectangleB, vec3) == false
        @test isapprox(ret[1], [0.5, 1.0, 1.5], rtol=1e-3)
        @test isapprox(ret[2], [2.5, 2.0, 1.5], rtol=1e-3)

        rectangleA = HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(1.0, 2.0, 3.0))
        rectangleB = HyperRectangle(Vec(1.0, 2.0, 0.0), Vec(1.0, 2.0, 3.0))
        ret = closest_points(rectangleA, rectangleB, vec3)
        @test collision_detection(rectangleA, rectangleB, vec3) == true

        rectangleA = HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(1.0, 2.0, 3.0))
        rectangleB = HyperRectangle(Vec(0.5, 0.0, 0.0), Vec(1.0, 2.0, 3.0))
        ret = closest_points(rectangleA, rectangleB, vec3)
        @test collision_detection(rectangleA, rectangleB, vec3) == true
    end # rectangle

end # gjk
