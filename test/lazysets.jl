using LinearAlgebra
using LazySets

@testset "lazysets" begin
    pt2 = Singleton([1.0,1.0])
    pt3 = Singleton([1.0,1.0,1.0])
    vec2 = SVector{2}([1.0, 0.0])
    vec3 = SVector{3}([1.0, 0.0, 0.0])

    @testset "Ball1" begin
        b1 = Ball1(zeros(2), 1.0)
        ret = minimum_distance(b1, pt2, vec2)
        @test ret ≈ 0.7071067811865476 rtol=0.001
        b1 = Ball1(zeros(3), 1.0)
        ret = minimum_distance(b1, pt3, vec3)
        @test ret ≈ 1.1547005383792515 rtol=0.001
    end

    @testset "Ball2" begin
        b2 = Ball2(zeros(2), 1.0)
        ret = minimum_distance(b2, pt2, vec2)
        @test ret ≈ 0.41421356237309515 rtol=0.001
        b2 = Ball2(zeros(3), 1.0)
        ret = minimum_distance(b2, pt3, vec3)
        @test ret ≈ 0.7320508075688773 rtol=0.001
    end

    @testset "BallInf" begin
        bi = BallInf(zeros(2), 1.0)
        ret = minimum_distance(bi, pt2, vec2)
        @test ret ≈ 0.0 rtol=0.001
        bi = BallInf(zeros(3), 1.0)
        ret = minimum_distance(bi, pt3, vec3)
        @test ret ≈ 0.0 rtol=0.001
    end

    @testset "Ballp" begin
        bp = Ballp(1.5, zeros(2), 1.0)
        ret = minimum_distance(bp, pt2, vec2)
        @test ret ≈ 0.5233148442327559 rtol=0.001
        bp = Ballp(1.5, zeros(3), 1.0)
        ret = minimum_distance(bp, pt3, vec3)
        @test ret ≈ 0.8993676299132733 rtol=0.001
        bp = Ballp(2.5, zeros(2), 1.0)
        ret = minimum_distance(bp, pt2, vec2)
        @test ret ≈ 0.3424400998368018 rtol=0.001
        bp = Ballp(2.5, zeros(3), 1.0)
        ret = minimum_distance(bp, pt3, vec3)
        @test ret ≈ 0.6159276335349732 rtol=0.001
    end

    #  Not testing anymore because rngs can be different
    #  on different versions of Base and StableRNGs is not
    #  fully supported by LazySets
    #
    #  @testset "Ellipsoid" begin
    #      e = rand(Ellipsoid, dim=2, rng=StableRNG(1))
    #      ret = minimum_distance(e, pt2, vec2)
    #      @test ret ≈ 0.5015857643287766  rtol=0.001
    #      e = rand(Ellipsoid, dim=3, seed=8)
    #      ret = minimum_distance(e, pt3, vec3)
    #      @test ret ≈ 0.17395550015075892 rtol=0.001
    #  end
    #
    #  @testset "HPolygon" begin
    #      hp = rand(HPolygon, seed=2)
    #      ret = minimum_distance(hp, pt2, vec2)
    #      @test ret ≈ 1.099023620263928 rtol=0.001
    #      hp = rand(HPolygonOpt, seed=3)
    #      ret = minimum_distance(hp, pt2, vec2)
    #      @test ret ≈ 0.19296447682913367 rtol=0.001
    #  end
    #
    #  @testset "Hyperrectangle" begin
    #      hr = rand(Hyperrectangle, dim=2, seed=3)
    #      ret = minimum_distance(hr, pt2, vec2)
    #      @test ret ≈ 2.546482967114228 rtol=0.001
    #      hr = rand(Hyperrectangle, dim=3, seed=4)
    #      ret = minimum_distance(hr, pt3, vec3)
    #      @test ret ≈ 2.4101143198871897 rtol=0.001
    #  end
    #
    #  @testset "LineSegment" begin
    #      ls = rand(LineSegment, dim=2, seed=5)
    #      ret = minimum_distance(ls, pt2, vec2)
    #      @test ret ≈ 0.44758226978293086 rtol=0.001
    #  end
    #
    #  @testset "VPolytope" begin
    #      vp = rand(VPolygon, dim=2, seed=10)
    #      ret = minimum_distance(vp, pt2, vec2)
    #      @test ret ≈ 1.2323399249078948 rtol=0.001
    #      vp = rand(VPolytope, dim=2, seed=12)
    #      ret = minimum_distance(vp, pt2, vec2)
    #      @test ret ≈ 0.5871018052878227 rtol=0.001
    #      vp = rand(VPolytope, dim=3, seed=14)
    #      ret = minimum_distance(vp, pt3, vec3)
    #      @test ret ≈ 0.9811157850634189 rtol=0.001
    #  end
    #
    #  @testset "Zonotope" begin
    #      z = rand(Zonotope, dim=2, seed=16)
    #      ret = minimum_distance(z, pt2, vec2)
    #      @test ret ≈ 0.9086568206474619 rtol=0.001
    #      z = rand(Zonotope, dim=3, seed=17)
    #      ret = minimum_distance(z, pt3, vec3)
    #      @test ret ≈ 0.1823104257339578 rtol=0.001
    #  end
end
