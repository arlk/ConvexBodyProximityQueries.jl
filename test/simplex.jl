using ConvexBodyProximityQueries: proj
@testset "projection operation" begin
    a⃗ = SVector{2}([1.0, 2.0])
    b⃗ = SVector{2}([-1.0, 1.0])
    @test all(proj(a⃗,b⃗) .≈ [0.2, 0.4])
    @test all(proj(b⃗,a⃗) .≈ [-0.5, 0.5])
    a⃗ = SVector{3}([1.0, 2.0, 2.0])
    b⃗ = SVector{3}([-1.0, 1.0, -2.0])
    @test all(proj(a⃗,b⃗) .≈ [-1.0, -2.0, -2.0]/3)
    @test all(proj(b⃗,a⃗) .≈ [0.5, -0.5, 1.0])
end # projection operation

#  TODO: fix tests
#  using GJK: findsimplex
#  @testset "find simplex" begin
#      @testset "find line simplex" begin
#          spx = SMatrix{2,2}([1.0 -1.0; 2.0 1.0])
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [1, 2])
#          @test all(ret[2] .≈ [0.6, -1.2])
#          @test ret[3] == false
#          spx = SMatrix{2,2}([1.0 0.5; 2.0 1.5])
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [2])
#          @test all(ret[2] .≈ [-0.5, -1.5])
#          @test ret[3] == false
#
#          spx = SMatrix{3,2}([1.0 -1.0; 2.0 1.0; 2.0 -2.0])
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [1, 2])
#          @test all(ret[2] .≈ [0.1428571428571429, -1.4285714285714286, 0.2857142857142858])
#          @test ret[3] == false
#          spx = SMatrix{3,2}([1.0 0.0; 2.0 1.5; 2.0 1.0])
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [2])
#          @test all(ret[2] .≈ [0.0, -1.5, -1.0])
#          @test ret[3] == false
#
#          @testset "degeneracies" begin
#              rvec = rand(2, 1) .+ 1.0
#              spx = [rvec -rvec]
#              ret = findsimplex(spx)
#              @test ret[3] .== true
#              rvec = rand(2, 1) .+ 1.0
#              spx = [rvec zero(rvec)]
#              ret = findsimplex(spx)
#              @test ret[3] .== true
#              rvec = rand(2, 1) .+ 1.0
#              spx = [zero(rvec) rvec]
#              ret = findsimplex(spx)
#              @test ret[3] .== true
#
#              rvec = rand(3, 1) .+ 1.0
#              spx = [rvec -rvec]
#              ret = findsimplex(spx)
#              @test ret[3] .== true
#              rvec = rand(3, 1) .+ 1.0
#              spx = [rvec zero(rvec)]
#              ret = findsimplex(spx)
#              @test ret[3] .== true
#              rvec = rand(3, 1) .+ 1.0
#              spx = [zero(rvec) rvec]
#              ret = findsimplex(spx)
#              @test ret[3] .== true
#          end # degneracies
#      end # find line simplex
#
#      @testset "find triangle simplex" begin
#          spx = [1.0 -1.0 0.0; 2.0 1.0 -1.0]
#          ret = findsimplex(spx)
#          @test ret[3] == true
#          spx = [0.0 1.0 -1.0; 3.0 2.0 1.0]
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [2, 3])
#          @test all(ret[2] .≈ [0.6, -1.2])
#          @test ret[3] == false
#          spx = [1.0 0.0 -1.0; 2.0 3.0 1.0]
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [1, 3])
#          @test all(ret[2] .≈ [0.6, -1.2])
#          @test ret[3] == false
#          spx = [1.0 -1.0 0.0; 2.0 1.0 1.0]
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [3])
#          @test all(ret[2] .≈ [0.0, -1.0])
#          @test ret[3] == false
#
#          spx = [1.0 -1.0 0.0; 2.0 1.0 -1.0; 2.0 1.0 1.0]
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [2, 1, 3])
#          @test all(ret[2] .≈ [0.4, 0.2, -1.0])
#          @test ret[3] == false
#          spx = [0.0 1.0 -1.0; 3.0 2.0 1.0; 2.0 1.0 1.0]
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [2, 3])
#          @test all(ret[2] .≈ [0.6, -1.2, -1.0])
#          @test ret[3] == false
#          spx = [1.0 -1.0 0.0; 2.0 1.0 1.0; 2.0 1.0 0.5]
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [3])
#          @test all(ret[2] .≈ [0.0, -1.0, -0.5])
#          @test ret[3] == false
#
#          @testset "degeneracies" begin
#              x = rand() + 1.0
#              y = rand() + 1.0
#              z = rand() + 1.0
#              spx = [x 0 -x; 0 y -y; z z -2z]
#              ret = findsimplex(spx)
#              @test ret[3] == true
#              spx = [x 0 0; 0 y 0; z -z 0]
#              ret = findsimplex(spx)
#              @test ret[3] == true
#
#              spx = [1.0 -1.0 0.0; -0.5 0.5 0.0]
#              ret = findsimplex(spx)
#              @test ret[3] == true
#              spx = [1.0 1.0 1.0; -3.0 2.0 -0.5]
#              ret = findsimplex(spx)
#              @test all(ret[1] .== [2, 3])
#              @test all(ret[2] .≈ [-1.0, 0.0])
#              @test ret[3] == false
#              spx = [1.0 1.0 1.0; -3.0 -2.0 -0.5]
#              ret = findsimplex(spx)
#              @test all(ret[1] .== [3])
#              @test all(ret[2] .≈ [-1.0, 0.5])
#              @test ret[3] == false
#
#              spx = [1.0 -1.0 0.0; -0.5 0.5 0.0; 1.0 -1.0 0.0]
#              ret = findsimplex(spx)
#              @test ret[3] == true
#              spx = [1.0 0.0 0.4; 0.0 1.0 0.6; 1.0 -1.0 -0.2]
#              ret = findsimplex(spx)
#              @test all(ret[1] .== [3])
#              @test all(ret[2] .≈ [-0.4, -0.6, 0.2])
#              @test ret[3] == false
#              spx = [1.0 0.0 0.5; 0.0 1.0 0.5; 1.0 0.0 0.5]
#              ret = findsimplex(spx)
#              @test all(ret[1] .== [2, 3])
#              @test all(ret[2] .≈ [-1.0, -2.0, -1.0]/3.0)
#              @test ret[3] == false
#          end # degneracies
#      end # find triangle simplex
#
#      @testset "find tetrahedron simplex" begin
#          spx = [1.0 0.0 -1.0 0.0; 2.0 -1.0 1.0 0.0; 2.0 1.0 1.0 -2.0]
#          ret = findsimplex(spx)
#          @test ret[3] == true
#          spx = [1.0 0.0 -1.0 0.0; 2.0 -1.0 1.0 0.0; 2.0 1.0 1.0 1.0]
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [2, 3, 4])
#          @test all(ret[2] .≈ [0.0, 0.0, -1.0])
#          @test ret[3] == false
#          spx = [-2.0 -2.0 2.0 1.0; -1.0 1.0 0.0 0.0; 1.0 1.0 1.0 0.0]
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [1, 2, 4])
#          @test all(ret[2] .≈ [-0.1, 0.0, -0.3])
#          @test ret[3] == false
#          spx = [-2.0 2.0 -2.0 1.0; -1.0 0.0 1.0 0.0; -1.0 -1.0 -1.0 0.0]
#          ret = findsimplex(spx)
#          @test all(ret[1] .== [3, 1, 4])
#          @test all(ret[2] .≈ [-0.1, 0.0, 0.3])
#          @test ret[3] == false
#
#          @testset "degeneracies" begin
#              x = rand() + 1.0
#              y = rand() + 1.0
#              z = rand() + 1.0
#              spx = [x 0 -x 0; 0 y -y 0; z z -z 0]
#              ret = findsimplex(spx)
#              @test ret[3] == true
#          end # degneracies
#      end # find tetrahedron simplex
#  end # find simplex

#  using GJK: findcombination
#  @testset "find combination" begin
#      @testset "find line combination" begin
#          spx = [1.0 1.0; -3.0 2.0]
#          vec = [1.0, 0.5]
#          @test all(findcombination(spx, vec) .≈ [0.3, 0.7])
#          vec = [1.0, 2.5]
#          @test all(findcombination(spx, vec) .≈ [0.0, 1.0])
#      end # find line combination
#
#      @testset "find triangle combination" begin
#          spx = [1.0 1.0 -1.0; 1.0 2.0 1.0; 2.0 2.0 1.0]
#          vec = [1.0, 4.0, 5.0]/3
#          @test all(findcombination(spx, vec) .≈ ones(3)/3.0)
#          vec = [-3.0, 0.0, 0.0]
#          @test all(findcombination(spx, vec) .≈ [0.0, 0.0, 1.0])
#          vec = [0.5, 0.0, 0.0]
#          @test all(findcombination(spx, vec) .≈ [0.4, 0.0, 0.6])
#          vec = [0.0, 2.0, 1.0]
#          @test all(findcombination(spx, vec) .≈ [0.0, 0.5, 0.5])
#      end # find triangle combination
#  end # find combination
