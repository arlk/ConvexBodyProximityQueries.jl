using ConvexBodyProximityQueries
using Test
using LinearAlgebra
using StaticArrays

@testset "ConvexBodyProximityQueries" begin
    include("support.jl")
    include("simplex.jl")
    include("proximity.jl")
end
