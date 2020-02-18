module ConvexBodyProximityQueries

export closest_points
export minimum_distance
export tolerance_verification
export collision_detection

export ConvexPolygon
export point, line, rect, square, randpoly
export @point, @line, @rect, @square


using LinearAlgebra
using StaticArrays
using Random: shuffle!

include("support.jl")
include("obstacles.jl")
include("simplex.jl")
include("proximity.jl")
include("plot.jl")

using Requires

function __init__()
    @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" include("lazysets.jl")
end

end # module
