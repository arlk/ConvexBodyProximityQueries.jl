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

end # module
