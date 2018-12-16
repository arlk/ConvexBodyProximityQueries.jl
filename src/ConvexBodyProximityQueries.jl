module ConvexBodyProximityQueries

export closest_points
export minimum_distance
export tolerance_verification
export collision_detection

using LinearAlgebra
using StaticArrays

include("support.jl")
include("simplex.jl")
include("proximity.jl")

end # module
