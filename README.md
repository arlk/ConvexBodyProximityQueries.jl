# ConvexBodyProximityQueries

[![Build Status](https://travis-ci.org/arlk/ConvexBodyProximityQueries.jl.svg?branch=master)](https://travis-ci.org/arlk/ConvexBodyProximityQueries.jl) [![codecov](https://codecov.io/gh/arlk/ConvexBodyProximityQueries.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/arlk/ConvexBodyProximityQueries.jl)

![](https://github.com/arlk/ConvexBodyProximityQueries.jl/raw/master/readme/logo.gif)

ConvexBodyProximityQueries.jl implements the Gilber-Johnson-Keerthi (GJK) Algorithm from their seminal paper on fast collision detection. The following query types are available for two convex objects:
 - Closest Points
 - Minimum Distance
 - Tolerance Verification
 - Collision Detection

## Usage

The package (by default) allows you to work with polytopes defined as an array of vertices, for example:
```julia
julia> using StaticArrays
julia> polyA = @SMatrix rand(2, 8)
2×8 SArray{Tuple{2,8},Float64,2,16}:
 0.732619   0.327745   0.0390878  0.477455  0.627223  0.502666  0.0529193  0.0523722
 0.0513408  0.0634308  0.892253   0.88009   0.100901  0.564782  0.789238   0.307813

julia> polyB = @SMatrix(rand(2, 8)) + 1.5
2×8 SArray{Tuple{2,8},Float64,2,16}:
 2.18993  1.75404  1.51373  1.60674  1.67257  2.14208  1.97779  2.24657
 2.32708  1.92212  2.32769  1.69457  1.85003  1.57441  1.93884  2.45361

julia> dir = @SVector(rand(2)) - 0.5
2-element SArray{Tuple{2},Float64,1,2}:
-0.4673435693835293
 0.4242237214159814
```

All the proximity queries can be performed simply by providing the polytope information and an initial searchdirection. In addition, `tolerance_verfication` requires an argument specifying the minimum tolerance of speration between two objects. :
```julia
julia> using BenchmarkTools
julia> @btime closest_points($polyA, $polyB, $dir)
  172.901 ns (0 allocations: 0 bytes)
([0.477455, 0.88009], [1.60674, 1.69457])

julia> @btime minimum_distance($polyA, $polyB, $dir)
  165.554 ns (0 allocations: 0 bytes)
1.3923553706117722

julia> @btime tolerance_verification($polyA, $polyB, $dir, $1.0)
  99.324 ns (0 allocations: 0 bytes)
true

julia> @btime collision_detection($polyA, $polyB, $dir)
  96.386 ns (0 allocations: 0 bytes)
false
```

If you want to use your custom convex objects, you can do so by extending the support function as:
```julia
import ConvexBodyProximityQueries.support
function ConvexBodyProximityQueries.support(obj::MyFancyShape, dir::SVector{N}) where {N}
  # do something
  return supporting_point::SVector{N}
end
```
_Note:_ This is how I intended the package to be used, the vanilla `support` function is quite naive and only works for a StaticArray of vertices. Here are some examples for some geometries found in [GeometryTypes.jl](https://github.com/JuliaGeometry/GeometryTypes.jl):
```julia
import ConvexBodyProximityQueries.support
using GeometryTypes: HyperSphere, HyperRectangle, HyperCube

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

@generated function ConvexBodyProximityQueries.support(cube::HyperCube{N, T}, dir::AbstractVector) where {N, T}
    exprs = Array{Expr}(undef, (N,))
    for i = 1:N
        exprs[i] = :(cube.width*(dir[$i] ≥ 0.0 ? 1.0 : -1.0)/2.0 + cube.origin[$i])
    end

    return quote
        Base.@_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return SVector{N, T}(elements)
    end
end
```

### Speed

As the core routines use StaticArrays, they are very well optimized and run quickly with no memory allocations. However, it is upto to the user to provide efficient code for the `support` and a good `init_dir` vector to squeeze the best performance from the functions.


## Examples

Minimum distance computation in 2D:

![](https://github.com/arlk/ConvexBodyProximityQueries.jl/raw/master/readme/collision2d.png)

Minimum distance computation in 3D:

![](https://github.com/arlk/ConvexBodyProximityQueries.jl/raw/master/readme/collision3d.png)

## Related Packages

[EnhancedGJK.jl](https://github.com/rdeits/EnhancedGJK.jl)

## References

Gilbert, E. G., D. W. Johnson, and S. S. Keerthi. “A Fast Procedure for Computing the Distance between Complex Objects in Three-Dimensional Space.” IEEE Journal on Robotics and Automation 4, no. 2 (April 1988): 193–203. https://doi.org/10.1109/56.2083.
