"""
    closest_points(p, q, dir; max_iter=100)

Compute the closest points between convex objects `p` and `q` if
they are not colliding. Provide an initial search direction
`dir` in the space of the problem.

Return result of type Tuple with StaticArrays of the two closest
points on each object.

# Examples
```julia-repl
julia> using StaticArrays
julia> p = SMatrix{2,3}([0.0 0.0 1.0; 0.0 2.0 1.0])
julia> q = SMatrix{2,3}([2.0 2.0 3.0; 0.0 2.0 1.0])
julia> dir = SVector{2}(1.0, 0.0)
julia> closest_points(p, q, dir)
([1.0, 1.0], [2.0, 1.0])
```
"""
function closest_points(p::Any, q::Any, init_dir::SVector; max_iter=100)
    collision, psimplex, qsimplex, dir, sz = gjk(p, q, init_dir, max_iter, minimum_distance_cond)

    if !collision
        pclose, qclose = nearestfromsimplex(psimplex, qsimplex, -dir, sz)
        return pclose, qclose
    end
end

"""
    minimum_distance(p, q, dir; max_iter=100)

Compute the minimum seperating distance between convex objects `p` and `q`.
Provide an initial search direction `dir` in the space of the problem.

# Examples
```julia-repl
julia> using StaticArrays
julia> p = SMatrix{2,3}([0.0 0.0 1.0; 0.0 2.0 1.0])
julia> q = SMatrix{2,3}([2.0 2.0 3.0; 0.0 2.0 1.0])
julia> dir = SVector{2}(1.0, 0.0)
julia> minimum_distance(p, q, dir)
1.0
```
"""
function minimum_distance(p::Any, q::Any, init_dir::SVector; max_iter=100)
    collision, psimplex, qsimplex, dir, sz = gjk(p, q, init_dir, max_iter, minimum_distance_cond)
    return collision ? 0.0 : norm(dir)
end

"""
    tolerance_verifcation(p, q, dir, τ; max_iter=100)

Compute if the convex objects `p` and `q` are at least
`τ` tolerance apart. Provide an initial search direction
`dir` in the space of the problem.

# Examples
```julia-repl
julia> using StaticArrays
julia> p = SMatrix{2,3}([0.0 0.0 1.0; 0.0 2.0 1.0])
julia> q = SMatrix{2,3}([2.0 2.0 3.0; 0.0 2.0 1.0])
julia> dir = SVector{2}(1.0, 0.0)
julia> tolerance_verification(p, q, dir, 0.25)
true
```
"""
function tolerance_verification(p::Any, q::Any, init_dir::SVector, τ::Real; max_iter=100)
    collision, = gjk(p, q, init_dir, max_iter, tolerance_verification_cond, τ)
    return !collision
end

"""
    collision_detection(p, q, dir; max_iter=100)

Compute if the convex objects `p` and `q` are colliding.
Provide an initial search direction `dir` in the space
of the problem.

# Examples
```julia-repl
julia> using StaticArrays
julia> p = SMatrix{2,3}([0.0 0.0 1.0; 0.0 2.0 1.0])
julia> q = SMatrix{2,3}([2.0 2.0 3.0; 0.0 2.0 1.0])
julia> dir = SVector{2}(1.0, 0.0)
julia> collision_detection(p, q, dir)
false
```
"""
function collision_detection(p::Any, q::Any, init_dir::SVector{N, T}; max_iter=100) where {N, T}
    collision, = gjk(p, q, init_dir, max_iter, tolerance_verification_cond, eps(T))
    return collision
end

# See Julia issue #28720 for why we have `where {F<:Function}` signature
function gjk(p::Any, q::Any, init_dir::SVector, max_iter::Int, distance::F, params...) where {F<:Function}
    sz = 1
    collision = false
    ps = support(p, init_dir)
    qs = support(q, -init_dir)
    psimplex = insertcolumn(ps)
    qsimplex = insertcolumn(qs)
    dir = qs - ps

    if sum(abs2, dir) ≤ 2*eps(Float64)
        return true, psimplex, qsimplex, dir, sz
    end

    ps = support(p, dir)
    qs = support(q, -dir)
    while check_gjk(psimplex - qsimplex, ps - qs, dir, sz, max_iter, distance, params...)
        sz += 1
        psimplex = insertcolumn(psimplex, ps, sz)
        qsimplex = insertcolumn(qsimplex, qs, sz)
        psimplex, qsimplex, dir, collision, sz = findsimplex(psimplex, qsimplex, sz)
        if collision
            break
        else
            ps = support(p, dir)
            qs = support(q, -dir)
        end

        max_iter -= 1
    end

    return collision, psimplex, qsimplex, dir, sz
end

function check_gjk(simplex::SMatrix, s::SVector, dir::SVector, sz::Int, iter::Int, distance::F, params...) where {F<:Function}
    if iter ≤ 0
        error("GJK failed: maximum number of iterations reached")
    elseif distance(s, dir, params...) || check_degeneracy(simplex, s, sz)
        return false
    else
        return true
    end
end

function minimum_distance_cond(s::SVector, dir::SVector)
    s⋅(-dir) ≥ sum(abs2, dir)
end

function tolerance_verification_cond(s::SVector, dir::SVector, tv::Real)
    s⋅(-dir) ≥ tv^2 || s⋅(-dir) ≥ sum(abs2, dir)
end

@generated function check_degeneracy(simplex::SMatrix{N, M, T}, s::SVector{N, T}, sz::Vararg{Int, 1}) where {N, M, T}
    exprs = :(false)
    for i = 1:M
        exprs = :($i > sz[1] ? $exprs : ($exprs || sum(abs2, simplex[:, $i] - s) ≤ 2*eps($T)))
    end
    return exprs
end
