"""
    closest_points(p, q, dir; max_iter=100, atol=0.0)

Compute the closest points between convex objects `p` and `q` if
they are not colliding within some tolerance. Provide an initial
search direction `dir` in the space of the problem.

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
function closest_points(p::Any, q::Any, init_dir::SVector{D, T}; max_iter=100, atol::T=sqrt(eps(T))*oneunit(T)) where {D, T}
    collision, dir, psimplex, qsimplex, sz = gjk(p, q, init_dir, max_iter, atol, minimum_distance_cond)

    if !collision
        pclose, qclose = nearestfromsimplex(psimplex, qsimplex, -dir, sz)
        return pclose, qclose
    end
end

"""
    minimum_distance(p, q, dir; max_iter=100, atol=0.0)

Compute the minimum seperating distance between convex objects
`p` and `q` within some tolerance. Provide an initial search
direction `dir` in the space of the problem.

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
function minimum_distance(p::Any, q::Any, init_dir::SVector{D, T}; max_iter=100, atol::T=sqrt(eps(T))*oneunit(T)) where {D, T}
    collision, dir, psimplex, qsimplex, sz = gjk(p, q, init_dir, max_iter, atol, minimum_distance_cond)
    return collision ? zero(eltype(dir)) : norm(dir)
end

"""
    tolerance_verification(p, q, dir, τ; max_iter=100, atol=0.0)

Compute if the convex objects `p` and `q` are at least `τ`
tolerance apart within some tolerance. Provide an initial
search direction `dir` in the space of the problem.

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
function tolerance_verification(p::Any, q::Any, init_dir::SVector{D, T}, τ::T; max_iter=100, atol::T=sqrt(eps(T))*oneunit(T)) where {D, T}
    collision, = gjk(p, q, init_dir, max_iter, atol, tolerance_verification_cond, τ)
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
function collision_detection(p::Any, q::Any, init_dir::SVector{D, T}; max_iter=100) where {D, T}
    collision, = gjk(p, q, init_dir, max_iter, Inf*oneunit(T), collision_detection_cond)
    return collision
end

# See Julia issue #28720 for why we have `where {F<:Function}` signature
function gjk(p::Any, q::Any, init_dir::SVector{N, T}, max_iter::Int, atol::T, query::F, params...) where {F<:Function, N, T}
    sz = 1
    collision = false
    distance_condition = false
    ps = support(p, init_dir)
    qs = support(q, -init_dir)
    psimplex = insertcolumn(ps)
    qsimplex = insertcolumn(qs)
    dir = qs - ps
    ntol2 = eps(T)*oneunit(T)^2 # numerical tolerance

    if sum(abs2, dir) ≤ ntol2
        return true, dir, psimplex, qsimplex, sz
    end

    ps = support(p, dir)
    qs = support(q, -dir)
    s = ps - qs
    while !distance_condition && !collision && !check_degeneracy(psimplex-qsimplex, s, sz)
        sz += 1
        psimplex = insertcolumn(psimplex, ps, sz)
        qsimplex = insertcolumn(qsimplex, qs, sz)
        psimplex, qsimplex, dir, collision, sz = findsimplex(psimplex, qsimplex, sz)
        ps = support(p, dir)
        qs = support(q, -dir)
        s = ps - qs
        collision = query(dir, collision, params...)
        distance_condition = (s⋅(-dir)) ≥ ntol2 && (sum(abs2, dir) - s⋅(-dir)) ≤ atol^2

        max_iter -= 1
        if max_iter ≤ 0
            error("GJK failed: maximum number of iterations reached")
        end
    end

    return collision, dir, psimplex, qsimplex, sz
end

function minimum_distance_cond(dir::SVector, collision::Bool)
    collision
end

function collision_detection_cond(dir::SVector, collision::Bool)
    collision
end

function tolerance_verification_cond(dir::SVector, collision::Bool, tv)
    tv > norm(dir) ? true : collision
end

@generated function check_degeneracy(simplex::SMatrix{N, M, T}, s::StaticVector{N, T}, sz::Vararg{Int, 1}) where {N, M, T}
    exprs = :(false)
    for i = 1:M
        exprs = :($i > sz[1] ? $exprs : ($exprs || norm(simplex[:, $i] - s) ≤ eps($T)*oneunit($T)))
    end
    return exprs
end
