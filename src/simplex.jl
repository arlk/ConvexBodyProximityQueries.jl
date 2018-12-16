"""
    proj(u, v)

Compute the vector projection of vector `v` onto vector `u`.
"""
function proj(u::AbstractVector, v::AbstractVector)
    (u ⋅ v)/(u ⋅ u)*u
end

@generated function insertcolumn(s::SVector{N, T}) where {N, T}
    S = (N, N+1)
    exprs = Array{Expr}(undef, S)
    itr = [1:n for n = S]
    for k = Base.product(itr...)
        exprs[k...] = :($(k[2]) > 1 ? zero($T) : s[$(k[1])])
    end

    return quote
        Base.@_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return SMatrix{N, N+1, T}(elements)
    end
end

@generated function insertcolumn(simplex::SMatrix{N, M, T}, s::SVector{N, T}, idx::Vararg{Int, 1}) where {N, M, T}
    S = (N, M)
    exprs = Array{Expr}(undef, S)
    itr = [1:n for n = S]
    for k = Base.product(itr...)
        exprs[k...] = :($(k[2]) > idx[1] ? zero($T) : ($(k[2]) == idx[1] ? s[$(k[1])] : simplex[$(k...)]))
    end

    return quote
        Base.@_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return SMatrix{N, M, T}(elements)
    end
end

@generated function pickcolumns(simplex::SMatrix{N, M, T}, idx::Vararg{Int, J}) where {N, M, T, J}
    S = (N, M)
    exprs = Array{Expr}(undef, S)
    itr = [1:n for n = S]
    for k = Base.product(itr...)
        exprs[k...] = :($(k[2]) > $J ? zero($T) : simplex[$(k[1]), idx[$(k[2])]])
    end

    return quote
        Base.@_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return SMatrix{N, M, T}(elements)
    end
end

@generated function reducecolumns(simplex::SMatrix{N, M, T}, ::Vararg{Any, J}) where {N, M, T, J}
    S = (N, J)
    exprs = Array{Expr}(undef, S)
    itr = [1:n for n = S]
    for k = Base.product(itr...)
        exprs[k...] = :(simplex[$(k...)])
    end

    return quote
        Base.@_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return SMatrix{N, J, T}(elements)
    end
end

"""
    findsimplex(simplex)

Compute the new simplices from a pair of given simplices.
Return the new search direction. Return a collision flag
if the origin was enclosed by the Minkowski simplex.
"""
function findsimplex(psimplex::SMatrix, qsimplex::SMatrix, sz::Int)
    if sz == 2
        findline(psimplex, qsimplex)
    elseif sz == 3
        findtriangle(psimplex, qsimplex)
    elseif sz == 4
        findtetrahedron(psimplex, qsimplex)
    end
end

function findline(psimplex::SMatrix{N}, qsimplex::SMatrix{N}) where {N}
    simplex = psimplex - qsimplex
    AB = simplex[:, 1] - simplex[:, 2]
    AO = -simplex[:, 2]
    if AB ⋅ AO > 0
        dir = AO - proj(AB, AO)
        collision = sum(abs2, dir) ≤ 2*eps(Float64)
        return psimplex, qsimplex, dir, collision, 2
    else
        dir = AO
        psimplex = pickcolumns(psimplex, 2)
        qsimplex = pickcolumns(qsimplex, 2)
        collision = sum(abs2, dir) ≤ 2*eps(Float64)
        return psimplex, qsimplex, dir, collision, 1
    end
end

function findtriangle(psimplex::SMatrix{N}, qsimplex::SMatrix{N}) where {N}
    simplex = psimplex - qsimplex
    AB = simplex[:, 2] - simplex[:, 3]
    AC = simplex[:, 1] - simplex[:, 3]
    BC = simplex[:, 1] - simplex[:, 2]
    AO = -simplex[:, 3]
    if (AC ⋅ AB * BC - AC ⋅ BC * AB) ⋅ AO > 0
        if AC ⋅ AO > 0
            psimplex = pickcolumns(psimplex, 1, 3)
            qsimplex = pickcolumns(qsimplex, 1, 3)
            dir = AO - proj(AC, AO)
            collision = sum(abs2, dir) ≤ 2*eps(Float64)
            return psimplex, qsimplex, dir, collision, 2
        else
            psimplex = pickcolumns(psimplex, 2, 3)
            qsimplex = pickcolumns(qsimplex, 2, 3)
            return findline(psimplex, qsimplex)
        end
    elseif (AB ⋅ BC * AB - AB ⋅ AB * BC) ⋅ AO > 0
        psimplex = pickcolumns(psimplex, 2, 3)
        qsimplex = pickcolumns(qsimplex, 2, 3)
        return findline(psimplex, qsimplex)
    else
        if sum(abs2, AO) ≤ 2*eps(Float64)
            return psimplex, qsimplex, AO, true, 3
        elseif sum(abs2, AC - proj(AB, AC)) ≤ 2*eps(Float64)
            psimplex = pickcolumns(psimplex, 2, 3)
            qsimplex = pickcolumns(qsimplex, 2, 3)
            return findline(psimplex, qsimplex)
        elseif N == 2
            return psimplex, qsimplex, AO, true, 3
        else
            ABC = AB × BC
            dir = proj(ABC, AO)
            if ABC ⋅ AO > 0
                psimplex = pickcolumns(psimplex, 2, 1, 3)
                qsimplex = pickcolumns(qsimplex, 2, 1, 3)
                collision = sum(abs2, dir) ≤ 2*eps(Float64)
                return psimplex, qsimplex, dir, collision, 3
            else
                collision = sum(abs2, dir) ≤ 2*eps(Float64)
                return psimplex, qsimplex, dir, collision, 3
            end
        end
    end
end

function findtetrahedron(psimplex::SMatrix{N}, qsimplex::SMatrix{N}) where {N}
    simplex = psimplex - qsimplex
    AB = simplex[:, 3] - simplex[:, 4]
    AC = simplex[:, 2] - simplex[:, 4]
    AD = simplex[:, 1] - simplex[:, 4]
    BC = simplex[:, 2] - simplex[:, 3]
    BD = simplex[:, 1] - simplex[:, 3]
    AO = -simplex[:, 4]
    if (AB × AC) ⋅ AO < 0
        psimplex = pickcolumns(psimplex, 2, 3, 4)
        qsimplex = pickcolumns(qsimplex, 2, 3, 4)
        return findtriangle(psimplex, qsimplex)
    elseif (AD × AB) ⋅ AO < 0
        psimplex = pickcolumns(psimplex, 2, 3, 4)
        qsimplex = pickcolumns(qsimplex, 2, 3, 4)
        return findtriangle(psimplex, qsimplex)
    elseif (AC × AD) ⋅ AO < 0
        psimplex = pickcolumns(psimplex, 2, 3, 4)
        qsimplex = pickcolumns(qsimplex, 2, 3, 4)
        return findtriangle(psimplex, qsimplex)
    else
        return psimplex, qsimplex, AO, true, 4
    end
end

"""
    nearestfromsimplex(psimplex, qsimplex, dir2origin)

Compute the nearest points between two simplexes given the
direction to in origin in the Minkowski difference space
"""
function nearestfromsimplex(psimplex::SMatrix{N, M}, qsimplex::SMatrix{N, M}, dir2origin::SVector{N}, sz::Int) where {N, M}
    if sz == 1
        return psimplex[:, 1], qsimplex[:, 1]
    elseif sz == 2
        λ = linecombination(psimplex - qsimplex, dir2origin)
        pclose = reducecolumns(psimplex, λ...)*SVector{2}(λ)
        qclose = reducecolumns(qsimplex, λ...)*SVector{2}(λ)
        return pclose, qclose
    elseif sz == 3
        λ = trianglecombination(psimplex - qsimplex, dir2origin)
        pclose = reducecolumns(psimplex, λ...)*SVector{3}(λ)
        qclose = reducecolumns(qsimplex, λ...)*SVector{3}(λ)
        return pclose, qclose
    end
end

function linecombination(simplex::SMatrix{N}, vec::SVector{N}) where {N}
    AV = vec - simplex[:, 2]
    AB = simplex[:, 1] - simplex[:, 2]
    λ = (AB ⋅ AV)/(AB ⋅ AB)
    if λ < 0
        return 0.0, 1.0 - λ
    else
        return λ, 1.0 - λ
    end
end

function trianglecombination(simplex::SMatrix{N}, vec::SVector{N}) where {N}
    AO = -simplex[:, 3]
    AV = vec - simplex[:, 3]
    AB = simplex[:, 2] - simplex[:, 3]
    AC = simplex[:, 1] - simplex[:, 3]
    BC = simplex[:, 1] - simplex[:, 2]

    if (AC ⋅ AB * BC - AC ⋅ BC * AB) ⋅ AV > 0
        if AC ⋅ AV > 0
            dir = AO - proj(AC, AV)
            idx = SMatrix{3, 2}(1, 0, 0, 0, 0, 1)
            simplex = pickcolumns(simplex, 1, 3)
            λAC1, λAC2 = linecombination(simplex, -dir)
            return λAC1, 0.0, λAC2
        else
            dir = AO - proj(AB, AV)
            simplex = pickcolumns(simplex, 2, 3)
            λAB1, λAB2 = linecombination(simplex, -dir)
            return 0.0, λAB1, λAB2
        end
    elseif (AB ⋅ BC * AB - AB ⋅ AB * BC) ⋅ AV > 0
        dir = AO - proj(AB, AV)
        simplex = pickcolumns(simplex, 2, 3)
        λAB1, λAB2 = linecombination(simplex, -dir)
        return 0.0, λAB1, λAB2
    else
        sABC = [AC AB]
        λ = (sABC' * sABC) \ (sABC' * AV)
        return λ[1], λ[2], 1-sum(λ)
    end
end
