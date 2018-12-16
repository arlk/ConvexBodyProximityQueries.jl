"""
    support(cvxpoly, dir)

Compute the point of contact between a convex object and its
supporting hyperplane defined by the given normal direction as
a StaticArray.
"""
function support(vertices::SMatrix, dir::SVector)
    @inbounds vertices[:, argmax(vertices'*dir)]
end
