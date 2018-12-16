"""
    support(cvxpoly, dir)

Compute the point of contact between a convex object and its
supporting hyperplane defined by the given normal direction as
a StaticArray.
"""
function support(vertices::SMatrix, dir::AbstractVector)
    @inbounds vertices[:, argmax(vertices'*dir)]
end

# TODO support for ordered vertices / simplices in a convex hull
#  struct CvxHull{A<:SMatrix}
#      pts::A
#  end
#  function support(vertices::CvxHull{SMatrix{2, M}}, dir::SVector{2, T}) where {M, T}
#      li = 1; ri = M
#      left = vertices[:, 1]'*dir
#      right = vertices[:, M]'*dir
#      cnt = 10
#      while ri - li > 2
#          mi = div(li + ri, 2)
#          mid = vertices[:, mi]'*dir
#          if left < right
#              left = mid
#              li = mi
#          else
#              right = mid
#              ri = mi
#          end
#          cnt -= 1
#      end
#      if left < right
#          return vertices[:, ri]
#      else
#          return vertices[:, li]
#      end
#  end
