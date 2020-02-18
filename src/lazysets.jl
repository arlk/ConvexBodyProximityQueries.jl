using .LazySets: LazySet, support_vector

function ConvexBodyProximityQueries.support(obj::LazySet{N}, dir::SVector{D, N}) where {D, N}
    return support_vector(dir, obj) |> SVector{D, N}
end
