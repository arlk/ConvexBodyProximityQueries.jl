struct ConvexPolygon{D, N, T}
    pts::SVector{N, SVector{D, T}}
end

function ConvexPolygon(pts::Array{Array{T, 1}}, ::Val{D}, ::Val{N}) where {D, N, T}
    ConvexPolygon(SVector{N, SVector{D, T}}(pts...))
end

function support(poly::ConvexPolygon{D, N, T}, dir::SVector{D, T}) where {D, N, T}
    @inbounds poly.pts[argmax(Ref(dir').*poly.pts)]
end

function point(pt::SVector{D, T}) where {D, T}
    ConvexPolygon(SVector{1, SVector{D, T}}(pt))
end

macro point(ex)
    :(point(@SVector $ex))
end

function line(ptA::SVector{D, T}, ptB::SVector{D, T}) where {D, T}
    ConvexPolygon(SVector{2, SVector{D, T}}(ptA, ptB))
end

macro line(ex)
    :(line(@SVector($(ex.args[1])), @SVector($(ex.args[2]))))
end

function rect(center::SVector{2, T}, widths::SVector{2, T}) where {T}
    pt1 = center + SVector{2}( 0.5,  0.5).*widths
    pt2 = center + SVector{2}(-0.5,  0.5).*widths
    pt3 = center + SVector{2}(-0.5, -0.5).*widths
    pt4 = center + SVector{2}( 0.5, -0.5).*widths
    ConvexPolygon(SVector{4, SVector{2, T}}(pt1, pt2, pt3, pt4))
end

function rect(center::SVector{3, T}, widths::SVector{3, T}) where {T}
    pt1 = center + SVector{3}( 0.5,  0.5,  0.5).*widths
    pt2 = center + SVector{3}(-0.5,  0.5,  0.5).*widths
    pt3 = center + SVector{3}(-0.5, -0.5,  0.5).*widths
    pt4 = center + SVector{3}( 0.5, -0.5,  0.5).*widths
    pt5 = center + SVector{3}( 0.5,  0.5, -0.5).*widths
    pt6 = center + SVector{3}(-0.5,  0.5, -0.5).*widths
    pt7 = center + SVector{3}(-0.5, -0.5, -0.5).*widths
    pt8 = center + SVector{3}( 0.5, -0.5, -0.5).*widths
    ConvexPolygon(SVector{8, SVector{3, T}}(pt1, pt2, pt3, pt4,
                                            pt5, pt6, pt7, pt8))
end

macro rect(ex)
    :(rect(@SVector($(ex.args[1])), @SVector($(ex.args[2]))))
end

function square(center::SVector{D, T}, width::Real) where {D, T}
    widths = width*@SVector(ones(D))
    rect(center, widths)
end

macro square(ex)
    :(square(@SVector($(ex.args[1])), $(ex.args[2])))
end

function randpoly(center=[0., 0.], θ=0.0; n::Int=10, scale=1.0)
    # from http://cglab.ca/~sander/misc/ConvexGeneration/convex.html
    xs = sort(rand(n))
    ys = sort(rand(n))
    minx = xs[1]; maxx = xs[end]
    miny = ys[1]; maxy = ys[end]
    xv = similar(xs)
    yv = similar(ys)
    ltop = minx; lbot = minx;
    llef = miny; lrig = miny;
    for i = 2:n-1
        if rand(Bool)
            xv[i-1] = xs[i] - ltop
            ltop = xs[i]
        else
            xv[i-1] = lbot - xs[i]
            lbot = xs[i]
        end
        if rand(Bool)
            yv[i-1] = ys[i] - llef
            llef = ys[i]
        else
            yv[i-1] = lrig - ys[i]
            lrig = ys[i]
        end
    end
    xv[n-1] = maxx - ltop
    xv[n] = lbot - maxx
    yv[n-1] = maxy - llef
    yv[n] = lrig - maxy
    shuffle!(xv)
    shuffle!(yv)
    vec = [[xv[i], yv[i]] for i = 1:n]
    sort!(vec, by = x -> atan(x[2], x[1]))
    v = [0.0, 0.0]
    minv = [0.0, 0.0]
    pts = similar(vec)
    for i = 1:n
        pts[i] = v
        v += vec[i]
        minv = min.(minv, v)
    end
    shift = [minx-0.5, miny-0.5] - minv
    pts .+= Ref(shift)
    rotate = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    pts = Ref(rotate).*pts
    pts *= scale
    pts .+= Ref(center)
    ConvexPolygon(pts, Val(2), Val(n))
end
