using Luxor
using Colors
using QHull: chull
using ConvexBodyProximityQueries: closest_points
using StaticArrays

function randdisk(sz)
    r = rand(sz)
    θ = 2*π*rand(sz)
    x = sqrt.(r).*cos.(θ)
    y = sqrt.(r).*sin.(θ)
    return hcat(x, y)
end

function randomcvx(sz::Int)
    pts = randdisk(sz)
    ch = chull(pts)
    sz = length(ch.vertices)
    return randomcvx(ch.points[ch.vertices, :], Val(sz))
end

function randomcvx(hull::Array, ::Val{N}) where {N}
    hulln = vcat(hull', ones(1, N))
    return SMatrix{3, N}(hulln)
end

logo = Movie(400, 283, "test")

function backdrop(scene, framenumber)
    background("white")
end

transform2(θ,c) = SMatrix{2,3}([cos(θ) -sin(θ) c[1]; sin(θ) cos(θ) c[2]])
polyA = randomcvx(35)
polyB = randomcvx(35)
polyC = randomcvx(35)

function drawpoly(polygon, θ, center)
    p = transform2(θ, center)*polygon
    plux = [Point(p[1:2,i]...) for i in 1:size(p, 2)]
    poly(plux, :fill)
    return p
end

function frame(scene, framenumber)
    global polyA, polyB, polyC
    centerA = SMatrix{2, 1}( 0.,-√2.)
    centerB = SMatrix{2, 1}( 2., √2.)
    centerC = SMatrix{2, 1}(-2., √2.)
    scale(50)

    # draw polygons
    sethue(Colors.RGB(96/255,173/255,81/255))
    pA = drawpoly(polyA, deg2rad(-framenumber/2), centerA)
    sethue(Colors.RGB(170/255,121/255,193/255))
    pB = drawpoly(polyB, deg2rad(framenumber), centerB)
    sethue(Colors.RGB(213/255,99/255,92/255))
    pC = drawpoly(polyC, deg2rad(framenumber), centerC)

    # get closest points between objects
    cp = map(p -> closest_points(p[1], p[2], SVector{2}(1., 1.)), [[pA, pB], [pB, pC], [pC, pA]])
    cp = map(p -> map(c -> Point(c.data...), p), cp)

    # draw line segments
    sethue(Colors.RGB(102/255,130/225,223/255))
    map(p -> line(p[1], p[2], :stroke), cp)

    # draw points
    sethue("blue")
    map(p -> map(c -> circle(c, 0.04, :fill), p), cp)
end

animate(logo, [
    Scene(logo, backdrop, 0:719),
    Scene(logo, frame, 0:719)
    ],
    creategif=true,
    framerate=30,
    pathname="logo.gif")
