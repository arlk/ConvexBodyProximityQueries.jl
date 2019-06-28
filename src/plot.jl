using RecipesBase

@recipe function f(a::ConvexPolygon)
    grid --> false
    axis --> false
    ticks --> false
    legend --> false
    aspect_ratio --> :equal
    linewidth --> 1.0
    linecolor --> :black
    fillcolor --> :black
    fillalpha --> 0.5
    fillrange := 0.0
    leg := false
    @series begin
        vertices = vcat(a.pts, a.pts[1:1])
        Tuple.(vertices)
    end
end
