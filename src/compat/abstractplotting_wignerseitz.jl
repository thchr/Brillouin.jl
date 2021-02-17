using AbstractPlotting
import AbstractPlotting: plot, plot!

const DEFAULT_PYPLOT_OPTS = (color=:mediumblue, linewidth=2)
function plot!(s::Scene, c::Cell; kwargs...)
    cam3d!(s)
    for (i,poly) in enumerate(c)
        lines!(s, push!(getindex.(poly, 1), poly[1][1]),
                  push!(getindex.(poly, 2), poly[1][2]),
                  push!(getindex.(poly, 3), poly[1][3]);
                  DEFAULT_PYPLOT_OPTS..., # default options
                  kwargs...               # possible keyword overrides
              )
    end
    return s
end
plot(c::Cell; kwargs...) = plot!(Scene(), c; kwargs...)