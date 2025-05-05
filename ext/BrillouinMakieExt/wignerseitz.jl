# ---------------------------------------------------------------------------------------- #
# `Cell`

@recipe(CellPlot, c) do scene
    Attributes(
        color = BZ_COL[],
        linewidth = 2,
        axis = NamedTuple(),
    )
end

function Makie.plot!(cp::CellPlot{Tuple{Cell{D}}}) where D
    oc = cp[1] # NB: this is an `Observable{Cell{D}}` not a `Cell{D}`, so needs [] to access
    segments = Makie.Observable(Makie.Point{D,Float64}[])
    function update_plot!(c)
        empty!(segments[])
        setting(c) !== CARTESIAN && (c = cartesianize(c))
        for (i, poly) in enumerate(c)
            append!(segments[], poly)
            push!(segments[], poly[1]) # end point
            i ≠ length(c) && push!(segments[], Point{D,Float64}(NaN)) # separator
        end
        segments[]=segments[] # trigger observable
    end
    Makie.Observables.onany(update_plot!, oc) # update plot if `oc` observable changes
    update_plot!(oc[]) # populate `segments` for initial call

    Makie.lines!(cp, segments; color = cp.color, linewidth = cp.linewidth)

    return cp
end

# ---------------------------------------------------------------------------------------- #

Makie.plottype(::Cell) = CellPlot # alias `cellplot` to `plot`

