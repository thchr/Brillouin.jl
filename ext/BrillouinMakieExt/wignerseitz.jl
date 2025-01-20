# ---------------------------------------------------------------------------------------- #
# `Cell`

@recipe(CellPlot) do scene
    Attributes(
        color = BZ_COL[],
        linewidth = 2
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

function Makie.plot!(ax::Makie.Block, c::Union{Observable{<:Cell}, <:Cell}; kws...)
    cellplot!(ax, c; kws...)
end
function Makie.plot!(c::Union{Observable{<:Cell}, <:Cell}; kws...)
    Makie.plot!(Makie.current_axis(), c; kws...)
end

@inline _Axis_by_dim(D)   = D == 3 ? Makie.Axis3   : Makie.Axis
@inline _aspect_by_dim(D) = D == 3 ? :data : Makie.DataAspect()
function Makie.plot(c::Union{Observable{Cell{D}}, Cell{D}};
                    hideaxis::Bool = true,
                    axis = NamedTuple(), figure = NamedTuple(), kws...) where D
    f = Makie.Figure(; figure...)
    ax = _default_bare_axis!(f, Val(D); hideaxis, axis)

    p = Makie.plot!(ax, c; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end

function _default_bare_axis!(f, ::Val{3}; hideaxis::Bool=true, axis=NamedTuple())
    f[1,1] = ax = Makie.Axis3(f;
        aspect=:data,
        viewmode=:fit,
        axis...)
    if hideaxis
        Makie.hidedecorations!(ax); ax.protrusions[] = 0 # cf. https://github.com/MakieOrg/Makie.jl/issues/2259
        Makie.hidespines!(ax)
    end
    return ax
end
function _default_bare_axis!(f, ::Val{2}; hideaxis::Bool=true, axis=NamedTuple())
    f[1,1] = ax = Makie.Axis(f;
        aspect=Makie.DataAspect(),
        axis...)
    hideaxis && (Makie.hidedecorations!(ax); Makie.hidespines!(ax))
    return ax
end