using LinearAlgebra
using .StaticArrays: normalize

# ---------------------------------------------------------------------------------------- #
# `KPath`
# NB: This file expects to be loaded *after* /src/requires/makie_wignerseitz.jl

@recipe(KPathPlot) do scene
    Attributes(
        linecolor = KPATH_COL[],
        linewidth = 3,
        linekws = (;),
        markercolor = KPATH_COL[],
        markersize = 20,
        markerkws = (;),
        textcolor = :black,
        textkws = (;)
    )
end

function Makie.plot!(kpp::KPathPlot{Tuple{KPath{D}}}) where D
    okp = kpp[1] # NB: an `Observable{KPath{D}}` not a `KPath{D}`, so needs [] to access

    segments  = Makie.Observable(Makie.Point{D,Float64}[])
    pts       = Makie.Observable(Makie.Point{D,Float64}[])
    labels    = Makie.Observable(String[])
    labelspos = Makie.Observable(Makie.Point{D,Float64}[])

    function update_plot!(kp)
        empty!(segments[])
        empty!(pts[])
        empty!(labels[])
        empty!(labelspos[])

        setting(kp) !== CARTESIAN && (kp = cartesianize(kp))
        for (i,path) in enumerate(paths(kp))
            kvs = [points(kp)[lab] for lab in path]
            append!(segments[], kvs)
            i ≠ length(paths(kp)) && push!(segments[], Point{D,Float64}(NaN)) # separator
        end
        pts[] = collect(values(points(kp))) # update and trigger observables
        labels[] = [string(klab) for klab in keys(points(kp))]
        labelspos[] = eval_label_positions(pts[], basis(kp))
    end
    Makie.Observables.onany(update_plot!, okp) # update plot if `oc` observable changes
    update_plot!(okp[]) # populate `segments`, `pts`, `labels`, `labelspos` for initial call

    Makie.lines!(kpp, segments;
                 color = kpp.linecolor, linewidth = kpp.linewidth, kpp.linekws...)
    Makie.scatter!(kpp, pts; 
                   color=kpp.markercolor, markersize=kpp.markersize, marker=:circle,
                   kpp.markerkws...)
    Makie.text!(kpp, labelspos[]; text=labels[],
                align = (0.5,0.5), fontsize = 22,
                color = kpp.textcolor, strokewidth = D==3 ? 2 : 0, strokecolor = :white,
                kpp.textkws...)

    return kpp
end

function Makie.plot!(ax::Makie.Block, kp::Union{Observable{<:KPath}, <:KPath}; kws...)
    kpathplot!(ax, kp; kws...)
end
function Makie.plot!(kp::Union{Observable{<:KPath}, <:KPath}; kws...)
    Makie.plot!(Makie.current_axis(), kp; kws...)
end

function Makie.plot(kp::Union{Observable{KPath{D}}, KPath{D}};
                    axis = NamedTuple(), figure = NamedTuple(), kws...) where D
    f = Makie.Figure(; figure...)
    ax = _default_bare_axis!(f, Val(D); axis)

    p = Makie.plot!(ax, kp; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end

# ---------------------------------------------------------------------------------------- #
# Utilities

function eval_label_positions(kvs, basis)
    maxnorm = maximum(norm, basis) # determine offset based on basis norms
    shiftlen = maxnorm * .065
    map(kvs) do kv
        if iszero(kv) # pull Γ _away_ from the mean of all other points
            kv + normalize(sum(kvs)) * (-shiftlen)
        else
            kv + normalize(kv) * shiftlen
        end
    end
end

# ---------------------------------------------------------------------------------------- #
# Hover interaction 
# (borrowed from https://github.com/JuliaPlots/GraphMakie.jl/blob/master/src/interaction.jl)

"""
    mutable struct HoverHandler{P, F}

Object to handle hovers on `plot::P`. Calls `fun` on hover.
"""
mutable struct HoverHandler{P,F}
    idx::Union{Nothing,Int}
    plot::Union{Nothing,P}
    fun::F
end

function process_interaction(handler::HoverHandler, event::MouseEvent, axis)
    if event.type === MouseEventTypes.over
        (element, idx) = convert_selection(mouse_selection(axis.scene)...)
        if element == handler.plot
            if handler.idx === nothing
                handler.idx = idx
                ret = handler.fun(true, handler.idx, event, axis)
                return ret isa Bool ? ret : false
            end
        else
            if handler.idx !== nothing
                ret = handler.fun(false, handler.idx, event, axis)
                handler.idx = nothing
                return ret isa Bool ? ret : false
            end
        end
    end
    return false
end

"""
    NodeHoverHandler(fun)

Initializes `HoverHandler` for markers. Calls function `fun(hoverstate, idx, event, axis)`
with `hoverstate=true` on hover and `false` at the end of hover. `idx` is the marker index.
Enabled by `register_interaction!(ax, :nodehover, NodeHoverHandler(action))`.
"""
NodeHoverHandler(fun::F) where {F} = HoverHandler{Scatter,F}(nothing, nothing, fun)

"""
    NodeHoverHeighlight(p::KPathPlot, factor=2)

Magnifies the `markersize` of node under cursor by `factor`.
Enabled by `register_interaction!(ax, :nodehover, NodeHoverHighlight(p))`
"""
function NodeHoverHighlight(p::KPathPlot, factor=2)
    if !(p.markersize[] isa Vector{<:Real})
        error("`markersize` must be a `Vector{<:Real}` for this interactivity to work")
    end
    action = (state, idx, _, _) -> begin
        old = p.markersize[][idx]
        p.markersize[][idx] = state ? old * factor : old / factor
        p.markersize[] = p.markersize[] #trigger observable
    end
    return NodeHoverHandler(action)
end