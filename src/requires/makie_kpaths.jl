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
                    hideaxis::Bool = true,
                    axis = NamedTuple(), figure = NamedTuple(), kws...) where D
    f = Makie.Figure(; figure...)
    ax = _default_bare_axis!(f, Val(D); hideaxis, axis)

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
