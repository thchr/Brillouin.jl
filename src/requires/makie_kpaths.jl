# ---------------------------------------------------------------------------------------- #
# `KPath`
# NB: This file expects to be loaded *after* /src/requires/makie_wignerseitz.jl

@recipe(KPathPlot) do scene
    Attributes(
        linecolor = KPATH_COL[],
        linewidth = 3,
        markercolor = KPATH_COL[],
        markersize = 20,
        marker = :circle
    )
end

function Makie.plot!(kpp::KPathPlot{Tuple{KPath{D}}}) where D
    okp = kpp[1] # NB: an `Observable{KPath{D}}` not a `KPath{D}`, so needs [] to access

    segments = Makie.Observable(Makie.Point{D,Float64}[])
    pts      = Makie.Observable(Makie.Point{D,Float64}[])
    function update_plot!(kp)
        empty!(segments[])
        setting(kp) !== CARTESIAN && (kp = cartesianize(kp))
        for (i,path) in enumerate(paths(kp))
            kvs = [points(kp)[lab] for lab in path]
            append!(segments[], kvs)
            i ≠ length(paths(kp)) && push!(segments[], Point{D,Float64}(NaN)) # separator
        end

        for (i,(lab, kv)) in enumerate(points(kp))
            push!(pts[], kv)
            # TODO: k-labels on hover?
        end
        pts[]=pts[] # trigger observables
    end
    Makie.Observables.onany(update_plot!, okp) # update plot if `oc` observable changes
    update_plot!(okp[]) # populate `segments` and `pts` for initial call

    Makie.lines!(kpp, segments; color = kpp.linecolor, linewidth = kpp.linewidth)
    Makie.scatter!(kpp, pts; 
                   color=kpp.markercolor, markersize=kpp.markersize, marker=kpp.marker)
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
    f[1,1] = ax = _Axis_by_dim(D)(f; aspect=_aspect_by_dim(D), axis...)
    Makie.hidedecorations!(ax); Makie.hidespines!(ax)

    p = Makie.plot!(ax, kp; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end