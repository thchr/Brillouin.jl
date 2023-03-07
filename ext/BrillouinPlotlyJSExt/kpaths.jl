# NB: This file should be loaded *after* /ext/BrillouinPlotlyExt/wignerseitz.jl

# ---------------------------------------------------------------------------------------- #

function plot(kp::KPath{3}, layout::Layout = Layout();
              config::PlotConfig = PlotConfig(responsive=true, displaylogo=false))
    setting(kp) !== CARTESIAN && (kp = cartesianize(kp))
    layout = merge(DEFAULT_PLOTLY_LAYOUT_3D, layout)

    tpaths = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(paths(kp)))
    for (i,path) in enumerate(paths(kp))
        kvs = [points(kp)[lab] for lab in path]
        tpaths[i] = PlotlyJS.scatter3d(
            x=getindex.(kvs, 1), y=getindex.(kvs, 2), z=getindex.(kvs, 3),
            hoverinfo = "none",
            mode="lines", line=attr(color=KPATH_COL[], width=8))
    end
    tpoints = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(points(kp)))
    for (i,(lab, kv)) in enumerate(points(kp))
        tpoints[i] = PlotlyJS.scatter3d(
            x=kv[1:1], y=kv[2:2], z=kv[3:3],
            mode="marker", hovertext=string(lab)*" = 2π"*string(round.(kv/(2π), digits=3)), 
            hoverinfo="text",
            marker=attr(color=KPATH_COL[], size=6, line=attr(color="white", width=1)))
    end
    return plot(vcat(tpoints, tpaths), layout; config=config)
end

# ---------------------------------------------------------------------------------------- #

function plot(kp::KPath{2}, layout::Layout = Layout();
              config::PlotConfig = PlotConfig(responsive=true, displaylogo=false))
    setting(kp) !== CARTESIAN && (kp = cartesianize(kp))
    layout = merge(DEFAULT_PLOTLY_LAYOUT_2D, layout)

    tpaths = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(paths(kp)))
    for (i,path) in enumerate(paths(kp))
        kvs = [points(kp)[lab] for lab in path]
        tpaths[i] = PlotlyJS.scatter(
            x=getindex.(kvs, 1), y=getindex.(kvs, 2),
            hoverinfo = "none",
            mode="lines", line=attr(color=KPATH_COL[], width=5))
    end
    tpoints = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(points(kp)))
    for (i,(lab, kv)) in enumerate(points(kp))
        tpoints[i] = PlotlyJS.scatter(
            x=kv[1:1], y=kv[2:2],
            mode="marker", hovertext=string(lab)*" = 2π"*string(round.(kv/(2π), digits=3)), 
            hoverinfo="text",
            marker=attr(color=KPATH_COL[], size=8, line=attr(color="white", width=1)))
    end
    return plot(vcat(tpoints, tpaths), layout; config=config)
end
# ---------------------------------------------------------------------------------------- #


function plot(c::Cell{D}, kp::KPath{D}, layout::Layout = Layout();
              config::PlotConfig = PlotConfig(responsive=true, displaylogo=false)
              ) where D

    D ∉ (2,3) && error("must be 2D or 3D Cell and KPath")
    # ::Cell
    Pᶜ  = plot(c, layout) # modifies `layout`; gets stored in resulting `Pᶜ`
    tsᶜ = Pᶜ.plot.data
    layout′ = Pᶜ.plot.layout # grab modified `layout`

    # get rid of "inside" axis lines; can be identified based on color
    filter!(tsᶜ) do t
        l = get(t, :line, nothing)
        l == nothing && return true
        c = get(l, :color, nothing)
        c == nothing && return true
        return !(c == BASIS_LIGHT_COL || c == AXIS_LIGHT_COL)
    end

    # ::KPath
    Pᵏᵖ  = plot(kp, layout′)
    tsᵏᵖ = Pᵏᵖ.plot.data

    # combine traces and plot
    ts = vcat(tsᶜ, tsᵏᵖ)
    return PlotlyJS.plot(ts, layout′; config=config)
end

# ---------------------------------------------------------------------------------------- #