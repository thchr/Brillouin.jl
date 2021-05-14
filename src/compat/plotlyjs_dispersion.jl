using .PlotlyJS
import .PlotlyJS: plot
                                          # from https://flatuicolors.com/palette/gb
const BAND_COL  = Ref("rgb(39,60,117)")   # "mazarine blue"
const KLINE_COL = Ref("rgb(220,221,225)") # "hint of pensive"

const DEFAULT_PLOTLY_LAYOUT_DISPERSION = Layout(
    showlegend=false,
    xaxis=attr(zeroline=false,
            showgrid=true, showbackground=false,
            gridcolor = KLINE_COL[], ticks="outside",
            ),
    yaxis=attr(zeroline=false,
            showgrid=false, showbackground=false,
            ticks="outside",
            ),
    hovermode = "closest",
    autosize = true,
    width = 480, height = 480,
    #margin=attr(l=50, r=5, b=15, t=10),

    plot_bgcolor="rgba(255, 255, 255, 1)", paper_bgcolor="rgba(255, 255, 255, 1)",
    )

"""
    plot(kpi::KPathInterpolant, bands, [layout]; kwargs...)

Plot a dispersion diagram for provided `bands` and **k**-path interpolant `kpi`.

`bands` must be an iterable of iterables of `<:Real`s (e.g., a `Vector{Vector{Float64}}`),
with the first iteration running over distinct energy bands, and the second running
over distinct **k**-points in `kpi`.
Note that the length of each iterant of `bands` must equal `length(kpi)`.

Alternatively, `bands` can be an `AbstractMatrix{<:Real}`, with columns interpreted as
distinct energy bands and rows as distinct **k**-points.

## Keyword arguments `kwargs`
- `ylims`: y-axis limits (default: quasi-tight around `bands`'s range)
- `ylabel`: y-axis label (default: "Energy")
"""
function plot(kpi::KPathInterpolant, bands,#::AVec{<:AVec{<:Real}},
              layout::Layout = DEFAULT_PLOTLY_LAYOUT_DISPERSION;
              ylims = nothing, ylabel = "Energy")
    # check input
    N = length(kpi)
    if !all(band -> length(band) == N, bands)
        throw(DimensionMismatch("mismatched dimensions of `kpi` and `bands`"))
    end
    # copy layout since we may want to mutate the layout and don't want to corrupt it
    layout = deepcopy(layout)

    # set default y-limits in layout, if not already set
    haskey(layout, :yaxis) || (layout[:yaxis] = attr())
    if isnothing(ylims)
        if !haskey(layout[:yaxis], :range)
            ylims = default_dispersion_ylims(bands)
            layout[:yaxis][:range] = ylims
        else
            ylims = layout[:yaxis][:range] # grab what was already in `layout`
        end
    else
        # overwrite if ylims was provided, regardless of what it is in `layout`
        layout[:yaxis][:range] = ylims
    end
    layout[:yaxis][:title] = ylabel

    # prepare to plot band diagram
    local_xs = cumdists.(kpi.kpaths)
    spacing = sum(last, local_xs) / 30

    # plot bands and k-lines/labels
    tbands = Vector{GenericTrace{Dict{Symbol,Any}}}()
    xticks = Vector{Float64}(undef, sum(length, kpi.labels))
    xlabs  = Vector{Symbol}(undef, sum(length, kpi.labels))
    tboxs  = Vector{GenericTrace{Dict{Symbol,Any}}}()
    start_idx = 1
    offset_x  = 0.0
    i = 0
    for (local_x, labels) in zip(local_xs, kpi.labels)
        stop_idx = start_idx+length(local_x)-1
        x = local_x .+ offset_x
        # plot bands
        for band in bands
            push!(tbands, 
                    PlotlyJS.scatter(x=x, y=band[start_idx:stop_idx],
                        hoverinfo="y", mode="lines", line=attr(color=BAND_COL[], width=3)))
        end
        # define xticks
        for (idx, lab) in labels
            xticks[i+=1] = x[idx]
            xlabs[i]     = lab
        end
        # draw boxes
        xlims = extrema(x)
        push!(tboxs, 
                PlotlyJS.scatter(x=xlims[[1, 2, 2, 1, 1]], y=ylims[[1,1,2,2,1]],
                    hoverinfo="none", mode="lines", line=attr(color="black", width=1)))
        # prepare for next iteration
        offset_x = last(x) + spacing
        start_idx = stop_idx + 1
    end
    
    # set xticks and x-range
    haskey(layout, :xaxis) || (layout[:xaxis] = attr())
    layout[:xaxis][:range]  = [0.0, offset_x - spacing]
    layout[:xaxis][:tickvals] = xticks
    layout[:xaxis][:ticktext] = xlabs

    return plot(vcat(tbands,tboxs), layout)
end
# `bands` can also be supplied as a matrix (w/ distinct bands in distinct columns)
function plot(kpi::KPathInterpolant, bands::AbstractMatrix{<:Real},
     layout::Layout = DEFAULT_PLOTLY_LAYOUT_DISPERSION; kwargs...)
     plot(kpi, eachcol(bands), layout; kwargs...)
end

function default_dispersion_ylims(bands)
    ylims = [mapfoldl(minimum, min, bands, init=Inf), 
             mapfoldl(maximum, max, bands, init=-Inf)]
    δ = (ylims[2]-ylims[1])/30
    if isapprox(ylims[1], 0, atol=1e-6)
        ylims[2] += δ
    else
        ylims .+= (-δ, δ)
    end
    return ylims
end