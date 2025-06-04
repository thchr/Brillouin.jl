# Implement `plot(kpi::KPathInterpolant, bands)` for Makie.jl: the key challenge to doing
# this being able to create multiple `Axis` objects (for disjoint parts of `kpi`), which
# is only possible via the SpecApi interface, which is slightly different but general than
# the `@recipe` interface

import Makie.SpecApi as S

## --------------------------------------------------------------------------------------- #

# keyword arguments must be explicitly marked in SpecApi :(
function Makie.used_attributes(::KPathInterpolant, ::AbstractVector{<:AbstractVector{<:Real}})
    (:color, :linewidth, :linestyle, :ylabel, :label, :annotations)
end

# the conversion method creates a grid of `Axis` objects with `Lines` plots inside
function Makie.convert_arguments(
    ::Type{<:AbstractPlot}, # makes `plot`, `lines`, `scatter` etc. dispatch here
    kpi::KPathInterpolant,
    bandsv::AbstractVector{<:AbstractVector{<:Real}}...;
    # keyword arguments below; must be explicit in SpecApi :(
    color     = BAND_COL[],
    linewidth = 3, # or, could have been set to `Makie.theme(:linewidth)` (=1.5, I think)
    linestyle = nothing,
    ylabel    = "Energy",
    label     = nothing, # legend labels
    annotations::Union{Dict{<:Union{String, Symbol}, 
                            <:Vector{<:Pair{<:Union{Int, UnitRange{Int}}, String}}},
                       Nothing} = nothing,
    )

    # check input
    Nk = length(kpi)
    for bands in bandsv
        if !all(band -> length(band) == Nk, bands)
            throw(DimensionMismatch("mismatched dimensions of `kpi` and `bandsv[i]`"))
        end
    end

    # preparation for `annotations`, if not nothing
    klab2idxs = if !isnothing(annotations)
        map(kpi.labels) do idx2klabs
            d = Dict{Symbol, Vector{Int}}() # NB: there can be multiple indices per k-point
            for (idx, klab) in idx2klabs
                push!(get!(()->Int[], d, klab), idx)
            end
            d
        end
    else
        nothing
    end

    # prepare to plot band diagram
    Npaths           = length(kpi.kpaths)
    Nbands           = sum(length, bandsv)
    Nbandsv          = length(bandsv)
    local_xs         = cumdists.(cartesianize(kpi).kpaths)
    local_xs_lengths = last.(local_xs)
    rel_xs_lengths   = local_xs_lengths./sum(local_xs_lengths)
    ylims = mapreduce(default_dispersion_ylims, extrema_tuplewise, bandsv; init=(Inf, -Inf))
    axs = Matrix{typeof(S.Axis())}(undef, 1, Npaths)
    start_idx = 1
    for (path_idx, (local_x, labels)) in enumerate(zip(local_xs, kpi.labels))
        # `lines` for energy bands, in the SpecApi style
        stop_idx = start_idx+length(local_x)-1
        plots = Vector{PlotSpec}(undef, Nbands)
        j = 0
        for (i, bands) in enumerate(bandsv)
            col = color isa AbstractVector ? color[i] : color
            lw = is_tup_or_vec(linewidth) ? linewidth[i] : linewidth
            ls = is_tup_or_vec(linestyle) ? linestyle[i] : linestyle
            lb = is_tup_or_vec(label)     ? label[i]     : label
            for band in bands
                plots[j+=1] = S.Lines(local_x, band[start_idx:stop_idx];
                                      color=col, linewidth=lw, linestyle=ls, label=lb)
            end
        end

        # add annotations, if not nothing
        if !isnothing(annotations)
            if Nbandsv ≠ 1
                error("annotations are not supported for multiple band arguments")
            end
            bands = bandsv[1]
            a_local_x = Vector{Float64}() # x-coordinates for annotations
            a_y       = Vector{Float64}() # y-coordinates for annotations
            text  = Vector{String}()
            for (klab, as) in annotations
                x_idxs = get(klab2idxs[path_idx], Symbol(klab), nothing)
                isnothing(x_idxs) && continue
                for x_idx in x_idxs
                    append!(a_local_x, Iterators.repeated(local_x[x_idx], length(as)))
                    for (bands_idxs, a_label) in as
                        if length(bands_idxs) == 1
                            push!(a_y, bands[bands_idxs[1]][start_idx - 1 + x_idx])
                        else
                            y_idx = (start_idx - 1) .+ x_idx
                            y = sum(bands[b][y_idx] for b in bands_idxs) / length(bands_idxs)
                            push!(a_y, y)
                        end
                        push!(text, x_idx == firstindex(local_x) ? " "*a_label  : 
                                    x_idx == lastindex(local_x)  ? a_label*" "  : a_label)
                    end
                end
            end
            if !isempty(a_local_x)
                # special-case alignment at labels at first/last x-points
                x_min, x_max = extrema(local_x)
                align = [x==x_min ? (:left,:bottom) :
                         x==x_max ? (:right,:bottom) : (:center,:bottom) for x in a_local_x]
                # TODO: replace by `annotate` when a new Makie version is released
                t = S.Text(a_local_x, a_y; text, align, color=color, offset=(0.0, 8.0))
                push!(plots, t)
            end
        end
        
        # add plots to an Axis object
        ax = S.Axis(; plots)

        # set axis properties (must be declaritive in SpecApi, not functional)
        ax.limits = (first(local_x), last(local_x), ylims...)
        ax.xticks = ([local_x[x_idx] for x_idx in keys(labels)], #= x-tick-coords =#
                     [string(x_lab) for x_lab in values(labels)] #= x-tick-labels =# )
        ax.xgridvisible = true
        ax.ygridvisible = false
        if path_idx == 1
            ax.ylabel = ylabel
        else
            ax.yticklabelsvisible = false
        end

        # assign axis & associated plot to position in grid of axes (shortly, a GridLayout)
        axs[path_idx] = ax

        # prepare for next iteration
        start_idx = stop_idx + 1
    end

    layout = S.GridLayout(axs; colsizes = Relative.(rel_xs_lengths), yaxislinks = vec(axs))
    return layout
end

is_tup_or_vec(::AbstractVector) = true
is_tup_or_vec(::Tuple)          = true
is_tup_or_vec(_)                = false

# ---------------------------------------------------------------------------------------- #
# make it possible to also provide bands as a `Matrix` (bands in columns)

function Makie.used_attributes(kpi::KPathInterpolant, bandsv::AbstractMatrix{<:Real}...)
    Makie.used_attributes(kpi, eachcol(bandsv[1]))
end
function Makie.convert_arguments(
    P::Type{<:AbstractPlot}, # makes `plot`, `lines`, `scatter` etc. dispatch here
    kpi::KPathInterpolant,
    bandsv::AbstractMatrix{<:Real}...;
    kws...
    )
    Makie.convert_arguments(P, kpi, (eachcol(bands) for bands in bandsv)...; kws...)
end
