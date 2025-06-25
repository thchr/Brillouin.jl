# Implement `plot(kpi::KPathInterpolant, bands)` for Makie.jl: the key challenge to doing
# this being able to create multiple `Axis` objects (for disjoint parts of `kpi`), which
# is only possible via the SpecApi interface, which is slightly different but general than
# the `@recipe` interface

import Makie.SpecApi as S

## --------------------------------------------------------------------------------------- #

# keyword arguments must be explicitly marked in SpecApi :(
function Makie.used_attributes(::KPathInterpolant, ::AbstractVector{<:AbstractVector{<:Real}})
    (:color, :linewidth, :linestyle, :ylabel, :label, :annotations, :ylims)
end

#=
Interim doc-string, until I figure out how to add one properly (TODO):

Plot a band structure over a `KPathInterpolant`, for vararg band arguments `bandsv`.

## Keyword arguments
- `annotations`: allows annotating text to high-symmetry points in the band structure.
  E.g., to label the 1st band at the M point, the 2nd and 3rd band at Y, and the
  average energy of the 4th and 5th band at X, `annotations` could be set to:
    Dict(:M => [1:1 => "Band 1 at X"], 
         :Y => [2:2 => "Band 2 at Y, 3 => "Band 3 at Y"] # scalar or a 1-element range is equivalent
         :X => [4:5 => "Average of bands 4 and 5 at X"])
- TODO: 
=#
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
    ylims::Union{Nothing, Tuple{Real, Real}} = nothing
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
    ylims = @something(
        ylims, 
        mapreduce(default_dispersion_ylims, extrema_tuplewise, bandsv; init=(Inf, -Inf))
    )
    # avoid subnormal zeros, cf. `ax.limits` NB comment below
    iszero(ylims[1]) && (ylims = (-1.0e-307, ylims[2]))
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
                        push!(text, a_label)
                    end
                end
            end
            if !isempty(a_local_x)
                # special-case alignment at labels at first/last x-points
                t = S.Annotation(a_local_x, a_y; text, color=color)
                m = S.Scatter(a_local_x, a_y; color=color, markersize=2.5*linewidth, 
                              strokecolor=:white, strokewidth=linewidth*.65)
                push!(plots, t, m)
            end
        end
        
        # add plots to an Axis object
        ax = S.Axis(; plots)

        # set axis properties (must be declaritive in SpecApi, not functional)
        ax.limits = (first(local_x) - 1.0e-307, last(local_x), ylims...)
        # NB: ↑ we subtract `1.0e-307` from the the first x-coordinate above (whose value is
        #     0.0) because there is an annoying bug, related to subnormal zeros, where the
        #     first x-tick will be omitted if its value and lowest `limits` value are both
        #     0.0, basically because 0.0 and -0.0 do not compare meaningfully under Julia's
        #     default choice of `set_zero_subnormals(true)`. The value 1.0e-307 is the
        #     smallest value x such that `0.0 - x ≠ 0.0` is true
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
