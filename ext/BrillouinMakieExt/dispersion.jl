# Implement `plot(kpi::KPathInterpolant, bands)` for Makie.jl: the key challenge to doing
# this being able to create multiple `Axis` objects (for disjoint parts of `kpi`), which
# is only possible via the SpecApi interface, which is slightly different but general than
# the `@recipe` interface

import Makie.SpecApi as S

## --------------------------------------------------------------------------------------- #

# keyword arguments must be explicitly marked in SpecApi :(
function Makie.used_attributes(::KPathInterpolant, ::AbstractVector{<:AbstractVector{<:Real}})
    (:color, :linewidth)
end

# the conversion method creates a grid of `Axis` objects with `Lines` plots inside
function Makie.convert_arguments(
    ::Type{<:AbstractPlot}, # makes `plot`, `lines`, `scatter` etc. dispatch here
    kpi::KPathInterpolant,
    bands::AbstractVector{<:AbstractVector{<:Real}};
    # keyword arguments below; must be explicit in SpecApi :(
    color     = BAND_COL[],
    linewidth = 3 # or, could have been set to `Makie.theme(:linewidth)` (=1.5, I think)
    )

    # check input
    Nk = length(kpi)
    if !all(band -> length(band) == Nk, bands)
        throw(DimensionMismatch("mismatched dimensions of `kpi` and `bands`"))
    end

    # prepare to plot band diagram
    Npaths           = length(kpi.kpaths)
    local_xs         = cumdists.(cartesianize(kpi).kpaths)
    local_xs_lengths = last.(local_xs)
    rel_xs_lengths   = local_xs_lengths./sum(local_xs_lengths)
    ylims = default_dispersion_ylims(bands)
    axs = Matrix{typeof(S.Axis())}(undef, 1, Npaths)
    start_idx = 1
    for (path_idx, (local_x, labels)) in enumerate(zip(local_xs, kpi.labels))
        # `lines` for energy bands, in the SpecApi style
        stop_idx = start_idx+length(local_x)-1
        plots = [S.Lines(local_x, band[start_idx:stop_idx]; color, linewidth) for band in bands]
        
        # add plots to an Axis object
        ax = S.Axis(; plots)

        # set axis properties (must be declaritive in SpecApi, not functional)
        ax.limits = (first(local_x), last(local_x), ylims...)
        ax.xticks = ([local_x[x_idx] for x_idx in keys(labels)], #= x-tick-coords =#
                     [string(x_lab) for x_lab in values(labels)] #= x-tick-labels =# )
        ax.xgridvisible = true
        ax.ygridvisible = false
        if path_idx == 1
            ax.ylabel = "Energy (arb. units)"
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

# ---------------------------------------------------------------------------------------- #
# make it possible to also provide bands as a `Matrix` (bands in columns)

function Makie.used_attributes(kpi::KPathInterpolant, bands::AbstractMatrix{<:Real})
    Makie.used_attributes(kpi, eachcol(bands))
end
function Makie.convert_arguments(
    P::Type{<:AbstractPlot}, # makes `plot`, `lines`, `scatter` etc. dispatch here
    kpi::KPathInterpolant,
    bands::AbstractMatrix{<:Real};
    color     = BAND_COL[],
    linewidth = 3
    )
    Makie.convert_arguments(P, kpi, eachcol(bands); color, linewidth)
end
