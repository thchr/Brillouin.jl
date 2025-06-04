## --------------------------------------------------------------------------------------- #
# Default colors

const TRANSPARENT_COL = Ref("rgba(255, 255, 255, 1)")

# --- wigner-seitz colors ---
# from https://flatuicolors.com/palette/gb
const BZ_COL          = Ref("rgb(47,54,64)")    # "electromagnetic"
const BASIS_COL       = Ref("rgb(39,60,117)")   # "pico void"
const BASIS_LIGHT_COL = Ref("rgb(212,216,227)") # 20% BASIS_COL, 80% white
const AXIS_COL        = Ref("rgb(194,54,22)")   # "harley davidson orange"
const AXIS_LIGHT_COL  = Ref("rgb(242,215,208)") # 20% AXIS_COL, 80% white

# --- k-path colors ---
# from https://flatuicolors.com/palette/ca
const KPATH_COL = Ref("rgb(95,39,205)")         # "nasu purple"

# --- dispersion colors ---
# from https://flatuicolors.com/palette/gb
const BAND_COL     = Ref("rgb(39,60,117)")      # "mazarine blue"
const KLINE_COL    = Ref("rgb(220,221,225)")    # "hint of pensive" (light gray)
const ANNOTATE_COL = Ref("rgb(53, 59, 72)")     # "blue nights" (dark gray)

## --------------------------------------------------------------------------------------- #
# Utilities

extrema_skipnan(x; kws...) = extrema(Iterators.filter(!isnan, x); kws...)
function extrema_tuplewise(x::Tuple{Real, Real}, y::Tuple{Real, Real})
    return (min(x[1], y[1]), max(x[2], y[2]))
end
function default_dispersion_ylims(bands::AbstractVector)
    ylims = mapreduce(band->extrema_skipnan(band; init=(Inf, -Inf)),
                      extrema_tuplewise,
                      bands; init=(Inf, -Inf))
    pad_ylims(ylims)
end
function default_dispersion_ylims(bands::AbstractMatrix{<:Real})
    return pad_ylims(extrema_skipnan(bands; init=(Inf, -Inf)))
end
function pad_ylims(ylims::Tuple{Real, Real})
    pad = (ylims[2]-ylims[1]) * 0.05
    if isapprox(ylims[1], zero(pad), atol=1e-4)
        return (ylims[1], ylims[2] + pad)
    else
        return (ylims[1] - pad, ylims[2] + pad)
    end
end