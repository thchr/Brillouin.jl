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

function default_dispersion_ylims(bands)
    ylims = mapreduce(extrema,
                      (x,y)->(min(x[1],y[1]), max(x[2],y[2])),
                      bands, init=(Inf, -Inf))
    δ = (ylims[2]-ylims[1])/20
    if isapprox(ylims[1], 0, atol=1e-5)
        return (ylims[1], ylims[2] + δ)
    else
        return (ylims[1] - δ, ylims[2] + δ)
    end
end
