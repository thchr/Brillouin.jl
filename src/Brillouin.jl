module Brillouin

using Reexport
using StaticArrays
# ---------------------------------------------------------------------------------------- #
const AVec = AbstractVector
const BasisLike{D} = AVec{<:SVector{D,<:Real}}
const SHOWDIGITS = 6
# ---------------------------------------------------------------------------------------- #
"""
    BasisEnum <: Enum{Int8}

Enum type with instances

```
CARTESIAN
LATTICE
```

Used to indicate whether the coordinates of a `KPath` instance are referred to a lattice
basis (`LATTICE`) or a cartesian basis (`CARTESIAN`).

See also [`setting`](@ref) and `set_setting!`.
"""
@enum BasisEnum CARTESIAN LATTICE

"""
    setting(x::Union{KPath, KPathInterpolant, Cell})

Return the basis setting of coordinates in `x`. The returned value is a member of the
[`BasisEnum`](@ref) enum with member values `LATTICE` (i.e. coordinates in the basis of the
lattice vectors) or `CARTESIAN` (i.e. coordinates in the Cartesian basis).
By default, methods in Brillouin return coordinates in the `LATTICE` setting.
"""
setting(x) = x.setting[]

set_setting!(x, new_setting) = (x.setting[] = new_setting)

# ---------------------------------------------------------------------------------------- #
"""
    basis(x::Union{KPath, KPathInterpolant, Cell})

Return the (reciprocal or direct) lattice basis associated with `x`, in Cartesian
coordinates.

Methods in Brillouin by default return points in the lattice basis, i.e., points are 
referred to `basis(x)`. This corresponds to the `setting(x) == LATTICE`. Coordinates may
instead be referred to a Cartesian basis, corresponding to `setting(x) == CARTESIAN` by
using [`cartesianize`](@ref). The result of  `basis(x)`, however, is invariant to this and
always refers to the lattice basis in Cartesian coordinates.
"""
basis(x) = x.basis

export basis, setting

# ---------------------------------------------------------------------------------------- #
# MAIN FUNCTIONALITY

include("KPaths.jl")
@reexport using .KPaths

include("WignerSeitz.jl")
@reexport using .WignerSeitz

# ---------------------------------------------------------------------------------------- #
# SNOOPPRECOMPILE

using PrecompileTools
include("precompile.jl")

# ---------------------------------------------------------------------------------------- #
# EXTENSIONS: loaded via Requires.jl on Julia versions <v1.9; otherwise via the Pkg
#             extension system

if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    function __init__()
        @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
            include("../ext/BrillouinMakieExt/BrillouinMakieExt.jl")
        end
        @require PlotlyJS="f0f68f2c-4968-5e81-91da-67840de0976a" begin
            include("../ext/BrillouinPlotlyJSExt/BrillouinPlotlyJSExt.jl")
        end
        @require Spglib="f761d5c5-86db-4880-b97f-9680a7cccfb5" begin
            include("../ext/BrillouinSpglibExt.jl")
        end
    end
end

# ---------------------------------------------------------------------------------------- #

end
