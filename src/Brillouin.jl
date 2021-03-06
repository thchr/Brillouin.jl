module Brillouin

using Reexport
using Requires
using StaticArrays
# ---------------------------------------------------------------------------------------- #
const AVec = AbstractVector
const BasisLike{D} = AVec{<:SVector{D,<:Real}}
const SHOWDIGITS = 6
# ---------------------------------------------------------------------------------------- #
# defines a a module `CrystallineBravaisVendor` with a single function `bravaistype`, which
# is just a hardcopy of output from Crystalline's equivalent function
include("CrystallineBravaisVendor.jl")
# ---------------------------------------------------------------------------------------- #
@enum BasisEnum CARTESIAN LATTICE
# ---------------------------------------------------------------------------------------- #
function latticize! end
latticize(v::AVec{<:Real}, basismatrix::AbstractMatrix{<:Real}) = basismatrix\v
latticize(v::AVec{<:Real}, basis::AVec{<:AVec{<:Real}}) = latticize(v, hcat(basis...))'
latticize(x) = (x.basisenum[] === CARTESIAN ? latticize!(deepcopy(x)) : deepcopy(x))

function cartesianize! end
cartesianize(v::AVec{<:Real}, basis::AVec{<:AVec{<:Real}}) = v'basis
cartesianize(x) = (x.basisenum[] === LATTICE ? cartesianize!(deepcopy(x)) : deepcopy(x))

"""
    basis(x::Union{KPath, KPathInterpolant, Cell})

Return the (reciprocal or direct) lattice basis associated with `x`, in Cartesian
coordinates.

Methods in Brillouin will by default return points in the lattice basis, i.e., points are
referred to `basis(x)`. This corresponds to the setting `x.basisenum[] == LATTICE`.
Coordinates may instead be referred to a Cartesian basis, corresponding to
`x.basisenum[] == CARTESIAN` by using [`cartesianize`](@ref). The result of `basis(x)`,
however, is invariant to this and always refers to the lattice basis in Cartesian
coordinates.
"""
basis(x) = x.basis

export latticize!, latticize,
    cartesianize!, cartesianize,
    basis
# ---------------------------------------------------------------------------------------- #

# MAIN FUNCTIONALITY

include("KPaths.jl")
@reexport using .KPaths

include("WignerSeitz.jl")
@reexport using .WignerSeitz

# ---------------------------------------------------------------------------------------- #

function __init__()
    # plotting extensions on GLMakie load
    @require AbstractPlotting="537997a7-5e4e-5d89-9595-2241ea00577e" begin
        include("requires/abstractplotting_wignerseitz.jl")
    end

    # plotting extensions on PlotlyJS load
    @require PlotlyJS="f0f68f2c-4968-5e81-91da-67840de0976a" begin
        include("requires/plotlyjs_wignerseitz.jl")
        include("requires/plotlyjs_kpaths.jl")
        include("requires/plotlyjs_dispersion.jl")
    end
end

# ---------------------------------------------------------------------------------------- #

end
