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
