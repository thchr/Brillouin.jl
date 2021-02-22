module Brillouin

using Reexport
using StaticArrays
# ---------------------------------------------------------------------------------------- #
const AVec = AbstractVector
const BasisLike{D} = AVec{<:SVector{D,<:Real}}
# ---------------------------------------------------------------------------------------- #
# defines a a module `CrystallineBravaisVendor` with a single function `bravaistype`, which
# is just a hardcopy of output from Crystalline's equivalent function
include("CrystallineBravaisVendor.jl")
# ---------------------------------------------------------------------------------------- #

include("KPath.jl")
@reexport using .KPath

include("WignerSeitz.jl")
@reexport using .WignerSeitz

end
