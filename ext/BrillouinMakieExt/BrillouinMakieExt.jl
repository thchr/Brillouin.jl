module BrillouinMakieExt
    
using Brillouin
using Brillouin: CARTESIAN
using LinearAlgebra: norm
using StaticArrays: normalize
if isdefined(Base, :get_extension)
    using Makie
else
    using ..Makie
end

include("../default_colors.jl")
include("wignerseitz.jl")
include("kpaths.jl")

end # BrillouinMakieExt