module BrillouinMakieExt
    
using Brillouin
using Brillouin: CARTESIAN, cumdists, cartesianize, KPathInterpolant
using LinearAlgebra: norm
using StaticArrays: normalize
if isdefined(Base, :get_extension)
    using Makie
else
    using ..Makie
end

include("../shared_plotting_utils.jl")
include("dispersion.jl")
include("wignerseitz.jl")
include("kpaths.jl")
include("plot_overload_hack.jl")

end # BrillouinMakieExt