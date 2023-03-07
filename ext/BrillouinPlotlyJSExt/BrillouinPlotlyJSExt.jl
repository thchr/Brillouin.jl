module BrillouinPlotlyJSExt

using Brillouin
using Brillouin: CARTESIAN
using LinearAlgebra: norm, dot
using StaticArrays
if isdefined(Base, :get_extension)
    using PlotlyJS
    import PlotlyJS: plot
else
    using ..PlotlyJS
    import ..PlotlyJS: plot
end

include("../default_colors.jl")
include("dispersion.jl")
include("wignerseitz.jl")
include("kpaths.jl")

end # BrillouinPlotlyJSExt