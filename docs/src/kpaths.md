# **k**-paths

To generate a **k**-path for, say, space group 202 (a cubic face-centered Bravais lattice), we can call [`irrfbz_path`](@ref), which will return a minimal path in the irreducible Brillouin zone:
```@example kpath
using Brillouin
sgnum = 202
Rs = [[1,0,0], [0,1,0], [0,0,1]] # conventional direct basis
kp = irrfbz_path(sgnum, Rs)
```
The path data is sourced from the [HPKOT paper](http://dx.doi.org/10.1016/j.commatsci.2016.10.015) (or, equivalently, the [SeeK-path](https://github.com/giovannipizzi/seekpath) Python package).

The resulting `KPath` structure initially gives the **k**-point coordinates in the basis of the *primitive* reciprocal basis. To convert to a Cartesian basis, we use [`cartesianize!`](@ref):
```@example kpath
pGs = 2π.*[[-1,1,1], [1,-1,1], [1,1,-1]] # primitive reciprocal basis
cartesianize!(kp, pGs)
```
We can visualize the **k**-path using [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl):
```@example kpath
using PlotlyJS
Pᵏ = plot(kp)
Main.HTMLPlot(Pᵏ) # hide
```

Usually though, it'll be more helpful to understand the path's geometry in the context of the associated Brillouin zone. To visualize this, we can simply plot the combination of a `Cell` (created via [`wignerseitz`](@ref)) and a `KPath`:
```@example kpath
c = wignerseitz(pGs)
Pᶜ⁺ᵏ = plot(c, kp)
Main.HTMLPlot(Pᶜ⁺ᵏ) # hide
```

## Interpolation
Interpolation of a `KPath` structure can be achieved using either [`interpolate(::KPath, ::Integer)`](@ref) or [`splice(::KPath, ::Integer)`](@ref).