# Wigner–Seitz cells

## Unit cells
The [Wigner–Seitz cell](https://en.wikipedia.org/wiki/Wigner–Seitz_cell) associated with an arbitrary lattice basis can be generated via [`wignerseitz`](@ref).
For example, to generate the unit cell of a (primitive) lattice with Bravais type cF, we might write:
```@example wignerseitz-cF
using Brillouin

Rs = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
cᴿ = wignerseitz(Rs)
```

We can plot the generated cells using e.g. [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) via `plot(cᴿ)` (or, alternatively, via a 3D-capable backend of [AbstractPlotting.jl](https://github.com/JuliaPlots/AbstractPlotting.jl)):
```@example wignerseitz-cF
using PlotlyJS
Pᴿ = plot(cᴿ)
Main.HTMLPlot(Pᴿ) # hide
```

## Brillouin zones
To generate Brillouin zones, we simply give the corresponding reciprocal lattice `Gs`:
```@example wignerseitz-cF
Gs = [[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]] # reciprocal basis of `Rs`
cᴳ = wignerseitz(Gs)
Pᴳ = plot(cᴳ)
Main.HTMLPlot(Pᴳ) # hide
```
