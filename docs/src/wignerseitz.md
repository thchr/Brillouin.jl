# Wigner–Seitz cells

## Unit cells
The [Wigner–Seitz cell](https://en.wikipedia.org/wiki/Wigner–Seitz_cell) associated with an arbitrary lattice basis can be generated via [`wignerseitz`](@ref).
For example, to generate the unit cell of a (primitive) lattice with Bravais type cF, we might write:
```@example wignerseitz-cF
using Brillouin

Rs = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
cᴿ = wignerseitz(Rs)
```
Note that the coordinates of the Wigner-Seitz cell vertices are referred to the basis `Rs`; to convert to Cartesian space, see [`cartesianize`](@ref) and [`cartesianize!`](@ref) (in-place).

We can plot the generated cells using e.g. [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) via `plot(cᴿ)` (or, alternatively, via a 3D-capable backend of [Makie.jl](https://github.com/JuliaPlots/Makie.jl)):
```@example wignerseitz-cF
using PlotlyJS
Pᴿ = plot(cᴿ)
Main.HTMLPlot(Pᴿ) # hide
```

## Brillouin zones
To generate Brillouin zones, we simply give the corresponding reciprocal lattice `Gs`:
```@example wignerseitz-cF
Gs = 2π.*[[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]] # reciprocal basis of `Rs`
cᴳ = wignerseitz(Gs)
Pᴳ = plot(cᴳ)
Main.HTMLPlot(Pᴳ) # hide
```

## Two dimensions

`wignerseitz` and `plot(::Cell)` works in two dimensions as well. As an example, we can illustrate the Wigner–Seitz unit cell of graphene (which has a hexagonal "hp" Bravais type):
```@example wignerseitz-2d
using Brillouin, PlotlyJS

Rs = [[1.0, 0.0], [-0.5, √3/2]]
cᴿ = wignerseitz(Rs)
Pᴿ = plot(cᴿ)
Main.HTMLPlot(Pᴿ) # hide
```
and its associated Brillouin zone:
```@example wignerseitz-2d
Gs = 2π.*[[1.0, 1/√3], [0.0, 2/√3]]
cᴳ = wignerseitz(Gs)
Pᴳ = plot(cᴳ)
Main.HTMLPlot(Pᴳ) # hide
```

## Reducing points to the interior of Wigner-Seitz cell

Points can be reduced to their equivalent image in the Wigner-Seitz cell via  [`reduce_to_wignerseitz`](@ref).
The resulting point has the least possible norm of all equivalent images.

As an example, consider the reduction of a point `v` initially outside the 2D Brillouin zone defined earlier:
```@example wignerseitz-2d
v = [0.8, 0.2]
v′ = reduce_to_wignerseitz(v, Gs)
```
We can visualize the reduction by adding the points to the previous plot of the Brillouin zone:
```@example wignerseitz-2d
vᶜ = cartesianize(v, Gs)
v′ᶜ = cartesianize(v′, Gs)

Pᴳ = deepcopy(Pᴳ)
addtraces!(Pᴳ, scatter(x=vᶜ[1:1],  y=vᶜ[2:2],  hovertext="v (unreduced)"))
addtraces!(Pᴳ, scatter(x=v′ᶜ[1:1], y=v′ᶜ[2:2], hovertext="v′ (reduced)"))
Main.HTMLPlot(Pᴳ) # hide
```
