# **k**-space paths

## Default paths

To generate a **k**-path for, say, the space group of diamond (space group 227; a cubic face-centered Bravais lattice), we can call [`irrfbz_path`](@ref), which will return a minimal path in the irreducible Brillouin zone:
```@example kpath
using Brillouin
sgnum = 227
Rs = [[1,0,0], [0,1,0], [0,0,1]] # conventional direct basis
kp = irrfbz_path(sgnum, Rs)
```
The path data is sourced from the [HPKOT paper](http://dx.doi.org/10.1016/j.commatsci.2016.10.015) (or, equivalently, the [SeeK-path](https://github.com/giovannipizzi/seekpath) Python package). 

The coordinates of the path are given with respect to the *primitive* reciprocal basis (here, `[[-2π,2π,2π], [2π,-2π,2π], [2π,2π,-2π]]`). To convert to a Cartesian basis, we can use [`cartesianize`](@ref) or [`cartesianize!`](@ref) (in-place):
```@example kpath
cartesianize(kp)
```
We can visualize the **k**-path using [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) (conversion to a Cartesian basis for plotting is automatic):
```@example kpath
using PlotlyJS
Pᵏ = plot(kp)
Main.HTMLPlot(Pᵏ) # hide
```

Usually, it'll be more helpful to understand the path's geometry in the context of the associated Brillouin zone. To visualize this, we can plot the combination of a `Cell` (created via [`wignerseitz`](@ref)) and a `KPath`:
```@example kpath
pGs = basis(kp)      # primitive reciprocal basis associated with k-path
c = wignerseitz(pGs) # associated Brillouin zone
Pᶜ⁺ᵏ = plot(c, kp)
Main.HTMLPlot(Pᶜ⁺ᵏ) # hide
```

## Interpolation
Interpolation of a `KPath` structure can be achieved using either [`interpolate(::KPath, ::Integer)`](@ref) or [`splice(::KPath, ::Integer)`](@ref), returning a `KPathInterpolant`.
As an example, `interpolate(kp, N)` returns an interpolation with a *target* of `N` interpolation points, distributed as equidistantly as possible (with the distance metric evaluated in Cartesian space):
```@example kpath
kpi = interpolate(kp, 100)
```
The returned `KPathInterpolant` implements the `AbstractVector` interface, with iterants returning `SVector{D, Float64}` elements.
To get a conventional "flat" vector, we can simply call `collect(kpi)`.

Internally, `KPathInterpolant` includes additional structure and information: namely, the high-symmetry points and associated labels along the path as well as a partitioning into connected vs. disconnected path segments.

## Band structure
The additional structure of `KPathInterpolation` enables convenient and clear visualizations of band structure diagrams in combination with PlotlyJS.

To illustrate this, suppose we are considering a tight-binding problem for an *s*-orbital situated at the 1a Wyckoff position. Such a problem has a single band with dispersion [^1] (assuming a cubic side length ``a = 1``):
```math
\epsilon(\mathbf{k})
=
4\gamma\Bigl(
    \cos \tfrac{1}{2}k_x \cos \tfrac{1}{2}k_y +
    \cos \tfrac{1}{2}k_y \cos \tfrac{1}{2}k_z +
    \cos \tfrac{1}{2}k_z \cos \tfrac{1}{2}k_x
    \Bigr)
```
with $k_{x,y,z}$ denoting coordinates in a Cartesian basis (which are related to the coordinates $k_{1,2,3}$ in a primitive reciprocal basis by $k_x = 2π(-k_1+k_2+k_3)$, $k_x = 2π(k_1-k_2+k_3)$, and $k_z = 2π(k_1+k_2-k_3)$).

We can calculate the energy band along our **k**-path using the interpolation object `kpi`. To do so, we define a function that implements $\epsilon(\mathbf{k})$ and broadcast it over the elements of `kpi`:
```@example kpath
function ϵ(k; γ::Real=1.0)
    kx = 2π*(-k[1]+k[2]+k[3])
    ky = 2π*(+k[1]-k[2]+k[3])
    kz = 2π*(+k[1]+k[2]-k[3])
    return 4γ * (cos(kx/2)*cos(ky/2) + cos(ky/2)*cos(kz/2) + cos(kz/2)*cos(kx/2))
end
band = ϵ.(kpi)
```

Finally, we can visualize the associated band using a Brillouin-overloaded PlotlyJS `plot`-call:
```@example kpath
P = plot(kpi, [band])
Main.HTMLPlot(P, 525) # hide
```

If we have multiple bands, say ``\epsilon_1(\mathbf{k}) = \epsilon(\mathbf{k})`` and ``\epsilon_2(\mathbf{k}) = 20 - \tfrac{1}{2}\epsilon(\mathbf{k})``, we can easily plot that by collecting the bands in a single vector (or concatenating into a matrix):
```@example kpath
band1 = ϵ.(kpi)
band2 = 20 .- (1/2).*band1
P¹² = plot(kpi, [band1, band2])
Main.HTMLPlot(P¹², 525) # hide
```

```@docs
plot(::KPathInterpolant, ::Any, ::Layout)
```

## K-paths for non-standard lattices
One can create a **k**-path for a non-standard unit cell by using Spglib.
```@example kpath
using Spglib
# Trigonal lattice, space group R-3m, 166
a = 1.0
c = 8.0
lattice_standard    = SVector(SVector(a*sqrt(3)/2,  a/2, c/3),
                              SVector(-a*sqrt(3)/2, a/2, c/3),
                              SVector(0.0, -a, c/3))
lattice_nonstandard = SVector(SVector(a*sqrt(3)/2, -a/2, c/3),
                              SVector(0.0, a, c/3),
                              SVector(-a*sqrt(3)/2, -a/2, c/3))

cell_standard = Spglib.Cell(lattice_standard, [[0, 0, 0]], [0])
cell_nonstandard = Spglib.Cell(lattice_nonstandard, [[0, 0, 0]], [0])

kp_standard = irrfbz_path(cell_standard)
kp_nonstandard = irrfbz_path(cell_nonstandard)
```

One can check that the generated **k**-paths for the non-standard and standard lattice vectors are equivalent by plotting the path and the Wigner-Seitz cell.
```@example kpath
using Bravais
PlotlyJS.plot(wignerseitz(Bravais.reciprocalbasis(Bravais.DirectBasis(lattice_standard))), kp_standard, Layout(title="standard cell"))
PlotlyJS.plot(wignerseitz(Bravais.reciprocalbasis(Bravais.DirectBasis(lattice_nonstandard))), kp_nonstandard, Layout(title="non-standard cell"))
```

[^1] See e.g. [http://www.physics.rutgers.edu/~eandrei/chengdu/reading/tight-binding.pdf](http://www.physics.rutgers.edu/~eandrei/chengdu/reading/tight-binding.pdf)