# **k**-space paths

## Default paths

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
Interpolation of a `KPath` structure can be achieved using either [`interpolate(::KPath, ::Integer)`](@ref) or [`splice(::KPath, ::Integer)`](@ref), returning a `KPathInterpolant`.
The `KPathInterpolant` implements the `AbstractVector` interface, with `SVector{D, Float64}` elements. Collecting a `KPathInterpolant` (via `collect`) therefore returns an explicit "flat" vector of interpolated **k**-points.

## Band structure
`KPathInterpolation`s enable convenient visualization of band structure diagrams by overloading PlotlyJS.
As an example, consider a structure in a face-centered cubic lattice (cF), as realized e.g. by diamond with space group 227.
First, we construct the **k**-space path:
```@example band-diagram
using Brillouin
Rs  = [[1,0,0], [0,1,0], [0,0,1]]          # conventional direct basis
pGs = 2π .* [[-1,1,1], [1,-1,1], [1,1,-1]] # primitive reciprocal basis
kp  = irrfbz_path(227, Rs)                 # k-space path
cartesianize!(kp, pGs)                     # convert k-space path to Cartesian
```

As before, we can visualize the path:
```@example band-diagram
using PlotlyJS
c = wignerseitz(pGs)
Pᶜ⁺ᵏ = plot(c, kp)
Main.HTMLPlot(Pᶜ⁺ᵏ) # hide
```

And we can construct an interpolant for it (here, with a target of 100 interpolation points):
```@example band-diagram
kpi = interpolate(kp, 100)
```

Now, suppose we are considering a tight-binding problem for an *s*-orbital situated at the 1a Wyckoff position. Such a problem has a single band with dispersion [^1] (assuming a cubic side length ``a = 1``):
```math
\epsilon(\mathbf{k}) =
4\gamma\Bigl(
    \cos \tfrac{1}{2}k_x \cos \tfrac{1}{2}k_y +
    \cos \tfrac{1}{2}k_y \cos \tfrac{1}{2}k_z +
    \cos \tfrac{1}{2}k_z \cos \tfrac{1}{2}k_x
    \Bigr)
```

We can calculate the associated energy band along our earlier **k**-path interpolation easily:
```@example band-diagram
function ε(k; γ::Real=1.0)
    4γ * (cos(k[1]/2)*cos(k[2]/2) + cos(k[2]/2)*cos(k[3]/2) + cos(k[3]/2)*cos(k[1]/2))
end
band = ε.(kpi)
```

And, finally, we can visualize the associated band (again, via PlotlyJS):
```@example band-diagram
P = plot(kpi, [band])
Main.HTMLPlot(P, 525) # hide
```

If we have multiple bands, say ``\epsilon_1(\mathbf{k}) = \epsilon(\mathbf{k})`` and ``\epsilon_2(\mathbf{k}) = 20 - \tfrac{1}{2}\epsilon(\mathbf{k})``, we can easily plot that by collecting the bands in a single vector (or concatenating into a matrix):
```@example band-diagram
band1 = ε.(kpi)
band2 = 20 .- (1/2).*band1
P¹² = plot(kpi, [band1, band2])
Main.HTMLPlot(P¹², 525) # hide
```

[^1] See e.g. http://www.physics.rutgers.edu/~eandrei/chengdu/reading/tight-binding.pdf