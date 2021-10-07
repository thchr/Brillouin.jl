# Brillouin

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://thchr.github.io/Brillouin.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://thchr.github.io/Brillouin.jl/dev)
[![Build Status](https://github.com/thchr/Brillouin.jl/workflows/CI/badge.svg)](https://github.com/thchr/Brillouin.jl/actions)
[![Coverage](https://codecov.io/gh/thchr/Brillouin.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/thchr/Brillouin.jl)


Brillouin.jl is a Julia package that provides tools to generate and visualize **k**-space paths and Brillouin zones for eigenvalue problems in crystals.
The **k**-path functionalities are inspired by the Python [SeeK-path package](https://github.com/giovannipizzi/seekpath) (and return equivalent paths in 3D).

## Examples

### Wigner-Seitz cells & Brillouin zones

To generate the Brillouin zone of a crystal in space group 147 (Hermann-Mauguin symbol, P-3; Bravais type, hP), we first define its reciprocal basis `Gs` (e.g., using [Bravais.jl](https://github.com/thchr/Crystalline.jl)) and then call Brillouin's `wignerseitz`:
```jl
julia> using Brillouin, 
julia> using Bravais: reciprocalbasis
julia> Rs = ([1.0, 0.0, 0.0], [-0.5, sqrt(3)/2, 0.0],   [0, 0, 1.25]) # direct basis for space group 147
julia> Gs = reciprocalbasis(Rs) # using Bravais to create the reciprocal basis
julia> cell = wignerseitz(Gs)   # construct associated Brillouin zone
Cell{3} (8 faces, 12 vertices):
  verts: [0.666667, -0.333333, -0.5]
         [0.333333, -0.666667, -0.5]
         [0.666667, -0.333333, 0.5]
         [0.333333, 0.333333, 0.5]
         [0.333333, 0.333333, -0.5]
         [0.333333, -0.666667, 0.5]
         [-0.333333, 0.666667, -0.5]
         [-0.666667, 0.333333, -0.5]
         [-0.333333, -0.333333, -0.5]
         [-0.333333, -0.333333, 0.5]
         [-0.666667, 0.333333, 0.5]
         [-0.333333, 0.666667, 0.5]
  faces: [5, 4, 3, 1]
         [8, 9, 10, 11]
         [2, 1, 3, 6]
         [2, 6, 10, 9]
         [7, 5, 1, 2, 9, 8]
         [4, 12, 11, 10, 6, 3]
         [4, 5, 7, 12]
         [11, 12, 7, 8]
  basis: [6.283185, 3.627599, -0.0]
         [0.0, 7.255197, 0.0]
         [0.0, -0.0, 5.026548]
```
The returned vertices are in the coordinates of the reciprocal basis (to convert, see `cartesianize(!)`); this is the default behavior in Brillouin. The basis is accessible with `basis(cell)`.

The Brillouin zone can be plotted using e.g. [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) (or 3D-capable backends of [Makie.jl](https://github.com/JuliaPlots/Makie.jl) such as [GLMakie.jl](https://github.com/JuliaPlots/GLMakie.jl)):
```jl
julia> using PlotlyJS
julia> plot(cell)
```
Examples of interactive visualizations are [included in the documentation](https://thchr.github.io/Brillouin.jl/stable/wignerseitz/).

### Minimal **k**-paths in the irreducible Brillouin zone

Given a symmetry setting and a lattice, specified by a space group number `sgnum` and a conventional direct basis `Rs` (respecting the conventions of the International Tables of Crystallography, Volume A), `irrfbz_path` will return a "minimal" **k**-path in the irreducible Brillouin zone. E.g.,
```jl
julia> sgnum = 147
julia> kp = irrfbz_path(sgnum, Rs)
KPath{3} (7 points, 3 paths, 13 points in paths):
 points: :M => [0.5, 0.0, 0.0]
         :A => [0.0, 0.0, 0.5]
         :H => [0.333333, 0.333333, 0.5]
         :K => [0.333333, 0.333333, 0.0]
         :Γ => [0.0, 0.0, 0.0]
         :L => [0.5, 0.0, 0.5]
         :H₂ => [0.333333, 0.333333, -0.5]
  paths: [:Γ, :M, :K, :Γ, :A, :L, :H, :A]
         [:L, :M]
         [:H, :K, :H₂]
  basis: [6.283185, 3.627599, -0.0]
         [0.0, 7.255197, 0.0]
         [0.0, -0.0, 5.026548]
```
Returned **k**-vector coordinates are referred to the basis of the primitive reciprocal cell (in the CDML setting). The associated transformation matrix between conventional and primitive bases can be obtained via [Bravais.jl's `primitivebasismatrix`].

The resulting object `kp` can be interpolated, using either `interpolate(kp, N)` or `splice(kp, N)` which return a `KPathInterpolant` iterable whose values interpolate the connected paths (and enable convenient plotting of [band structure diagrams](https://thchr.github.io/Brillouin.jl/stable/kpaths/#Band-structure)). 
See also visualization [examples in documentation](https://thchr.github.io/Brillouin.jl/stable/kpaths/).