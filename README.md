# Brillouin

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://thchr.github.io/Brillouin.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://thchr.github.io/Brillouin.jl/dev)
[![Build Status](https://github.com/thchr/Brillouin.jl/workflows/CI/badge.svg)](https://github.com/thchr/Brillouin.jl/actions)
[![Coverage](https://codecov.io/gh/thchr/Brillouin.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/thchr/Brillouin.jl)


Brillouin.jl is a Julia package that provides tools to generate and visualize **k**-space paths and Brillouin zones that arise for crystalline eigenvalue problems.
The functionalities are inspired by the Python [SeeK-path package](https://github.com/giovannipizzi/seekpath).

## Examples

Generate the Brillouin zone of a crystal in space group 147 (Bravais type hP; using [Crystalline.jl](https://github.com/thchr/Crystalline.jl) to easily obtain an associated reciprocal basis):
```jl
julia> using Brillouin, Crystalline
julia> Rs   = ([1.0, 0.0, 0.0], [-0.5, sqrt(3)/2, 0.0],   [0, 0, 1.25]) # or `Crystalline.directbasis(147, 3)`
julia> Gs   = reciprocalbasis(Rs)
julia> cell = wignerseitz(Gs)
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
The returned vertices are in the coordiantes of the provided reciprocal basis (to convert, see `cartesianize(!)`); this is the default behavior in Brillouin.

The Brillouin zone can be plotted using e.g. [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) (or 3D-capable backends of [AbstractPlotting.jl](https://github.com/JuliaPlots/AbstractPlotting.jl) such as [GLMakie.jl](https://github.com/JuliaPlots/GLMakie.jl)):
```jl
julia> using PlotlyJS
julia> plot(cell)
```
Examples of interactive visualizations are [included in the documentation](https://thchr.github.io/Brillouin.jl/stable/wignerseitz/). Visualizations will automatically be displayed in Cartesian coordinates.

Irreducible **k**-paths are returned by `irrfbz_path`, and can similarly be visualized (see [examples in documentation](https://thchr.github.io/Brillouin.jl/stable/kpaths/)):
```jl
julia> kp = irrfbz_path(147, Rs)
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
The resulting object can be interpolated, using either `interpolate(kp, N)` or `splice(kp, N)`, which produces an `KPathInterpolant` iterable whose elements interpolate the connected paths (and enable convenient plotting of [band structure diagrams](https://thchr.github.io/Brillouin.jl/stable/kpaths/#Band-structure)).