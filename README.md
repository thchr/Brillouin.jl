# Brillouin

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://thchr.github.io/Brillouin.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://thchr.github.io/Brillouin.jl/dev)
[![Build Status](https://github.com/thchr/Brillouin.jl/workflows/CI/badge.svg)](https://github.com/thchr/Brillouin.jl/actions)
[![Coverage](https://codecov.io/gh/thchr/Brillouin.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/thchr/Brillouin.jl)


Brillouin.jl is a Julia package that provides tools to generate and visualize **k**-space paths and Brillouin zones that arise for crystalline eigenvalue problems.
The functionalities are inspired by the Python [SeeK-path package](https://github.com/giovannipizzi/seekpath).

## Examples

Generate the Brillouin zone of a crystal in space group 147 (Bravais type hP):
```jl
julia> using Brillouin, Crystalline
julia> Rs   = directbasis(147)
julia> Gs   = reciprocalbasis(Rs)
julia> cell = wignerseitz(Gs)
Cell{3} (8 faces, 12 vertices):
 verts: [-0.0, -0.577, -0.866]
        [0.0, 0.577, -0.866]
        [0.5, 0.289, -0.866]
        [0.5, -0.289, -0.866]
        [-0.0, -0.577, 0.866]
        [0.0, 0.577, 0.866]
        [-0.5, 0.289, 0.866]
        [-0.5, -0.289, 0.866]
        [-0.5, -0.289, -0.866]
        [-0.5, 0.289, -0.866]
        [0.5, -0.289, 0.866]
        [0.5, 0.289, 0.866]
 faces: [4, 3, 12, 11]
        [2, 6, 12, 3]
        [11, 12, 6, 7, 8, 5]
        [10, 9, 8, 7]
        [7, 6, 2, 10]
        [2, 3, 4, 1, 9, 10]
        [5, 8, 9, 1]
        [4, 11, 5, 1]
```
The resulting Brillouin zone can be plotted using e.g. [`GLMakie.jl`](https://github.com/JuliaPlots/GLMakie.jl) or other 3D-capable backends of [`AbstractPlotting.jl`](https://github.com/JuliaPlots/AbstractPlotting.jl):
```jl
julia> using GLMakie
julia> plot(cell)
```