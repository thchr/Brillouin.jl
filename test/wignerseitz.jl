using Brillouin, Test
using LinearAlgebra
using StaticArrays

@testset "WignerSeitz" begin
    # --- `wignerseitz` ---
    # hexagonal lattice (example from space group 147)
    Rs = convert(SVector{3,SVector{3,Float64}}, [[1.0,0,0], [-.5,√3/2,0], [0,0,1.25]])      # direct basis
    Gs = convert(SVector{3,SVector{3,Float64}}, [[2π,2π/√3,0], [0,4π/√3,0], [0,0,2π/1.25]]) # reciprocal basis
    cell = wignerseitz(Gs)
    @test length(cell) == 8
    @test length(vertices(cell)) == 12
    @test basis(cell) == Gs

    # test that everything works the same if we use ordinary vectors instead of SVectors
    @test wignerseitz(convert(Vector{SVector{3,Float64}}, Gs)) == cell
    @test wignerseitz(convert(Vector{Vector{Float64}}, Gs))    == cell

    # ------------------------------------------------------------------------------------ #
    # for CI and testing generally, we cannot depend on a specific sorting/ordering
    # from `wignerseitz` (and ultimately, `PySpatial.Voronoi` or `PySpatial.ConvexHull`),
    # so to be robust to this, we compare output after sorting to a canonical order
    function canonicalize!(c::Cell)
        vs = vertices(c)
        fs = faces(c)
        # sort the vertices in some canonical order (we also use `round` to be robust to 
        # floating point differences across architectures)
        I = sortperm(vs, by=v->round.(v, digits=5))
        permute!(vs, I)
        # relabel the face-indices to match the new vertex sorting
        replace!.(fs, (I .=> eachindex(I))...)
        # sort the faces in some canonical order; note that each face
        # must themselves also be sorted (because they are cyclical)
        map!(circshift_minimal_first, fs, fs)
        sort!(fs)

        return c
    end

    function circshift_minimal_first(v::Vector)
        _, idx = findmin(v)
        return circshift(v, length(v)-idx+1)
    end

    # actually canonicalize `cell`
    canonicalize!(cell)

    # ------------------------------------------------------------------------------------ #

    # test show/display (canonicalization is crucial here)
    cell_show_reference = """
    Cell{3} (8 faces, 12 vertices):
      verts: [-0.666667, 0.333333, -0.5]
             [-0.666667, 0.333333, 0.5]
             [-0.333333, -0.333333, -0.5]
             [-0.333333, -0.333333, 0.5]
             [-0.333333, 0.666667, -0.5]
             [-0.333333, 0.666667, 0.5]
             [0.333333, -0.666667, -0.5]
             [0.333333, -0.666667, 0.5]
             [0.333333, 0.333333, -0.5]
             [0.333333, 0.333333, 0.5]
             [0.666667, -0.333333, -0.5]
             [0.666667, -0.333333, 0.5]
      faces: [1, 2, 6, 5]
             [1, 3, 4, 2]
             [1, 5, 9, 11, 7, 3]
             [2, 4, 8, 12, 10, 6]
             [3, 7, 8, 4]
             [5, 6, 10, 9]
             [7, 11, 12, 8]
             [9, 10, 12, 11]
      basis: [6.283185, 3.627599, 0.0]
             [0.0, 7.255197, 0.0]
             [0.0, 0.0, 5.026548]"""

    cell_show_referenceᶜ = """
    Cell{3} (8 faces, 12 vertices):
      verts: [-4.18879, 0.0, -2.513274]
             [-4.18879, 0.0, 2.513274]
             [-2.094395, -3.627599, -2.513274]
             [-2.094395, -3.627599, 2.513274]
             [-2.094395, 3.627599, -2.513274]
             [-2.094395, 3.627599, 2.513274]
             [2.094395, -3.627599, -2.513274]
             [2.094395, -3.627599, 2.513274]
             [2.094395, 3.627599, -2.513274]
             [2.094395, 3.627599, 2.513274]
             [4.18879, 0.0, -2.513274]
             [4.18879, 0.0, 2.513274]
      faces: [1, 2, 6, 5]
             [1, 3, 4, 2]
             [1, 5, 9, 11, 7, 3]
             [2, 4, 8, 12, 10, 6]
             [3, 7, 8, 4]
             [5, 6, 10, 9]
             [7, 11, 12, 8]
             [9, 10, 12, 11]
      basis: [6.283185, 3.627599, 0.0]
             [0.0, 7.255197, 0.0]
             [0.0, 0.0, 5.026548]""" # cartesian representation
    test_show(cell, cell_show_reference)
    test_show(Brillouin.cartesianize(cell), cell_show_referenceᶜ)

    # cartesianize/latticize(!)
    @test Brillouin.latticize!(Brillouin.cartesianize!(deepcopy(cell))) ≈ cell
    cellᶜ = Brillouin.cartesianize!(deepcopy(cell))
    @test Brillouin.cartesianize!(Brillouin.latticize(cellᶜ, basis(cell)[[3,2,1]])) ≈ cellᶜ

    # test indexing/iteration of Cell struct
    @test cell[1] ≈ [[-2/3,  1/3, -1/2],
                     [-2/3,  1/3, 1/2],
                     [-1/3,  2/3, 1/2],
                     [-1/3,  2/3, -1/2]]
    @test cell[4] ≈ [[-2/3,  1/3, 1/2],
                     [-1/3, -1/3, 1/2],
                     [ 1/3, -2/3, 1/2],
                     [ 2/3, -1/3, 1/2],
                     [ 1/3,  1/3, 1/2],
                     [-1/3,  2/3, 1/2]]
    @test collect(cell) == [cell[i] for i in eachindex(cell)] # iteration vs. indexing

    # test error on out-of-bounds indexing error
    @test_throws BoundsError cell[9]

    # test :triangles output
    cell′ = wignerseitz(Gs; merge = true)
    @test all(v -> length(v) == 3, vertices(cell′))   

    # ------------------------------------------------------------------------------------ #
    # 2D test case
    Rs_2d = convert(SVector{2,SVector{2,Float64}}, [[1.0,0], [-.5,√3/2]])
    cell_2d = wignerseitz(Rs_2d)
    @test length(vertices(cell_2d)) == 6
    @test length(only(faces(cell_2d))) == 6
    @test basis(cell_2d) == Rs_2d

    # 1D test case
    Rs_1d = convert(SVector{1,SVector{1,Float64}}, [[1.5]])
    cell_1d = wignerseitz(Rs_1d)
    @test vertices(cell_1d) == [[-0.5], [0.5]]
    @test only(faces(cell_1d)) == [1,2]
    @test basis(cell_1d) == Rs_1d
    @test_throws DimensionMismatch wignerseitz([SVector(1.0), SVector(1.0)])
end

@testset "reduce_to_wignerseitz" begin
    Rs = convert(SVector{3,SVector{3,Float64}}, [[1.0,0,0], [-.5,√3/2,0], [0,0,1.25]])
    cell = wignerseitz(Rs)
    vs = vertices(cell)

    # we don't guarantee anything about points exactly on the border of the Wigner-Seitz
    # both anything else should be solid: hence, check that points slightly inside stay
    # put after reduction, while points outside definitely change
    @test all(r->reduce_to_wignerseitz(r, Rs) ≈ r, .999*vs)
    @test all(r->norm(reduce_to_wignerseitz(r, Rs) - r) > norm(r)/1.1, 1.001*vs)

    # pick a point inside the cell, then shuffle it out by a few lattice vectors; check
    # that it is shuffled back in
    r = SVector(0.12, 0.32, -0.24)
    r′ = r + SVector(3, -4, 8)
    @test reduce_to_wignerseitz(r, Rs) ≈ r
    @test reduce_to_wignerseitz(r′, Rs) ≈ r

    # test a point in the rectangular unit cell that isn't in the Wigner-Seitz cell
    ruc = SVector(-0.49, 0.49, 0.49)
    @test !(reduce_to_wignerseitz(ruc, Rs) ≈ ruc)

    # test utility accessors
    r′′ = convert(Vector{Float64}, r)
    @test reduce_to_wignerseitz(r′′, Rs) ≈ r
    @test_throws DimensionMismatch reduce_to_wignerseitz(r′′[1:2], Rs)
    @test_throws DimensionMismatch reduce_to_wignerseitz(r′′, [Rs[1][1:3], Rs[2][1:3], Rs[3][1:2]])
end