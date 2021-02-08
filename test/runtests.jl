using Brillouin
using Test
using StaticArrays

@testset "KPath" begin
    # --- `cumdist` ---
    kvs_2d = [[0,0], [0,1], [1,1], [1,0], [0,0], [1,1], [-1,-2], [2,-1]]
    @test cumdists(kvs_2d) ≈ [0,1,2,3,4,4+√2,4+√2+√13, 4+√2+√13+√10]
    kvs_3d = [[2,1.5,3], [4,5,6], [-1.1, 2, -2.3]]
    @test cumdists(kvs_3d) ≈ [0, √25.25, √25.25+√103.9]
end


@testset "WignerSeitz" begin
    # --- `wignerseitz` ---
    # hexagonal lattice (example from space group 147)
    Rs = SVector{3,Float64}.([[1,0,0], [-.5,√3/2,0], [0,0,1.25]])        # direct basis
    Gs = SVector{3,Float64}.([[2π,2π/√3,0], [0,4π/√3,0], [0,0,2π/1.25]]) # reciprocal basis
    cell = wignerseitz(Gs)
    @test length(cell) == 8
    @test length(vertices(cell)) == 12

    # for CI and testing generally, we cannot depend on a specific sorting/ordering
    # from `wignerseitz` (and ultimately, `PySpatial.Voronoi` or `PySpatial.ConvexHull`),
    # so to be robust to this, we compare output after sorting to a canonical order
    function canonicalize!(c::WignerSeitz.Cell)
        vs = vertices(c)
        fs = faces(c)
        # sort the vertices in some canonical order (we also use `round` to be robust to 
        # floating point differences across architectures)
        I = sortperm(vs, by=v->round.(v, digits=5))
        permute!(vs, I)
        # relabel the face-indices to match the new vertex sorting
        replace!.(fs, (I .=> eachindex(I))...)
        # sort the faces in some canonical order
        sort!(fs)

        return c
    end
    let io = IOBuffer()
        show(io, MIME"text/plain"(), canonicalize!(cell))
        show_str = """
        Cell{3} (8 faces, 12 vertices):
         verts: [-4.189, 0.0, -2.513]
                [-4.189, 0.0, 2.513]
                [-2.094, -3.628, -2.513]
                [-2.094, -3.628, 2.513]
                [-2.094, 3.628, -2.513]
                [-2.094, 3.628, 2.513]
                [2.094, -3.628, -2.513]
                [2.094, -3.628, 2.513]
                [2.094, 3.628, -2.513]
                [2.094, 3.628, 2.513]
                [4.189, 0.0, -2.513]
                [4.189, 0.0, 2.513]
         faces: [1, 3, 4, 2]
                [2, 6, 5, 1]
                [5, 9, 11, 7, 3, 1]
                [7, 8, 4, 3]
                [7, 11, 12, 8]
                [9, 10, 12, 11]
                [10, 6, 2, 4, 8, 12]
                [10, 9, 5, 6]
        """
        
        @test String(take!(io)) == show_str
    end
end
