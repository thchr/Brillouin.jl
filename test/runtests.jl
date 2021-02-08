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
    Rs = SVector{3,Float64}.([1,0,0], [-.5,√3/2,0], [0,0,1.25])         # direct basis
    Gs = SVector{3,Float64}.([[2π,0,0], [2π/√3,4π/√3,0], [0,0,2π/1.25]]) # reciprocal basis
    cell = wignerseitz(Gs)
    @test length(cell) == 8
    @test length(vertices(cell)) == 12

    let io = IOBuffer()
        show(io, MIME"text/plain"(), cell)
        show_str = """
        Cell{3} (8 faces, 12 vertices):
         verts: [0.486, 4.291, -2.513]
                [-3.142, 2.964, 2.513]
                [-3.142, 2.964, -2.513]
                [0.486, 4.291, 2.513]
                [3.142, 2.964, -2.513]
                [3.142, -2.964, -2.513]
                [-0.486, -4.291, -2.513]
                [-3.142, -2.964, -2.513]
                [-3.142, -2.964, 2.513]
                [-0.486, -4.291, 2.513]
                [3.142, -2.964, 2.513]
                [3.142, 2.964, 2.513]
         faces: [5, 12, 11, 6]
                [7, 6, 11, 10]
                [8, 3, 1, 5, 6, 7]
                [8, 7, 10, 9]
                [9, 2, 3, 8]
                [10, 11, 12, 4, 2, 9]
                [1, 3, 2, 4]
                [12, 5, 1, 4]
        """
        
        @test String(take!(io)) == show_str
    end
end
