using Brillouin, Test
using Brillouin.KPaths: Bravais
using StaticArrays
using Spglib

@testset "KPath for Spglib cell" begin
    # Trigonal lattice, space group R-3m, 166
    # Example: Bi2Se3 (https://materialsproject.org/materials/mp-541837/)
    a = 1.0
    c = 8.0
    lattice_nonstandard = SVector(SVector(a*sqrt(3)/2, -a/2, c/3),
                                  SVector(0.0, a, c/3),
                                  SVector(-a*sqrt(3)/2, -a/2, c/3))
    conv_lattice = SVector(SVector(a*sqrt(3), 0, 0),
                           SVector(-a*sqrt(3)/2, a*3/2, 0),
                           SVector(0, 0, c))
    sgnum = 166
    lattice_standard = SVector(Bravais.primitivize(Bravais.DirectBasis(collect(conv_lattice)), Bravais.centering(sgnum, 3)))

    θ = 10 / 180 * pi
    rotation = [[cos(θ) sin(θ) 0]; [-sin(θ) cos(θ) 0]; [0 0 1]]
    lattice_rotated = Ref(rotation) .* lattice_standard

    cell_standard = Spglib.Cell(lattice_standard, [[0, 0, 0]], [0])
    cell_nonstandard = Spglib.Cell(lattice_nonstandard, [[0, 0, 0]], [0])
    cell_rotated = Spglib.Cell(lattice_rotated, [[0, 0, 0]], [0])

    kp_standard = irrfbz_path_for_cell(cell_standard)
    kp_nonstandard = irrfbz_path_for_cell(cell_nonstandard)
    kp_rotated = irrfbz_path_for_cell(cell_rotated)

    @test Bravais.reciprocalbasis(kp_standard.basis) ≈ lattice_standard
    @test Bravais.reciprocalbasis(kp_nonstandard.basis) ≈ lattice_nonstandard
    @test Bravais.reciprocalbasis(kp_rotated.basis) ≈ lattice_rotated

    @test kp_nonstandard.paths == kp_standard.paths
    @test kp_rotated.paths == kp_standard.paths

    for (lab, kv) in kp_standard.points
        @test kp_nonstandard.points[lab] ≈ kv
        @test kp_rotated.points[lab] ≈ kv
    end
end
