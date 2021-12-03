using Brillouin, Test
using Brillouin.KPaths: Bravais
using Spglib

@testset "KPath for Spglib cell" begin
    # Trigonal lattice, space group R-3m, 166
    # Example: Bi2Se3 (https://materialsproject.org/materials/mp-541837/)
    a = 1.0
    c = 8.0
    pRs_standard    = [[a*√3/2,  a/2, c/3],
                       [-a*√3/2, a/2, c/3],
                       [0, -a, c/3]]
    pRs_nonstandard = [[a*√3/2, -a/2, c/3],
                       [0, a, c/3],
                       [-a*√3/2, -a/2, c/3]]

    θ = 10 / 180 * pi
    rotation = [[cos(θ) sin(θ) 0]; [-sin(θ) cos(θ) 0]; [0 0 1]]
    pRs_rotated = Ref(rotation) .* pRs_standard

    cell_standard = Spglib.Cell(pRs_standard, [[0, 0, 0]], [0])
    cell_nonstandard = Spglib.Cell(pRs_nonstandard, [[0, 0, 0]], [0])
    cell_rotated = Spglib.Cell(pRs_rotated, [[0, 0, 0]], [0])

    kp_standard = irrfbz_path(cell_standard)
    kp_nonstandard = irrfbz_path(cell_nonstandard)
    kp_rotated = irrfbz_path(cell_rotated)

    @test Bravais.reciprocalbasis(kp_standard.basis) ≈ pRs_standard
    @test Bravais.reciprocalbasis(kp_nonstandard.basis) ≈ pRs_nonstandard
    @test Bravais.reciprocalbasis(kp_rotated.basis) ≈ pRs_rotated

    @test kp_nonstandard.paths == kp_standard.paths
    @test kp_rotated.paths == kp_standard.paths

    for (lab, kv) in kp_standard.points
        @test kp_nonstandard.points[lab] ≈ kv
        @test kp_rotated.points[lab] ≈ kv
    end
end
