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

    @test Bravais.reciprocalbasis(basis(kp_standard)) ≈ pRs_standard
    @test Bravais.reciprocalbasis(basis(kp_nonstandard)) ≈ pRs_nonstandard
    @test Bravais.reciprocalbasis(basis(kp_rotated)) ≈ pRs_rotated

    @test paths(kp_nonstandard) == paths(kp_standard)
    @test paths(kp_rotated) == paths(kp_standard)

    for (lab, kv) in points(kp_standard)
        @test points(kp_nonstandard)[lab] ≈ kv
        @test points(kp_rotated)[lab] ≈ kv
    end

    # The case where the input cell is a supercell of a smaller primitive cell
    pRs_super = [pRs_standard[1], pRs_standard[2], 2 .* pRs_standard[3]]
    cell_super = Spglib.Cell(pRs_super, [[0, 0, 0], [0, 0, 0.5]], [0, 0])
    kp_super = irrfbz_path(cell_super)

    @test Bravais.reciprocalbasis(basis(kp_super)) ≈ pRs_super
    @test paths(kp_super) == paths(kp_standard)

    # Check that the Cartesian coordinates of the high-symmetry k points are identical
    kp_standard_cart = cartesianize(kp_standard)
    kp_super_cart = cartesianize(kp_super)
    for (lab, kv) in points(kp_standard_cart)
        @test points(kp_super_cart)[lab] ≈ kv
    end
end
