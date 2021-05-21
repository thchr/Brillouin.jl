import Crystalline

@testset "KPaths" begin
    # --- `cumdist` ---
    kvs_2d = [[0,0], [0,1], [1,1], [1,0], [0,0], [1,1], [-1,-2], [2,-1]]
    @test cumdists(kvs_2d) ≈ [0,1,2,3,4,4+√2,4+√2+√13, 4+√2+√13+√10]
    kvs_3d = [[2,1.5,3], [4,5,6], [-1.1, 2, -2.3]]
    @test cumdists(kvs_3d) ≈ [0, √25.25, √25.25+√103.9]

    # --- `KPath` ---
    sgnum = 227
    Rs = [[1,0,0], [0,1,0], [0,0,1]] # conventional direct basis
    kp = irrfbz_path(sgnum, Rs)

    # test that `Rs` has no impact in sgnum 227 (cubic face-centered, FCC)
    @test kp == irrfbz_path(sgnum, nothing)

    # `getindex`
    @test collect(kp) == [:Γ => [0.0, 0.0, 0.0]
                          :X => [0.5, 0.0, 0.5]
                          :U => [0.625, 0.25, 0.625]
                          :K => [0.375, 0.375, 0.75]
                          :Γ => [0.0, 0.0, 0.0]
                          :L => [0.5, 0.5, 0.5]
                          :W => [0.5, 0.25, 0.75]
                          :X => [0.5, 0.0, 0.5]]

    # `paths`
    @test paths(kp) == [[:Γ, :X, :U], [:K, :Γ, :L, :W, :X]]

    # `show(::IO, ::MIME"text/plain", ::KPath)`
    kp_show_reference = """
    KPath{3} (6 points, 2 paths, 8 points in paths):
     points: :U => [0.625, 0.25, 0.625]
             :W => [0.5, 0.25, 0.75]
             :K => [0.375, 0.375, 0.75]
             :Γ => [0.0, 0.0, 0.0]
             :L => [0.5, 0.5, 0.5]
             :X => [0.5, 0.0, 0.5]
      paths: [:Γ, :X, :U]
             [:K, :Γ, :L, :W, :X]"""
    test_show(kp, kp_show_reference)

    # `cartesianize`
    pGs = 2π.*[[-1,1,1], [1,-1,1], [1,1,-1]] # primitive reciprocal basis
    kp′ = Brillouin.KPaths.cartesianize(kp, pGs)
    @test all(zip(points(kp), points(kp′))) do ((klab, kv), (klab′, kv′))
        klab == klab′ && kv'pGs ≈ kv′
    end

    # vendored `bravaistype` vs upstream Crystalline; check that they are in sync
    for D in 1:3, sgnum in 1:Crystalline.MAX_SGNUM[D]
        @test Crystalline.bravaistype(sgnum, D, normalize=false) == Brillouin.CrystallineBravaisVendor.bravaistype(sgnum, D)
    end
    
    # `extended_bravais` and `irrfbz_path`
    for sgnum in 1:230
        Rs = begin
            # the triclinic case (space group 1 & 2) needs a bit of care, since we Brillouin
            # only allows all-obtuse or all-acute basis systems; just hardcode two examples
            if sgnum == 1     # check all-obtuse triclinic `Rs`
                Crystalline.DirectBasis([1.0, 0.0, 0.0], [1.18, 1.54, 0.0], [0.95, 0.44, 1.32])
            elseif sgnum == 2 # check all-acute triclinic `Rs`
                Crystalline.DirectBasis([1.0, 0.0, 0.0], [-0.53, 0.49, 0.0], [0.13, -0.30, 0.74])
            else              # just go with `directbasis` for all other space groups
                Crystalline.directbasis(sgnum, 3)
            end
        end
        bt = Crystalline.bravaistype(sgnum, 3, normalize=false)
        # extended Bravais types
        ebt = Brillouin.KPaths.extended_bravais(sgnum, bt, Rs)
        @test contains(string(ebt), bt)

        # hard to test output of `irrfbz_path` systematically; just test that it returns a
        # `KPath`, and matches `get_points` and `pathsd`
        kp = irrfbz_path(sgnum, Rs)
        @test kp isa KPath
        @test points(kp) == Brillouin.KPaths.get_points(ebt, Rs)
        @test paths(kp)  == Brillouin.KPaths.pathsd[ebt]
    end
    @test_throws DomainError Brillouin.KPaths.extended_bravais(110, "Q", nothing)  # "undefined bravais type"
    @test_throws DomainError Brillouin.KPaths.extended_bravais(194, "cP", nothing) # `_throw_conflicting_sgnum_and_bravais`
    @test_throws DomainError Brillouin.KPaths.extended_bravais(194, "cF", nothing) # `_throw_conflicting_sgnum_and_bravais`
    @test_throws DomainError Brillouin.KPaths.extended_bravais(38, "tI", nothing)  # `_throw_basis_required`
    Rs′ = Crystalline.DirectBasis([1, 0, 0], [0.3, 0.8, 0], [-1.6, 0.8, 0.9]) # neither all-obtuse nor all-acute
    @test_throws DomainError Brillouin.KPaths.extended_bravais(1, "aP", Rs′)       # `_throw_basis_required`

    # --- `KPathInterpolant` ---
    # `interpolate`
    N = 100
    kpi = interpolate(kp, N)
    @test abs(length(kpi) - N) ≤ 1 # allow deviation by 1 point in total (too tight?)
    
    # all points in `kp` must be explicitly contained in `kpi`
    @test all(∈(kpi), values(points(kp)))

    # `splice`
    kps = splice(kp, N)
    highsym_points = sum(path->length(path), paths(kp))
    segments       = sum(path->length(path)-1, paths(kp))
    @test length(kps) == highsym_points + segments*N
end