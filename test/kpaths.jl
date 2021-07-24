import Crystalline
using Brillouin
using Brillouin.KPaths: cartesianize, latticize

@testset "KPath & KPathInterpolant" begin
    # --- `cumdist` ---
    kvs_2d = [[0,0], [0,1], [1,1], [1,0], [0,0], [1,1], [-1,-2], [2,-1]]
    @test cumdists(kvs_2d) ≈ [0,1,2,3,4,4+√2,4+√2+√13, 4+√2+√13+√10]
    kvs_3d = [[2,1.5,3], [4,5,6], [-1.1, 2, -2.3]]
    @test cumdists(kvs_3d) ≈ [0, √25.25, √25.25+√103.9]

    # --- `KPath` ---
    sgnum = 227 # face-centered cubic (cF)
    Rsᶜᶸᵇⁱᶜ = [[1,0,0], [0,1,0], [0,0,1]] # conventional direct basis
    kp = irrfbz_path(sgnum, Rsᶜᶸᵇⁱᶜ)

    # test that `Rs` has no impact in sgnum 227 (cubic face-centered, FCC)
    @test kp == irrfbz_path(sgnum, [rand(3) for _ in 1:3])

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
             [:K, :Γ, :L, :W, :X]
      basis: [-6.283185, 6.283185, 6.283185]
             [6.283185, -6.283185, 6.283185]
             [6.283185, 6.283185, -6.283185]"""
    test_show(kp, kp_show_reference)

    # `cartesianize` & `latticize`
    pGs = 2π.*[[-1,1,1], [1,-1,1], [1,1,-1]] # primitive reciprocal basis
    kp′ = Brillouin.cartesianize(kp)
    @test all(zip(points(kp), points(kp′))) do ((klab, kv), (klab′, kv′))
        klab == klab′ && kv'pGs ≈ kv′
    end
    kp′′ = Brillouin.latticize(kp′)
    @test keys(points(kp′′)) == keys(points(kp)) # point labels agree
    @test paths(kp′′)        == paths(kp)        # segments agree
    @test all(values(points(kp′′)) .≈ values(points(kp))) # point coordinates approx. agree

    # vendored `bravaistype` vs upstream Crystalline; check that they are in sync
    for D in 1:3, sgnum in 1:Crystalline.MAX_SGNUM[D]
        @test Crystalline.bravaistype(sgnum, D, normalize=false) == Brillouin.CrystallineBravaisVendor.bravaistype(sgnum, D)
    end
    
    # `extended_bravais` and `irrfbz_path`
    sgnums = (1:2, 1:17, 1:230)
    for D in (1, 2, 3)
        Dᵛ = Val(D)
        for sgnum in sgnums[D]
            Rs = begin
                # the 3D triclinic case (space group 1 & 2) needs more care, since Brillouin
                # only allows all-obtuse or all-acute bases; just hardcode two examples
                if sgnum == 1 && D == 3    # check all-obtuse triclinic `Rs`
                    Crystalline.DirectBasis([1.0, 0.0, 0.0], [1.18, 1.54, 0.0], [0.95, 0.44, 1.32])
                elseif sgnum == 2 && D == 3 # check all-acute triclinic `Rs`
                    Crystalline.DirectBasis([1.0, 0.0, 0.0], [-0.53, 0.49, 0.0], [0.13, -0.30, 0.74])
                else                        # go with `directbasis` for everything else
                    Crystalline.directbasis(sgnum, Dᵛ)
                end
            end
            bt = Crystalline.bravaistype(sgnum, D, normalize=false)
            # extended Bravais types
            ebt = Brillouin.KPaths.extended_bravais(sgnum, bt, Rs, Dᵛ)
            @test contains(string(ebt), bt)

            # hard to test output of `irrfbz_path` systematically; just test that it returns
            # a `KPath`, and matches `get_points` and `get_paths`
            kp = irrfbz_path(sgnum, Rs)
            @test kp isa KPath{D}
            @test points(kp) == Brillouin.KPaths.get_points(ebt, Rs, Dᵛ)
            @test paths(kp)  == Brillouin.KPaths.get_paths(ebt, Dᵛ)
        end
    end
    @test_throws DomainError Brillouin.KPaths.extended_bravais(110, "Q", nothing, Val(3))  # "undefined bravais type"
    @test_throws DomainError Brillouin.KPaths.extended_bravais(194, "cP", nothing, Val(3)) # `_throw_conflicting_sgnum_and_bravais`
    @test_throws DomainError Brillouin.KPaths.extended_bravais(194, "cF", nothing, Val(3)) # `_throw_conflicting_sgnum_and_bravais`
    @test_throws DomainError Brillouin.KPaths.extended_bravais(38, "tI", nothing, Val(3))  # `_throw_basis_required`
    Rs′ = Crystalline.DirectBasis([1, 0, 0], [0.3, 0.8, 0], [-1.6, 0.8, 0.9]) # neither all-obtuse nor all-acute
    @test_throws DomainError Brillouin.KPaths.extended_bravais(1, "aP", Rs′, Val(3))       # `_throw_basis_required`


    # --- `KPathInterpolant` ---
    # `interpolate`
    N = 100
    kpi = interpolate(kp, N)
    @test abs(length(kpi) - N) ≤ 1 # allow deviation by 1 point in total (too tight?)
    
    # all points in `kp` must be explicitly contained in `kpi`
    @test all(∈(kpi), values(points(kp)))
    @test kpi == interpolate(kp, length=N)

    # `splice`
    kps = splice(kp, N)
    highsym_points = sum(path->length(path), paths(kp))
    segments       = sum(path->length(path)-1, paths(kp))
    @test length(kps) == highsym_points + segments*N

    # `cartesianize` 
    # test that `cartesianize` commutes with `interpolate` when applied to a `kp`/`kpi`
    cntr = Crystalline.centering(5, 2)         # centering type 'c' in plane group 5
    Rs   = Crystalline.directbasis(5, Val(2))
    Gs   = Crystalline.reciprocalbasis(Rs)
    pGs  = Crystalline.primitivize(Gs, cntr)

    kp   = irrfbz_path(5, Rs, Val(2))
    kpi  = interpolate(Brillouin.cartesianize(kp), 100)
    kpi′ = Brillouin.cartesianize!(interpolate(kp, 100))

    @test typeof(kpi) === typeof(kpi′) === KPathInterpolant{2}
    @test kpi′ ≈ kpi # test that `interpolate` commutes with `latticize`/`cartesianize`
    @test Brillouin.cartesianize!(latticize!(kpi)) ≈ kpi

    # interpolate via `density` kwarg
    kp = irrfbz_path(227, Rsᶜᶸᵇⁱᶜ)
    for ρ in [1, 10, 20, 21, 100, 3.5]
        kpi_ρ = interpolate(kp, density=ρ)
        
        # test a lower bound on how many points must be in kpi
        segments = sum(p -> length(p) - 1, kp.paths)
        tot_dist = sum(last.(cumdists.(Brillouin.cartesianize(kpi_ρ).kpaths)))
        min_points = tot_dist*ρ - segments
        @test length(kpi_ρ) > min_points

        # TODO: test that the density of k-points is _at least_ equal to ρ at every segment
    end

    # interpolate with wrong kwargs
    @test_throws ArgumentError interpolate(kp)
    @test_throws ArgumentError interpolate(kp; length=10, density=2)
end # testset