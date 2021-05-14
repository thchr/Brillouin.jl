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
    str = sprint((io,v) -> show(io, MIME"text/plain"(), v), kp)
    str_reference = """
    KPath{3} (6 points, 2 paths, 8 points in paths):
     points: :U => [0.625, 0.25, 0.625]
             :W => [0.5, 0.25, 0.75]
             :K => [0.375, 0.375, 0.75]
             :Γ => [0.0, 0.0, 0.0]
             :L => [0.5, 0.5, 0.5]
             :X => [0.5, 0.0, 0.5]
      paths: [:Γ, :X, :U]
             [:K, :Γ, :L, :W, :X]"""
    @test str == str_reference

    # `cartesianize`
    pGs = 2π.*[[-1,1,1], [1,-1,1], [1,1,-1]] # primitive reciprocal basis
    kp′ = Brillouin.KPaths.cartesianize(kp, pGs)
    @test all(zip(points(kp), points(kp′))) do ((klab, kv), (klab′, kv′))
        klab == klab′ && kv'pGs ≈ kv′
    end

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