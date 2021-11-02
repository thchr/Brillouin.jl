using Brillouin, PlotlyJS, Test
using Brillouin.KPaths: Bravais
using StaticArrays

# NB: For the plot functionality, there's not that many meaningful things to test against
#     in the output since it's a graphical thing. So, instead, we just check that things 
#     run without erroring by checking that the return type is as expected

@testset "PlotlyJS" begin
    # ------------------------------------------------------------------------------------ #
    # Wigner-Seitz cells and k-space path visualization
    # 3D
    sgnum = 147
    Rs = convert(SVector{3,SVector{3,Float64}}, [[1.0,0,0], [-.5,√3/2,0], [0,0,1.25]])
    Gs = Bravais.reciprocalbasis(Rs)

    cell = wignerseitz(Gs)
    kp = irrfbz_path(sgnum, Rs)
    @test plot(cell) isa PlotlyJS.SyncPlot
    @test plot(kp) isa PlotlyJS.SyncPlot
    @test plot(cell, kp) isa PlotlyJS.SyncPlot

    # 2D
    sgnum = 5
    cntr = Bravais.centering(sgnum, 2) # c-type centered
    Rs = convert(SVector{2,SVector{2,Float64}}, [[1.0,0], [0,.8]])
    pGs = Bravais.primitivize(Bravais.reciprocalbasis(Rs), cntr)

    cell = wignerseitz(pGs)
    kp = irrfbz_path(sgnum, Rs)
    @test plot(cell) isa PlotlyJS.SyncPlot
    @test plot(kp) isa PlotlyJS.SyncPlot
    @test plot(cell, kp) isa PlotlyJS.SyncPlot

    # ------------------------------------------------------------------------------------ #
    # Band structure plots (example from docs)
    sgnum = 227
    Rs = [[1,0,0], [0,1,0], [0,0,1]]
    kp = irrfbz_path(sgnum, Rs)
    pGs = basis(kp)
    kpi = interpolate(kp, 100)

    function ϵ(k; γ::Real=1.0)
        kx = 2π*(-k[1]+k[2]+k[3])
        ky = 2π*(+k[1]-k[2]+k[3])
        kz = 2π*(+k[1]+k[2]-k[3])
        return 4γ * (cos(kx/2)*cos(ky/2) + cos(ky/2)*cos(kz/2) + cos(kz/2)*cos(kx/2))
    end
    
    band1 = ϵ.(kpi)
    band2 = 20 .- (1/2).*band1

    # basic functionality
    @test plot(kpi, [band1]) isa PlotlyJS.SyncPlot
    @test plot(kpi, reduce(hcat, [band1, band2])) isa PlotlyJS.SyncPlot # matrix input

    # annotations
    offset = 2maximum(band2)
    band3 = band1 .+ offset
    annotations = Dict(:X => [1:2 => "X, bands 1:2", 30.0 => "X, some frequency"], 
                       :Γ => [1 => "Γ, band 1", 3 => "Γ, band 3"])
    @test plot(kpi, [band1, band2, band3]; annotations) isa PlotlyJS.SyncPlot

    # band highlights
    band_highlights = Dict(2:3 => attr(color=:black, width=3))
    @test plot(kpi, [band1, band2, band3]; band_highlights) isa PlotlyJS.SyncPlot

    # title/y-label
    @test plot(kpi, [band1]; title="hello world", ylabel="ωa/2πc") isa PlotlyJS.SyncPlot
end