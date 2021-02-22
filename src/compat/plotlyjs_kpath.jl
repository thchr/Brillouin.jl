sgnum = 122
cntr = centering(sgnum, 3)
Rs = directbasis(sgnum, 3)
pRs = primitivize(Rs, cntr)
Gs = reciprocalbasis(Rs)
pGs = reciprocalbasis(pRs)
paths_kvs, paths_labs, lab2kv = irrfbz_path(sgnum, 5, Rs, splice=true)
c = wignerseitz(pGs)

tpaths = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(paths_kvs))
for (i,path) in enumerate(paths_kvs)
    kvs = [kv'pGs for kv in path]
    tpaths[i] = PlotlyJS.scatter3d(
        x=getindex.(kvs, 1), y=getindex.(kvs, 2), z=getindex.(kvs, 3),
        mode="lines", line=attr(color="purple", width=8))
end
tpoints = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(lab2kv))
for (i,(lab, kv)) in enumerate(lab2kv)
    kv = kv'pGs
    tpoints[i] =  PlotlyJS.scatter3d(
        x=kv[1:1], y=kv[2:2], z=kv[3:3],
        mode="marker", hovertext=lab, hoverinfo="text",
        marker=attr(color="purple", size=6, line=attr(color="white", width=1)))
end
P = plot(c);
ts = P.plot.data
layout = P.plot.layout
#plot(Plot(tpaths, layout))
plot(Plot(vcat(ts, tpaths, tpoints), layout))