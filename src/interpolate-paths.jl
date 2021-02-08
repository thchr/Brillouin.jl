# Get the cumulative k-path length
function cumdists(kvs::AVec{<:AVec{<:Real}})
    N     = length(kvs)
    D     = length(first(kvs))
    dists = Vector{Float64}(undef, N)
    
    dists[1] = 0.0
    @inbounds for i in 2:N
        s = 0.0
        for j in 1:D # use a loop to avoid allocations
            s += abs2(kvs[i][j] - kvs[i-1][j])
        end
        dists[i] = dists[i-1] + sqrt(s)
    end
    return dists
end

"""
    interp_paths_from_labs(lab2kv, paths_labs, Nk::Integer; 
        [splice::Bool=false, legacy::Bool=false])           --> 

Return interpolated k-paths through vertices in `lab2kv` via labels in `paths_labs`.

## Keyword arguments
- `splice` (default, `false`): if `false` approximately `Nk` sampling points are used 
  *in total* across all paths. If `true`, `Nk` points are used per linear segment of each
  path.
- `legacy` (default, `false`): if `true`, an old semi-broken implementation is used to 
  distribute the frequency of interpolation points. This is only available to support
  old legacy datasets. Do not use.
"""
function interp_paths_from_labs(lab2kv, paths_labs, Nk::Integer;
            splice::Bool=false, 
            legacy::Bool=false)
    paths_verts = map(labs->[lab2kv[lab] for lab in labs], paths_labs)

    if legacy
        # old stupid legacy approach; see __legacy_get_kvpath
        splice && throw("splice and legacy cannot be simultaneously true")
        return __legacy_get_kvpath.(paths_verts, Nk)
    end
    
    if !splice
        # calculate length of each path
        cumdist_per_path = zeros(Float64, length(paths_verts))
        for (p,path_verts) in enumerate(paths_verts)
            for i in Base.OneTo(length(path_verts)-1)
                kᵢ, kᵢ₊₁ = path_verts[i], path_verts[i+1]
                cumdist_per_path[p] += norm(kᵢ₊₁ .- kᵢ)
            end
        end
        cumdist = sum(cumdist_per_path)

        # number of points per path
        Nksᵖ = round.(Int64, (Nk/cumdist).*cumdist_per_path, RoundUp)

        paths_kvs = map(zip(paths_verts, Nksᵖ)) do args # this map-form seems needed for
            interpolate_kvpath(args[1], args[2])[1]     # inference to succeed...
        end                                             # (a generator does e.g. not infer)
        return paths_kvs
        
    else
        return splice_kvpath.(paths_verts, Nk)
    end

end

# This is an old an pretty broken implementation, but we keep it around because some of the
# datasets (e.g. sgs 68, 86, & 230) were created with this implementation, and we need it to
# be able to refresh the kpaths associated with those calculations.
function __legacy_get_kvpath(kvs::AVec{<:AVec{<:Real}}, Ninterp::Integer)
    Nkpairs = length(kvs)-1
    dists = Vector{Float64}(undef, Nkpairs)
    @inbounds for i in Base.OneTo(Nkpairs)
        dists[i] = norm(kvs[i] .- kvs[i+1])
    end
    meandist = sum(dists)/Nkpairs

    kvpath = [float.(kvs[1])]
    @inbounds for i in Base.OneTo(Nkpairs)
        # try to maintain an even distribution of k-points along path
        Ninterp_i = round(Int64, dists[i]./meandist*Ninterp)
        # new k-points
        newkvs = range(kvs[i],kvs[i+1],length=Ninterp_i)
        # append new kvecs to kpath
        append!(kvpath, (@view newkvs[2:end]))
    end
    return kvpath
end