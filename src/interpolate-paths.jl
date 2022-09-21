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

# ---------------------------------------------------------------------------------------- #
# INTERPOLATIONS OF ::KPath

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct KPathInterpolant{D} <: AbstractPath{SVector{D, Float64}}
    kpaths  :: Vector{Vector{SVector{D, Float64}}}
    labels  :: Vector{Dict{Int, Symbol}}
    basis   :: ReciprocalBasis{D}
    setting :: Base.RefValue{BasisEnum}
end

size(kpi::KPathInterpolant) = (sum(length, kpi.kpaths),)
Base.@propagate_inbounds function composite_index(kpi::KPathInterpolant, i::Int)
    # index into the `i`th point in the "flattened" `kpi.kpaths`
    @boundscheck i < 0 && throw(BoundsError(kpi, i))
    j = 1
    i′ = i′′ = i
    while (i′′ -= length(kpi.kpaths[j])) > 0
        i′ = i′′
        j += 1
        @boundscheck j > length(kpi.kpaths) && throw(BoundsError(kpi, i))
    end
    @boundscheck i′ > length(kpi.kpaths[j]) && throw(BoundsError(kpi, i))
    return j, i′
end
Base.@propagate_inbounds function getindex(kpi::KPathInterpolant, i::Int)
    j, i′ = composite_index(kpi, i)
    return @inbounds kpi.kpaths[j][i′]
end
IndexStyle(::Type{<:KPathInterpolant}) = IndexLinear()
Base.@propagate_inbounds function setindex!(kpi::KPathInterpolant{D}, 
                                            val::SVector{D, Float64}, i::Int) where D
    j, i′ = composite_index(kpi, i)
    @inbounds kpi.kpaths[j][i′] = val
    return kpi
end

"""
    cartesianize!(kpi::KPathInterpolant)

In-place transform an interpolated **k**-path `kpi` in a lattice basis to a Cartesian basis
with (primitive) reciprocal lattice vectors `basis(kpi)`.
"""
function cartesianize!(kpi::KPathInterpolant)
    setting(kpi) === CARTESIAN && return kpi
    for i in eachindex(kpi)
        @inbounds kpi[i] = cartesianize(kpi[i], kpi.basis)
    end
    set_setting!(kpi, CARTESIAN)
    return kpi
end

"""
    latticize!(kpi::KPathInterpolant{D})

In-place transform an interpolated **k**-path `kpi` in a Cartesian basis to a lattice basis
with (primitive) reciprocal lattice vectors `basis(kpi)`.
"""
function latticize!(kpi::KPathInterpolant)
    setting(kpi) === LATTICE && return kpi
    basismatrix = reduce(hcat, basis(kpi))
    for i in eachindex(kpi)
        @inbounds kpi[i] = latticize(kpi[i], basismatrix)
    end
    set_setting!(kpi, LATTICE)
    return kpi
end

"""
    latticize(kpi::KPathInterpolant{D}, basis::AbstractVector{<:AbstractVector{<:Real}})

Transform an interpolated **k**-path `kpi` in a Cartesian basis to a lattice basis with
(primitive) reciprocal lattice vectors `basis`.

If `kpi` is not in a Cartesian basis (i.e., if `setting(kpi) == LATTICE`), `kpi` is returned
as-is.
"""
function latticize(kpi::KPathInterpolant{D}, basis::AVec{<:AVec{<:Real}}) where D
    setting(kpi) === LATTICE && return kpi
    basismatrix = convert(SMatrix{D,D,Float64,D*D}, reduce(hcat, basis))
    kpi′ = typeof(kpi)([copy(kpathsⱼ) for kpathsⱼ in kpi.kpaths], # ⇐ mutate this below
                       kpi.labels, basis, Ref(LATTICE))
    for i in eachindex(kpi′)
        @inbounds kpi′[i] = latticize(kpi′[i], basismatrix)
    end
    return kpi′
end

"""
    interpolate(kp::KPath, N::Integer)
    interpolate(kp::KPath; N::Integer, density::Real) --> KPathInterpolant

Return an interpolant of `kp` with `N` points distributed approximately equidistantly
across the full **k**-path (equidistance is measured in a Cartesian metric).

Note that the interpolant may contain slightly fewer or more points than `N` (typically
fewer) in order to improve equidistance. `N` can also be provided as a keyword argument.

As an alternative to specifying the desired total number of interpolate points via `N`,
a desired density per unit (reciprocal) length can be specified via the keyword argument
`density`.
"""
function interpolate(kp::KPath{D}, N::Integer) where D
    kpᶜ = setting(kp) === CARTESIAN ? kp : cartesianize(kp)
    distss = map(paths(kpᶜ)) do path
        map(1:length(path)-1) do i
            norm(points(kpᶜ)[path[i]] - points(kpᶜ)[path[i+1]])
        end
    end
    totaldist = sum(dists->sum(dists), distss)

    kipaths = [Vector{SVector{D, Float64}}() for _ in 1:length(paths(kp))]
    labels = [Dict{Int, Symbol}() for _ in 1:length(paths(kp))]
    for (j, (path, dists)) in enumerate(zip(paths(kp), distss))
        push!(labels[j], 1 => first(path))
        for i in 1:length(path)-1
            Nᵢ = max(convert(Int, div(N*dists[i], totaldist, RoundUp)), 2)
            append!(kipaths[j], range(points(kp)[path[i]], points(kp)[path[i+1]];
                                      length=Nᵢ)[1:end-1])
            push!(labels[j], length(kipaths[j])+1 => path[i+1])
        end
        push!(kipaths[j], points(kp)[last(path)])
    end

    return KPathInterpolant(kipaths, labels, basis(kp), Ref(setting(kp)))
end

function interpolate(kp::KPath;
            length::Union{Integer, Nothing}=nothing,
            density::Union{Real, Nothing}=nothing)
    if length === nothing && density === nothing
        throw(ArgumentError("must set either `length` or `density` keyword arguments explicitly or use positional argument for `length`"))
    elseif length !== nothing && density !== nothing
        throw(ArgumentError("cannot define `length` and `density` simultaneously"))
    elseif length !== nothing && density === nothing
        return interpolate(kp, length)
    elseif length === nothing && density !== nothing
        return interpolate_via_density(kp, density)
    end
    
    error("unreachable reached; please file an issue")
end

function interpolate_via_density(kp::KPath{D}, density::Real) where D
    kpᶜ = setting(kp) === CARTESIAN ? kp : cartesianize(kp)

    kipaths = [Vector{SVector{D, Float64}}() for _ in 1:length(paths(kp))]
    labels = [Dict{Int, Symbol}() for _ in 1:length(paths(kp))]
    for (j, path) in enumerate(paths(kpᶜ)) # over connected paths
        push!(labels[j], 1 => first(path))
        for i in 1:length(path)-1 # over linear segments
            klabᵢ   = path[i]
            klabᵢ₊₁ = path[i+1]

            # compute required number of points for current segment
            kᶜᵢ   = points(kpᶜ)[klabᵢ]   # start
            kᶜᵢ₊₁ = points(kpᶜ)[klabᵢ₊₁] # end
            distᵢ = norm(kᶜᵢ₊₁ - kᶜᵢ)
            Nᵢ    = max(ceil(Int, distᵢ*density), 2)

            # compute interpolation of current segment (excluding end point)
            kᵢ    = points(kp)[klabᵢ]
            kᵢ₊₁  = points(kp)[klabᵢ₊₁]
            ks    = range(kᵢ, kᵢ₊₁, length=Nᵢ)[1:end-1]

            append!(kipaths[j], ks)
            push!(labels[j], length(kipaths[j])+1 => klabᵢ₊₁)
        end
        push!(kipaths[j], points(kp)[last(path)]) # add previously omitted end point
    end

    return KPathInterpolant(kipaths, labels, basis(kp), Ref(setting(kp)))
end


"""
    splice(kp::KPath, N::Integer) --> KPathInterpolant

Return an interpolant of `kp` with `N` points inserted into each **k**-path segment of `kp`.
"""
function splice(kp::KPath{D}, N::Integer) where D
    kipaths = [Vector{SVector{D, Float64}}() for _ in 1:length(paths(kp))]
    labels = [Dict{Int, Symbol}() for _ in 1:length(paths(kp))]
    for (j, path) in enumerate(paths(kp))
        push!(labels[j], 1 => first(path))
        for i in 1:length(path)-1
            append!(kipaths[j], (@view range(points(kp)[path[i]], 
                                             points(kp)[path[i+1]], 
                                             length=N+2)[1:end-1]))
            push!(labels[j], length(kipaths[j])+1 => path[i+1])
        end
        push!(kipaths[j], points(kp)[last(path)])
    end

    return KPathInterpolant(kipaths, labels, kp.basis, Ref(setting(kp)))
end

# ---------------------------------------------------------------------------------------- #
# INTERPOLATIONS OF ::AVec{<:AVec{<:Real}}

# TODO: Merge the implementations of `interpolate` and `splice` for `AVec{<:AVec{<:Real}}`
#       vs. `KPath`.

"""
    interpolate(kvs::AbstractVector{<:AbstractVector{<:Real}}, N::Integer)
                                                    --> Vector{Vector{<:Real}}

Return an interpolated **k**-path between discrete **k**-points in `kvs`, with approximately
`N` interpolation points in total (typically fewer).

Note that, in general, it is not possible to do this so that all interpolated **k**-points
are equidistant; samples are however *exactly* equidistant across each linear segment
defined by points in `kvs` and *approximately* equidistant across all segments.

See also [`interpolate(::KPath, ::Integer)`](@ref) and [`splice`](@ref).

!!! warning "Future deprecation"
    This method signature is likely to be deprecated in future versions of Brillouin.jl.
"""
function interpolate(kvs::AVec{<:AVec{<:Real}}, N::Integer)
    Nkpairs = length(kvs)-1
    elT     = typeof(first(first(kvs))/1)
    dists   = Vector{elT}(undef, Nkpairs)
    @inbounds for i in 1:Nkpairs
        dists[i] = norm(kvs[i] .- kvs[i+1])
    end
    totaldist  = sum(dists)
    N_per_dist = N/totaldist

    # TODO: preallocate
    kvpath = [float.(kvs[1])]
    @inbounds for i in 1:Nkpairs
        # try to maintain an even distribution of k-points along path
        Nᵢ = round(Int64, dists[i]*N_per_dist, RoundUp) # points in current segment
        Nᵢ = max(Nᵢ, 2) # must be at least two points
        new_kvs = range(kvs[i],kvs[i+1],length=Nᵢ)
        append!(kvpath, (@view new_kvs[2:end]))         # append `new_kvs` to `kvpath`
    end
    return kvpath
end

"""
    splice(kvs::AbstractVector{<:AbstractVector{<:Real}}, N::Integer)
                                                    --> Vector{Vector{<:Real}}

Return an interpolated **k**-path between the discrete **k**-points in `kvs`, with `N`
interpolation points inserted in each segment defined by pairs of adjacent **k**-points.

See also [`splice(::KPath, ::Integer)`](@ref) and [`interpolate`](@ref).

!!! warning "Future deprecation"
    This method signature is likely to be deprecated in future versions of Brillouin.jl.
"""
function splice(kvs::AVec{<:AVec{<:Real}}, N::Integer)
    D       = length(first(kvs))
    Nkpairs = length(kvs)-1
    N⁺²     = N+2

    elT    = typeof(first(first(kvs))/1)
    kvpath = [Vector{elT}(undef, D) for _ in Base.OneTo(Nkpairs+1 + Nkpairs*N)]
    kvpath[1] = kvs[1]
    start, stop = 2, N⁺²
    @inbounds for i in 1:Nkpairs
        new_kvs = range(kvs[i],kvs[i+1],length=N⁺²)
        @views kvpath[start:stop] .= new_kvs[2:end] # insert `new_kvs` in `kvpath`
        start  = stop+1
        stop  += N⁺²-1
    end
    return kvpath
end