# Paths in the representation domain of the BZ, obtained from the [SeeK] publication. 
# Methods return **k**-vector coordinates in a primitive reciprocal basis ("reciprocal 
# crystallographic primitive cell" in the [SeeK]'s nomenclature).
#
# ^[Seek] http://dx.doi.org/10.1016/j.commatsci.2016.10.015; see also online interface at 
#         https://www.materialscloud.org/work/tools/seekpath

module KPaths

# ---------------------------------------------------------------------------------------- #
export irrfbz_path, KPath, points, paths, cartesianize!, KPathInterpolant
# ---------------------------------------------------------------------------------------- #
using ..CrystallineBravaisVendor: bravaistype
using ..Brillouin: AVec, BasisLike, SHOWDIGITS
using LinearAlgebra: norm, dot, ×
using StaticArrays
using DocStringExtensions

import Base: show, summary, getindex, IndexStyle, size
# ---------------------------------------------------------------------------------------- #

include("codegen_kpoints.jl")  # defines `get_points` and `pathsd` (& other subfunctions)
include("bravais-branches.jl") # defines `extended_bravais`

# ---------------------------------------------------------------------------------------- #

abstract type AbstractPath{T} <: AbstractVector{T} end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct KPath{D} <: AbstractPath{Pair{Symbol, SVector{D, Float64}}}
    points :: Dict{Symbol, SVector{D,Float64}}
    paths  :: Vector{Vector{Symbol}}
end

"""
    points(kp::KPath) -> Dict{Symbol, SVector{D,Float64}}

Return a dictionary of the **k**-points (values) and associated **k**-labels (keys)
referenced in `kp`.
"""
points(kp::KPath) = kp.points

"""
    paths(kp::KPath) -> Dict{Symbol, SVector{D,Float64}}

Return a vector of vectors, with each vector describing a connected path between between
**k**-points referenced in `kp` (see also [`points(::KPath)`](@ref)).
"""
paths(kp::KPath) = kp.paths

Base.@propagate_inbounds function getindex(kp::KPath, i::Int)
    # index into the `i`th point in the "flattened" `paths(kp)`
    @boundscheck i < 0 && throw(BoundsError(kp, i))
    j = 1
    i′ = i′′ = i
    while (i′′ -= length(paths(kp)[j])) > 0
        i′ = i′′
        j += 1
        @boundscheck j > length(paths(kp)) && throw(BoundsError(kp, i))
    end
    @boundscheck i′ > length(paths(kp)[j]) && throw(BoundsError(kp, i))

    lab = @inbounds paths(kp)[j][i′]
    return lab => points(kp)[lab]
end
IndexStyle(::Type{<:KPath}) = IndexLinear()
size(kp::KPath) = (sum(length, paths(kp)),)

function summary(io::IO, kp::KPath{D}) where D
    print(io, "KPath{", D, "} (", length(points(kp)), " points, ", 
              length(paths(kp)), " paths, ", length(kp), " points in paths)")
end
function show(io::IO, ::MIME"text/plain", kp::KPath)
    summary(io, kp)
    println(io, ":")
    # points
    print(io, " points: ")
    foreach(enumerate(points(kp))) do (i,(lab,kv))
        i ≠ 1 && print(io, "         ")
        println(io, ":", lab, " => ", round.(kv, digits=SHOWDIGITS))
    end
    # paths
    print(io, "  paths: ")
    foreach(enumerate(paths(kp))) do (i,path)
        i ≠ 1 && print(io, "         ")
        if i ≠ length(paths(kp))
            println(io, path)
        else
            print(io, path)
        end
    end
end
# ---------------------------------------------------------------------------------------- #

"""
    irrfbz_path(sgnum::Integer, Rs::Union{Nothing, $(BasisLike)}=nothing)   -->  ::KPath

Returns a **k**-path (`::KPath`) in the (primitive) irreducible Brillouin zone that includes
all distinct high-symmetry lines and points as well as parts of the Brillouin zone boundary.

`Rs` refers to the direct basis of the conventional unit cell. For some space groups, it
is needed to disambiguate the "extended" Bravais types that may differ depending on the
lengths of the lattice vectors (because the Brillouin zone may depend on these lengths).
If the requested space group is known to not fall under this case, `Rs` can be supplied
as `nothing` (default).

Note that the returned **k**-points are given in the basis of the **primitive** reciprocal
basis (see [`cartesianize!`](@ref)).

To interpolate the resulting `KPath`, see [`interpolate(::KPath, ::Integer)`](@ref)
and [`splice(::KPath, ::Integer)`](@ref).

## Data and referencing
All data is sourced from the SeeK HPKOT publication: please cite the original work [^1].

All paths currently assume time-reversal symmetry (or, equivalently, inversion symmetry), 
corresponding to the SeeK's `[with inversion]` setting. If neither inversion nor
time-reversal, include the "inverted" -**k** paths as well manually.

[^1] Hinuma, Pizzi, Kumagai, Oba, & Tanaka, *Band structure diagram paths based on
     crystallography*, 
     [Comp. Mat. Sci. **128**, 140 (2017)](http://dx.doi.org/10.1016/j.commatsci.2016.10.015)
"""
function irrfbz_path(sgnum::Integer, Rs::Union{Nothing, AVec{<:AVec{<:Real}}}=nothing)
    if Rs isa AVec
        all(R->length(R)==3, Rs) || throw(DomainError(Rs, "Currently only implemented for 3D space groups"))
        Rs = convert.(SVector{3,Float64}, Rs)
    end
    bt = bravaistype(sgnum, 3)
    ext_bt = extended_bravais(sgnum, bt, Rs)

    # get data about path
    lab2kvs = get_points(ext_bt, Rs)
    paths = pathsd[ext_bt]

    return KPath(lab2kvs, paths)
end

"""
    cartesianize!(kp::KPath, Gs::BasisLike)

Transform a **k**-path `kp` to a Cartesian coordinate system using a primitive basis `Gs`.
Modifies the underlying dictionary in `kp` in-place.
"""
function cartesianize!(kp::KPath{D}, Gs::Union{BasisLike{D}, AVec{<:AVec{<:Real}}}) where D
    for (lab, kv) in points(kp)
        points(kp)[lab] = kv'Gs
    end
    return kp
end
function cartesianize(kp::KPath{D}, Gs::Union{BasisLike{D}, AVec{<:AVec{<:Real}}}) where D
    return cartesianize!(deepcopy(kp), Gs)
end

# ---------------------------------------------------------------------------------------- #

include("interpolate-paths.jl")
export interpolate, splice, cumdists

# ---------------------------------------------------------------------------------------- #
end # module