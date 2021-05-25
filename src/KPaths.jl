# Paths in the representation domain of the BZ, obtained from the [SeeK] publication. 
# Methods return **k**-vector coordinates in a primitive reciprocal basis ("reciprocal 
# crystallographic primitive cell" in the [SeeK]'s nomenclature).
#
# ^[Seek] http://dx.doi.org/10.1016/j.commatsci.2016.10.015; see also online interface at 
#         https://www.materialscloud.org/work/tools/seekpath

module KPaths

# ---------------------------------------------------------------------------------------- #
export irrfbz_path, KPath, points, paths, cartesianize!, latticize!, KPathInterpolant
# ---------------------------------------------------------------------------------------- #
using ..CrystallineBravaisVendor: bravaistype, boundscheck_sgnum
using ..Brillouin: AVec, BasisLike, SHOWDIGITS
using LinearAlgebra: norm, dot, ×
using StaticArrays
using DocStringExtensions

import Base: show, summary, getindex, setindex!, IndexStyle, size
# ---------------------------------------------------------------------------------------- #

include("codegen-kpoints.jl")  # defines `get_points_3d`, `pathsd_3d` (& other subfunctions)
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
    # TODO: basis :: SVector{D, SVector{D, Float64}}
end

"""
    points(kp::KPath{D}) -> Dict{Symbol, SVector{D,Float64}}

Return a dictionary of the **k**-points (values) and associated **k**-labels (keys)
referenced in `kp`.
"""
points(kp::KPath) = kp.points

"""
    paths(kp::KPath) -> Vector{Vector{Symbol}}

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
    irrfbz_path(sgnum::Integer,
                Rs::Union{Nothing, AVec{<:AVec{<:Real}}}=nothing,
                [::Union{Val(D), Integer},]=Val(3))                     -->  ::KPath{D}

Returns a **k**-path (`::KPath`) in the (primitive) irreducible Brillouin zone that includes
all distinct high-symmetry lines and points as well as relevant parts of the Brillouin zone
boundary.

The dimension `D` (1, 2, or 3) is specified as the third input argument, preferably as a
static `Val{D}` type parameter (or, type-unstably, as an `<:Integer`). Defaults to `Val(3)`.

`Rs` refers to the direct basis of the conventional unit cell. For some space groups, it
is needed to disambiguate the "extended" Bravais types that may differ depending on the
lengths of the lattice vectors (because the Brillouin zone may depend on these lengths).
If the requested space group is known to not fall under this case, `Rs` can be supplied
as `nothing` (default).
If `Rs` is a subtype of a [`StaticVector`](@ref), the dimension can be inferred from its
(static) size; no explicit input of dimension is needed in this case.

## Notes
- The returned **k**-points are given in the basis of the **primitive** reciprocal basis
  (see [`cartesianize!`](@ref)).
- To interpolate a `KPath`, see [`interpolate(::KPath, ::Integer)`](@ref) and
  [`splice(::KPath, ::Integer)`](@ref).
- All paths currently assume time-reversal symmetry (or, equivalently, inversion symmetry).
  If neither are present, include the "inverted" -**k** paths manually.

## Data and referencing
3D paths are sourced from the SeeK HPKOT publication: please cite the original work [^1].

[^1] Hinuma, Pizzi, Kumagai, Oba, & Tanaka, *Band structure diagram paths based on
     crystallography*, 
     [Comp. Mat. Sci. **128**, 140 (2017)](http://dx.doi.org/10.1016/j.commatsci.2016.10.015)
"""
function irrfbz_path(sgnum::Integer,
                     Rs::Union{Nothing, AVec{<:AVec{<:Real}}}=nothing,
                     Dᵛ::Val{D}=Val(3)) where D
                     
    if !isnothing(Rs)
        D′ = length(Rs)
        if D′ ≠ D || any(R -> length(R) ≠ D, Rs)
            throw(DimensionMismatch("inconsistent input dimensions in `Rs` or `D`"))
        end
        Rs = convert(SVector{D, SVector{D, Float64}}, Rs)
    end
    bt = bravaistype(sgnum, D)
    ext_bt = extended_bravais(sgnum, bt, Rs, Dᵛ)

    # get data about path
    lab2kvs = get_points(ext_bt, Rs, Dᵛ)
    paths = get_paths(ext_bt, Dᵛ)

    return KPath(lab2kvs, paths)
end
irrfbz_path(sgnum::Integer, Rs::StaticVector{D, <:AVec{<:Real}}) where D          = irrfbz_path(sgnum, Rs, Val(D))
irrfbz_path(sgnum::Integer, Rs::Union{Nothing, AVec{<:AVec{<:Real}}}, D::Integer) = irrfbz_path(sgnum, Rs, Val(D))

function get_points(ext_bt::Symbol,
                    Rs::Union{Nothing, AVec{<:AVec{<:Real}}},
                    Dᵛ::Val{D}) where D
    if Dᵛ === Val(3)
        return get_points_3d(ext_bt, Rs)
    elseif Dᵛ === Val(2)
        return get_points_2d(ext_bt, Rs)
    elseif Dᵛ === Val(1)
        ext_bt === :lp || throw(DomainError(ext_bt, "invalid extended Bravais type"))
        return Dict(:Γ => SVector(0.0), :X => SVector(0.5))
    else
        throw(DomainError(D, "unsupported input dimension"))
    end
end

function get_paths(ext_bt::Symbol, Dᵛ::Val{D}) where D
    if Dᵛ === Val(3)
        paths = get(pathsd_3d, ext_bt, nothing)
    elseif Dᵛ === Val(2)
        paths = get(pathsd_2d, ext_bt, nothing)
    elseif Dᵛ === Val(1)
        paths = ext_bt === :lp ? [[:Γ, :X]] : nothing
    else
        throw(DomainError(D, "unsupported input dimension"))
    end

    paths === nothing && throw(DomainError(ext_bt, "invalid extended Bravais type"))
    return paths
end

# ---------------------------------------------------------------------------------------- #

include("interpolate-paths.jl")
export interpolate, splice, cumdists

# ---------------------------------------------------------------------------------------- #

"""
    cartesianize!(kp::KPath{D}, Gs::BasisLike{D})

Transform a **k**-path `kp` in a lattice basis to a Cartesian basis with (primitive)
basis vectors `Gs`.
Modifies `kp` in-place.
"""
function cartesianize!(kp::KPath{D}, Gs::Union{BasisLike{D}, AVec{<:AVec{<:Real}}}) where D
    for (lab, kv) in points(kp)
        points(kp)[lab] = kv'Gs
    end
    return kp
end
cartesianize(kp::Union{KPath{D}, KPathInterpolant{D}},
             Gs::Union{BasisLike{D}, AVec{<:AVec{<:Real}}}) where D = cartesianize!(deepcopy(kp), Gs)

"""
    latticize!(kp::KPath{D}, Gs::BasisLike{D})

Transform a **k**-path `kpi` in a Cartesian basis to a lattice basis with (primitive)
reciprocal lattice vectors `Gs`.
Modifies `kp` in-place.
"""
function latticize!(kp::KPath{D}, Gs::Union{BasisLike{D}, AVec{<:AVec{<:Real}}}) where D
    Gm = hcat(Gs...)
    for (lab, kv) in points(kp)
        points(kp)[lab] = Gm\kv
    end
    return kp
end
latticize(kp::Union{KPath{D}, KPathInterpolant{D}},
          Gs::Union{BasisLike{D}, AVec{<:AVec{<:Real}}}) where D = latticize!(deepcopy(kp), Gs)

# ---------------------------------------------------------------------------------------- #
end # module