# Paths in the representation domain of the BZ, obtained from the [SeeK] publication. 
# Methods return **k**-vector coordinates in a primitive reciprocal basis ("reciprocal 
# crystallographic primitive cell" in the [SeeK]'s nomenclature).
#
# ^[Seek] http://dx.doi.org/10.1016/j.commatsci.2016.10.015; see also online interface at 
#         https://www.materialscloud.org/work/tools/seekpath

module KPaths

# ---------------------------------------------------------------------------------------- #

using Bravais:
    Bravais, # so it is possible to call Brillouin.Bravais from the "outside"
    bravaistype,
    boundscheck_sgnum,
    reciprocalbasis,
    primitivize,
    transform,
    ReciprocalBasis,
    DirectBasis,
    primitivebasismatrix
import Bravais:
    cartesianize,
    cartesianize!,
    latticize,
    latticize!
using ..Brillouin:
    AVec,
    BasisLike,
    SHOWDIGITS,
    BasisEnum,
    CARTESIAN,
    LATTICE,
    setting,
    set_setting!
import ..Brillouin:
    basis
using LinearAlgebra:
    norm,
    dot,
    ×
using StaticArrays
using DocStringExtensions

import Base: show, summary, getindex, setindex!, IndexStyle, size

# ---------------------------------------------------------------------------------------- #

export irrfbz_path, KPath, points, paths, KPathInterpolant
export cartesianize, cartesianize!, latticize, latticize!

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
    points  :: Dict{Symbol, SVector{D,Float64}}
    # TODO: Make values of `points` of type `ReciprocalPoint{D}`?
    paths   :: Vector{Vector{Symbol}}
    basis   :: ReciprocalBasis{D}
    setting :: Base.RefValue{BasisEnum}
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
        println(io, path)
    end
    # basis
    print(io, "  basis: ")
    foreach(enumerate(basis(kp))) do (i,v)
        i ≠ 1 && print(io, "         ")
        if i ≠ length(basis(kp))
            println(io, round.(v, digits=SHOWDIGITS))
        else
            print(io, round.(v, digits=SHOWDIGITS))
        end
    end
end
# ---------------------------------------------------------------------------------------- #

"""
    irrfbz_path(sgnum::Integer, Rs, [::Union{Val(D), Integer},]=Val(3))  -->  ::KPath{D}

Returns a **k**-path (`::KPath`) in the (primitive) irreducible Brillouin zone for a space
group with number `sgnum`, (conventional) direct lattice vectors `Rs`, and dimension `D`.
The path includes all distinct high-symmetry lines and points as well as relevant parts of
the Brillouin zone boundary.

The dimension `D` (1, 2, or 3) is specified as the third input argument, preferably as a
static `Val{D}` type parameter (or, type-unstably, as an `<:Integer`). Defaults to `Val(3)`.

`Rs` refers to the direct basis of the conventional unit cell, i.e., not the primitive 
direct basis vectors. The setting of `Rs` must agree with the conventional setting choices
in the International Tables of Crystallography, Volume A (the "ITA conventional setting").
If `Rs` is a subtype of a `StaticVector` or `NTuple`, the dimension can be inferred from its
(static) size; in this case, this dimension will take precedence (i.e. override, if
different) over any dimension specified in the third input argument.

## Notes
- The returned **k**-points are given in the basis of the **primitive** reciprocal basis in
  the CDML setting. To obtain the associated transformation matrices between the
  ITA conventional setting and the CDML primitive setting, see `primitivebasismatrix` of
  Bravais.jl](https://thchr.github.io/Crystalline.jl/stable/bravais/) (or, equivalently, the
  relations defined Table 2 of [^1]).
  To transform to a Cartesian basis, see [`cartesianize!`](@ref).
- To interpolate a `KPath`, see [`interpolate(::KPath, ::Integer)`](@ref) and
  [`splice(::KPath, ::Integer)`](@ref).
- All paths currently assume time-reversal symmetry (or, equivalently, inversion symmetry).
  If neither are present, include the "inverted" -**k** paths manually.

## Data and referencing
3D paths are sourced from the SeeK-path publication: please cite the original work
[^2].

## References
[^1]: Aroyo et al., [Acta Cryst. A70, 126 (2014)](https://doi.org/10.1107/S205327331303091X).
[^2]: Hinuma, Pizzi, Kumagai, Oba, & Tanaka, *Band structure diagram paths based on
      crystallography*, [Comp. Mat. Sci. **128**, 140 (2017)]
      (http://dx.doi.org/10.1016/j.commatsci.2016.10.015).
"""
function irrfbz_path(sgnum::Integer, Rs, Dᵛ::Val{D}=Val(3)) where D
    D′ = length(Rs)::Int
    D′ ≠ D && throw(DimensionMismatch("inconsistent dimensions of `Rs` and `Dᵛ`"))
    any(R -> length(R)::Int ≠ D, Rs) && throw(DimensionMismatch("inconsistent element dimensions in `Rs`"))
    Rs = convert(DirectBasis{D}, Rs)
    @static if VERSION ≥ v"1.9.0-DEV.1433"
        # for earlier versions of Julia, this could fail due to an internal compiler bug,
        # cf. issue #21 and https://github.com/JuliaLang/julia/issues/46871; this was fixed
        # in https://github.com/JuliaLang/julia/pull/46882; so only typeassert on Julia
        # versions where this is fixed.
        # the typeassert is here in order to improve inference for callers of `irrfbz_path`
        # where `Rs` is poorly inferred or entirely uninferred (e.g., `typeof(Rs) === Any`)
        Rs::DirectBasis{D}
    end

    # (extended) bravais type
    bt = bravaistype(sgnum, D; normalize=false)
    ext_bt = extended_bravais(sgnum, bt, Rs, Dᵛ)

    # get data about path
    lab2kvs = get_points(ext_bt, Rs, Dᵛ)
    paths = get_paths(ext_bt, Dᵛ)

    # unshuffle HPKOT-standard primitive setting-difference in oA and mC settings
    unshuffle_hpkot_setting!(lab2kvs, bt, D)
    
    # compute (primitive) reciprocal basis
    cntr = last(bt)
    pGs = primitivize(reciprocalbasis(Rs), cntr)

    return KPath(lab2kvs, paths, pGs, Ref(LATTICE))
end
function irrfbz_path(sgnum::Integer, 
            Rs::Union{StaticVector{D, <:AVec{<:Real}}, NTuple{D, <:AVec{<:Real}}}) where D
    return irrfbz_path(sgnum, Rs, Val(D))
end
irrfbz_path(sgnum::Integer, Rs::AVec{<:AVec{<:Real}}, D::Integer) = irrfbz_path(sgnum, Rs, Val(D))

# Note: output k-points are returned in the primitive setting defined in Table 3 of the 
# HPKOT paper. This setting differs from the standard crystallographic primitive setting
# for 'oA' and 'mC' Bravais types; "unshuffling" this is done in `unshuffle_hpkot_setting!`
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

function unshuffle_hpkot_setting!(lab2kvs, bt, D)
    # the points in `labs2kvs` are returned in the reciprocal basis of the primitive HPKOT
    # setting (Table 3, HPKOT paper). That setting differs from the primitive CDML setting
    # that we follow (Table 2, https://doi.org/10.1107/S205327331303091X) for Bravais types
    # 'oA' and 'mC'. Hence, we need to shuffle it back to CDML primitive settings.
    # Specifically, HPKOT has centering matrix `P′` which differs from the CDML centering
    # matrix `P`. To cast the points back to the CDML primitive setting we first transform
    # the points to the conventional setting (by multiplying by `((P′)⁻¹)ᵀ`) and then
    # transform back to the CDML primitive setting (by  multiplying with `Pᵀ`):
    if D == 3 && (bt == "oA" || bt == "mC")
        cntr = last(bt)
        if bt == "oA"
            P′ = @SMatrix [0 0 1; 0.5 0.5 0; -0.5 0.5 0] # HPKOT centering matrix
        else # bt == "mC"
            P′ = @SMatrix [0.5 -0.5 0; 0.5 0.5 0; 0 0 1]
        end
        P = primitivebasismatrix(cntr, Val(3))           # CDML centering matrix
        PᵀP′⁻¹ᵀ = P'*inv(P′)'
        for (lab, kv′) in lab2kvs
            kv = PᵀP′⁻¹ᵀ*kv′  # k-vector coordinates transform w/ transpose of matrix (as
            lab2kvs[lab] = kv # explained in e.g. Bravais.jl's src/transform.jl)
        end
    end
    return lab2kvs
end

# ---------------------------------------------------------------------------------------- #

include("interpolate-paths.jl")
export interpolate, splice, cumdists

# ---------------------------------------------------------------------------------------- #

"""
    cartesianize(kp::KPath{D})

Transform a **k**-path `kp` in a lattice basis to a Cartesian basis with (primitive)
reciprocal basis vectors `basis(kp)`.
"""
function cartesianize(x::Union{KPath, KPathInterpolant})
    return setting(x) === LATTICE ? cartesianize!(deepcopy(x)) : deepcopy(x)
end

"""
    latticize(kp::KPath)

Transform a **k**-path `kp` in a Cartesian basis to a lattice basis with (primitive)
reciprocal lattice vectors `basis(kp)`.
"""
function latticize(x::Union{KPath, KPathInterpolant})
    return setting(x) === CARTESIAN ? latticize!(deepcopy(x)) : deepcopy(x)
end

"""
    cartesianize!(kp::KPath{D})

Transform a **k**-path `kp` in a lattice basis to a Cartesian basis with (primitive)
reciprocal basis vectors `basis(kp)`. Modifies `kp` in-place.
"""
function cartesianize!(kp::KPath{D}) where D
    setting(kp) === CARTESIAN && return kp
    for (lab, kv) in points(kp)
        @inbounds points(kp)[lab] = cartesianize(kv, basis(kp))
    end
    set_setting!(kp, CARTESIAN)
    return kp
end

"""
    latticize!(kp::KPath)

Transform a **k**-path `kp` in a Cartesian basis to a lattice basis with (primitive)
reciprocal lattice vectors `basis(kp)`. Modifies `kp` in-place.
"""
function latticize!(kp::KPath)
    setting(kp) === LATTICE && return kp
    basismatrix = reduce(hcat, basis(kp))
    for (lab, kv) in points(kp)
        @inbounds points(kp)[lab] = latticize(kv, basismatrix)
    end
    set_setting!(kp, LATTICE)
    return kp
end

"""
    latticize(kp::KPath{D}, basis::AbstractVector{<:AbstractVector{<:Real}})

Transform a **k**-path `kp` in a Cartesian basis to a lattice basis with (primitive)
reciprocal lattice vectors `basis`.

If `kp` is not in a Cartesian basis (i.e., if `setting(kp) == LATTICE`), `kp` is returned
as-is.
"""
function latticize(kp::KPath{D}, basis::AVec{<:AVec{<:Real}}) where D
    setting(kp) === LATTICE && return kp
    basismatrix = convert(SMatrix{D,D,Float64,D*D}, reduce(hcat, basis))
    kvs′ = Dict{Symbol, SVector{D,Float64}}()
    for (lab, kv) in points(kp)
        @inbounds kvs′[lab] = latticize(kv, basismatrix)
    end
    return typeof(kp)(kvs′, paths(kp), basis, Ref(LATTICE))
end

# ---------------------------------------------------------------------------------------- #
end # module