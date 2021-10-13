module WignerSeitz

# ---------------------------------------------------------------------------------------- #

using ..Brillouin: 
    AVec,
    SHOWDIGITS,
    BasisLike,
    BasisEnum, CARTESIAN, LATTICE, setting, set_setting!,
    cartesianize, latticize
import ..Brillouin:
    latticize!,
    cartesianize!,
    basis
using StaticArrays
using LinearAlgebra:
    norm,
    dot,
    ×
using PyCall
using DocStringExtensions

import Base: getindex, size, IndexStyle, show, summary

# ---------------------------------------------------------------------------------------- #

export Cell, wignerseitz, faces, vertices, reduce_to_wignerseitz

# ---------------------------------------------------------------------------------------- #

const PySpatial = PyNULL()
function __init__()
    # bringing in SciPy's Spatial module (for `Voronoi` and `ConvexHull`)
    copy!(PySpatial, pyimport_conda("scipy.spatial", "scipy"))
end

# ---------------------------------------------------------------------------------------- #
# STRUCTURES

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Cell{D} <: AVec{Vector{SVector{D, Float64}}}
    verts   :: Vector{SVector{D, Float64}}
    faces   :: Vector{Vector{Int}}
    basis   :: SVector{D, SVector{D, Float64}}
    setting :: Ref{BasisEnum} # internal field
end

faces(c::Cell)    = c.faces
vertices(c::Cell) = c.verts

# abstract array interface
Base.@propagate_inbounds getindex(c::Cell, i::Int) = vertices(c)[faces(c)[i]]
size(c::Cell) = size(c.faces)
IndexStyle(::Type{<:Cell}) = IndexLinear()

function summary(io::IO, c::Cell{D}) where D
    print(io, "Cell{", D, "} (", 
              length(c), " faces, ", 
              length(vertices(c)), " vertices)")
end
function show(io::IO, ::MIME"text/plain", c::Cell{D}) where D
    summary(io, c)
    println(io, ":")

    println(io, "  verts: ", round.(first(vertices(c)), digits=SHOWDIGITS))
    foreach(v -> println(io, "         ", round.(v, digits=SHOWDIGITS)), @view vertices(c)[2:end])
     
    println(io, "  faces: ", first(faces(c)))
    foreach(f -> println(io, "         ", f), @view faces(c)[2:end])

    print(io, "  basis: ")
    foreach(enumerate(basis(c))) do (i,v) # dance to ensure final line isn't `\n`-terminated
        i ≠ 1 && print(io, "         ")
        if i ≠ length(basis(c))
            println(io, round.(v, digits=SHOWDIGITS))
        else
            print(io, round.(v, digits=SHOWDIGITS))
        end
    end
end
# ---------------------------------------------------------------------------------------- #
# MAIN FUNCTION

"""
    wignerseitz(basis::AbstractVector{<:SVector{D}}; merge::Bool = true, Nmax = 3)
    wignerseitz(basis::AbstractVector{<:AbstractVector}; merge::Bool = true, Nmax = 3)
                                                                --> Cell{D}

Given a provided basis, return a `Cell{D}` containing the vertices and associated (outward
oriented) faces of the Wigner-Seitz cell defined by `basis` in `D` dimensions. The returned
vertices are given in the the basis of `basis` (see [`cartesianize!`](@ref) for conversion).

## Keyword arguments
- `merge` (default, `true`): if `true`, co-planar faces are merged to form polygonal
  planar faces (e.g., triangles, quadrilaterals, and ngons generally). If `false`, raw
  "unprocessed" triangles (`D=3`) and segments (`D=2`) are returned instead. `merge` has no
  impact for `D=1`.
- `Nmax` (default, `3`): includes `-Nmax:Nmax` points in the initial lattice used to
  generate the underlying Voronoi tesselation. It is unwise to set this to anything lower
  than 3 without explicitly testing convergence; and probably unnecessary to increase it
  beyond 3 as well.
"""
function wignerseitz(basis::AVec{<:SVector{D,<:Real}};
            merge::Bool = true,
            Nmax::Integer = 3) where D
    # "supercell" lattice of G-vectors
    Ns = -Nmax:Nmax
    lattice = Vector{SVector{D,Float64}}(undef, length(Ns)^D)
    idx_cntr = 0
    for (idx, I) in enumerate(CartesianIndices(ntuple(_->Ns, Val(D))))
                   V =  basis[1]*I[1]
        D >= 2 && (V += basis[2]*I[2])
        D == 3 && (V += basis[3]*I[3])
        lattice[idx] = V
        iszero(I) && (idx_cntr = idx)
    end

    vor = PySpatial.Voronoi(lattice) # voronoi tesselation of lattice

    # grab all the vertices of the central voronoi region enclosing origo
    verts_cntr =  # NB: offsets by 1 due to Julia 1-based vs. Python 0-based indexing
        [vor.vertices[idx+1,:] for idx in vor.regions[vor.point_region[idx_cntr]+1]]

    # get convex hull of central vertices
    hull = PySpatial.ConvexHull(verts_cntr)
    c    = convert_to_cell(hull, basis)
    c    = reorient_normals!(c)

    # return either raw simplices or "merged" polygons
    if merge
        return latticize!(merge_coplanar!(c))
    else
        return latticize!(c)
    end
end
# overload for input as ordinary vectors; type-unstable obviously, but convenient sometimes
function wignerseitz(basis::AVec{<:AVec{<:Real}}; kwargs...)
    D = length(first(basis))
    all(V->length(V) == D, @view basis[2:end]) || error(DomainError(basis, "provided `basis` must have identical dimensions"))

    return wignerseitz(SVector{D,Float64}.(basis); kwargs...)
end

function wignerseitz(basis::AVec{<:SVector{1,<:Real}}; kwargs...)
    a = only(only(basis))
    Cell{1}([(@SVector [-a/2]), (@SVector [a/2])], [1,2], basis)
end
# ---------------------------------------------------------------------------------------- #
# UTILITIES & SUBFUNCTIONS

function convert_to_cell(hull, basis::AVec{<:SVector{D,<:Real}}) where D
    vs′ = hull.points         # vertices
    fs′ = hull.simplices .+ 1 # faces

    vs = SVector{D, Float64}.(eachrow(vs′))
    fs = Vector{Int}.(eachrow(fs′))
    return Cell(vs, fs, SVector{D, SVector{D, Float64}}(basis), Ref(CARTESIAN))
end


"""
    is_outward_oriented(c, n)
    
Return whether a face, specified by its center `c` and normal `n`, is outward pointing,
assuming the face is part of a convex hull centered around origo.
"""
is_outward_oriented(c::SVector{D}, n::SVector{D}) where D = dot(c,n) > zero(eltype(c))

function reorient_normals!(c::Union{D}) where D
    D == 1 && return c # no sense of "outward" in 1D
    fs = faces(c)
    for (i,f) in enumerate(fs)
        cntr = face_center(c, i)
        n    = face_normal(c, i)
        if !is_outward_oriented(cntr, n)
            reverse!(fs[i])
        end
    end
    return c
end

function merge_coplanar!(c::Cell{3})
    fs = faces(c)
    ns = face_normals(c)

    # merge the all co-planar face; step forward in a "conquering" fashion
    i = 1
    while i ≤ length(fs)
        fᵢ = fs[i]
        maybe_more_to_merge = true
        while maybe_more_to_merge 
            # this outer loop is here to ensure we do the j-loop until convergence (each
            # j-loop can create new merge-options if two faces were indeed merged)
            j  = i+1
            maybe_more_to_merge = false
            while j ≤ length(fs)
                fⱼ = fs[j]
                # can be co-planar faces if they share two or more vertices; we only 
                # support merging with two vertices currently, so just check for that
                if sum(∈(fᵢ), fⱼ) == 2 && is_coplanar(ns[i], ns[j])
                    fs[i] = fᵢ = merge_coplanar(fᵢ, fⱼ)
                    deleteat!(fs, j), deleteat!(ns, j)
                    was_merged = true
                    maybe_more_to_merge = true
                else
                    was_merged = false
                end
                was_merged || (j += 1)
            end
        end
        i += 1
    end

    return c
end

# compute an (oriented) face-normal from a face `f` and a list of associated vertices `vs`
function face_normal(vs::AVec{<:SVector{3,<:Real}}, f::Vector{Int})
    v₁₂ = vs[f[2]] - vs[f[1]]
    v₂₃ = vs[f[3]] - vs[f[2]]
    n = v₁₂ × v₂₃

    return n/norm(n)
end
function face_normal(vs::AVec{<:SVector{2,<:Real}}, f::Vector{Int})
    # define "canonical direction" as 1->2
    v₁₂ = vs[f[2]] - vs[f[1]]
    n = SVector{2, eltype(v₁₂)}(v₁₂[2], -v₁₂[1]) # = n = [y, -x]
    
    return n/norm(n)
end
face_normal(c::Cell, i::Int) = face_normal(vertices(c), faces(c)[i])
face_normals(c::Cell)        = face_normal.(Ref(vertices(c)), faces(c))

# compute a *rough* center position of a polygon
function face_center(vs::Vector{SVector{D,Float64}}, f::Vector{Int}) where D
    return sum(i -> vs[i], f)/D
end
face_center(c::Cell, i::Int) = face_center(vertices(c), faces(c)[i])

# assumes normal vectors to be normalized and identically oriented
is_coplanar(n1, n2, atol::Real=1e-8) = isapprox(dot(n1, n2), 1.0, atol=atol)

function merge_coplanar(fᵢ::Vector{Int}, fⱼ::Vector{Int})
    idxs_sharedᵢ = sort!(findall(∈(fⱼ), fᵢ))
    idxs_sharedⱼ = [findfirst(==(fᵢ[i]), fⱼ)::Int for i in idxs_sharedᵢ]

    fᵢ[idxs_sharedᵢ] == fⱼ[idxs_sharedⱼ] || error("failed face identification")
    length(idxs_sharedᵢ) != 2 && throw(DomainError(length(idxs_sharedᵢ), "unexpected: faces share fewer or more than two vertices"))

    # now we do a semi-terrifying dance to "follow-along" the edges of face i, then 
    # transition into face j, and finally follow any remaining edges in face j; the trick
    # is to jump at from face i to j and vice versa when we hit a shared vertex; it's made
    # tedious by the fact that the face is not "cyclical" and we have to emulate that 
    # ourselves
    Dᵢ, Dⱼ = length(fᵢ), length(fⱼ)
    f = Vector{Int}(undef, Dᵢ+Dⱼ-2)
    n = 1
    # start in ngon i
    idxsᵢ, neighbourⱼ = if idxs_sharedᵢ[1] == 1 && idxs_sharedᵢ[2] == Dᵢ
        2:Dᵢ, 2
    else
        1:idxs_sharedᵢ[1], 1
    end
    for idxᵢ in idxsᵢ
        f[n] = fᵢ[idxᵢ]
        n += 1
    end

    # go to ngon j
    if idxs_sharedⱼ[neighbourⱼ] == Dⱼ
        idxsⱼ = Iterators.flatten((1:Dⱼ-1,))
    else
        if idxs_sharedⱼ[neighbourⱼ] == 1
            idxsⱼ = Iterators.flatten((idxs_sharedⱼ[neighbourⱼ]+1:Dⱼ,))
        else
            other_neighborⱼ = neighbourⱼ == 2 ? 1 : 2
            idxsⱼ = Iterators.flatten((idxs_sharedⱼ[neighbourⱼ]+1:Dⱼ, 1:idxs_sharedⱼ[other_neighborⱼ]))
        end
    end
    for idxⱼ in idxsⱼ
        f[n] = fⱼ[idxⱼ]
        n += 1
    end

    # done with ngon j; go back to ngon i
    remaining_idxsᵢ = idxs_sharedᵢ[2]+1:Dᵢ
    for idxᵢ in remaining_idxsᵢ
        f[n] = fᵢ[idxᵢ]
        n += 1
    end

    return f
end

function merge_coplanar!(c::Cell{2})
    fs′ = segments2polygon(faces(c))
    empty!(faces(c))
    push!(faces(c), fs′)
    return c
end
function segments2polygon(fs::Vector{Vector{Int}})
    # assumes there is a *single* well-defined polygon defined by `fs`
    N   = size(fs,1) # number of segments
    fs′ = Vector{Int}(undef, N)
    fs′[1:2] .= fs[1] # starting point
    for i in 3:N
        j = findfirst(f -> f[1]==(fs′[i-1]), fs)
        fs′[i] = fs[j][2]
    end
    return fs′
end

# ---------------------------------------------------------------------------------------- #

function latticize!(c::Cell)
    setting(c) === LATTICE && return c
    basismatrix = hcat(c.basis...)
    c.verts .= latticize.(c.verts, Ref(basismatrix))
    set_setting!(c, LATTICE)
    return c
end
function cartesianize!(c::Cell)
    setting(c) === CARTESIAN && return c
    c.verts .= cartesianize.(c.verts, Ref(c.basis))
    set_setting!(c, CARTESIAN)
    return c
end

# ---------------------------------------------------------------------------------------- #

reduce_to_symmetric_unitcell(x::Real) = mod(x + 1//2, 1) - 1//2
reduce_to_symmetric_unitcell(v::AbstractVector{<:Real}) = reduce_to_symmetric_unitcell.(v)

"""
    reduce_to_wignerseitz(v::StaticVector, Vs::BasisLike)  -->  v′

Return the periodic image `v′` of the point `v` in the basis `Vs`.

`v` is assumed to be provided in the lattice basis (i.e., relative to `Vs`) and `v′` is
returned in similar fashion.

The returned point `v′` lies in the Wigner-Seitz cell (or its boundary) defined by `Vs`,
has the least possible norm among all equivalent images of `v`, and differs from `v` at
most by integer lattice translations such that `mod(v, 1) ≈ mod(v′, 1)`.
"""
function reduce_to_wignerseitz(v::StaticVector{D, <:Real}, Vs::BasisLike{D}) where D
    # reduce `v` to "rectangular" unit cell with coordinates in [-1/2, 1/2)
    v₀ = reduce_to_symmetric_unitcell(v)
    v₀ᶜ = cartesianize(v₀, Vs)
    d₀ = norm(v₀ᶜ)
    # check whether `v₀` has the smallest norm of all adjacent equivalent points; if not,
    # we initiate a recursive search the new possible minimum until convergence
    return search_shortest_norm_among_neighbors_recursive(v₀, d₀, Vs)
end
function reduce_to_wignerseitz(v::AVec{<:Real}, Vs::AVec{<:AVec{<:Real}})
    D = length(v)
    D == length(Vs) || throw(DimensionMismatch("dimensions of `v` and `Vs` must match"))
    all(V -> length(V) == D, Vs) || throw(DimensionMismatch("internal dimensions of `Vs` do not match"))

    # NB: type-unstable; but only really exists as a convenience accessor...
    v_static  = convert(SVector{D, eltype(v)}, v)
    Vs_static = convert(SVector{D, SVector{D, eltype(first(Vs))}}, Vs)

    return reduce_to_wignerseitz(v_static, Vs_static)
end

function search_shortest_norm_among_neighbors_recursive(
            v₀::StaticVector{D, <:Real}, d₀::Real, Vs::BasisLike{D}) where D

    for I in CartesianIndices(ntuple(_->-1:1, Val(D)))
        iszero(I) && continue
        v′ = v₀ .+ Tuple(I)
        v′ᶜ = cartesianize(v′ᶜ, Vs)
        d′ = norm(v′ᶜ)
        if d′ < d₀
            return search_shortest_norm_among_neighbors_recursive(v′, d′, Vs)
        end
    end
    return v₀
end
end # module