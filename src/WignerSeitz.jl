module WignerSeitz

# ---------------------------------------------------------------------------------------- #

using ..Brillouin: AVec
using Requires
using StaticArrays
using LinearAlgebra: dot, cross, norm
using PyCall

import Base: getindex, size, IndexStyle, show, summary

# ---------------------------------------------------------------------------------------- #

export Cell, wignerseitz, faces, vertices, basis

# ---------------------------------------------------------------------------------------- #

const PySpatial = PyNULL()
function __init__()
    # bringing in SciPy's Spatial module (for `Voronoi` and `ConvexHull`)
    copy!(PySpatial, pyimport_conda("scipy.spatial", "scipy"))

    # plotting extensions on GLMakie load
    @require AbstractPlotting="537997a7-5e4e-5d89-9595-2241ea00577e" begin
        include("compat/abstractplotting_wignerseitz.jl")
    end

    # plotting extensions on PlotlyJS load
    @require PlotlyJS="f0f68f2c-4968-5e81-91da-67840de0976a" begin
        include("compat/plotlyjs_wignerseitz.jl")
    end
end

# ---------------------------------------------------------------------------------------- #
# STRUCTURES

struct Cell{D} <: AVec{Vector{SVector{D, Float64}}} # over of polygons
    verts :: Vector{SVector{D, Float64}}
    faces :: Vector{Vector{Int}}
    basis :: SVector{D, SVector{D, Float64}}
end
faces(c::Cell)    = c.faces
vertices(c::Cell) = c.verts
basis(c::Cell)    = c.basis

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

    println(io, " verts: ", round.(first(vertices(c)), digits=3))
    foreach(v -> println(io, "        ", round.(v, digits=3)), @view vertices(c)[2:end])
     
    println(io, " faces: ", first(faces(c)))
    foreach(f -> println(io, "        ", f), @view faces(c)[2:end])

    println(io, " basis: ", round.(first(basis(c)), digits=3))
    foreach(v -> println(io, "        ", round.(v, digits=3)), @view basis(c)[2:D])
end
# ---------------------------------------------------------------------------------------- #
# MAIN FUNCTION

"""
    wignerseitz(Vs::AbstractVector{<:SVector{D}}; merge::Bool = true, Nmax = 3)
    wignerseitz(Vs::AbstractVector{AbstractVector}; merge::Bool = true, Nmax = 3)
                                                                --> Cell{D}

Return a `Cell{D}` structure, containing the vertices and associated (outward oriented)
faces, of the Wigner-Seitz cell defined by a basis `Vs` in `D` dimensions.

## Keyword arguments
- `merge` (default, `true`): if `:true`, co-planar faces are merged to form polygonal
  planar faces (e.g., triangles, quadrilaterals, and ngons generally). If `false`, raw
  "unprocessed" triangles (`D=3`) and segments (`D=2`) are returned instead. `merge` has no
  impact for `D=1`.
- `Nmax` (default, `3`): includes `-Nmax:Nmax` points in the initial lattice used to
  generate the underlying Voronoi tesselation. It is unwise to set this to anything lower
  than 3 without explicitly testing convergence; and probably unnecessary to increase it
  beyond 3 as well.
"""
function wignerseitz(Vs::AVec{<:SVector{D,<:Real}};
            merge::Bool = true,
            Nmax::Integer = 3) where D
    # "supercell" lattice of G-vectors
    Ns = -Nmax:Nmax
    lattice = Vector{SVector{D,Float64}}(undef, length(Ns)^D)
    idx_cntr = 0
    for (idx, I) in enumerate(CartesianIndices(ntuple(_->Ns, Val(D))))
                   V =  Vs[1]*I[1]
        D >= 2 && (V += Vs[2]*I[2])
        D == 3 && (V += Vs[3]*I[3])
        lattice[idx] = V
        iszero(I) && (idx_cntr = idx)
    end

    vor = PySpatial.Voronoi(lattice) # voronoi tesselation of lattice

    # grab all the vertices of the central voronoi region enclosing origo
    verts_cntr =  # NB: offsets by 1 due to Julia 1-based vs. Python 0-based indexing
        [vor.vertices[idx+1,:] for idx in vor.regions[vor.point_region[idx_cntr]+1]]

    # get convex hull of central vertices
    hull = PySpatial.ConvexHull(verts_cntr)
    c    = convert_to_cell(hull, Vs)
    c    = reorient_normals!(c)

    # return either raw simplices or "merged" polygons
    if merge
        return merge_coplanar!(c)
    else
        return c
    end
end
# overload for input as ordinary vectors; type-unstable obviously, but convenient sometimes
function wignerseitz(Vs::AVec{<:AVec{<:Real}}; kwargs...)
    D = length(first(Vs))
    all(V->length(V) == D, @view Vs[2:end]) || error(DomainError(Vs, "provided vectors `Vs` must have identical dimension"))

    return wignerseitz(SVector{D,Float64}.(Vs); kwargs...)
end

function wignerseitz(Vs::AVec{<:SVector{1,<:Real}}; kwargs...)
    a = only(only(Vs))
    Cell{1}([(@SVector [-a/2]), (@SVector [a/2])], [1,2], Vs)
end
# ---------------------------------------------------------------------------------------- #
# UTILITIES & SUBFUNCTIONS

function convert_to_cell(hull, Vs::AVec{<:SVector{D,<:Real}}) where D
    vs′ = hull.points         # vertices
    fs′ = hull.simplices .+ 1 # faces

    vs = SVector{D, Float64}.(eachrow(vs′))
    fs = Vector{Int}.(eachrow(fs′))
    return Cell(vs, fs, SVector{D, SVector{D, Float64}}(Vs))
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
    n = cross(v₁₂, v₂₃)

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

end # module