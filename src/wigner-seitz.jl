module WignerSeitz

using StaticArrays, Crystalline
using GeometryBasics
using PyCall
using GLMakie
using Crystalline: Basis
using LinearAlgebra: dot, cross
using Colors

const PySpatial = PyNULL()
function __init__()
    copy!(PySpatial, pyimport_conda("scipy.spatial", "scipy"))
end

export wignerseitz

"""
    wignerseitz(Vs::Basis, output::Val = Val(:polygons); Nmax::Integer = 3)

Return the vertices and associated (outward oriented) faces of the Wigner-Seitz cell defined
by a basis `Vs`.

If the `output = Val(:polygons)` (default), co-planar faces are merged to form polygonal
planar faces of arbitrary order. If instead `output = Val(:triangles)` "unprocessed"
triangles (face of order 3) are returned instead.
"""
function wignerseitz(Vs::Basis{D}, ::Val{O} = Val(:polygons); Nmax::Integer = 3) where {D, O}
    # "supercell" lattice of G-vectors
    Ns = -Nmax:Nmax
    lattice = Vector{Point{D,Float64}}(undef, length(Ns)^D)
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

    # get convex hull of central verts
    hull = PySpatial.ConvexHull(verts_cntr)
    m    = convert_to_mesh(hull, Val(D))
    m    = reorient_normals!(m)

    # return either triangles or polygons
    if O == :polygons
        return merge_coplanar(m)
    elseif O == :triangles
        return m
    else
        error(DomainError(O, "invalid output type"))
    end
end


function convert_to_mesh(hull, ::Val{D}) where D
    vs′ = hull.points         # vertices
    fs′ = hull.simplices .+ 1 # faces

    vs = Point{D}.(eachrow(vs′))
    fs = NgonFace{D}.(eachrow(fs′))
    return Mesh(vs, fs)
end


"""
    is_outward_pointing(c, n)
    
Return whether a face, specified by its center `c` and normal `n`, is outward pointing,
assuming the face is part of a convex hull centered around origo.
"""
function is_outward_pointing(c, n)
    return dot(c,n) > 0
end

function reorient_normals!(m::Mesh)
    fs = faces(m)
    for (i,f) in enumerate(fs)
        c, n = face_normal(m[i]), face_center(m[i])
        if !is_outward_pointing(c, n)
            fs[i] = typeof(f)(reverse(f)...)
        end
    end
    return m
end

function merge_coplanar(m::Mesh{D}) where D
    fs = faces(m)
    ns = face_normal.(m)

    # merge the all co-planar face; step forward in a "conquering" fashion
    mfs = Vector{NgonFace}(fs) # merged faces
    i = 1
    while i ≤ length(mfs)
        fᵢ = mfs[i]
        j  = i+1
        while j ≤ length(mfs)
            fⱼ = mfs[j]
            # can be co-planar faces if they share two or more vertices; we only 
            # support merging with two vertices currently, so just check for that
            if sum(∈(fᵢ), fⱼ) == D-1 && is_coplanar(ns[i], ns[j])
                mfs[i] = fᵢ = merge_coplanar(fᵢ, fⱼ)
                deleteat!(mfs, j), deleteat!(ns, j)
                was_merged = true
            else
                was_merged = false
            end
            was_merged || (j += 1)
        end
        i += 1
    end

    return coordinates(m), mfs
end

# compute an (oriented) face-normal
function face_normal(p::GeometryBasics.Ngon{D}) where D
    D == 2 && error("cannot compute a face normal for a 2D Ngon")

    v₁₂ = p[2] - p[1]
    v₂₃ = p[3] - p[2]
    n = cross(v₁₂, v₂₃)
    #n = GeometryBasics.orthogonal_vector(p[1], p[2], p[3])
    return GeometryBasics.normalize(n)
end

face_center(p::GeometryBasics.Ngon) = sum(p)/length(p)

# assumes normal vectors to be normalized and identically oriented
is_coplanar(n1, n2, atol::Real=1e-8) = isapprox(dot(n1, n2), 1.0, atol=atol)

function merge_coplanar(fᵢ::NgonFace{Dᵢ,T}, fⱼ::NgonFace{Dⱼ,T}) where {Dᵢ,Dⱼ,T}
    idxs_sharedᵢ = sort!(findall(∈(fⱼ), fᵢ))
    idxs_sharedⱼ = [findfirst(==(fᵢ[i]), fⱼ)::Int for i in idxs_sharedᵢ]

    fᵢ[idxs_sharedᵢ] == fⱼ[idxs_sharedⱼ] || error("failed face identification")
    length(idxs_sharedᵢ) != 2 && throw(DomainError(length(idxs_sharedᵢ), "unexpected: faces share fewer or more than two vertices"))

    # now we do a semi-terrifying dance to "follow-along" the edges of face i, then 
    # transition into face j, and finally follow any remaining edges in face j; the trick
    # is to jump at from face i to j and vice versa when we hit a shared vertex; it's made
    # tedious by the fact that the face is not "cyclical" and we have to emulate that 
    # ourselves
    f = MVector{Dᵢ+Dⱼ-2, Int}(undef)
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
        other_neighborⱼ = neighbourⱼ == 2 ? 1 : 2
        if idxs_sharedⱼ[neighbourⱼ] == 1
            idxsⱼ = Iterators.flatten((idxs_sharedⱼ[neighbourⱼ]+1:Dⱼ,))
        else
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

    return NgonFace{Dᵢ+Dⱼ-2}(f...)
end

function visualize_mesh(m::Mesh, s::Scene=Scene())
    #mesh!(s, m, color=:blue, transparent=true)
    wireframe!(s, m, color=:black, linewidth=4)
    ns = face_normal.(m)
    cs = face_center.(m)
    #for (n, c) in zip(ns, cs)
    @show ns
        arrows!(s,
                getindex.(cs, 1), getindex.(cs, 2), getindex.(cs, 3),
                getindex.(ns, 1), getindex.(ns, 2), getindex.(ns, 3),)
    #end
    display(s)
end

function visualize_ngons(vs, fs, s::Scene=Scene())
    cam3d!(s)
    for (i,f) in enumerate(fs)
        path = vs[f]
        lines!(s, vcat(collect(getindex.(path, 1)), path[1][1]), 
                  vcat(collect(getindex.(path, 2)), path[1][2]), 
                  vcat(collect(getindex.(path, 3)), path[1][3]),
                  color=RGB(.15,.25,.8), linewidth=2
                  )
    end
    display(s)
end

end