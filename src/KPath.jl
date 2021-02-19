# Paths in the representation domain of the BZ, obtained from the [SeeK] publication. 
# Methods return *k*-vector coordinates in a primitive reciprocal basis ("reciprocal 
# crystallographic primitive cell" in the [SeeK]'s nomenclature).
#
# ^[Seek] http://dx.doi.org/10.1016/j.commatsci.2016.10.015; see also online interface at 
#         https://www.materialscloud.org/work/tools/seekpath

module KPath

# ---------------------------------------------------------------------------------------- #
export irrfbz_path
# ---------------------------------------------------------------------------------------- #
using ..CrystallineBravaisVendor: bravaistype
using ..Brillouin: AVec, BasisLike
using LinearAlgebra: norm
using StaticArrays
# ---------------------------------------------------------------------------------------- #

include("interpolate-paths.jl")
export splice_path, interpolate_path, cumdists

const interpolate_kpath = interpolate_path
const splice_kpath      = splice_path
@deprecate interpolate_kpath interpolate_path true
@deprecate splice_kpath splice_path true

# ---------------------------------------------------------------------------------------- #

include("codegen_kpoints.jl")
include("bravais-branches.jl")

# ---------------------------------------------------------------------------------------- #

"""
    irrfbz_path(sgnum::Integer, Nk::Integer, Rs::Union{Nothing, $(BasisLike)}=nothing; 
        pathtype::String="SeeK", has_inversion_or_tr::Bool=true,
        splice::Bool=false, legacy::Bool=false)
                                                --> paths_kvs, paths_labs, lab2kv

Returns a **k**-path in the (primitive) irreducible Brillouin zone that includes all 
distinct high-symmetry lines and points as well as parts of the Brillouin zone boundary.

`Rs` refers to the direct basis of the conventional unitcell. For some space groups, it
is needed to disambiguate the "extended" Bravais types that may differ depending on the
lengths of the lattice vectors (because the Brillouin zone may depend on these lengths).
If the requested space group is known to not fall under this case, `Rs` can be supplied
as `nothing` (default).

## Data and referencing
All data is sourced from the SeeK publication[^1]: please cite the original work.

All paths currently assume time-reversal symmetry (or, equivalently, inversion symmetry), 
corresponding to the SeeK's `[with inversion]` setting. If neither inversion nor
time-reversal symmetry is present, additional paths may be required (SeeK's `[no inversion]`
setting).

[1] Hinuma, Pizzi, Kumagai, Oba, & Tanaka, *Band structure diagram paths based on
    crystallography*, Comp. Mat. Sci. **128**, 140 (2017)](http://dx.doi.org/10.1016/j.commatsci.2016.10.015)
    (see also online interface at https://www.materialscloud.org/work/tools/seekpath).
"""
function irrfbz_path(sgnum::Integer, Nk::Integer, 
                Rs::Union{Nothing, BasisLike{D}}=nothing;
                pathtype::String="SeeK", has_inversion_or_tr::Bool=true,
                splice::Bool=false, legacy::Bool=false) where D

    if Rs isa AbstractVector
        D ≠ 3 && throw(DomainError(D, "Currently only implemented for 3D space groups"))
    end
    bt = bravaistype(sgnum, 3)

    # get high-symmetry points and labels and associated paths
    if has_inversion_or_tr
        lab2kv, paths_labs = _irrfbz_path(bt, Rs, sgnum, pathtype)
    else
        lab2kv, paths_labs = _irrfbz_path_without_inversion_or_tr(bt, Rs, sgnum, pathtype)
    end

    # interpolate k-paths
    paths_kvs = interp_paths_from_labs(lab2kv, paths_labs, Nk, splice, legacy) 

    return paths_kvs, paths_labs, lab2kv
end

# always returns `lab2kv, paths_labs`
function _irrfbz_path(bt::String, Rs::Union{Nothing, BasisLike{D}}, sgnum::Integer,
                      pathtype::String) where D
    if bt == "tP"       # ⇒ extended Bravais type tP1
        # sgs 75:78, 81, 83:86, 89:106, 111:118, 123:138
        return irrfbz_data_tP1(; pathtype=pathtype)
    elseif bt == "tI"  # ⇒ extended Bravais type tI1 and tI2
        # sgs 79, 80, 82, 87, 88, 97, 98, 107, 108:110, 119:122, 139:142
        Rs === nothing && _throw_requires_direct_basis()
        a, c = norm(Rs[1]), norm(Rs[3])
        if c ≤ a     # tI1
            return irrfbz_data_tI1(a, c; pathtype=pathtype)
        else #c > a  # tI2
            return irrfbz_data_tI2(a, c; pathtype=pathtype)
        end

    elseif bt == "cI"  # ⇒ extended Bravais type cI1
        # sgs 197, 199, 204, 206, 211, 214, 217, 220, 229, 230
        return irrfbz_data_cI1(; pathtype=pathtype)

    elseif bt == "oC"   # ⇒ extended Bravais types oC1 and oC2
        # sgs 20, 21, 35:41, 63:68
        Rs === nothing && _throw_requires_direct_basis()
        a,b,_ = norm.(Rs)
        if a < b     # oC1
            return irrfbz_data_oC1(a, b; pathtype=pathtype)
        else #a > b  # oC2
            return irrfbz_data_oC2(a, b; pathtype=pathtype)
        end
    elseif bt == "hP"  # ⇒ extended Bravais types hP1 and hP2
        if sgnum ∈ (143:149..., 151, 153, 157, 159:163...)  # hP1 ("long path")
            return irrfbz_data_hP1(; pathtype=pathtype)
        else                                                # hP2 ("short path")
            return irrfbz_data_hP2(; pathtype=pathtype)
        end

    else
        throw(DomainError(bt, "The requested Bravais type has not yet been implemented"))
    end
end

# As relevant in systems that have neither inversion nor time-reversal symmetry
function _irrfbz_path_without_inversion_or_tr(bt::String, Rs::BasisLike, sgnum::Integer,
                                              pathtype::String)
    if bt == "tI"  # ⇒ extended Bravais type tI1 and tI2
        # sgs 79, 80, 82, 87, 88, 97, 98, 107, 108:110, 119:122, 139:142
        a, c = norm(Rs[1]), norm(Rs[3])
        if c ≤ a     # tI1
            return irrfbz_data_tI1_without_inversion_or_tr(a, c; pathtype=pathtype)
        else #c > a  # tI2
            throw("Not yet implemented")
            #return irrfbz_data_tI2_without_inversion_or_tr(a, c; pathtype=pathtype)
        end

    else
        throw(DomainError(bt, "The requested Bravais type has not been implemented with "*
                              "broken inversion/TR symmetry yet"))
    end
end

# ---------------------------------------------------------------------------------------- #

# Table 82 of [SeeK]: extended Bravais type oC1 (requiring a<b) (with inversion)
function irrfbz_data_oC1(a::Real, b::Real; pathtype::String="SeeK")
    a<b || throw(DomainError((a,b), "requires a<b"))

    ζ = (1.0 + (a/b)^2)/4
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :Γ  => [0.0, 0.0, 0.0], :Y  => [-0.5, 0.5, 0.0], :T  => [-0.5, 0.5, 0.5],
        :Z  => [0.0, 0.0, 0.5], :S  => [0.0, 0.5, 0.0],  :R  => [0.0, 0.5, 0.5],
        :Σ₀ => [ζ, ζ, 0.0],     :C₀ => [-ζ, 1.0-ζ, 0.0], :A₀ => [ζ, ζ, 0.5], 
        :E₀ => [-ζ, 1.0-ζ, 0.5]
        )
    if pathtype == "SeeK"       # suggested path from [SeeK]
        paths_labs = [[:Γ, :Y, :C₀],  [:Σ₀, :Γ, :Z, :A₀],  [:E₀, :T, :Y], [:Γ, :S, :R, :Z, :T]]
    elseif pathtype == "simple" # simpler but slightly longer (and maybe redundant) path
        paths_labs = [[:Γ, :S, :R, :Z, :T, :Y, :Γ, :Σ₀, :A₀, :Z, :Γ], [:T, :E₀, :C₀, :Y]]
    else
        _throw_pathtype(pathtype)
    end

    return lab2kv, paths_labs
end

# Table 83 of [SeeK]: extended Bravais type oC2 (requiring a>b) (with inversion)
function irrfbz_data_oC2(a::Real, b::Real; pathtype::String="SeeK")
    a>b || throw(DomainError((a,b), "requires a>b"))

    ζ = (1.0 + (b/a)^2)/4
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :Γ  => [0.0, 0.0, 0.0],   :Y  => [0.5, 0.0, 0.0],   :T  => [0.5, 0.0, 0.5],
        :Z  => [0.0, 0.0, 0.5],   :S  => [0.25, 0.25, 0.0], :R  => [0.25, 0.25, 0.5],
        :Δ₀ => [0.0, ζ, 0.0],     :F₀ => [0.5, 0.5-ζ, 0.0], :B₀ => [0.0, ζ, 0.5],
        :G₀ => [0.5, 0.5-ζ, 0.5]
        #:T₂ => [0.5, 0.0, -0.5],  :Z₂ => [0.0, 0.0, -0.5],  :R₂ => [0.25, 0.25, -0.5],
        #:B₂ => [0.0, ζ, -0.5],    :G₂ => [0.5, 0.5-ζ, -0.5]
        )
    # Note: the (...)₂ points are not really needed, but are listed in SeeK anyway (in case,
    #       I think, of oC2 space groups without inversion). Here we assume inversion, so
    #       we don't need them, but we keep them as out-commented in case we ever will...

    pathtype == "SeeK" || _throw_pathtype(pathtype)
    paths_labs = [[:Γ, :Y, :F₀], [:Δ₀, :Γ, :Z, :B₀], [:G₀, :T, :Y], [:Γ, :S, :R, :Z, :T]]

    return lab2kv, paths_labs
end

# Table 71 of [SeeK]: extended Bravais type cI1 (with inversion)
function irrfbz_data_cI1(; pathtype::String="SeeK")
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :Γ => [0.0, 0.0, 0.0],   :H′ => [0.5, -0.5, 0.5], :N => [0.0, 0.0, 0.5], 
        :P => [0.25, 0.25, 0.25], :H => [0.5, 0.5, 0.5]
        # NB: There is some clashes between what the H-point is in SG 230. SeeK defines the
        #     H-point as [1,-1,1]/2, but bandreps and ISOTROPY defines it as [1,1,1]/2. 
        #     Notably, Bilbao itself is inconsistent in its naming; e.g. BANDREPS follows
        #     ISOTROPY's convention, but KVEC agrees with SeeK's convention. Here, we 
        #     follow ISOTROPY's convention, and so have renamed SeeK's H ⇒ H′.
        )

    if pathtype == "SeeK"
        paths_labs = [[:Γ, :H′, :N, :Γ, :P, :H′], [:P, :N]]
    elseif pathtype == "bandreps"
        paths_labs = [[:H′, :Γ, :P, :H, :N, :Γ], [:P, :N]]
    else
        _throw_pathtype(pathtype)
    end

    return lab2kv, paths_labs
end

# Table 72 of [Seek]: extended Bravais type tP1 (with inversion)
function irrfbz_data_tP1(; pathtype::String="SeeK")
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :A => [0.5, 0.5, 0.5], :Γ => [0.0, 0.0, 0.0], :M => [0.5, 0.5, 0.0], 
        :R => [0.0, 0.5, 0.5], :X => [0.0, 0.5, 0.0], :Z => [0.0, 0.0, 0.5]
        )

    pathtype == "SeeK" || _throw_pathtype(pathtype)
    paths_labs = [[:Γ, :X, :M, :Γ, :Z, :R, :A, :Z], [:X, :R], [:M, :A]]

    return lab2kv, paths_labs
end

# Table 73 (c < a): extended Bravais type tI1 (with inversion)
function irrfbz_data_tI1(a::Real, c::Real; pathtype::String="SeeK")
    c ≤ a || throw(DomainError((a,c), "requires c≤a"))

    η = (1.0 + (c/a)^2)/4
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :Γ => [0, 0, 0],          :M => [-0.5, 0.5, 0.5], :N => [0, 0.5, 0], 
        :P => [0.25, 0.25, 0.25], :X => [0, 0, 0.5],      :Z => [η, η, -η],
        :Z₀ => [-η, 1-η, η]
        )

    pathtype == "SeeK" || _throw_pathtype(pathtype)
    paths_labs = [[:Γ, :X, :M, :Γ, :Z], [:Z₀, :M], [:X, :P, :N, :Γ]]

    return lab2kv, paths_labs
end

# Table 74 (c > a): extended Bravais type tI2 (with inversion)
function irrfbz_data_tI2(a::Real, c::Real; pathtype::String="SeeK")
    c ≥ a || throw(DomainError((a,c), "requires c≥a"))

    η = (1.0 + (a/c)^2)/4
    ζ = ((a/c)^2)/2
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :Γ => [0, 0, 0],          :M => [0.5, 0.5, -0.5], :X => [0, 0, 0.5], 
        :P => [0.25, 0.25, 0.25], :N => [0, 0.5, 0],      :S₀ => [-η, η, η], 
        :S => [η, 1-η, -η],       :R => [-ζ, ζ, 0.5],     :G => [0.5, 0.5, -ζ]
        )

    pathtype == "SeeK" || _throw_pathtype(pathtype)
    paths_labs = [[:Γ, :X, :P, :N, :C, :M, :S], [:S₀, :Γ], [:X, :R], [:G, :M]]

    return lab2kv, paths_labs
end

# Table 84 ("long path"): extended Bravais type hP1
function irrfbz_data_hP1(; pathtype::String="SeeK")
    pathtype == "SeeK" || _throw_pathtype(pathtype)
    paths_labs = [[:Γ, :M, :K, :Γ, :A, :L, :H, :A], [:L, :M], [:H, :K, :H₂]]

    return _irrfbz_data_hP(), paths_labs
end

# Table 84 ("short path"): extended Bravais type hP2
function irrfbz_data_hP2(; pathtype::String="SeeK")
    pathtype == "SeeK" || _throw_pathtype(pathtype)
    paths_labs = [[:Γ, :M, :K, :Γ, :A, :L, :H, :A], [:L, :M], [:H, :K]] # doesn't include K-H₂ path as in hP2

    return _irrfbz_data_hP(), paths_labs
end

# Table 84 k-points only
function _irrfbz_data_hP()
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :Γ => [0, 0, 0],       :A => [0, 0, 1/2],      :K => [1/3, 1/3, 0],
        :H => [1/3, 1/3, 1/2], :H₂=> [1/3, 1/3, -1/2], :M => [1/2, 0, 0],
        :L => [1/2, 0, 1/2]
        ) # H₂ only needed for "long path" (hP1)
end

function irrfbz_data_tI1_without_inversion_or_tr(a::Real, c::Real; pathtype::String="SeeK")
    c ≤ a || throw(DomainError((a,c), "requires c≤a"))

    η = (1.0 + (c/a)^2)/4
    lab2kv, _ = irrfbz_data_tI1(a, c, pathtype="SeeK")
    # push negative-k (i.e. inverted) partners to `lab2kv`
    for (lab, kv) in copy(lab2kv)
        lab == :Γ && continue
        lab2kv[Symbol(lab, '′')] = -kv
    end

    if pathtype == "SeeK"
        paths_labs = [[:Γ, :X, :M, :Γ, :Z], [:Z₀, :M], [:X, :P, :N, :Γ, :X′, :M′, :Γ, :Z′], 
                      [:Z₀′, :M′], [:X′, :P′, :N′, :Γ]]
    elseif pathtype == "simple" # ±symmetric path (same length as "SeeK" path though)
        paths_labs = [[:Γ, :X, :M, :Γ, :Z], [:Z₀, :M], [:X, :P, :N, :Γ, :N′, :P′, :X′], 
                      [:M′, :Z₀′], [:Z′, :Γ, :M′, :X′, :Γ]]
    else
        _throw_pathtype(pathtype)
    end

    return lab2kv, paths_labs
end

# ---------------------------------------------------------------------------------------- #
_throw_pathtype(pathtype::String) = throw(DomainError(pathtype, "undefined pathtype"))
_throw_requires_direct_basis()    = throw(DomainError(nothing, "the k-path of the requested space group requires information about the direct basis; `Rs` must be of type `$(BasisLike{3})` here"))
# ---------------------------------------------------------------------------------------- #
end # module