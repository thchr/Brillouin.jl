module CrystallineBravaisVendor

using ..Brillouin: BasisLike
using StaticArrays

# ---------------------------------------------------------------------------------------- #
# this is just a raw output of `bravaistype` from Crystalline, to avoid depending on all
# of Crystalline just for `bravaistype`; basically just hardcoded "memoization"; note that
# we do *not* normalize "oA" to "oC" (i.e., this is `bravaistype(_, _, normalize=false)`)

@inline function boundscheck_sgnum(sgnum::Integer, D::Integer)
    if D == 3
        if sgnum ∉ 1:230
            throw(DomainError(sgnum, "sgnum must be between 1 and 230 in dimension 3"))
        end
    elseif D == 2
        if sgnum ∉ 1:17
            throw(DomainError(sgnum, "sgnum must be between 1 and 17 in dimension 2"))
        end
    elseif D == 1
        if sgnum ∉ 1:2
            throw(DomainError(sgnum, "sgnum must be between 1 and 2 in dimension 1"))
        end
    else
        throw(DomainError(D, "dimensions must be 1, 2, or 3"))
    end
end

const BRAVAISTYPE_2D = (
    "mp", "mp", "op", "op", "oc", "op", "op", "op", "oc", "tp", "tp", "tp", "hp", "hp",
    "hp", "hp", "hp"
)
const BRAVAISTYPE_3D = (
    "aP", "aP", "mP", "mP", "mC", "mP", "mP", "mC", "mC", "mP", "mP", "mC", "mP", "mP",
    "mC", "oP", "oP", "oP", "oP", "oC", "oC", "oF", "oI", "oI", "oP", "oP", "oP", "oP",
    "oP", "oP", "oP", "oP", "oP", "oP", "oC", "oC", "oC", "oA", "oA", "oA", "oA", "oF",
    "oF", "oI", "oI", "oI", "oP", "oP", "oP", "oP", "oP", "oP", "oP", "oP", "oP", "oP",
    "oP", "oP", "oP", "oP", "oP", "oP", "oC", "oC", "oC", "oC", "oC", "oC", "oF", "oF",
    "oI", "oI", "oI", "oI", "tP", "tP", "tP", "tP", "tI", "tI", "tP", "tI", "tP", "tP",
    "tP", "tP", "tI", "tI", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tI", "tI",
    "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tI", "tI", "tI", "tI", "tP", "tP",
    "tP", "tP", "tP", "tP", "tP", "tP", "tI", "tI", "tI", "tI", "tP", "tP", "tP", "tP",
    "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tI", "tI",
    "tI", "tI", "hP", "hP", "hP", "hR", "hP", "hR", "hP", "hP", "hP", "hP", "hP", "hP",
    "hR", "hP", "hP", "hP", "hP", "hR", "hR", "hP", "hP", "hP", "hP", "hR", "hR", "hP",
    "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP",
    "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "cP", "cF",
    "cI", "cP", "cI", "cP", "cP", "cF", "cF", "cI", "cP", "cI", "cP", "cP", "cF", "cF",
    "cI", "cP", "cP", "cI", "cP", "cF", "cI", "cP", "cF", "cI", "cP", "cP", "cP", "cP",
    "cF", "cF", "cF", "cF", "cI", "cI"
)
function bravaistype(sgnum::Integer, Dᵛ::Val{D}=Val(3)) where D
    @boundscheck boundscheck_sgnum(sgnum, D)
    if D == 3
        return BRAVAISTYPE_3D[sgnum]
    elseif D == 2
        return BRAVAISTYPE_2D[sgnum]
    else
        return "lp"
    end
end
bravaistype(sgnum::Integer, D::Integer) = bravaistype(sgnum, Val(D))

# ---------------------------------------------------------------------------------------- #
# vendor the tools from Crystalline needed to convert back and forth between
# conventional/reciprocal lattice bases

# ... relies on `bravaistype` *not* normalizing "oA" to "oC"
centering(sgnum::Integer, Dᵛ::Val{D}) where D = last(bravaistype(sgnum, Dᵛ))

const PRIMITIVE_BASIS_MATRICES = (
    # 1D
    Base.ImmutableDict('p'=>SMatrix{1,1,Float64}(1)),                # primitive
    # 2D
    Base.ImmutableDict('p'=>SMatrix{2,2,Float64}([1 0; 0 1]),        # primitive/simple
                  'c'=>SMatrix{2,2,Float64}([1 1; -1 1]./2)),   # centered      
    # 3D
    Base.ImmutableDict(
        'P'=>SMatrix{3,3,Float64}([1 0 0; 0 1 0; 0 0 1]),       # primitive/simple
        'F'=>SMatrix{3,3,Float64}([0 1 1; 1 0 1; 1 1 0]./2),    # face-centered
        'I'=>SMatrix{3,3,Float64}([-1 1 1; 1 -1 1; 1 1 -1]./2), # body-centered
        'R'=>SMatrix{3,3,Float64}([2 -1 -1; 1 1 -2; 1 1 1]./3), # rhombohedrally-centered
        'A'=>SMatrix{3,3,Float64}([2 0 0; 0 1 -1; 0 1 1]./2),   # base-centered (along x)
        'C'=>SMatrix{3,3,Float64}([1 1 0; -1 1 0; 0 0 2]./2))   # base-centered (along z)
    )
@inline function primitivebasismatrix(cntr::Char, ::Val{D}=Val(3)) where D
    D ∉ 1:3 && throw(DomainError(D, "dimension must be 1, 2, or 3"))
    return PRIMITIVE_BASIS_MATRICES[D][cntr]
end

function reciprocalbasis(Rs::BasisLike{D}) where D
    Rm = hcat(Rs...)
    Gm = 2π.*inv(transpose(Rm))
    return SVector{D}(ntuple(i->SVector{D,Float64}(Gm[:,i]), Val(D)))
end

function transform_reciprocal(Gs::BasisLike{D}, P::AbstractMatrix{<:Real}) where D
    Gm′ = hcat(Gs...)/P'
    return SVector{D}(ntuple(i->SVector{D}(Gm′[:,i]), Val(D)))
end

for (f, op_P) in ((:primitivize_reciprocal, :identity), (:conventionalize_reciprocal, :inv))
    @eval function $f(basis::BasisLike{D}, sgnum::Integer) where D
        cntr = centering(sgnum, Val(D))
        return $f(basis, cntr)
    end
    @eval function $f(Gs::BasisLike{D}, cntr::Char) where D
        if cntr == 'P' || cntr == 'p' # the conventional and primitive bases coincide
            return Gs
        else         
            P = primitivebasismatrix(cntr, Val(D))        
            return transform_reciprocal(Gs, $op_P(P))
        end
    end
end

end # module