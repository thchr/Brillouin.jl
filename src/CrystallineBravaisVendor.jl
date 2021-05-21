module CrystallineBravaisVendor
# this is just a raw output of `bravaistype` from Crystalline, to avoid depending on all
# of Crystalline just for `bravaistype`; basically just hardcoded "memoization"; note that
# we do not normalize "oA" to "oC" (i.e., this is `bravaistype(_, _, normalize=false)`)

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

end # module