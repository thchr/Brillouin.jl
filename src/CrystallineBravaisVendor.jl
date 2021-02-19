module CrystallineBravaisVendor
# this is just a raw output of `bravaistype` from Crystalline, to avoid depending on all
# of Crystalline just for `bravaistype`; basically just hardcoded "memoization"

const BRAVAISTYPE_2D = (
    "mp", "mp", "op", "op", "oc", "op", "op", "op", "oc", "tp", "tp", "tp", "hp", "hp",
    "hp", "hp", "hp"
)
const BRAVAISTYPE_3D = (
    "aP", "aP", "mP", "mP", "mC", "mP", "mP", "mC", "mC", "mP", "mP", "mC", "mP", "mP",
    "mC", "oP", "oP", "oP", "oP", "oC", "oC", "oF", "oI", "oI", "oP", "oP", "oP", "oP",
    "oP", "oP", "oP", "oP", "oP", "oP", "oC", "oC", "oC", "oC", "oC", "oC", "oC", "oF",
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
function bravaistype(sgnum::Integer, D::Integer=3)
    D == 3 && return BRAVAISTYPE_3D[sgnum]
    D == 2 && return BRAVAISTYPE_2D[sgnum]
    D == 1 && return "lp"
end

end # module