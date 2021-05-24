@noinline function _throw_conflicting_sgnum_and_bravais(bt, sgnum)
    throw(DomainError((bt, sgnum), 
            "provided bravais type and space group number are mutually inconsistent"))
end
@noinline function _throw_triclinic_angles(Rs)
    throw(DomainError(Rs, "Triclinic Bravais lattice must be specificied with a "*
                          "lattice that is either all-acute or all-obtuse; provided "*
                          "lattice system is neither"))
end
@noinline function _throw_basis_required(Rs)
    throw(DomainError(Rs, "`Rs` must be supplied for the considered Bravais type "))
end

# --- 3D ---
# this is a translation and simplification (omitting several warnings / edge case handling)
# of SeeK's `get_path` extended bravais type branch table. It should be pretty much the
# same as Table 94 of the HPKOT paper (i.e. paper "behind" SeeK). It simply returns the
# extended Bravais symbol as a `Symbol`. Note the input direct basis is expected to
# be a conventional cell in ITA settings.
function extended_bravais(sgnum::Integer,
                          bt::String,
                          Rs::Union{Nothing, StaticVector{3, <:StaticVector{3, <:Real}}},
                          Dᵛ::Val{3})

    @boundscheck boundscheck_sgnum(sgnum, 3)

    if bt == "cP"
        if 195 ≤ sgnum ≤ 206
            return :cP1
        elseif 207 ≤ sgnum ≤ 230
            return :cP2
        else
            _throw_conflicting_sgnum_and_bravais(bt, sgnum)
        end

    elseif bt == "cF"
        if 195 ≤ sgnum ≤ 206
            return :cF1
        elseif 207 ≤ sgnum ≤ 230
            return :cF2
        else
            _throw_conflicting_sgnum_and_bravais(bt, sgnum)
        end

    elseif bt == "cI"
        return :cI1

    elseif bt == "tP"
        return :tP1

    elseif bt == "tI"
        a, _, c = basisnorms(Rs)
        if c ≤ a
            return :tI1
        else
            return :tI2
        end

    elseif bt == "oP"
        return :oP1

    elseif bt == "oF"
        a, b, c = basisnorms(Rs)
        if inv(a^2) > inv(b^2) + inv(c^2)
            return :oF1
        elseif inv(c^2) > inv(a^2) + inv(b^2)
            return :oF2
        else
            return :oF3
        end

    elseif bt == "oI"
        a, b, c = basisnorms(Rs)
        if c>a && c>b     # (c largest)
            return :oI1
        elseif a>b && a>c # (a largest)
            return :oI2
        else # b≥a && b≥c   (b largest, or equal)
            return :oI3
        end

    elseif bt == "oC"
        a, b, _ = basisnorms(Rs)
        if a ≤ b
            return :oC1
        else
            return :oC2
        end

    elseif bt == "oA"
        _, b, c = basisnorms(Rs)
        if b ≤ c
            return :oA1
        else
            return :oA2
        end

    elseif bt == "hP"
        if sgnum in (143:149..., 151, 153, 157, 159:163...)
            return :hP1
        else
            return :hP2
        end

    elseif bt == "hR"
        a, _, c = basisnorms(Rs)
        if sqrt(3)a ≤ sqrt(2)c
            return :hR1
        else
            return :hR2
        end

    elseif bt == "mP"
        return :mP1

    elseif bt == "mC"
        a, b, c = basisnorms(Rs)
        cosβ = dot(Rs[3], Rs[1])/(c*a)
        if b < a * sqrt(1 - cosβ^2)
            return :mC1
        else
            if -a*cosβ/c + a^2*(1 - cosβ^2) / b^2 ≤ 1
                return :mC2
            else
                return :mC3
            end
        end

    elseif bt == "aP"      
        # SeeK-path only provides paths for all-acute or all-obtuse reciprocal aP cells
        # (more precisely, SeeK-path will transform to this setting automatically) so check:

        # get reciprocal basis inter-angles angles
        Gs = (Rs[2]×Rs[3], Rs[3]×Rs[1], Rs[1]×Rs[2]) # NB: "scale" prefactor omitted here
        aᴳ, bᴳ, cᴳ = basisnorms(Gs)
        cosαᴳ = dot(Gs[2], Gs[3]) / (bᴳ * cᴳ)
        cosβᴳ = dot(Gs[1], Gs[3]) / (aᴳ * cᴳ)
        cosγᴳ = dot(Gs[1], Gs[2]) / (aᴳ * bᴳ)

        # detect whether all-obtuse or all-acute reciprocal basis
        if cosαᴳ ≤ 0 && cosβᴳ ≤ 0 && cosγᴳ ≤ 0     # all-obtuse
            return :aP2
        elseif cosαᴳ ≥ 0 && cosβᴳ ≥ 0 && cosγᴳ ≥ 0 # all-acute
            return :aP3
        else
            _throw_triclinic_angles(Rs)
        end

    else
        throw(DomainError(bt, "undefined bravais type"))

    end
end # function

# --- 2D ---
function extended_bravais(sgnum::Integer,
                          bt::String,
                          Rs::Union{Nothing, StaticVector{2, <:StaticVector{2, <:Real}}},
                          Dᵛ::Val{2})
    @boundscheck boundscheck_sgnum(sgnum, 2)

    if bt == "hp"
        return :hp1
    elseif bt == "tp"
        return :tp1
    elseif bt == "op"
        if sgnum ∈ (3,4)
            return :op1
        elseif sgnum ∈ (6,7,8)
            return :op2
        else
            _throw_conflicting_sgnum_and_bravais(bt, sgnum)
        end
    elseif bt == "oc"
        a, b = basisnorms(Rs)
        if a ≤ b
            return :oc1
        else # a>b
            return :oc2
        end
    elseif bt == "mp"
        a, b = basisnorms(Rs)
        cosα = dot(Rs[1], Rs[2]) / (a * b)
        if cosα ≥ 0 # ∠(R₁, R₂) ≤ 90
            return :mp1
        else        # ∠(R₁, R₂) > 90
            return :mp2
        end
    else
        throw(MethodError("dimension 2 not currently implemented"))
    end
end

# --- 1D ---
function extended_bravais(sgnum::Integer,
                          bt::String,
                          Rs::Union{Nothing, StaticVector{1, <:StaticVector{1, <:Real}}},
                          Dᵛ::Val{1})

    @boundscheck boundscheck_sgnum(sgnum, 1)
    return :lp # trivial case; no "extended" bravais types
end

function basisnorms(Rs::Union{AVec{V}, NTuple{3, V}}) where V <: StaticVector{3, <:Real}
    a = norm(Rs[1])
    b = norm(Rs[2])
    c = norm(Rs[3])
    return a, b, c
end
function basisnorms(Rs::Union{AVec{V}, NTuple{2, V}}) where V <: StaticVector{2, <:Real}
    a = norm(Rs[1])
    b = norm(Rs[2])
    return a, b
end
basisnorms(Rs::Nothing) = _throw_basis_required(Rs)