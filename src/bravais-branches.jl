@noinline function _throw_conflicting_sgnum_and_bravais(bt, sgnum)
    throw(DomainError((bt, sgnum), 
            "provided bravais type and space group number are mutually inconsistent"))
end
@noinline function _warn_monoclinic()
    @warn """
          SeeK has extra conditions on the parameters of monoclinic systems, beyond
          those of ITA; please check that you fulfil these conditions by comparing
          with the HPKOT paper
          """
end

# this is a translation and simplification (omitting several warnings / edge case handling)
# of SeeK's `get_path` extended bravais type branch table. It should be pretty much the
# same as Table 94 of the HPKOT paper (i.e. paper "behind" SeeK). It simply returns the
# extended Bravais symbol as a `Symbol`. Note the input direct basis is expected to
# be a conventional cell in ITA settings.
function extended_bravais(sgnum::Integer,
                          bt::String,
                          Rs::AbstractVector{<:SVector{3,<:Real}})
    a = norm(Rs[1])
    b = norm(Rs[2]) 
    c = norm(Rs[3])

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
        if c ≤ a
            return :tI1
        else
            return :tI2
        end

    elseif bt == "oP"
        return :oP1

    elseif bt == "oF"
        if inv(a^2) > inv(b^2) + inv(c^2)
            return :oF1
        elseif inv(c^2) > inv(a^2) + inv(b^2)
            return :oF2
        else
            return :oF3
        end

    elseif bt == "oI"
        if c>a && c>b     # (c largest)
            return :oI1
        elseif a>b && a>c # (a largest)
            return :oI2
        else # b≥a && b≥c   (b largest, or equal)
            return :oI3
        end

    elseif bt == "oC"
        if a ≤ b
            return :oC1
        else
            return :oC2
        end

    elseif bt == "oA"
        if b <= c
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
        if sqrt(3)a ≤ sqrt(2)c
            return :hR1
        else
            return :hR2
        end

    elseif bt == "mP"
        # TODO
        _warn_monoclinic()
        return :mP1

    elseif bt == "mC"
        # TODO
        _warn_monoclinic()
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
        # TODO
        error(DomainError("aP", "bravais type is not presently supported; try SeeK-path"))

    else
        throw(DomainError(bt, "undefined bravais type"))

    end
end # function