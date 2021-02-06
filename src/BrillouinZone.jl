
module BrillouinZone
using Crystalline: DirectBasis, ReciprocalBasis, reciprocalbasis, bravaistype
using LinearAlgebra: norm

export brillouin_zone

"""
    brillouin_zone(sgnum::Integer, Rs::DirectBasis{D}) where D

Takes a **conventional** unit cell basis; returns the facets of the associated Brillouin 
zone in a **primitive** basis. ITA conventions are assumed.

If `sgnum == 0`, a cubic Brillouin zone is returned.
"""
function brillouin_zone(sgnum::Integer, Rs::DirectBasis{D}) where D
    D ≠ 3 && throw(DomainError(D, "Currently only implemented for 3D space groups"))

    sgnum == 0 && return bz_data_cube()

    bt = bravaistype(sgnum, D)

    # get the polygonal BZ "faces" from manual tables...
    # TODO: Really should be doing this in a general fashion, with Voronoi cells etc, as in
    #       Seek's implementation (see https://github.com/giovannipizzi/seekpath/blob/3fe68e7884eafa1f693bcfdd3cb4715b9765174e/seekpath/brillouinzone/brillouinzone.py#L37)
    #       The current "write-down-the-bz-faces-manually" approach does _not_ scale well...
    #       An ideal situation would be to use a wrapper of e.g. Voro++.
    if bt == "aP"
        # sgs 1 and 2
        lab2kv, cntcs = bz_data_aP1()

    elseif bt == "tP"
        # sgs 75:78, 81, 83:86, 89:106, 111:118, 123:138
        lab2kv, cncts = bz_data_tP1()

    elseif bt == "tI"
        # sgs 79, 80, 82, 87, 88, 97, 98, 107, 108:110, 119:122, 139:142
        a, c = norm(Rs[1]), norm(Rs[3])
        if c ≤ a     # tI1
            lab2kv, cncts = bz_data_tI1(a, c)
        else #c > a  # tI2
            throw(DomainError(bt, "The requested Bravais type has not yet been implemented"))
    #        lab2kv, cncts = bz_data_tI2(a, c)
        end

    elseif bt == "cI"
        # sgs 197, 199, 204, 206, 211, 214, 217, 220, 229, 230
        lab2kv, cncts = bz_data_cI();

    elseif bt == "oC"
        # sgs 20, 21, 35:41, 63:68
        a, b = norm(Rs[1]), norm(Rs[2])
        if a ≤ b     # oC1
            lab2kv, cncts = bz_data_oC1(a, b)
        else #a > b  # oC2
            lab2kv, cncts = bz_data_oC2(a, b)
        end

    elseif bt == "hP"
        # sgs 143:145, 147, 149:154, 156:159, 162:165, 168:194
        lab2kv, cncts = bz_data_hP()

    else
        throw(DomainError(bt, "The requested Bravais type has not yet been implemented"))
    end

    return lab2kv, cncts
end

# NOTE: All points with subscripts not equal to 0 or 2 or with primes after, are not 
#       following any conventional labelling but are just ad-hoc labels developed to keep 
#       sane here: they are derived from their equivalent 0 or 2 subscripted versions.
function bz_data_oC1(a::Real, b::Real)
    a≤b || throw(DomainError((a,b), "requires a≤b"))

    ζ = (1.0 + (a/b)^2)/4
    ξ = 1/2 - ζ
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :A₀  => [ζ, ζ, 1/2],                :A₂  => [ζ, ζ, -1/2],
        :A₀′ => [-ζ, -ζ, 1/2],              :A₂′ => [-ζ, -ζ, -1/2],
        :E₀  => [-1/2+ξ, 1/2+ξ, 1/2],       :E₂  => [-1/2+ξ, 1/2+ξ, -1/2],
        :E₁  => [-1/2-ξ, 1/2-ξ, 1/2],       :E₃  => [-1/2-ξ, 1/2-ξ, -1/2],
        :E₀′ => [-(-1/2+ξ), -(1/2+ξ), 1/2], :E₂′ => [-(-1/2+ξ), -(1/2+ξ), -1/2],
        :E₁′ => [-(-1/2-ξ), -(1/2-ξ), 1/2], :E₃′ => [-(-1/2-ξ), -(1/2-ξ), -1/2]
        )
    cncts = [[:A₀, :E₀, :E₁, :A₀′, :E₀′, :E₁′, :A₀], # tops
             [:A₂, :E₃′, :E₂′, :A₂′, :E₃, :E₂, :A₂],
             [:E₀, :E₁, :E₃, :E₂],                   # sliver sides
             [:E₀′, :E₂′, :E₃′, :E₁′],
             [:A₀, :A₂, :E₂, :E₀,    :A₀],           # big sides
             [:A₀′, :E₁, :E₃, :A₂′,  :A₀′],
             [:A₀′,  :A₂′, :E₂′, :E₀′, :A₀′],
             [:A₀, :E₁′, :E₃′, :A₂,  :A₀]]

    return lab2kv, cncts
end

function bz_data_oC2(a::Real, b::Real)
    a≥b || throw(DomainError((a,b), "requires a≥b"))

    ζ = (1.0 + (b/a)^2)/4
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :T  => [1/2, 0.0, 1/2],   :T₂ => [1/2, 0.0, -1/2],
        :B₀ => [0.0, ζ, 1/2],     :B₂ => [0.0, ζ, -1/2],
        :G₀ => [1/2, 1/2-ζ, 1/2], :G₂ => [1/2, 1/2-ζ, -1/2]
        # TODO: Remaining points
        )
    throw(DomainError("oC2", "The requested Bravais type has not yet been implemented"))
    # TODO: cncts
end

function bz_data_tP1()
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :A  => [1/2, 1/2, 1/2],     :A₁ => [-1/2, 1/2, 1/2],   # positive "z"
        :A₂ => [-1/2, -1/2, 1/2],   :A₃ => [1/2, -1/2, 1/2],   
        :Aᶻ  => [1/2, 1/2, -1/2],   :A₁ᶻ => [-1/2, 1/2, -1/2], # negative "z"
        :A₂ᶻ => [-1/2, -1/2, -1/2], :A₃ᶻ => [1/2, -1/2, -1/2])
        
    cncts = [[:A,  :A₁,  :A₂,  :A₃,  :A],   # tops
             [:Aᶻ, :A₃ᶻ, :A₂ᶻ, :A₁ᶻ, :Aᶻ],
             [:A,  :Aᶻ,  :A₁ᶻ, :A₁,  :A],   # sides
             [:A₁, :A₁ᶻ, :A₂ᶻ, :A₂,  :A₁],
             [:A₂, :A₂ᶻ, :A₃ᶻ, :A₃,  :A₂],
             [:A₃, :A₃ᶻ, :Aᶻ,  :A,   :A₃]]

    return lab2kv, cncts
end

function bz_data_tI1(a, c)
    c ≤ a || throw(DomainError((a,c), "requires c≤a"))

    η = (1.0 + (c/a)^2)/4
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :P   => [1/4, 1/4, 1/4],    :P⁻   => [-1/4, -1/4, 3/4],   # "up G₃"
        :Z₀  => [-η, 1-η, η],       :Z₀⁻  => [-1+η, η, 1-η],
        :Z₁  => [1-η, -η, η],       :Z₁⁻  => [η, -1+η, 1-η],
        :Pᶻ  => [1/4, 1/4, 1/4-1],  :P⁻ᶻ  => [-1/4, -1/4, 3/4-1], # "down G₃"
        :Z₀ᶻ => [-η, 1-η, η-1],     :Z₀⁻ᶻ => [-1+η, η, 1-η-1],
        :Z₁ᶻ => [1-η, -η, η-1],     :Z₁⁻ᶻ => [η, -1+η, 1-η-1],
        :Z   => [η, η, -η],         :Z⁻   => [-η, -η, η],         # "z-tips"
        :NP₀  => [-1/4, 3/4, -1/4], :NP₁  => [3/4, -1/4, -1/4],   # weird corner points
        :NP₀⁻  => [-3/4, 1/4, 1/4], :NP₁⁻  => [1/4, -3/4, 1/4]    
        )
    # :X => [0, 0, 1/2],  :M => [-1/2, 1/2, 1/2], 
    cncts = [[:P, :Z₀, :Z₀⁻, :P⁻, :Z₁⁻, :Z₁, :P],
             [:Pᶻ, :Z₁ᶻ, :Z₁⁻ᶻ, :P⁻ᶻ, :Z₀⁻ᶻ, :Z₀ᶻ, :Pᶻ],
             [:Z₀, :NP₀, :Z₀ᶻ, :Z₀⁻ᶻ, :NP₀⁻, :Z₀⁻, :Z₀],
             [:Z₁, :Z₁⁻, :NP₁⁻, :Z₁⁻ᶻ, :Z₁ᶻ, :NP₁, :Z₁],
             [:Z, :NP₀, :Z₀, :P, :Z], [:Z, :Pᶻ, :Z₀ᶻ, :NP₀, :Z],
             [:Z, :P, :Z₁, :NP₁, :Z], [:Z, :NP₁, :Z₁ᶻ, :Pᶻ, :Z],
             [:Z⁻, :P⁻, :Z₀⁻, :NP₀⁻, :Z⁻], [:Z⁻, :NP₀⁻, :Z₀⁻ᶻ, :P⁻ᶻ, :Z⁻],
             [:Z⁻, :NP₁⁻, :Z₁⁻, :P⁻, :Z⁻], [:Z⁻, :P⁻ᶻ, :Z₁⁻ᶻ, :NP₁⁻, :Z⁻],
            ]
    return lab2kv, cncts
end

function bz_data_hP()
    lab2kv = Dict{Symbol, Vector{Float64}}( # NB: our subscript-labeling is nonconventional; does not match that kpath notation
        :H₁  => [1/3, 1/3, 1/2],   :H₂  => [2/3, -1/3, 1/2],  # positive "z"
        :H₃  => [1/3, -2/3, 1/2],  :H₄  => [-1/3, -1/3, 1/2],  
        :H₅  => [-2/3, 1/3, 1/2],  :H₆  => [-1/3, 2/3, 1/2],
        :H₁ᶻ => [1/3, 1/3, -1/2],  :H₂ᶻ => [2/3, -1/3, -1/2], # negative "z"
        :H₃ᶻ => [1/3, -2/3, -1/2], :H₄ᶻ => [-1/3, -1/3, -1/2],  
        :H₅ᶻ => [-2/3, 1/3, -1/2], :H₆ᶻ => [-1/3, 2/3, -1/2]
        )
    cncts = [[:H₁,  :H₂,  :H₃,  :H₄,  :H₅,  :H₆,  :H₁],  # top
             [:H₁ᶻ, :H₂ᶻ, :H₃ᶻ, :H₄ᶻ, :H₅ᶻ, :H₆ᶻ, :H₁ᶻ], # bottom
             [:H₁,  :H₂,  :H₂ᶻ, :H₁ᶻ, :H₁],              # sides
             [:H₂,  :H₃,  :H₃ᶻ, :H₂ᶻ, :H₂],
             [:H₃,  :H₄,  :H₄ᶻ, :H₃ᶻ, :H₃],
             [:H₄,  :H₅,  :H₅ᶻ, :H₄ᶻ, :H₄],
             [:H₅,  :H₆,  :H₆ᶻ, :H₅ᶻ, :H₅],
             [:H₆,  :H₁,  :H₁ᶻ, :H₆ᶻ, :H₆]]

    return lab2kv, cncts
end

function bz_data_cI()
    lab2kv = Dict(
        :H′ => [1/2, -1/2, 1/2],   :H′ʸ => [-1/2, 1/2, -1/2],  
        :P₁ => [1/4, 1/4, 1/4],    :P₁ʸ => [-1/4, 3/4, -1/4],   :C₁ => [-1/2, 1/2, 1/2],
        :P₂ => [-1/4, -1/4, 3/4],  :P₂ʸ => [-3/4, 1/4, 1/4],    :C₂ => [-1/2, -1/2, 1/2],
        :P₃ => [1/4, -3/4, 1/4],   :P₃ʸ => [-1/4, -1/4, -1/4],  :C₃ => [1/2, -1/2, -1/2],
        :P₄ => [3/4, -1/4, -1/4],  :P₄ʸ => [1/4, 1/4, -3/4],    :C₄ => [1/2, 1/2, -1/2],
    )

    cncts = [[:P₁, :H′, :P₄, :C₄], [:P₂, :H′, :P₁, :C₁],        # tops (along cartesian y-axis)
             [:P₃, :H′, :P₂, :C₂], [:P₄, :H′, :P₃, :C₃],
             [:P₁ʸ, :C₄, :P₄ʸ, :H′ʸ], [:P₂ʸ, :C₁, :P₁ʸ, :H′ʸ],  # bottoms
             [:P₃ʸ, :C₂, :P₂ʸ, :H′ʸ], [:P₄ʸ, :C₃, :P₃ʸ, :H′ʸ],
             [:P₁, :C₄, :P₁ʸ, :C₁], [:P₂, :C₁, :P₂ʸ, :C₂],      # sides
             [:P₃, :C₂, :P₃ʸ, :C₃], [:P₄, :C₃, :P₄ʸ, :C₄]
            ]
    
    return lab2kv, cncts
end


# a perfect cube, useful for plotting in the "paralleliped BZ"
function bz_data_cube()
    lab2kv = Dict{Symbol, Vector{Float64}}(
        :A => [-1/2, -1/2, 1/2],  :B => [1/2, -1/2, 1/2],   # positive "z"
        :C => [1/2, 1/2, 1/2],    :D => [-1/2, 1/2, 1/2],   
        :E => [-1/2, -1/2, -1/2], :F => [-1/2, 1/2, -1/2],  # negative "z"
        :G => [1/2, 1/2, -1/2],   :H => [1/2, -1/2, -1/2])
        
    cncts = [[:A, :B, :C, :D, :A],   # z-sides (first z>0, then z<0)
             [:E, :F, :G, :H, :E],
             [:C, :D, :F, :G, :C],   # y-sides (first y>0, then y<0)
             [:A, :B, :H, :E, :A],
             [:B, :C, :G, :H, :B],   # x-sides (first x>0, then x<0)
             [:A, :D, :F, :E, :A]]

    return lab2kv, cncts
end

end

##

# --- visualization example ---
# using PyPlot
# using3D()
#
# a, b = 1.0, 1.1; Rs = DirectBasis([a,0,0], [0,b,0], [0,0,1.0])
# Gs = reciprocalbasis(Rs)
# lab2kv, cncts=Main.BrillouinZone.brillouin_zone(68, Rs)
#
# close("all")
# ax=figure().gca(projection="3d")
#
# for (key,val) in lab2kv
#     kv = sum(val.*Gs)
#     ax.plot(kv[1:1], kv[2:2], kv[3:3], "o")
#     println(string(key))
#     ax.text(kv[1], kv[2], kv[3], string(key))
# end
# 
# for cnct in cncts
#     kvpath = [sum(lab2kv[lab].*Gs) for lab in cnct]
#     plot(getindex.(kvpath, 1), getindex.(kvpath, 2), getindex.(kvpath, 3), "-k")
# end
# xlim(-1/2,1/2), ylim(-1/2,1/2), zlim(-1/2,1/2)