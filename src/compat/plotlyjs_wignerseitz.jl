using .PlotlyJS
import .PlotlyJS: plot

# ---------------------------------------------------------------------------------------- #
# CONSTANTS

# default layout
const DEFAULT_PLOTLY_LAYOUT  = Layout(
    showlegend=false,
    scene=attr(
        xaxis=attr(tickvals=[], zeroline=false,
                showgrid=false, showbackground=false,
                title=attr(text="")#"k<sub>x</sub>")
                ),
        yaxis=attr(tickvals=[], zeroline=false,
                showgrid=false, showbackground=false,
                title=attr(text="")#"k<sub>y</sub>")
                ),
        zaxis=attr(tickvals=[], zeroline=false,
                showgrid=false, showbackground=false,
                title=attr(text=""),#"k<sub>z</sub>")   
                ),
        aspectmode = "data",
        camera=attr(up = attr(x=0, z=1, y=0),),
        ),
    margin=attr(l=0, r=0, b=0, t=0),
    autosize=false,
    #width=200, height=200,
    plot_bgcolor="rgba(0, 0, 0, 0)", paper_bgcolor="rgba(0, 0, 0, 0)",
    )

# colors from the "british" flatcolors color palette: https://flatuicolors.com/palette/gb
const BZ_COL          = "rgb(47,54,64)"    # "electromagnetic"
const BASIS_COL       = "rgb(39,60,117)"   # "pico void"
const BASIS_LIGHT_COL = "rgb(212,216,227)" # 20% BASIS_COL, 80% white
const AXIS_COL        = "rgb(194,54,22)"   # "harley davidson orange"
const AXIS_LIGHT_COL  = "rgb(242,215,208)" # 20% AXIS_COL, 80% white

# ---------------------------------------------------------------------------------------- #

function plot(c::Cell{3}, layout::Layout=DEFAULT_PLOTLY_LAYOUT)
    scale = maximum(norm, basis(c))

    # BZ
    tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(c))
    for (i,poly) in enumerate(c)
        tbz[i] = PlotlyJS.scatter3d(
            x=push!(getindex.(poly, 1), poly[1][1]),
            y=push!(getindex.(poly, 2), poly[1][2]),
            z=push!(getindex.(poly, 3), poly[1][3]); 
            mode="lines", hovertext="BZ", hoverinfo="text+x+y+z",
            line=attr(color=BZ_COL, width=3)
            )
    end

    # lattice vectors
    tgs    = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 6)
    tgtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 3)
    intersects = axis_intersections(c, basis(c))
    for (i,V) in enumerate(basis(c))
        V′ = V./norm(V)
        name = "<b>v</b><sub>$(i)</sub>"
        for j in (1,2) # in/outside BZ
            start = j == 1 ? 0.0           : intersects[i]
            stop  = j == 1 ? intersects[i] : 1.0
            V₀ = start*V .+ (j == 1 ? 0.0   : 0.025)*scale*V′
            V₁ = stop*V  .- (j == 1 ? 0.025 : 0.0)  *scale*V′
            tgs[i+(j-1)*3] = PlotlyJS.scatter3d(
                x=[V₀[1],V₁[1]], y=[V₀[2],V₁[2]], z=[V₀[3],V₁[3]];
                mode="lines", hovertext=name, hoverinfo="text",
                line=attr(color=ifelse(j==1, BASIS_LIGHT_COL, BASIS_COL),
                          width=ifelse(j==1, 4, 6))
                )
        end
        tgtips[i] = PlotlyJS.cone(
            x=[V[1]], u=[V′[1]], y=[V[2]], v=[V′[2]], z=[V[3]], w=[V′[3]],
            sizeref=.1*scale, showscale=false, anchor="tail", 
            colorscale=[[0, BASIS_LIGHT_COL], [1, BASIS_COL]],
            hovertext=name, hoverinfo="text+x+y+z")
    end

    # Cartesian axes
    cartVs = cartesian_axes(Val(3))
    name_axs = ("x", "y", "z") 
    taxs    = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 6)
    taxtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 3)
    intersects = axis_intersections(c, cartVs)
    for (i, V) in enumerate(cartVs)
        name = "<b>"*string(name_axs[i])*"</b>"
        V′ = V./norm(V)
        for j in (1,2) # in/outside BZ
            start = j == 1 ? 0.0           : intersects[i]
            stop  = j == 1 ? intersects[i] : intersects[i]
            V₀ = start*V .+ (j == 1 ? 0.0   : 0.025)*scale*V′
            V₁ = stop*V  .- (j == 1 ? 0.025 : -0.2) *scale*V′
            taxs[i+(j-1)*3] = PlotlyJS.scatter3d(
                x=[V₀[1],V₁[1]], y=[V₀[2],V₁[2]], z=[V₀[3],V₁[3]];
                mode="lines", hovertext=name, hoverinfo="text",
                line=attr(color=ifelse(j==1, AXIS_LIGHT_COL, AXIS_COL),
                          width=ifelse(j==1, 4, 6))
                )
            if j == 2
                taxtips[i] = PlotlyJS.cone(
                    x=[V₁[1]], u=[V[1]],
                    y=[V₁[2]], v=[V[2]],
                    z=[V₁[3]], w=[V[3]],
                    sizeref=.1*scale, showscale=false, anchor="tail", 
                    colorscale=[[0, AXIS_LIGHT_COL], [1, AXIS_COL]],
                    hovertext=name, hoverinfo="text")
            end
        end
    end

    # combine traces and plot
    ts = vcat(tbz, tgs, tgtips, taxs, taxtips)
    p = PlotlyJS.Plot(ts, layout)
    P = PlotlyJS.plot(p)
    PlotlyJS.display_blink(P)

    return P, ts
end

# ---------------------------------------------------------------------------------------- #
# UTILITIES 

# the "outward" intersections of of lines with direction `Vs` though origo with a cell
function axis_intersections(c::Cell{3},
            Vs::AbstractVector{<:AbstractVector}=cartesian_axes(Val(D)))

    intersects = MVector{3,Float64}(undef)
    fill!(intersects, Inf)
    # intersection between plane (r-rₒ) ⋅ n = 0 and line r(t) = l₀ + lt:
    #   t = (r₀ - l₀) ⋅ n / (l ⋅ n)     [https://wikipedia.org/wiki/Line–plane_intersection]
    for (i,poly) in enumerate(c)
        r₀ = sum(poly)/length(poly)
        n  = face_normal(c, i)
        for (j,V) in enumerate(Vs)
            t  = dot(r₀, n)/dot(V, n)
            if t > 0.0 && t < intersects[j]
                intersects[j] = t
            end
        end
    end
    return intersects
end

# a set of cartesian basis vectors in D dimensions
cartesian_axes(Dᵛ::Val{D}) where D = SVector(ntuple(i->SVector(ntuple(j->i==j ? 1.0 : 0.0, Dᵛ)), Dᵛ))

# ---------------------------------------------------------------------------------------- #