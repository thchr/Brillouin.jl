using .PlotlyJS
import .PlotlyJS: plot
using LinearAlgebra: norm, dot
using StaticArrays

using .WignerSeitz: face_normal
# ---------------------------------------------------------------------------------------- #
# CONSTANTS

# default layout
const DEFAULT_PLOTLY_LAYOUT_3D  = Layout(
    showlegend=false,
    scene=attr(
        xaxis=attr(tickvals=[], zeroline=false,
                showgrid=false, showbackground=false,
                title=attr(text=""),
                ),
        yaxis=attr(tickvals=[], zeroline=false,
                showgrid=false, showbackground=false,
                title=attr(text=""),
                ),
        zaxis=attr(tickvals=[], zeroline=false,
                showgrid=false, showbackground=false,
                title=attr(text=""),
                ),
        aspectmode = "data",
        camera=attr(up = attr(x=0, z=1, y=0), center = attr(x=0, y=0, z=0)),
        dragmode = "turntable",
        ),
    margin=attr(l=0, r=0, b=0, t=0),
    autosize=false,
    #width=200, height=200,
    plot_bgcolor="rgba(255, 255, 255, 1)", paper_bgcolor="rgba(255, 255, 255, 1)",
    )

# colors from the "british" flatcolors color palette: https://flatuicolors.com/palette/gb
const BZ_COL          = Ref("rgb(47,54,64)")    # "electromagnetic"
const BASIS_COL       = Ref("rgb(39,60,117)")   # "pico void"
const BASIS_LIGHT_COL = Ref("rgb(212,216,227)") # 20% BASIS_COL, 80% white
const AXIS_COL        = Ref("rgb(194,54,22)")   # "harley davidson orange"
const AXIS_LIGHT_COL  = Ref("rgb(242,215,208)") # 20% AXIS_COL, 80% white

# ---------------------------------------------------------------------------------------- #
# 3D

function plot(c::Cell{3}, layout::Layout=DEFAULT_PLOTLY_LAYOUT_3D;
              config::PlotConfig = PlotConfig(responsive=true, displaylogo=false))

    c.basisenum[] !== CARTESIAN && (c = cartesianize(c))
    scale = maximum(norm, basis(c))

    # BZ
    tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(c))
    for (i,poly) in enumerate(c)
        tbz[i] = PlotlyJS.scatter3d(
            x=push!(getindex.(poly, 1), poly[1][1]),
            y=push!(getindex.(poly, 2), poly[1][2]),
            z=push!(getindex.(poly, 3), poly[1][3]); 
            mode="lines", hovertext="Cell", hoverinfo="text+x+y+z",
            line=attr(color=BZ_COL[], width=3)
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
                line=attr(color=ifelse(j==1, BASIS_LIGHT_COL[], BASIS_COL[]),
                          width=ifelse(j==1, 5, 6))
                )
        end
        tgtips[i] = PlotlyJS.cone(
            x=[V[1]], u=[V′[1]], y=[V[2]], v=[V′[2]], z=[V[3]], w=[V′[3]],
            sizeref=.1*scale, showscale=false, anchor="tail",
            colorscale=[[0, BASIS_COL[]], [1, BASIS_COL[]]],
            hovertext=name, hoverinfo="text+x+y+z")
    end

    # Cartesian axes
    cart_basis = cartesian_axes(Val(3))
    name_axs = ("x", "y", "z")
    taxs    = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 6)
    taxtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 3)
    intersects = axis_intersections(c, cart_basis)
    for (i, V) in enumerate(cart_basis)
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
                line=attr(color=ifelse(j==1, AXIS_LIGHT_COL[], AXIS_COL[]),
                          width=ifelse(j==1, 5, 6))
                )
            if j == 2
                taxtips[i] = PlotlyJS.cone(
                    x=[V₁[1]], u=[V[1]],
                    y=[V₁[2]], v=[V[2]],
                    z=[V₁[3]], w=[V[3]],
                    sizeref=.1*scale, showscale=false, anchor="tail",
                    colorscale=[[0, AXIS_COL[]], [1, AXIS_COL[]]],
                    hovertext=name, hoverinfo="text")
            end
        end
    end

    # combine traces and plot
    ts = vcat(tbz, tgs, tgtips, taxs, taxtips)
    return PlotlyJS.plot(ts, layout; config=config)
end

# ---------------------------------------------------------------------------------------- #
# 2D

# default layout
const DEFAULT_PLOTLY_LAYOUT_2D  = Layout(
    showlegend=false,
    xaxis=attr(tickvals=[], zeroline=false,
            showgrid=false, showbackground=false,
            title=attr(text=""),
            ),
    yaxis=attr(tickvals=[], zeroline=false,
            showgrid=false, showbackground=false,
            title=attr(text=""),
            scaleanchor="x", scaleratio=1
            ),
    aspectmode = "data",
    hovermode="closest",
    margin=attr(l=0, r=0, b=0, t=0),
    autosize=false,
    plot_bgcolor="rgba(255, 255, 255, 1)", paper_bgcolor="rgba(255, 255, 255, 1)",
    annotations=PlotlyBase.PlotlyAttribute[]
    )

function plot(c::Cell{2}, layout::Layout=DEFAULT_PLOTLY_LAYOUT_2D;
              config::PlotConfig = PlotConfig(responsive=true, displaylogo=false))
    layout = deepcopy(layout) # because we have to mutate to get arrows...
    
    c.basisenum[] !== CARTESIAN && (c = cartesianize(c))
    
    scale = maximum(norm, basis(c))
    max_x, max_y = maximum(v->abs(v[1]), basis(c)), maximum(v->abs(v[2]), basis(c))
    get!(layout[:xaxis], :range, [-max_x-scale/15, max_x+scale/15])
    get!(layout[:yaxis], :range, [-max_y-scale/15, max_y+scale/15])

    # Cell boundaries
    tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(c))
    for (i,poly) in enumerate(c)
        tbz[i] = PlotlyJS.scatter(
            x=push!(getindex.(poly, 1), poly[1][1]),
            y=push!(getindex.(poly, 2), poly[1][2]);
            mode="lines", hovertext="Cell", hoverinfo="text+x+y",
            line=attr(color=BZ_COL[], width=3)
            )
    end

    # lattice vectors
    tgs    = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 4)
    tgtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 2)
    intersects = axis_intersections(c, basis(c))
    for (i,V) in enumerate(basis(c))
        V′ = V./norm(V)
        name = "<b>v</b><sub>$(i)</sub>"
        for j in (1,2) # in/outside BZ
            start = j == 1 ? 0.0           : intersects[i]
            stop  = j == 1 ? intersects[i] : 1.0
            V₀ = start*V .+ (j == 1 ? 0.0   : 0.025)*scale*V′
            V₁ = stop*V  .- (j == 1 ? 0.025 : 0.0)  *scale*V′
            tgs[i+(j-1)*2] = PlotlyJS.scatter(
                x=[V₀[1],V₁[1]], y=[V₀[2],V₁[2]];
                mode="lines", hovertext=name, hoverinfo=ifelse(j==1,"text","text+x+y"),
                line=attr(color=ifelse(j==1, BASIS_LIGHT_COL[], BASIS_COL[]),
                          width=ifelse(j==1, 5, 6))
                )
        end
        # arrow heads have to be added as annotations to layout in 2D :/
        haskey(layout, :annotations) || (layout[:annotations] = PlotlyBase.PlotlyAttribute[])
        push!(layout[:annotations],
            attr(x=V[1]+.05V′[1]*scale, y=V[2]+.05V′[2]*scale,   # awful fidgeting; plotly's
                 ax=V[1]-.05V′[1]*scale, ay=V[2]-.05V′[2]*scale, # arrows are stupid
                 xref="ax", yref="ay", axref="x", ayref="y",
                 showarrow=true, arrowhead=2, arrowwidth=6, arrowsize=.5,
                 arrowcolor=BASIS_COL[]))
    end

    # Cartesian axes
    cart_basis = cartesian_axes(Val(2))
    name_axs = ("x", "y") 
    taxs    = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 4)
    taxtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 2)
    intersects = axis_intersections(c, cart_basis)
    for (i, V) in enumerate(cart_basis)
        name = "<b>"*string(name_axs[i])*"</b>"
        V′ = V./norm(V)
        for j in (1,2) # in/outside BZ
            start = j == 1 ? 0.0           : intersects[i]
            stop  = j == 1 ? intersects[i] : intersects[i]
            V₀ = start*V .+ (j == 1 ? 0.0   : 0.025)*scale*V′
            V₁ = stop*V  .- (j == 1 ? 0.025 : -0.2) *scale*V′

            taxs[i+2(j-1)] = PlotlyJS.scatter(
                x=[V₀[1],V₁[1]], y=[V₀[2],V₁[2]];
                mode="lines", hovertext=name, hoverinfo="text",
                line=attr(color=ifelse(j==1, AXIS_LIGHT_COL[], AXIS_COL[]),
                            width=ifelse(j==1, 5, 6))
                )
            if j == 2
                # arrow heads have to be added as annotations to layout in 2D :/
                haskey(layout, :annotations) || (layout[:annotations] = PlotlyBase.PlotlyAttribute[])
                push!(layout[:annotations],
                    attr(x=V₁[1]+.05V′[1]*scale, y=V₁[2]+.05V′[2]*scale,   # awful fidgeting; plotly's
                         ax=V₁[1]-.05V′[1]*scale, ay=V₁[2]-.05V′[2]*scale, # arrows are stupid
                         xref="ax", yref="ay", axref="x", ayref="y",
                         showarrow=true, arrowhead=2, arrowwidth=6, arrowsize=.5,
                         arrowcolor=AXIS_COL[]))
            end
        end
    end

    # combine traces and plot
    ts = vcat(tbz, tgs, taxs)
    return PlotlyJS.plot(ts, layout; config=config)
end

# ---------------------------------------------------------------------------------------- #
# UTILITIES 

# the "outward" intersections of of lines with direction `Vs` though origo with a cell face
function axis_intersections(c::Cell{3},
            Vs::AbstractVector{<:AbstractVector}=cartesian_axes(Val(3)))

    intersects = MVector{3,Float64}(undef)
    fill!(intersects, Inf)
    # intersection between plane (r-rₒ) ⋅ n = 0 and line r(t) = l₀ + lt:
    #   t = (r₀ - l₀) ⋅ n / (l ⋅ n)     [https://wikipedia.org/wiki/Line–plane_intersection]
    for (i,poly) in enumerate(c)
        r₀ = sum(poly)/length(poly)
        n  = face_normal(c, i)
        for (j,V) in enumerate(Vs)
            t = dot(r₀, n)/dot(V, n)
            if t > 0.0 && t < intersects[j]
                intersects[j] = t
            end
        end
    end
    return intersects
end

# the "outward" intersections of of lines with direction `Vs` though origo with a cell segment
function axis_intersections(c::Cell{2},
    Vs::AbstractVector{<:AbstractVector}=cartesian_axes(Val(2)))

    intersects = MVector{2,Float64}(undef)
    fill!(intersects, Inf)
    # intersection between line rₐ(t) a + nₐt and line rᵦ(t) = β + nᵦt can be gotten by 
    # linear algebra: [-nₐ|nᵦ][tₐ,tᵦ] = a-β. Here we, pick rᵦ for the axes `Vs` (β=0) and
    # let rₐ refer to cell segments. We're then only interested in tᵦ:
    fs = faces(c)
    vs = vertices(c)
    for i in eachindex(vs) # number of vertices and segments are equal in 2D
        if length(fs) == 1 # `merge = true`
            idxs = i ≠ length(vs) ? (i,i+1) : (i,1)
        else
            idxs = (fs[i][1], fs[i][2])
        end
        a  = vs[idxs[1]]
        nₐ = vs[idxs[2]] - a
        for (j,V) in enumerate(Vs)
            nᵦ = V
            tᵦ = (nₐ[2]*a[1]-nₐ[1]*a[2])/(-nₐ[1]*nᵦ[2] + nₐ[2]*nᵦ[1]) # ~ inv([-nₐ|nᵦ])
            if tᵦ > 0.0 && tᵦ < intersects[j]
                intersects[j] = tᵦ
            end
        end
    end
    return intersects
end

# a set of cartesian basis vectors in D dimensions
cartesian_axes(Dᵛ::Val{D}) where D = SVector(ntuple(i->SVector(ntuple(j->i==j ? 1.0 : 0.0, Dᵛ)), Dᵛ))
