# We want to make `plot(kp)`, `plot(c)`, `plot(c, kp)` work out of the box, and we want to
# apply the axis-settings from `_default_bare_axis!` to the `KPath` and `Cell` plots;
# as far as I have been able to tell, there is currently no good way to do this in Makie -
# it is possible to do with `SpecApi` and `Makie.convert_arguments` - but apparently not if
# we simultaneously use a `@recipe`.
# If we hadn't used a `@recipe`, we should have been able to do something like this:
#       ```
#       function Makie.convert_arguments(::Type{Plot{plot}}, c::Cell{D}) where D
#           ax = S.Axis(
#               plots = [S.CellPlot(c)];
#               xticklabelsvisible = false,
#               xticksvisible=false, xgridvisible=false, yticklabelsvisible=false, 
#               yticksvisible=false, ygridvisible=false, bottomspinevisible=false, 
#               topspinevisible=false, leftspinevisible=false, rightspinevisible=false)
#           return S.GridLayout([ax])
#       end
#       ```
#       (the big chunk of `x/ytick, x/ygrid, bottom/top/left/rightspine` is manually doing
#       `hidedecorations!(ax)` and `hidespines!(ax)` but in declarative style)
#
# For now, we do the following, which is not great, since it overload `plot` which is
# against the design philosophy of Makie, but it works.

# -----------------------------
# plot(:KPath)

function Makie.plot(
    kp::Union{Observable{KPath{D}}, KPath{D}};
    axis = NamedTuple(),
    figure = NamedTuple(),
    kws...
) where D
    f = Makie.Figure(; figure...)
    ax = _default_bare_axis!(f, Val(D); axis)
    p = Makie.plot!(ax, kp; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end

# -----------------------------
# plot(::Cell)

function Makie.plot(
    c::Union{Observable{Cell{D}}, Cell{D}};
    axis = NamedTuple(),
    figure = NamedTuple(),
    kws...
) where D
    f = Makie.Figure(; figure...)
    ax = _default_bare_axis!(f, Val(D); axis)
    p = Makie.plot!(ax, c; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end

# -----------------------------
# plot(::Cell, :KPath)

function Makie.plot(
    c::Union{Observable{Cell{D}}, Cell{D}},
    kp::Union{Observable{KPath{D}}, KPath{D}};
    axis = NamedTuple(),
    figure = NamedTuple(),
    kws...
) where D
    f = Makie.Figure(; figure...)
    ax = _default_bare_axis!(f, Val(D); axis)
    p = Makie.plot!(ax, c, kp; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end

## --------------------------------------------------------------------------------------- #
# default axis settings in 3D and 1D+2D

function _default_bare_axis!(f, ::Val{3}; axis=NamedTuple())
    f[1,1] = ax = Makie.Axis3(f;
        aspect=:data,
        viewmode=:fit,
        axis...)
    Makie.hidedecorations!(ax)
    ax.protrusions[] = 0 # cf. https://github.com/MakieOrg/Makie.jl/issues/2259
    Makie.hidespines!(ax)

    return ax
end

function _default_bare_axis!(f, ::Union{Val{1}, Val{2}}; axis=NamedTuple())
    f[1,1] = ax = Makie.Axis(f;
        aspect=Makie.DataAspect(),
        axis...)
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    return ax
end