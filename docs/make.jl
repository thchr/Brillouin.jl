using Documenter
using Brillouin

# ---------------------------------------------------------------------------------------- #
# make PlotlyJS plots showable in ```@example ``` blocks, following the approach suggested
# in https://github.com/fredrikekre/Literate.jl/issues/126
using PlotlyJS
struct HTMLPlot
    p
    h::Int # desired display height in pixels
end
HTMLPlot(p) = HTMLPlot(p, 400)
const ROOT_DIR = joinpath(@__DIR__, "build")
const PLOT_DIR = joinpath(ROOT_DIR, "plots")
function Base.show(io::IO, ::MIME"text/html", p::HTMLPlot)
    mkpath(PLOT_DIR)
    path = joinpath(PLOT_DIR, string(hash(p) % UInt32, ".html"))
    PlotlyJS.savehtml(p.p, path, :remote)
    print(io, "<object type=\"text/html\" data=\"../$(relpath(path, ROOT_DIR))\" style=\"width:100%;height:$(p.h)px;\"></object>")
end
# ---------------------------------------------------------------------------------------- #

makedocs(;
    modules=[Brillouin],
    authors="Thomas Christensen <tchr@mit.edu> and contributors",
    repo="https://github.com/thchr/Brillouin.jl/blob/{commit}{path}#L{line}",
    sitename="Brillouin.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://thchr.github.io/Brillouin.jl",
        assets = ["assets/custom.css"], # increase logo size
    ),
    pages=[
        "Home" => "index.md",
        "Wigner–Seitz cells" => "wignerseitz.md",
        "k-space paths" => "kpaths.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/thchr/Brillouin.jl",
)

# make sure we actually do `using Brillouin` before calling doctests (as described in
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Module-level-metadata)
DocMeta.setdocmeta!(Brillouin, :DocTestSetup, :(using Brillouin); recursive=true)