using Brillouin
using Documenter

makedocs(;
    modules=[Brillouin],
    authors="Thomas Christensen <tchr@mit.edu> and contributors",
    repo="https://github.com/thchr/Brillouin.jl/blob/{commit}{path}#L{line}",
    sitename="Brillouin.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://thchr.github.io/Brillouin.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/thchr/Brillouin.jl",
)
