using LaVIDA
using Documenter

makedocs(;
    modules=[LaVIDA],
    authors="Paul <ptiede91@gmail.com> and contributors",
    repo="https://github.com/ptiede/LaVIDA.jl/blob/{commit}{path}#L{line}",
    sitename="LaVIDA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ptiede.github.io/LaVIDA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "getting_started.md",
        "function_index.md",
    ],
)

deploydocs(;
    repo="github.com/ptiede/LaVIDA.jl",
)
