using VIDA
using Documenter

makedocs(;
    modules=[VIDA],
    authors="Paul <ptiede91@gmail.com> and contributors",
    repo="https://github.com/ptiede/VIDA.jl/blob/{commit}{path}#L{line}",
    sitename="VIDA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ptiede.github.io/VIDA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "getting_started.md",
        #"function_index.md",
    ],
)

deploydocs(;
    repo="github.com/ptiede/VIDA.jl",
)
