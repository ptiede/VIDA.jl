using Documenter

# if haskey(ENV, "GITHUB_ACTIONS")
#     ENV["JULIA_DEBUG"] = "Documenter"
# end


using VIDA
using Plots
using Literate

EXAMPLE1 = joinpath(@__DIR__, "..", "example", "introduction.jl")
EXAMPLE2 = joinpath(@__DIR__, "..", "example", "custom_template.jl")
OUTPUT = joinpath(@__DIR__, "src/generated")

Literate.markdown(EXAMPLE1, OUTPUT, documenter=true)
Literate.markdown(EXAMPLE2, OUTPUT, documenter=true)

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
        "interface.md",
        "generated/introduction.md",
        "generated/custom_template.md",
        "api/function_index.md",
    ],
)

deploydocs(;
    repo="github.com/ptiede/VIDA.jl",
    push_preview=true
)
