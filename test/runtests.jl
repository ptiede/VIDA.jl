using VIDA
using Test
using Plots

@testset "VIDA.jl" begin
    # Write your own tests here.
    tests = [
        "filters",
        "images",
        "divergences",
        "visualizations",
        "extractor"
    ]

    res = map(tests) do t
        @eval module $(Symbol("Test_", t))
            include($t*".jl")
        end
        return
    end
end
