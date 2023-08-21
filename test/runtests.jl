using Plots
using VIDA
using Test

@testset "VIDA.jl" begin
    # Write your own tests here.
    tests = [
        "templates",
        "images",
        "movies",
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
