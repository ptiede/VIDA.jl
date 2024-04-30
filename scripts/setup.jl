using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.update()

@info "Testing whether you can import VIDA"
using VIDA
@info "If there are no errors then VIDA imported correctly"
