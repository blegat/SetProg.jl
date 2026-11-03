using Aqua
using SetProg
using Test

@testset "Aqua" begin
    # SetProg intentionally adds JuMP membership constraints for halfspaces.
    Aqua.test_all(SetProg; piracies = (treat_as_own = [Polyhedra.HalfSpace],))
end
