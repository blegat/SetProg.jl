using Aqua
using SetProg
using Test

@testset "Aqua" begin
    # Also catches precompilation failures such as overwritten methods.
    Aqua.test_persistent_tasks(SetProg)
end
