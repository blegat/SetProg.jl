using Test
using SetProg
const MOIT = MOI.Test

function bridged_mock(mock_optimize!::Function...;
                      # We use a UF to add support for SOSCone, ... so that
                      # we don't have to set the variables created by the SOS
                      # bridges
                      model = MOIU.UniversalFallback(MOI.Utilities.Model{Float64}()))
    mock = MOI.Utilities.MockOptimizer(model)
    bridged = MOI.Bridges.full_bridge_optimizer(mock, Float64)
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!...)
    return bridged
end
