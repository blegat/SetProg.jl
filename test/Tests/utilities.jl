using Test, JuMP
const MOIT = MOI.Test
const MOIB = MOI.Bridges

function _model(optimizer::MOI.AbstractOptimizer)
    MOI.empty!(optimizer)
    return direct_model(optimizer)
end

function _model(factory::OptimizerFactory)
    return Model(factory)
end

#"""
#    @test_suite setname subsets
#
#Defines a function `setname_test(model, config, exclude)` that runs the tests
#defined in the dictionary `setname_tests` with the model `model` and config
#`config` except the tests whose dictionary key is in `exclude`. If `subsets` is
#`true` then each test runs in fact multiple tests hence the `exclude` argument
#is passed as it can also contains test to be excluded from these subsets of
#tests.
#"""
macro test_suite(setname, subsets=false)
    testname = Symbol(string(setname) * "_test")
    testdict = Symbol(string(testname) * "s")
    if subsets
        runtest = :( f(model, config, exclude) )
    else
        runtest = :( f(model, config) )
    end
    esc(:(
      function $testname(model::Union{$MOI.ModelLike, OptimizerFactory},
                         config::$MOI.Test.TestConfig,
                         exclude::Vector{String} = String[])
            for (name,f) in $testdict
                if name in exclude
                    continue
                end
                @testset "$name" begin
                    $runtest
                end
            end
        end
    ))
end

# Utilities for building the mock `optimize!` from the solution of a solver
# Variables primal values for inner bridged model
function print_value(v, atol)
    i = round(v)
    if isapprox(v, i, atol=atol)
        print(float(i))
    else
        print(v)
    end
end
function inner_variable_value(model, atol=1e-4)
    inner = backend(model)
    println("optimize!(mock) = MOIU.mock_optimize!(mock,")
    println(JuMP.termination_status(model))
    if JuMP.primal_status(model) != MOI.NO_SOLUTION
        values = MOI.get(inner, MOI.VariablePrimal(),
                         MOI.get(inner, MOI.ListOfVariableIndices()))
        print("(MOI.FEASIBLE_POINT, [")
        for (i, v) in enumerate(values)
            if i > 1
                print(", ")
            end
            print_value(v, atol)
        end
        print("])")
    else
        print("tuple()")
    end
    println(",")
    if JuMP.dual_status(model) != MOI.NO_SOLUTION
        for (F, S) in MOI.get(inner, MOI.ListOfConstraints())
            print("($F, $S) => [")
            for ci in MOI.get(inner, MOI.ListOfConstraintIndices{F, S}())
                print(MOI.get(inner, MOI.ConstraintDual(), ci))
                print(", ")
            end
            println("])")
        end
    end
    println(")")
end
# Constraint dual values for inner bridged model
