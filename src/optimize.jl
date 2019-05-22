function load(model::JuMP.Model, d::Data)
    d.state = Loading
    for variable in d.variables
        load(model, variable)
    end
    for (index, constraint) in d.constraints
        cref = load(model, constraint)
        if cref !== nothing
            d.transformed_constraints[index] = cref
        end
    end
    if d.objective !== nothing
        load(model, d.objective)
    end
    d.state = Modeling
end

# set as `optimize_hook` to JuMP Model in `data` so it is called at
# `JuMP.optimize!` if at least one set variable is created
function optimize_hook(model::JuMP.AbstractModel)
    d = data(model)
    set_space(d, model)
    # In case `optimize!` is called then the problem is modified and then it is
    # called again we need to clear first the space that might be wrong
    synchronize_perspective(d)
    clear_spaces(d)
    create_spaces(d)
    load(model, d)
    JuMP.optimize!(model, ignore_optimize_hook = true)
end
