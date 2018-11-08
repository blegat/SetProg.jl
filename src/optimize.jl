function load(model::JuMP.Model, d::Data)
    for variable in d.variables
        load(d, variable)
    end
    for (index, constraint) in d.constraints
        load(d, constraint)
    end
    if d.objective !== nothing
        load(model, d.objective)
    end
end

# set as `optimize_hook` to JuMP Model in `data` so it is called at
# `JuMP.optimize!` if at least one set variable is created
function optimize_hook(model::JuMP.AbstractModel)
    d = data(model)
    set_space(d, model)
    load(model, d)
    JuMP.optimize!(model, ignore_optimize_hook = true)
end
