# Numerical continuation infrastructure. Dispatch objectives, temporary-row
# policy and physical acceptance belong to the caller. A residual check never
# replaces the caller's raw solver gate.
function _inverter_start_snapshot(m)
    variables=JuMP.all_variables(m)
    constraints=JuMP.has_duals(m) ?
        JuMP.all_constraints(m;include_variable_in_set_constraints=true) : JuMP.ConstraintRef[]
    (variables=variables,primals=JuMP.value.(variables),
        constraints=constraints,duals=JuMP.dual.(constraints))
end

function _restore_inverter_start!(m,start)
    start===nothing && return
    for (variable,value) in zip(start.variables,start.primals)
        JuMP.is_valid(m,variable) && isfinite(value) && JuMP.set_start_value(variable,value)
    end
    for (constraint,value) in zip(start.constraints,start.duals)
        JuMP.is_valid(m,constraint) && isfinite(value) && JuMP.set_dual_start_value(constraint,value)
    end
    if JuMP.solver_name(m)=="Ipopt"
        JuMP.set_attribute(m,"warm_start_init_point","yes")
        JuMP.set_attribute(m,"warm_start_bound_push",1e-9)
        JuMP.set_attribute(m,"warm_start_bound_frac",1e-9)
        JuMP.set_attribute(m,"warm_start_mult_bound_push",1e-9)
    end
end

function _inverter_solver_continuation!(m,groups,weights;solve_stage!,accepted)
    !isempty(weights) && iszero(last(weights)) ||
        throw(ArgumentError("final continuation weight must be zero"))
    all(w->isfinite(w)&&w>=0,weights) || throw(ArgumentError("invalid continuation weights"))
    attempts=Any[];start=nothing
    function run_stage!(name,weight)
        result=solve_stage!(name,weight)
        # Snapshot all primals/duals before any model/start mutation. Failed
        # stages must not seed continuation, even when they have finite values.
        accepted(result) && (start=_inverter_start_snapshot(m))
        push!(attempts,result)
    end
    cuts=[cut for (_,group) in groups for cut in group]
    try
        if !isempty(groups)
            run_stage!(:neutral_initialization,first(weights))
            for (index,(name,group)) in enumerate(groups)
                _restore_inverter_start!(m,start)
                for cut in group
                    JuMP.delete(m,cut)
                end
                index<length(groups) && run_stage!(Symbol("neutral_release_"*string(name)),first(weights))
            end
        end
        all(c->!JuMP.is_valid(m,c),cuts) || error("temporary initialization constraint survived")
        for weight in weights
            _restore_inverter_start!(m,start)
            run_stage!(:full_capability,weight)
        end
    finally
        # An exception must not leave a caller's model silently restricted.
        for cut in cuts
            JuMP.is_valid(m,cut) && JuMP.delete(m,cut)
        end
    end
    (attempts=attempts,temporary_constraints_remaining=count(c->JuMP.is_valid(m,c),cuts))
end
