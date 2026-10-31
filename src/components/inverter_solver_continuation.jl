# Numerical continuation infrastructure. Dispatch objectives, temporary-row
# policy and physical acceptance belong to the caller. A residual check never
# replaces the caller's raw solver gate.
function _validate_inverter_initialization(initialization)
    initialization in (:automatic,:spread,:consistent) ||
        throw(ArgumentError("initialization must be :automatic, :spread or :consistent"))
end

function _inverter_positive_iteration_budget(options)
    for (name,value) in (options isa NamedTuple ? pairs(options) : options)
        string(name)=="max_iter" && value isa Real && value<=0 && return false
    end
    true
end

function _inverter_numerical_restart(status,options)
    _inverter_positive_iteration_budget(options) || return false
    status.termination_status in ("OTHER_ERROR","NUMERICAL_ERROR",
        "ALMOST_LOCALLY_SOLVED","SLOW_PROGRESS","ITERATION_LIMIT")
end

function _inverter_numerical_attempt(m,status,initialization,objective,physical;stage=:final)
    (initialization=initialization===:consistent ? :consistent : :spread,
        stage=stage,selection_objective=objective,raw=status,
        raw_solver_status=JuMP.raw_status(m),physical_accepted=physical,
        nlp_scaling_method=JuMP.solver_name(m)=="Ipopt" ? JuMP.get_attribute(m,"nlp_scaling_method") : nothing)
end

function _inverter_zero_objective_seed!(ctx,plants,loss,accepted;rail_refinement=true)
    m=_opf_model(ctx)
    JuMP.@objective(m,Min,loss)
    _optimize_inverter_rails!(ctx,plants;enabled=rail_refinement)
    status=SolveStatus(_solve_outcome(m))
    physical=status.has_primal && accepted(status)
    attempt=_inverter_numerical_attempt(m,status,:consistent,:loss,physical;stage=:loss_initialization)
    # Time limits end this solve path. Other failed seeds supply no starts.
    status.termination_status=="TIME_LIMIT" && return (attempts=(attempt,),continue_final=false)
    start=status.publishable && physical ? _inverter_start_snapshot(m) : nothing
    JuMP.@objective(m,Min,0.)
    _restore_inverter_start!(m,start)
    (attempts=(attempt,),continue_final=true)
end

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
