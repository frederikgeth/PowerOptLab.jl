# Solver-based local control reference. Each stage retains the physical plant;
# no weighted sum is used to trade reactive support against active export.

"""Common result interface for solver-based local inverter service allocation."""
abstract type AbstractInverterServiceAllocationResult <: AbstractSolveResult end

"""
    InverterServiceAllocationResult

Result of [`allocate_inverter_services`](@ref). `stages` preserves every local
solve in Q-tracking, P-export, then I₂-support order. `allocation` is the last
stage, including a failed candidate. `current_command` is `nothing` unless all
stages publish. `plant_feasibility` is `:feasible` only with an accepted stage or
an independent diagnostic witness; otherwise it is `:unresolved`.
`policy_feasibility` is `:feasible` only when the final stage publishes. An NLP
failure is inconclusive, even if the unrestricted plant diagnostic succeeds.
These are local nonconvex results, not global lexicographic optimality proofs.
"""
struct InverterServiceAllocationResult <: AbstractInverterServiceAllocationResult
    reactive_request_var::Float64
    negative_sequence_request_A::ComplexF64
    target::Symbol
    reactive_mode::Symbol
    stages::Vector{Pair{Symbol,ConstrainedInverterResult}}
    allocation::ConstrainedInverterResult
    negative_scale::Float64
    current_command::Union{Nothing,NTuple{3,ComplexF64}}
    plant_feasibility::Symbol
    policy_feasibility::Symbol
    plant_diagnostic::Union{Nothing,InverterCurrentProjectionResult}
end

solve_status(r::InverterServiceAllocationResult)=solve_status(r.allocation)
physical_feasibility(r::InverterServiceAllocationResult)=physical_feasibility(r.allocation)
solve_diagnostics(r::InverterServiceAllocationResult)=(
    reactive_request_var=r.reactive_request_var,target=r.target,reactive_mode=r.reactive_mode,
    negative_sequence_request_A=r.negative_sequence_request_A,negative_scale=r.negative_scale,
    plant_feasibility=r.plant_feasibility,policy_feasibility=r.policy_feasibility,
    stages=[(name=name,solve=solve_status(stage),operation=solve_diagnostics(stage)) for (name,stage) in r.stages])

"""
    allocate_inverter_services(inv, request; grid_voltage, reactive_request,
        active_fraction=1., volt_watt_basis=:rated, negative_sequence_request=0im,
        target=:converter, reactive_mode=:within_request,
        dc_source=DCSourceBounds(p_min=0.), diagnose_plant=true,
        s_base=1e5, optimizer=Ipopt.Optimizer, verbose=false, solver_options=())

Allocate a three-leg command at measured, fixed POC RMS voltage phasors using
the full averaged plant and continuous rail refinement. Stage 1 minimizes total
POC Q tracking error; stage 2 holds its achieved Q and maximizes POC P; stage 3,
when requested, holds both powers and maximizes the scale of the requested I₂.
There is no finite objective weight that exchanges higher-priority service for
lower-priority service. Every stage independently checks physical limits, fixed
voltage, sequence-command and preserved-power residuals before continuing.

`:within_request` restricts Q between zero and `reactive_request`. Explicit
`:protective_override` permits Q outside that interval, seeking the closest
locally feasible Q when voltage protection requires it. This changes the policy;
it is never activated as a silent fallback. I₂ follows the requested complex ray
with scale in [0,1] on the selected converter/grid arm. Positive-sequence current
delivers the chosen total POC powers, including filter/shunt effects.

The volt-watt ceiling uses `p_rated*active_fraction` for `:rated` and
`active_fraction*(available_dc - (p_dc-p_poc))` for `:available`; DC availability
also remains a hard bound. Default source bounds prevent DC charging and allow
AC loss imports at zero PV. `grid_voltage` refers to the plant neutral or earth
when absent, in phase-terminal order. Supports THREE_LEG, non-grid-forming
plants, including centered-PWM thermal allocation with an open carrier-frequency
source. PWM/nonlinear-energy stages preserve service within a 0.0001 var/W SI
budget to avoid a singular tangent lock, without relaxing hardware constraints.
Four-leg phase service uses `allocate_phase_inverter_services`.

This solver-based local reference can be re-evaluated at successive measurements.
Its solve time, feeder equilibrium, stability and response time are not certified.
It does not replace the fixed-structure `SequenceController` stamp. Local solver
success does not establish globally best Q/P/I₂ allocation. On a failed first
stage, an optional unrestricted physical projection supplies a plant-feasibility
witness without issuing that alternative as the controller command.
"""
function allocate_inverter_services(inv::AdvancedInverter,request::InverterControlRequest;
        grid_voltage,reactive_request::Real,active_fraction::Real=1.,
        volt_watt_basis::Symbol=:rated,negative_sequence_request::Number=0im,
        target::Symbol=:converter,reactive_mode::Symbol=:within_request,
        dc_source::DCSourceBounds=DCSourceBounds(p_min=0.),diagnose_plant::Bool=true,
        s_base::Real=1e5,optimizer=Ipopt.Optimizer,verbose::Bool=false,solver_options=())
    inv.topology===:THREE_LEG || throw(ArgumentError("service allocation currently requires THREE_LEG"))
    target in (:converter,:grid) || throw(ArgumentError("target must be :converter or :grid"))
    reactive_mode in (:within_request,:protective_override) || throw(ArgumentError("unsupported reactive_mode"))
    volt_watt_basis in (:rated,:available) || throw(ArgumentError("volt_watt_basis must be :rated or :available"))
    isfinite(s_base) && s_base>0 || throw(ArgumentError("s_base must be positive and finite"))
    isfinite(reactive_request) || throw(ArgumentError("reactive_request must be finite"))
    isfinite(active_fraction) && 0<=active_fraction<=1 || throw(ArgumentError("active_fraction must lie in [0,1]"))
    i2=ComplexF64(negative_sequence_request)
    isfinite(i2) || throw(ArgumentError("negative_sequence_request must be finite"))
    q=Float64(reactive_request); fraction=Float64(active_fraction)
    v=_projection_phasors(grid_voltage,"grid_voltage")
    policy = reactive_mode===:within_request && iszero(q) ?
        InverterOperatingConstraints(reactive=FixedReactivePower(0.)) : InverterOperatingConstraints()
    device=ConstrainedInverter(inv,policy;current_allocation=:free,dc_source)
    net=_projection_network(inv,v)
    validate_device(device,(net,);periods=1)
    handles=Ref{_ConstrainedInverterHandles}()
    negative=Ref{Any}(nothing)
    reactive_scale=Ref{Any}(nothing)
    error_variable=Ref{Any}(nothing)
    error_rows=Any[]
    hook! = ctx -> begin
        handles[]=stamp_device!(ctx,device;request)
        h=handles[].plant; m=_opf_model(ctx); scale=inv.s_max/h.sb
        cap=volt_watt_basis===:rated ? request.p_rated/h.sb*fraction :
            fraction*(_available_dc(device,request)/h.sb-(h.p_dc-h.p_poc))
        JuMP.@constraint(m,(h.p_poc-cap)/scale<=0)
        re,im=target===:converter ? (h.cri,h.cii) : (h.gri,h.gii)
        seq=_sequence_components(re,im)[3]
        current_scale=InverterControlRatings(inv,_current_target(target)).i_scale/h.ib
        alpha=iszero(i2) ? 1. : JuMP.@variable(m,lower_bound=0.,upper_bound=1.,start=.5)
        negative[]=alpha
        JuMP.@constraint(m,(seq[1]-alpha*real(i2)/h.ib)/current_scale==0)
        JuMP.@constraint(m,(seq[2]-alpha*imag(i2)/h.ib)/current_scale==0)
        if reactive_mode===:within_request
            if iszero(q)
                JuMP.@objective(m,Min,0.)
            else
                qs=JuMP.@variable(m,lower_bound=0.,upper_bound=1.,start=.5)
                reactive_scale[]=qs
                JuMP.@constraint(m,(h.q_poc-qs*q/h.sb)/scale==0)
                JuMP.@objective(m,Max,qs)
            end
        else
            error=JuMP.@variable(m,lower_bound=0.,start=.1)
            error_variable[]=error
            push!(error_rows,JuMP.@constraint(m,(h.q_poc-q/h.sb)/scale<=error))
            push!(error_rows,JuMP.@constraint(m,-(h.q_poc-q/h.sb)/scale<=error))
            JuMP.@objective(m,Min,error)
        end
    end
    ctx=build_opf_model(net;per_unit=true,s_base=Float64(s_base),add_objective=false,
        model_hook! = hook!,optimizer,verbose)
    m=_opf_model(ctx);h=handles[].plant
    effective_options=Dict{String,Any}()
    if JuMP.solver_name(m)=="Ipopt"
        merge!(effective_options,Dict("tol"=>1e-9,"max_iter"=>1000,
            "mu_strategy"=>"adaptive","nlp_scaling_method"=>"none"))
    end
    for (name,value) in (solver_options isa NamedTuple ? pairs(solver_options) : solver_options)
        effective_options[string(name)]=value
    end
    _set_inverter_solver_options!(m,effective_options)
    stages=Pair{Symbol,ConstrainedInverterResult}[]
    qlock=plock=nothing
    alpha_value=NaN
    function capture(name)
        outcome=_solve_outcome(m);status=SolveStatus(outcome)
        e=extract_device(device,handles[],status)
        report=e.plant.physical_feasibility
        if status.has_primal
            actual_v=complex.(JuMP.value.(h.vrpoc),JuMP.value.(h.vipoc)).*h.vb
            report=_projection_boundary_report(report,actual_v,v,0.,0.)
            ire,iim=target===:converter ? (h.cri,h.cii) : (h.gri,h.gii)
            actual_i=complex.(JuMP.value.(ire),JuMP.value.(iim)).*h.ib
            actual2=_phasor_sequences(Tuple(actual_i))[3]
            actual_grid=complex.(JuMP.value.(h.gri),JuMP.value.(h.gii)).*h.ib
            actual_power=sum(actual_v.*conj.(actual_grid))+im*inv.b_filter_shunt*sum(abs2,actual_v)
            alpha_value=iszero(i2) ? 1. : JuMP.value(negative[])
            checks=InverterLimitCheck[
                _inverter_limit_check("service_negative_sequence",:A,abs(actual2-alpha_value*i2),0.,0.)]
            if !iszero(i2)
                push!(checks,_inverter_limit_check("service_negative_ray",:A,alpha_value*abs(i2),0.,abs(i2)))
            end
            for (label,value,locked) in (("Q",imag(actual_power),qlock),("P",real(actual_power),plock))
                locked===nothing || push!(checks,_inverter_limit_check("service_preserved_$label",:VA,
                    abs(value-locked),0.,0.))
            end
            reactive_mode===:within_request && push!(checks,_inverter_limit_check(
                "service_reactive_ray",:VA,imag(actual_power),min(0.,q),max(0.,q)))
            # For the available basis, P <= f*(A-Pdc+P) is equivalent to
            # (1-f)*P+f*Pdc <= f*A. Coefficients are nonnegative, so the SI
            # auditor's reported-Pdc residual supplies a conservative enclosure.
            dc_residual=only(c.value for c in report.checks if c.name=="reported_p_dc_residual")
            pdc_upper=JuMP.value(h.p_dc)*h.sb+dc_residual
            protected_power=volt_watt_basis===:rated ? real(actual_power) :
                (1-fraction)*real(actual_power)+fraction*pdc_upper
            ceiling=fraction*(volt_watt_basis===:rated ? request.p_rated : _available_dc(device,request))
            push!(checks,_inverter_limit_check("service_volt_watt",:VA,protected_power,-Inf,ceiling))
            report=_extend_inverter_feasibility(report,checks)
        end
        plant=merge(e.plant,(physical_feasibility=report,))
        operation=merge(e.operation,(allocation_stage=name,raw_solver_status=JuMP.raw_status(m),
            physical_feasible=report.passed,physical_failed_checks=[c.name for c in report.checks if !c.passed]))
        network=_extract_result(ctx,outcome)
        r=ConstrainedInverterResult(string(outcome.termination_status),plant,operation,
            e.converter_terminal,e.grid_phase_current,network["bus"],_physical_solve_status(status,report))
        push!(stages,name=>r)
        r
    end
    _optimize_inverter_rails!(ctx,((inv,h),))
    result=capture(:reactive_tracking)
    if solve_status(result).publishable
        actual_power=sum(result.plant.poc_phase_voltage.*conj.(result.grid_phase_current))+
            im*inv.b_filter_shunt*sum(abs2,result.plant.poc_phase_voltage)
        qlock=imag(actual_power)
        starts=[(x,JuMP.value(x)) for x in JuMP.all_variables(m)]
        qscale=reactive_scale[]===nothing ? nothing : JuMP.value(reactive_scale[])
        JuMP.@objective(m,Max,h.p_poc/(inv.s_max/h.sb))
        if reactive_mode===:within_request
            # Fix the achieved service coordinate instead of adding a Q equality
            # almost coincident with an active Q interval bound. The latter
            # produces a spurious infeasibility at ordinary Australian requests.
            if reactive_scale[]!==nothing
                inv.pwm_strategy===:NONE && inv.dc_energy_model===:linearized ? JuMP.fix(reactive_scale[],qscale;force=true) :
                    _lock_pwm_service!(m,reactive_scale[],qscale,abs(q))
            end
        else
            JuMP.@constraint(m,h.q_poc/(inv.s_max/h.sb)==qlock/inv.s_max)
            for row in error_rows;JuMP.delete(m,row);end
            JuMP.delete(m,error_variable[])
        end
        for (x,value) in starts
            JuMP.is_valid(m,x) && JuMP.set_start_value(x,value)
        end
        _optimize_inverter_rails!(ctx,((inv,h),))
        result=capture(:active_export)
    end
    if solve_status(result).publishable && !iszero(i2)
        plock=real(sum(result.plant.poc_phase_voltage.*conj.(result.grid_phase_current)))
        starts=[(x,JuMP.value(x)) for x in JuMP.all_variables(m)]
        if inv.pwm_strategy===:NONE && inv.dc_energy_model===:linearized
            JuMP.@constraint(m,h.p_poc/(inv.s_max/h.sb)==plock/inv.s_max)
        else
            JuMP.@constraint(m,(plock-_PWM_SERVICE_PRESERVATION_VA)/inv.s_max<=
                h.p_poc/(inv.s_max/h.sb)<=(plock+_PWM_SERVICE_PRESERVATION_VA)/inv.s_max)
        end
        JuMP.@objective(m,Max,negative[])
        for (x,value) in starts; JuMP.set_start_value(x,value); end
        _optimize_inverter_rails!(ctx,((inv,h),))
        result=capture(:negative_sequence_support)
    end
    accepted=solve_status(result).publishable
    witnessed=any(solve_status(r).publishable for (_,r) in stages)
    diagnostic=nothing
    if !witnessed && diagnose_plant
        u1=_phasor_sequences(Tuple(v))[2]
        desired1=abs(u1)>1e-9 ? conj(complex(min(request.p_available,request.p_rated*fraction),q)/(3u1)) : 0.0im
        diagnostic=project_inverter_current(inv,request;grid_voltage=v,
            requested_current=_phase_from_sequences(desired1,i2),target,dc_source,s_base,
            optimizer,verbose,solver_options=effective_options)
        witnessed=diagnostic.plant_feasibility===:feasible
    end
    command=accepted ? (target===:converter ? result.converter_terminal.phase_current : result.grid_phase_current) : nothing
    InverterServiceAllocationResult(q,i2,target,reactive_mode,stages,result,alpha_value,command,
        witnessed ? :feasible : :unresolved,accepted ? :feasible : :unresolved,diagnostic)
end

"""
    evaluate_voltage_aware(controlled, measurement, request; grid_voltage, kwargs...)

Evaluate a solver-based Q→P→I₂ allocation from the exact mean/worst/sequence
VV+VW voltage laws. `measurement` contains sensed phasors; `grid_voltage` must
separately declare the physical POC phasors because sensing and power references
can differ. Requires reactive priority and power-first sequence priority. Supports
NoUnbalanceControl or voltage-admittance droop with `ripple_blend=0`; a
current-dependent ripple-compensation request needs a separate feedback model.
The result is a local fixed-measurement command, not a network equilibrium or
dynamic-stability certificate. No existing controller stamp is changed.
"""
function evaluate_voltage_aware(controlled::ControlledDevice{<:AdvancedInverter,<:SequenceController},
        measurement::InverterControlMeasurement,request::InverterControlRequest;grid_voltage,kwargs...)
    c=controlled.controller
    c.limiter.pq_priority===:var || throw(ArgumentError("voltage-aware service allocation requires reactive priority"))
    c.limiter isa SequencePriorityLimiter && c.limiter.priority!==:power && throw(ArgumentError("Q→P→I2 requires sequence priority=:power"))
    c.unbalance isa NegativeSequenceAdmittanceDroop && !iszero(c.unbalance.ripple_blend) &&
        throw(ArgumentError("current-dependent ripple compensation is not supported by this fixed-measurement allocator"))
    _,u1,u2=_phasor_sequences(measurement.phase_voltage)
    abs(u1)>=c.power_voltage_floor || throw(ArgumentError("measured |U1| is below the controller power_voltage_floor"))
    fraction,q,_,_=_positive_exact(c.positive,u1,abs.(measurement.phase_voltage),request)
    _,_,i2,_=_unbalance_exact(c.unbalance,u1,u2,0.0im)
    target=c.current_target isa ConverterCurrentTarget ? :converter : :grid
    allocate_inverter_services(controlled.device,request;grid_voltage,reactive_request=q,
        active_fraction=fraction,volt_watt_basis=_volt_watt_basis(c.positive),
        negative_sequence_request=i2,target,kwargs...)
end
