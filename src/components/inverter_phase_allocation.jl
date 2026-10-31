# Four-leg phase service with one shared DC source and explicit neutral physics.

"""
    PerPhaseVoltVarWatt(; volt_var=nothing, volt_watt=nothing,
        volt_watt_basis=:rated, current_target=:converter,
        power_voltage_floor=1., sensing_reference=:device)

Four-leg voltage-aware local policy. Each phase uses its own sensed RMS voltage;
Q and rated P bases are one third of the total request bases. Reactive priority
maximizes one common [0,1] service coordinate for the three Q requests, then holds
it while maximizing total P (100 µvar preservation budget in PWM/nonlinear mode).
POC active powers are equal across phases; that is
not balanced current under unbalanced voltage. Available-basis VW uses one third
of shared AC availability, including total solved losses. No separate phase PV
sources or per-phase VA ratings are introduced. Zero/negative sequence and the
neutral current follow from the phase powers and full physical filter.

Use [`evaluate_voltage_aware`](@ref); fixed-graph implicit controller stamping is
not provided. This is an exact-curve, solver-based quasi-static reference, not a
firmware timing, global-optimality or dynamic-stability guarantee. Requires an
explicit neutral, FOUR_LEG and non-GFM plant. Centered-PWM capacitor allocation
is supported with an open carrier-frequency DC source. Omitted current limits
remain omitted; in particular no neutral rating is inferred from a phase rating.
"""
struct PerPhaseVoltVarWatt <: AbstractInverterControlLaw
    volt_var::Union{Nothing,PiecewiseLinearLaw}
    volt_watt::Union{Nothing,PiecewiseLinearLaw}
    volt_watt_basis::Symbol
    current_target::AbstractCurrentTarget
    power_voltage_floor::Float64
    sensing_reference::Union{Symbol,String}
    function PerPhaseVoltVarWatt(;volt_var=nothing,volt_watt=nothing,
            volt_watt_basis::Symbol=:rated,current_target=:converter,
            power_voltage_floor::Real=1.,sensing_reference=:device)
        _validate_positive_curves(volt_watt,volt_var)
        volt_watt_basis in (:rated,:available) || throw(ArgumentError("unsupported volt-watt basis"))
        isfinite(power_voltage_floor) && power_voltage_floor>0 || throw(ArgumentError("phase voltage floor must be positive and finite"))
        new(volt_var,volt_watt,volt_watt_basis,_current_target(current_target),
            Float64(power_voltage_floor),_sensing_reference(sensing_reference))
    end
end

sensing_terminal(c::PerPhaseVoltVarWatt,inv::AdvancedInverter)=
    c.sensing_reference===:device ? inv.neutral : c.sensing_reference===:ground ? nothing : c.sensing_reference

function validate_device(c::ControlledDevice{<:AdvancedInverter,<:PerPhaseVoltVarWatt},nets=();periods::Int=length(nets))
    inv=c.device
    inv.topology===:FOUR_LEG && inv.neutral!==nothing || throw(ArgumentError("per-phase control requires FOUR_LEG with an explicit neutral"))
    validate_device(ConstrainedInverter(inv,InverterOperatingConstraints();current_allocation=:free),nets;periods)
    reference=sensing_terminal(c.controller,inv)
    if reference!==nothing
        reference in inv.phase_terminals && throw(ArgumentError("phase terminal cannot be a sensing reference"))
        for net in nets
            bus=get(get(net,"bus",Dict()),inv.bus,nothing)
            bus===nothing && continue
            reference in bus["terminal_names"] || throw(ArgumentError("sensing reference '$reference' is absent"))
        end
    end
    c
end

stamp_device!(ctx,c::ControlledDevice{<:AdvancedInverter,<:PerPhaseVoltVarWatt};request::InverterControlRequest)=
    throw(ArgumentError("PerPhaseVoltVarWatt uses voltage-aware allocation; implicit controller stamping is not implemented"))

solve_controlled_inverter(net::Dict{String,Any},c::ControlledDevice{<:AdvancedInverter,<:PerPhaseVoltVarWatt},
    request::InverterControlRequest;kwargs...)=throw(ArgumentError(
        "PerPhaseVoltVarWatt uses solve_voltage_aware_inverter_fleet or fixed-voltage evaluate_voltage_aware"))

"""
    PhaseInverterServiceAllocationResult

Result of [`allocate_phase_inverter_services`](@ref). `reactive_scale` is the
common fraction of `phase_reactive_request_var`; `stages` retains raw solver
and physical evidence. `current_command` is present only when every policy
stage passed. Plant/policy feasibility is `:feasible` when witnessed, otherwise
`:unresolved`; a failed local solve does not prove infeasibility. A plant
diagnostic, when requested, never replaces the failed controller command.
"""
struct PhaseInverterServiceAllocationResult <: AbstractInverterServiceAllocationResult
    phase_reactive_request_var::NTuple{3,Float64}
    phase_active_fraction::NTuple{3,Float64}
    reactive_scale::Float64
    target::Symbol
    stages::Vector{Pair{Symbol,ConstrainedInverterResult}}
    allocation::ConstrainedInverterResult
    current_command::Union{Nothing,NTuple{3,ComplexF64}}
    plant_feasibility::Symbol
    policy_feasibility::Symbol
    plant_diagnostic::Union{Nothing,InverterCurrentProjectionResult}
end
solve_status(r::PhaseInverterServiceAllocationResult)=solve_status(r.allocation)
physical_feasibility(r::PhaseInverterServiceAllocationResult)=physical_feasibility(r.allocation)
solve_diagnostics(r::PhaseInverterServiceAllocationResult)=(
    phase_reactive_request_var=r.phase_reactive_request_var,phase_active_fraction=r.phase_active_fraction,
    reactive_scale=r.reactive_scale,target=r.target,plant_feasibility=r.plant_feasibility,
    policy_feasibility=r.policy_feasibility,stages=[(name=n,solve=solve_status(s)) for (n,s) in r.stages])

function _phase_service_checks(inv,request,voltage,current,pdc,q,fraction,lambda,basis,available)
    power=voltage.*conj.(current).+im*inv.b_filter_shunt.*abs2.(voltage)
    checks=InverterLimitCheck[]
    for k in 1:3
        push!(checks,_inverter_limit_check("phase_Q_tracking_$k",:VA,abs(imag(power[k])-lambda*q[k]),0.,0.))
        push!(checks,_inverter_limit_check("phase_Q_ray_$k",:VA,imag(power[k]),min(0.,q[k]),max(0.,q[k])))
        push!(checks,_inverter_limit_check("phase_active_balance_$k",:VA,abs(real(power[k])-real(sum(power))/3),0.,0.))
        protected=basis===:rated ? real(power[k]) : (1-fraction[k])*real(power[k])+fraction[k]*pdc/3
        ceiling=fraction[k]*(basis===:rated ? request.p_rated : available)/3
        push!(checks,_inverter_limit_check("phase_volt_watt_$k",:VA,protected,-Inf,ceiling))
    end
    checks
end

"""
    allocate_phase_inverter_services(inv, request; grid_voltage,
        reactive_request, active_fraction=(1.,1.,1.), volt_watt_basis=:rated,
        target=:converter, dc_source=DCSourceBounds(p_min=0.), kwargs...)

Four-leg local Q-vector-ray → equal-phase-P allocation at fixed physical POC
phase-neutral phasors. Hardware/source/continuous rails and phase-policy checks
gate every stage. A failed stage issues no command and proves no global plant or
policy infeasibility. A separate unrestricted projection may supply a physical
witness without replacing the failed policy. Current limits are RMS. The
optional centered-PWM term includes an all-duty AC-ripple-induced DC thermal
reserve; full switching dynamics still require representative EMT validation.
"""
function allocate_phase_inverter_services(inv::AdvancedInverter,request::InverterControlRequest;
        grid_voltage,reactive_request,active_fraction=(1.,1.,1.),volt_watt_basis::Symbol=:rated,
        target::Symbol=:converter,reactive_mode::Symbol=:within_request,
        dc_source::DCSourceBounds=DCSourceBounds(p_min=0.),diagnose_plant::Bool=true,
        s_base::Real=1e5,optimizer=Ipopt.Optimizer,verbose::Bool=false,solver_options=())
    inv.topology===:FOUR_LEG && inv.neutral!==nothing || throw(ArgumentError("phase allocation requires FOUR_LEG with explicit neutral"))
    reactive_mode===:within_request || throw(ArgumentError("phase allocation currently supports within-request Q rays only"))
    target in (:converter,:grid) || throw(ArgumentError("unsupported current target"))
    volt_watt_basis in (:rated,:available) || throw(ArgumentError("unsupported volt-watt basis"))
    isfinite(s_base) && s_base>0 || throw(ArgumentError("s_base must be positive and finite"))
    length(reactive_request)==length(active_fraction)==3 || throw(DimensionMismatch("provide three Q requests and watt fractions"))
    q=Float64.(collect(reactive_request));fraction=Float64.(collect(active_fraction))
    all(isfinite,q) && all(x->isfinite(x) && 0<=x<=1,fraction) || throw(ArgumentError("invalid phase requests"))
    v=_projection_phasors(grid_voltage,"grid_voltage")
    device=ConstrainedInverter(inv,InverterOperatingConstraints();current_allocation=:free,dc_source)
    net=_projection_network(inv,v);validate_device(device,(net,);periods=1)
    handles=Ref{_ConstrainedInverterHandles}();coordinate=Ref{Any}(nothing)
    hook! = ctx -> begin
        handles[]=stamp_device!(ctx,device;request);h=handles[].plant;m=_opf_model(ctx)
        scale=inv.s_max/(3h.sb)
        p=[JuMP.@expression(m,h.vrpoc[k]*h.gri[k]+h.vipoc[k]*h.gii[k]) for k in 1:3]
        reactive=[JuMP.@expression(m,h.vipoc[k]*h.gri[k]-h.vrpoc[k]*h.gii[k]+
            inv.b_filter_shunt*(h.vb/h.ib)*(h.vrpoc[k]^2+h.vipoc[k]^2)) for k in 1:3]
        lambda=all(iszero,q) ? 1. : JuMP.@variable(m,lower_bound=0.,upper_bound=1.,
            start=inv.pwm_dc_voltage_envelope ? 0. : .5)
        coordinate[]=lambda
        for k in 1:3
            JuMP.@constraint(m,(reactive[k]-lambda*q[k]/h.sb)/scale==0)
            k==1 || JuMP.@constraint(m,(p[k]-p[1])/scale==0)
            cap=volt_watt_basis===:rated ? request.p_rated/(3h.sb)*fraction[k] :
                fraction[k]*(_available_dc(device,request)/(3h.sb)-(h.p_dc/3-p[k]))
            JuMP.@constraint(m,(p[k]-cap)/scale<=0)
        end
        JuMP.@objective(m,Max,lambda)
    end
    ctx=build_opf_model(net;per_unit=true,s_base=Float64(s_base),add_objective=false,model_hook! = hook!,optimizer,verbose)
    m=_opf_model(ctx);h=handles[].plant
    options=Dict{String,Any}()
    JuMP.solver_name(m)=="Ipopt" && merge!(options,Dict("tol"=>1e-9,"max_iter"=>1000,
        "mu_strategy"=>(inv.pwm_dc_voltage_envelope ? "monotone" : "adaptive"),"nlp_scaling_method"=>"none"))
    for (name,value) in (solver_options isa NamedTuple ? pairs(solver_options) : solver_options)
        options[string(name)]=value
    end
    _set_inverter_solver_options!(m,options)
    stages=Pair{Symbol,ConstrainedInverterResult}[];lambda_value=NaN
    function capture(name)
        outcome=_solve_outcome(m);status=SolveStatus(outcome)
        e=extract_device(device,handles[],status);report=e.plant.physical_feasibility
        if status.has_primal && report.checked
            actual_v=complex.(JuMP.value.(h.vrpoc),JuMP.value.(h.vipoc)).*h.vb
            actual_i=complex.(JuMP.value.(h.gri),JuMP.value.(h.gii)).*h.ib
            lambda_value=coordinate[] isa Real ? coordinate[] : JuMP.value(coordinate[])
            pdc=only(c.value for c in report.checks if c.name=="dc_source_power_max")
            report=_projection_boundary_report(report,actual_v,v,0.,0.)
            report=_extend_inverter_feasibility(report,_phase_service_checks(inv,request,
                actual_v,actual_i,pdc,q,fraction,lambda_value,volt_watt_basis,_available_dc(device,request)))
        end
        plant=merge(e.plant,(physical_feasibility=report,))
        operation=merge(e.operation,(allocation_stage=name,reactive_service_scale=lambda_value,
            raw_solver_status=JuMP.raw_status(m),physical_failed_checks=[c.name for c in report.checks if !c.passed]))
        network=_extract_result(ctx,outcome)
        result=ConstrainedInverterResult(string(outcome.termination_status),plant,operation,
            e.converter_terminal,e.grid_phase_current,network["bus"],_physical_solve_status(status,report))
        push!(stages,name=>result);result
    end
    _optimize_inverter_rails!(ctx,((inv,h),));result=capture(:phase_reactive_tracking)
    if solve_status(result).publishable
        starts=[(x,JuMP.value(x)) for x in JuMP.all_variables(m)]
        if !(coordinate[] isa Real)
            if inv.pwm_strategy===:NONE && inv.dc_energy_model===:linearized
                JuMP.fix(coordinate[],lambda_value;force=true)
            else
                _lock_pwm_service!(m,coordinate[],lambda_value,sum(abs,q))
            end
        end
        JuMP.@objective(m,Max,h.p_poc/(inv.s_max/h.sb))
        for (x,value) in starts;JuMP.set_start_value(x,value);end
        _optimize_inverter_rails!(ctx,((inv,h),));result=capture(:phase_active_export)
    end
    accepted=solve_status(result).publishable;witnessed=any(solve_status(s).publishable for (_,s) in stages)
    diagnostic=nothing
    if !witnessed && diagnose_plant
        desired=[abs(v[k])>1e-9 ? conj(complex(request.p_available/3,q[k])/v[k]) : 0im for k in 1:3]
        diagnostic=project_inverter_current(inv,request;grid_voltage=v,requested_current=desired,
            target=:grid,dc_source,s_base,optimizer,verbose,solver_options=options)
        witnessed=diagnostic.plant_feasibility===:feasible
    end
    command=accepted ? (target===:converter ? result.converter_terminal.phase_current : result.grid_phase_current) : nothing
    PhaseInverterServiceAllocationResult(Tuple(q),Tuple(fraction),lambda_value,target,stages,result,command,
        witnessed ? :feasible : :unresolved,accepted ? :feasible : :unresolved,diagnostic)
end

function _phase_voltage_requests(c::PerPhaseVoltVarWatt,measurement,request)
    magnitudes=abs.(measurement.phase_voltage)
    all(v->v>=c.power_voltage_floor,magnitudes) || throw(ArgumentError("a sensed phase is below power_voltage_floor"))
    q=ntuple(k->c.volt_var===nothing ? 0. : request.q_scale/3*_pwl_exact(c.volt_var,magnitudes[k]),3)
    f=ntuple(k->c.volt_watt===nothing ? 1. : _pwl_exact(c.volt_watt,magnitudes[k]),3)
    q,f
end

function evaluate_voltage_aware(controlled::ControlledDevice{<:AdvancedInverter,<:PerPhaseVoltVarWatt},
        measurement::InverterControlMeasurement,request::InverterControlRequest;grid_voltage,kwargs...)
    validate_device(controlled)
    q,f=_phase_voltage_requests(controlled.controller,measurement,request)
    target=controlled.controller.current_target isa ConverterCurrentTarget ? :converter : :grid
    allocate_phase_inverter_services(controlled.device,request;grid_voltage,reactive_request=q,
        active_fraction=f,volt_watt_basis=controlled.controller.volt_watt_basis,target,kwargs...)
end
