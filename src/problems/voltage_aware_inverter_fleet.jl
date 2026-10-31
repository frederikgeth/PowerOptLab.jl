# A quasi-static numerical fixed point of the physical local allocation map and
# the native network equations. Virtual feeder probes are not transient dispatch.

"""Result of [`solve_voltage_aware_inverter_fleet`](@ref).

`devices` contains local Q→P→I₂ allocations at the final verification voltages.
`replayed` and `physical_feasibility(result)` describe the imposed grid currents
at the returned **network** voltages. `network_equation_checks` evaluates native
zero-equality rows in SI; it is a model residual audit, not an independent circuit
implementation. `reason` describes the iteration outcome; `solve` retains the
last raw solver evidence and additionally gates publication on convergence,
physical reports and policy/current closure. Disconnected devices have no local
allocation or inverter report because their entire AC model is absent.
"""
struct VoltageAwareInverterFleetResult <: AbstractSolveResult
    termination_status::String
    reason::Symbol
    converged::Bool
    iterations::Int
    devices::Dict{String,AbstractInverterServiceAllocationResult}
    spec::ControlledInverterFleetSpec
    operating_modes::Dict{String,Symbol}
    grid_currents::Dict{String,NTuple{3,ComplexF64}}
    replayed::Dict{String,NamedTuple}
    physical_reports::Dict{String,InverterFeasibilityReport}
    network::Dict{String,Any}
    network_equation_checks::Vector{InverterLimitCheck}
    history::Vector{NamedTuple}
    solve::SolveStatus
end

solve_status(r::VoltageAwareInverterFleetResult)=r.solve
physical_feasibility(r::VoltageAwareInverterFleetResult)=r.physical_reports
solve_diagnostics(r::VoltageAwareInverterFleetResult)=(reason=r.reason,converged=r.converged,
    iterations=r.iterations,operating_modes=r.operating_modes,
    failed_physical_devices=[id for (id,p) in r.physical_reports if !p.passed],
    failed_network_equations=[c.name for c in r.network_equation_checks if !c.passed],history=r.history)

function _voltage_aware_modes(ids,modes)
    if modes===nothing || modes isa Symbol
        value=modes===nothing ? :energized : modes
        result=Dict(id=>value for id in ids)
    elseif modes isa AbstractDict
        result=Dict(String(id)=>mode for (id,mode) in modes)
        length(result)==length(modes) || throw(ArgumentError("duplicate operating-mode ids after string conversion"))
        Set(keys(result))==Set(ids) || throw(ArgumentError("operating_modes keys must match fleet ids"))
    else
        throw(ArgumentError("operating_modes must be a mode symbol or a per-device dictionary"))
    end
    all(mode in (:energized,:disconnected) for mode in values(result)) ||
        throw(ArgumentError("operating modes must be :energized or :disconnected"))
    Dict{String,Symbol}(result)
end

# Native equality rows are audited independently of MOI's success flags. The
# largest coordinate base is conservative for blocks spanning several buses.
# This evaluates registered model rows; circuit replay remains separate evidence.
function _voltage_aware_network_checks(ctx,status)
    checks=InverterLimitCheck[]
    status.has_primal || return checks,false
    schema=BMOPFTools.opf_diagnostic_schema(ctx)
    get(schema.capabilities,"semantic_blocks_registered",false) || return checks,false
    ac=schema.coordinate_bases["ac"]
    allbases=vcat(collect(values(ac)),collect(values(schema.coordinate_bases["dc"])))
    scales=Dict(:current=>maximum(b.current for b in allbases),
        :voltage=>maximum(b.voltage for b in allbases),:power=>maximum(b.power for b in allbases))
    kcl_count=0
    for block in schema.semantic_blocks
        block.kind===:constraint && block.set_contract===:zero_equality || continue
        haskey(scales,block.quantity) || continue
        constraints=[JuMP.constraint_object(BMOPFTools.opf_object(ctx,key)) for key in block.members]
        all(c.set isa JuMP.MOI.EqualTo for c in constraints) || continue
        residual=block.model_to_canonical*[JuMP.value(c.func)-c.set.value for c in constraints]
        scale=scales[block.quantity]
        if all(k.family in (:kcl_r,:kcl_i) for k in block.members)
            bus=String(first(block.members).index[1])
            scale=ac[bus].current
            kcl_count+=1
        end
        for (component,value) in zip(block.components,residual)
            push!(checks,_inverter_limit_check("network/$(block.id)/$component",block.physical_unit,
                abs(value)*scale,0.,0.))
        end
    end
    checks,kcl_count>0 && !isempty(checks)
end

function _solve_voltage_aware_current_network(net,spec,currents,modes;
        s_base,optimizer,verbose,solver_options)
    builder=BMOPFTools.OpfDeviceBuilder(:PowerOptLab,function(ctx,ids)
        for id in ids
            _disable_native_ibr_port!(ctx,id)
            modes[id]===:disconnected && continue
            inv=inverter_spec(spec.devices[id]);bases=BMOPFTools.opf_ac_coordinate_bases(ctx,inv.bus)
            for (k,phase) in enumerate(inv.phase_terminals)
                cr,ci=reim(currents[id][k]/bases.current)
                vr,vi=_dv(ctx,inv.bus,phase,inv.neutral)
                b=inv.b_filter_shunt*bases.voltage/bases.current
                ir=cr+b*vi;ii=ci-b*vr
                BMOPFTools.add_terminal_injection!(ctx,inv.bus,phase,ir,ii)
                inv.neutral===nothing || BMOPFTools.add_terminal_injection!(ctx,inv.bus,inv.neutral,-ir,-ii)
            end
        end
        ctx
    end)
    build_spec=BMOPFTools.OpfBuildSpec(component_builders=Dict((:ibr,id)=>builder for id in keys(spec.devices)))
    ctx=build_opf_model(net;per_unit=true,s_base=Float64(s_base),add_objective=false,
        build_spec,optimizer,verbose,model_hook! = ctx -> JuMP.@objective(_opf_model(ctx),Min,0.))
    m=_opf_model(ctx)
    if JuMP.solver_name(m)=="Ipopt"
        _set_solver_options!(m,("tol"=>1e-9,"max_iter"=>1000,"mu_strategy"=>"adaptive","nlp_scaling_method"=>"none"))
    end
    _set_inverter_solver_options!(m,solver_options)
    BMOPFTools.optimize_opf!(ctx)
    outcome=_solve_outcome(m);status=SolveStatus(outcome)
    checks,checked=_voltage_aware_network_checks(ctx,status)
    accepted=SolveStatus(status.termination_status,status.primal_status,status.has_primal,
        status.feasible,status.optimal,status.publishable && checked && all(c.passed for c in checks))
    network=_network_without_replaced_ibrs(_extract_result(ctx,outcome),keys(spec.devices))
    network["raw_solver_status"]=JuMP.raw_status(m)
    network["equation_audit_scope"]="registered native zero-equality rows in SI; not independent circuit replay"
    (status=accepted,network=network,checks=checks,checked=checked)
end

function _voltage_aware_replayed(inv,voltage,grid_current,source,request)
    state=_projection_requested_state(inv,collect(voltage),collect(grid_current),:grid,0.)
    state===nothing && return nothing
    available=min(request.p_available,something(source.p_max,request.p_available))
    report=inverter_physical_feasibility(inv;state...,dc_min=source.p_min,dc_max=available)
    power=sum(state.grid_voltage.*conj.(state.grid_current))+im*inv.b_filter_shunt*sum(abs2,state.grid_voltage)
    conv=sum(state.converter_voltage.*conj.(state.converter_current))
    pdc=only(c.value for c in report.checks if c.name=="dc_source_power_max")
    (;state,physical_feasibility=report,p_poc=real(power),q_poc=imag(power),p_conv=real(conv),
        q_conv=imag(conv),p_dc=pdc,converter_terminal=_converter_terminal_result(
            Tuple(state.converter_voltage),Tuple(state.converter_current)))
end

function _voltage_aware_coupling_report(point,controlled,request,allocation,sensed,mode)
    inv=controlled.device;c=controlled.controller
    actual_current=c.current_target isa ConverterCurrentTarget ?
        point.state.converter_current : point.state.grid_current
    local_current=allocation.current_command
    local_power=sum(allocation.allocation.plant.poc_phase_voltage.*conj.(allocation.allocation.grid_phase_current))+
        im*inv.b_filter_shunt*sum(abs2,allocation.allocation.plant.poc_phase_voltage)
    checks=InverterLimitCheck[]
    for k in 1:3
        push!(checks,_inverter_limit_check("equilibrium_current_$k",:A,abs(actual_current[k]-local_current[k]),0.,0.))
    end
    for (label,value,reference) in (("P",point.p_poc,real(local_power)),("Q",point.q_poc,imag(local_power)))
        push!(checks,_inverter_limit_check("equilibrium_$label",:VA,abs(value-reference),0.,0.))
    end
    _,u1,_=_phasor_sequences(sensed)
    fraction,q,_,_=_positive_exact(c.positive,u1,abs.(sensed),request)
    mode===:within_request && push!(checks,_inverter_limit_check("equilibrium_Q_ray",:VA,point.q_poc,min(0.,q),max(0.,q)))
    power=_volt_watt_basis(c.positive)===:rated ? point.p_poc :
        (1-fraction)*point.p_poc+fraction*point.p_dc
    ceiling=fraction*(_volt_watt_basis(c.positive)===:rated ? request.p_rated :
        allocation.allocation.operation.effective_dc_availability_W)
    push!(checks,_inverter_limit_check("equilibrium_volt_watt",:VA,power,-Inf,ceiling))
    target2=_phasor_sequences(Tuple(actual_current))[3]
    push!(checks,_inverter_limit_check("equilibrium_I2",:A,
        abs(target2-allocation.negative_scale*allocation.negative_sequence_request_A),0.,0.))
    _extend_inverter_feasibility(point.physical_feasibility,checks)
end

function _voltage_aware_coupling_report(point,
        controlled::ControlledDevice{<:AdvancedInverter,<:PerPhaseVoltVarWatt},request,
        allocation::PhaseInverterServiceAllocationResult,sensed,mode)
    inv=controlled.device;c=controlled.controller
    actual=c.current_target isa ConverterCurrentTarget ? point.state.converter_current : point.state.grid_current
    checks=InverterLimitCheck[_inverter_limit_check("equilibrium_current_$k",:A,
        abs(actual[k]-allocation.current_command[k]),0.,0.) for k in 1:3]
    local_v=collect(allocation.allocation.plant.poc_phase_voltage)
    local_i=collect(allocation.allocation.grid_phase_current)
    reference=local_v.*conj.(local_i).+im*inv.b_filter_shunt.*abs2.(local_v)
    power=point.state.grid_voltage.*conj.(point.state.grid_current).+
        im*inv.b_filter_shunt.*abs2.(point.state.grid_voltage)
    for k in 1:3,(label,value) in (("P",real(power[k]-reference[k])),("Q",imag(power[k]-reference[k])))
        push!(checks,_inverter_limit_check("equilibrium_phase_$(label)_$k",:VA,abs(value),0.,0.))
    end
    q,f=_phase_voltage_requests(c,InverterControlMeasurement(sensed),request)
    append!(checks,_phase_service_checks(inv,request,point.state.grid_voltage,point.state.grid_current,
        point.p_dc,q,f,allocation.reactive_scale,c.volt_watt_basis,
        allocation.allocation.operation.effective_dc_availability_W))
    _extend_inverter_feasibility(point.physical_feasibility,checks)
end

function _voltage_aware_voltage_residual(old,new,spec,ids)
    residual=0.
    for id in ids
        inv=inverter_spec(spec.devices[id])
        for terminal in keys(old[inv.bus])
            left=old[inv.bus][terminal];right=new[inv.bus][terminal]
            residual=max(residual,abs(complex(left["vr"],left["vi"])-complex(right["vr"],right["vi"])))
        end
        for device in (inv,spec.devices[id])
            residual=max(residual,maximum(abs,_network_phase_state(old,device).-_network_phase_state(new,device)))
        end
    end
    residual
end

function _voltage_aware_blend(old,new,alpha)
    mixed=deepcopy(new)
    for (bus,terminals) in mixed,(terminal,data) in terminals
        left=old[bus][terminal]
        v=(1-alpha)*complex(left["vr"],left["vi"])+alpha*complex(data["vr"],data["vi"])
        data["vr"]=real(v);data["vi"]=imag(v);data["vm"]=abs(v);data["va"]=angle(v)
    end
    mixed
end

"""
    solve_voltage_aware_inverter_fleet(net, spec; operating_modes=:energized,
        initial_currents=nothing, relaxation=1., max_iterations=30,
        atol_voltage=1e-7, atol_current=1e-6, reactive_mode=:within_request,
        dc_source=DCSourceBounds(p_min=0.), s_base=1e5,
        optimizer=Ipopt.Optimizer, verbose=false, solver_options=(), progress=nothing)

Couple [`evaluate_voltage_aware`](@ref) to full native network power flows.
Each local allocation uses separate physical and sensed voltages. The feeder map
imposes allocated **grid-arm** currents plus the declared POC shunt. Intermediate
network probes omit inverter hardware constraints and are explicitly virtual;
otherwise a saturated rail could prevent evaluation of a perturbed iterate.
Every returned candidate is reconstructed and audited at the actual feeder voltage.
Acceptance requires voltage and current-map closure, local stage success, native
SI equation checks and full inverter/source/policy checks. Solver failures,
iteration limits and controller voltage-domain failures remain unpublishable.

`operating_modes` is a common symbol or a complete per-device dictionary.
`:energized` retains active controls and loss imports at zero PV. `:disconnected`
means an open grid-side contactor: the entire converter, filter and POC shunt are
absent from the AC model. DC standby losses are outside this disconnected model.
All-disconnected operation is a pure network power flow with zero iterations.
Unselected native IBRs remain native. Inputs are not mutated; inverter sizing is
unchanged. Supports three-leg sequence control and four-leg per-phase Q→P
allocation; the latter has equal POC phase P and a shared DC source.

This is a quasi-static numerical fixed point, not an EMT/time-domain controller
simulation, global uniqueness proof or latency/stability certificate. Registered
native equation residuals are not an independent implementation of the circuit.
`progress`, when supplied, receives each stored iteration diagnostics tuple.
"""
function solve_voltage_aware_inverter_fleet(net::Dict{String,Any},spec::ControlledInverterFleetSpec;
        operating_modes=:energized,initial_currents=nothing,relaxation::Real=1.,max_iterations::Integer=30,
        atol_voltage::Real=1e-7,atol_current::Real=1e-6,reactive_mode::Symbol=:within_request,
        dc_source::DCSourceBounds=DCSourceBounds(p_min=0.),s_base::Real=1e5,
        optimizer=Ipopt.Optimizer,verbose::Bool=false,solver_options=(),progress=nothing)
    _validate_controlled_fleet(net,spec)
    0<relaxation<=1 && isfinite(relaxation) || throw(ArgumentError("relaxation must lie in (0,1]"))
    max_iterations>=1 || throw(ArgumentError("max_iterations must be positive"))
    all(isfinite(x) && x>0 for x in (atol_voltage,atol_current,s_base)) ||
        throw(ArgumentError("voltage/current tolerances and s_base must be positive and finite"))
    reactive_mode in (:within_request,:protective_override) || throw(ArgumentError("unsupported reactive_mode"))
    ids=sort!(collect(keys(spec.devices)));modes=_voltage_aware_modes(ids,operating_modes)
    active=[id for id in ids if modes[id]===:energized]
    for id in active
        controlled=spec.devices[id];c=controlled.controller
        if c isa SequenceController
            c.limiter.pq_priority===:var || throw(ArgumentError("voltage-aware coupling requires reactive priority"))
            c.limiter isa SequencePriorityLimiter && c.limiter.priority!==:power &&
                throw(ArgumentError("voltage-aware coupling requires Q-P-I2 order"))
            c.unbalance isa NegativeSequenceAdmittanceDroop && !iszero(c.unbalance.ripple_blend) &&
                throw(ArgumentError("current-dependent ripple compensation is not supported"))
        elseif c isa PerPhaseVoltVarWatt
            reactive_mode===:within_request || throw(ArgumentError("phase allocation supports within-request Q only"))
        else
            throw(ArgumentError("unsupported voltage-aware controller"))
        end
        validate_device(ConstrainedInverter(controlled.device,InverterOperatingConstraints();dc_source),(net,);periods=1)
    end
    zero=ntuple(_->0.0im,3)
    currents=initial_currents===nothing ? Dict(id=>zero for id in ids) : first(_validated_fleet_currents(spec,initial_currents))
    for id in ids
        modes[id]===:disconnected && (currents[id]=zero)
        native=net["ibr"][id]
        if uppercase(String(get(native,"topology","FOUR_LEG")))=="THREE_LEG" && abs(sum(currents[id]))>1e-5
            throw(ArgumentError("initial grid current for three-wire '$id' must sum to zero"))
        end
    end
    solved=_solve_voltage_aware_current_network(net,spec,currents,modes;s_base,optimizer,verbose,solver_options)
    devices=Dict{String,AbstractInverterServiceAllocationResult}();points=Dict{String,NamedTuple}()
    reports=Dict(id=>_unchecked_inverter_feasibility() for id in active)
    history=NamedTuple[];converged=false;iterations=0;reason=:NETWORK_SOLVE_FAILED
    if solved.status.publishable
        for id in active
            point=_voltage_aware_replayed(spec.devices[id].device,
                _network_phase_state(solved.network["bus"],spec.devices[id].device),currents[id],dc_source,spec.requests[id])
            point===nothing && continue
            points[id]=point;reports[id]=point.physical_feasibility
        end
    end
    function allocations(bus)
        result=Dict{String,AbstractInverterServiceAllocationResult}()
        for id in active
            controlled=spec.devices[id]
            physical=_network_phase_state(bus,controlled.device);sensed=_network_phase_state(bus,controlled)
            c=controlled.controller
            voltage=c isa PerPhaseVoltVarWatt ? minimum(abs,sensed) : abs(_phasor_sequences(sensed)[2])
            voltage>=c.power_voltage_floor || return result,:CONTROLLER_VOLTAGE_DOMAIN
            allocation=evaluate_voltage_aware(controlled,InverterControlMeasurement(sensed),spec.requests[id];
                grid_voltage=physical,reactive_mode,dc_source,s_base,optimizer,verbose,solver_options)
            result[id]=allocation
            solve_status(allocation).publishable || return result,:ALLOCATION_FAILED
        end
        result,:OK
    end
    if solved.status.publishable && isempty(active)
        converged=true;reason=:DISCONNECTED_POWER_FLOW
    elseif solved.status.publishable
        measured=deepcopy(solved.network["bus"])
        reason=:ITERATION_LIMIT
        for iteration in 1:max_iterations
            iterations=iteration
            devices,failure=allocations(measured)
            if failure!==:OK
                reason=failure
                failed=first(id for id in active if !haskey(devices,id) || !solve_status(devices[id]).publishable)
                row=(iteration=iteration,reason=reason,failed_allocation_ids=[failed],
                    voltage_residual_inf_V=NaN,current_residual_inf_A=NaN,
                    failed_physical_devices=count(r->!r.passed,values(reports)),network_accepted=solved.status.publishable)
                push!(history,row);progress===nothing || progress(row)
                break
            end
            applied=Dict(id=>modes[id]===:disconnected ? zero : devices[id].allocation.grid_phase_current for id in ids)
            next=_solve_voltage_aware_current_network(net,spec,applied,modes;s_base,optimizer,verbose,solver_options)
            solved=next;currents=applied
            if !solved.status.publishable
                reason=:NETWORK_SOLVE_FAILED;empty!(points)
                for id in active;reports[id]=_unchecked_inverter_feasibility();end
                row=(iteration=iteration,reason=reason,failed_allocation_ids=String[],
                    voltage_residual_inf_V=NaN,current_residual_inf_A=NaN,
                    failed_physical_devices=length(active),network_accepted=false)
                push!(history,row);progress===nothing || progress(row)
                break
            end
            bus=solved.network["bus"]
            voltage_residual=_voltage_aware_voltage_residual(measured,bus,spec,active)
            empty!(points)
            for id in active
                reports[id]=_unchecked_inverter_feasibility()
                point=_voltage_aware_replayed(spec.devices[id].device,_network_phase_state(bus,spec.devices[id].device),
                    currents[id],dc_source,spec.requests[id])
                point===nothing && continue
                points[id]=point;reports[id]=point.physical_feasibility
            end
            current_residual=Inf
            if voltage_residual<=atol_voltage
                verification,failure=allocations(bus)
                if failure!==:OK
                    devices=verification;reason=failure
                    failed=first(id for id in active if !haskey(devices,id) || !solve_status(devices[id]).publishable)
                    row=(iteration=iteration,reason=reason,failed_allocation_ids=[failed],
                        voltage_residual_inf_V=voltage_residual,current_residual_inf_A=NaN,
                        failed_physical_devices=count(r->!r.passed,values(reports)),network_accepted=true)
                    push!(history,row);progress===nothing || progress(row)
                    break
                end
                current_residual=maximum(maximum(abs,verification[id].allocation.grid_phase_current.-currents[id]) for id in active)
                for id in active
                    haskey(points,id) || continue
                    sensed=_network_phase_state(bus,spec.devices[id])
                    reports[id]=_voltage_aware_coupling_report(points[id],spec.devices[id],spec.requests[id],verification[id],
                        sensed,reactive_mode)
                    points[id]=merge(points[id],(physical_feasibility=reports[id],))
                end
                if current_residual<=atol_current && length(points)==length(active) && all(r.passed for r in values(reports))
                    devices=verification;converged=true;reason=:CONVERGED
                end
            end
            row=(iteration=iteration,reason=converged ? :CONVERGED : :ITERATING,failed_allocation_ids=String[],
                voltage_residual_inf_V=voltage_residual,current_residual_inf_A=current_residual,
                failed_physical_devices=count(r->!r.passed,values(reports)),network_accepted=solved.status.publishable)
            push!(history,row);progress===nothing || progress(row)
            converged && break
            measured=_voltage_aware_blend(measured,bus,Float64(relaxation))
        end
    end
    raw=solved.status
    if reason===:ALLOCATION_FAILED && !isempty(devices)
        failed=findfirst(id->!solve_status(devices[id]).publishable,sort!(collect(keys(devices))))
        failed===nothing || (raw=solve_status(devices[sort!(collect(keys(devices)))[failed]]))
    end
    status=SolveStatus(raw.termination_status,raw.primal_status,raw.has_primal,raw.feasible,raw.optimal,
        raw.publishable && converged && all(r.passed for r in values(reports)))
    VoltageAwareInverterFleetResult(raw.termination_status,reason,converged,iterations,devices,deepcopy(spec),
        modes,currents,points,reports,solved.network,solved.checks,history,status)
end

"""Deterministic SI rows for a voltage-aware coupled fleet.

Only `accepted=true` rows belong in scientific equilibrium comparisons.
Disconnected AC powers/currents are zero; DC power and inverter certificate fields
are `nothing` because that device model is absent. Voltages remain customer-bus
observables, not a claim that disconnected electronics are measuring them.
"""
function voltage_aware_inverter_rows(result::VoltageAwareInverterFleetResult)
    rows=NamedTuple[]
    for id in sort!(collect(keys(result.spec.devices)))
        controlled=result.spec.devices[id];inv=controlled.device;request=result.spec.requests[id]
        active=result.operating_modes[id]===:energized
        point=get(result.replayed,id,nothing)
        report=get(result.physical_reports,id,nothing)
        allocation=get(result.devices,id,nothing)
        current=point===nothing ? nothing : point.state.converter_current
        grid=point===nothing ? nothing : point.state.grid_current
        check_value(name)=report===nothing ? nothing :
            (c=findfirst(c->c.name==name,report.checks);c===nothing ? nothing : report.checks[c].value)
        sensed=try
            _network_phase_state(result.network["bus"],controlled)
        catch err
            err isa ArgumentError || rethrow()
            nothing
        end
        vuf=sensed===nothing ? NaN : (u=_phasor_sequences(sensed);iszero(abs(u[2])) ? NaN : 100abs(u[3])/abs(u[2]))
        power(name)=!active ? 0. : point===nothing ? NaN : getproperty(point,name)
        push!(rows,(device_id=id,bus=inv.bus,operating_mode=result.operating_modes[id],model_present=active,
            accepted=result.solve.publishable,reason=result.reason,iterations=result.iterations,
            physical_checked=report===nothing ? nothing : report.checked,
            physical_feasible=report===nothing ? nothing : report.passed,
            p_available_W=request.p_available,p_poc_W=power(:p_poc),q_poc_var=power(:q_poc),
            p_dc_W=!active ? nothing : power(:p_dc),
            converter_loss_W=!active ? nothing : point===nothing ? NaN : point.p_dc-point.p_conv,
            filter_loss_W=!active ? nothing : point===nothing ? NaN : point.p_conv-point.p_poc,
            maximum_converter_current_A=!active ? 0. : current===nothing ? NaN : maximum(abs,current),
            maximum_grid_current_A=!active ? 0. : grid===nothing ? NaN : maximum(abs,grid),
            neutral_total_rms_current_A=!active ? 0. : current===nothing ? NaN :
                hypot(abs(sum(current)),inv.pwm_ac_neutral_reserve),
            neutral_current_rating_A=inv.In_max,
            converter_negative_sequence_A=!active ? 0. : current===nothing ? NaN : abs(_phasor_sequences(Tuple(current))[3]),
            capacitor_thermal_current_A=check_value("capacitor_current"),
            dc_ripple_voltage_V=check_value("dc_voltage_ripple"),
            continuous_rail_lower_V=report===nothing ? nothing : report.rail.lower_V,
            vmin_sensed_V=sensed===nothing ? NaN : minimum(abs,sensed),
            vmax_sensed_V=sensed===nothing ? NaN : maximum(abs,sensed),vuf_percent=vuf,
            reactive_request_var=allocation===nothing ? nothing : allocation isa PhaseInverterServiceAllocationResult ?
                sum(allocation.phase_reactive_request_var) : allocation.reactive_request_var,
            negative_sequence_scale=allocation isa InverterServiceAllocationResult ? allocation.negative_scale : nothing,
            phase_reactive_service_scale=allocation isa PhaseInverterServiceAllocationResult ? allocation.reactive_scale : nothing,
            local_termination=allocation===nothing ? "NOT_RUN" : solve_status(allocation).termination_status))
    end
    rows
end
