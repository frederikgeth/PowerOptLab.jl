# Saved-phasor reconstruction for the matched benchmark; no solver auxiliaries.
include("coordinated_engineering.jl")

engineering_hardware_fields(inv)=[[string(k),_engineering_json_safe(getfield(inv,k))] for k in fieldnames(AdvancedInverter)]
engineering_hardware_sha256(inv)=bytes2hex(sha256(JSON3.write(engineering_hardware_fields(inv))))

const ENGINEERING_LIMIT_NAMES=Set(["capacitor_current","dc_voltage_min","dc_voltage_max","dc_voltage_ripple",
    "dc_source_power_min","dc_source_power_max","poc_apparent_power","converter_apparent_power",
    "neutral_current","neutral_peak_current"])

"""Names of hardware/source limits within their activity threshold in a physical report."""
function engineering_active_limits(report)
    bound(c)=c.name in ENGINEERING_LIMIT_NAMES || any(startswith(c.name,prefix) for prefix in
        ("converter_current_","grid_current_","converter_peak_current_","grid_peak_current_"))
    active=String[]
    for c in report.checks
        bound(c) || continue
        threshold=c.unit===:A ? .001 : c.unit===:V ? .01 : .05
        c.margin<=threshold && push!(active,c.name)
    end
    report.rail.lower_V<=.01 && push!(active,"modulation")
    active
end

# Local control and the equal-service OPF must both satisfy the local VV/VW policy.
function engineering_record_metrics(point,mixed;mode,validate_local=mode in ("local","opf_equal_service"))
    phasors(x)=ComplexF64[complex(z...) for z in x]
    case_id=point["case"] isa NamedTuple ? point["case"].case_id : point["case"]["case_id"]
    rows=NamedTuple[];reports=Dict{String,Any}();policies=Dict{String,Bool}()
    for id in sort!(collect(keys(mixed.spec.devices)))
        controlled=mixed.spec.devices[id];inv=controlled.device;request=mixed.spec.requests[id]
        record=point["devices"][id];ts=point["bus"][inv.bus]
        current=phasors(record["converter_current"]);grid=phasors(record["grid_current"])
        u=phasors(record["converter_voltage"])
        v=ComplexF64[phasor(ts[p])-(inv.neutral===nothing ? 0im : phasor(ts[inv.neutral])) for p in inv.phase_terminals]
        sensed=PowerOptLab._network_phase_state(point["bus"],controlled)
        report=inverter_physical_feasibility(inv;converter_voltage=u,converter_current=current,
            grid_voltage=v,grid_current=grid,dc_min=0.,dc_max=request.p_available,
            pwm_thermal_headroom=PowerOptLab._PWM_ALLOCATION_HEADROOM)
        reports[id]=report
        value(name)=begin
            found=findfirst(c->c.name==name,report.checks)
            found===nothing ? nothing : report.checks[found].value
        end
        power=sum(v.*conj.(grid));converted=sum(u.*conj.(current));pdc=value("dc_source_power_max")
        sequence=PowerOptLab._phasor_sequences(Tuple(sensed))
        phase_power=sensed.*conj.(grid)
        if controlled.controller isa PerPhaseVoltVarWatt
            q,f=PowerOptLab._phase_voltage_requests(controlled.controller,InverterControlMeasurement(sensed),request)
            qactual=imag.(phase_power);den=sum(abs2,q)
            lambda=iszero(den) ? 1. : sum(qactual.*q)/den
            qerror=sum(abs.(qactual.-q));qrequested=sum(q);qabs=sum(abs,q)
            policy=-1e-6<=lambda<=1+1e-6 && maximum(abs.(qactual.-lambda.*q))<=.01 &&
                all(real(phase_power[k])<=request.p_rated*f[k]/3+.001 for k in 1:3)
            capability=maximum(abs.(real.(phase_power).-real(phase_power[1])))<=.001
        else
            fraction,qrequested,_,_=PowerOptLab._positive_exact(controlled.controller.positive,
                sequence[2],abs.(sensed),request)
            f=(fraction,fraction,fraction);qabs=abs(qrequested);qerror=abs(imag(power)-qrequested)
            lambda=iszero(qabs) ? 1. : imag(power)/qrequested
            policy=min(0.,qrequested)-.001<=imag(power)<=max(0.,qrequested)+.001 &&
                real(power)<=request.p_rated*fraction+.001
            capability=abs(PowerOptLab._phasor_sequences(Tuple(current))[3])<=1.01e-5
        end
        policies[id]=capability && (!validate_local || policy)
        envelope=PowerOptLab._inverter_current_envelope(inv)
        phase_rms=maximum(hypot.(abs.(current),envelope.converter_rms))
        phase_peak=maximum(sqrt(2).*abs.(current).+envelope.phase_peak)
        neutral_rms=hypot(abs(sum(current)),envelope.neutral_rms)
        neutral_peak=sqrt(2)*abs(sum(current))+envelope.neutral_peak
        active=engineering_active_limits(report)
        planning=point["planning_voltage_max_V"]
        planning!==nothing && maximum(abs,sensed)>=planning-.01 && push!(active,"planning_voltage")
        push!(rows,(case_id=case_id,mode=mode,device_id=id,bus=inv.bus,
            local_controller=mixed.labels[id],topology=inv.topology,
            hardware_sha256=engineering_hardware_sha256(inv),
            p_available_W=request.p_available,p_poc_W=real(power),q_poc_var=imag(power),
            p_dc_W=pdc,curtailment_W=iszero(request.p_available) ? 0. : max(0.,request.p_available-pdc),
            converter_loss_W=pdc-real(converted),filter_loss_W=real(converted-power),
            q_request_var=qrequested,q_request_absolute_var=qabs,reactive_service_error_var=qerror,
            reactive_request_fraction=iszero(qabs) ? nothing : clamp(lambda,0.,1.),
            volt_watt_fraction_min=minimum(f),vmax_sensed_V=maximum(abs,sensed),
            vuf_percent=100abs(sequence[3])/abs(sequence[2]),
            phase_rms_bound_A=phase_rms,phase_peak_bound_A=phase_peak,
            neutral_rms_bound_A=neutral_rms,neutral_peak_bound_A=neutral_peak,
            capacitor_thermal_A=value("capacitor_current"),
            dc_minimum_bound_V=value("dc_voltage_min"),dc_maximum_bound_V=value("dc_voltage_max"),
            rail_margin_V=report.rail.lower_V,binding_constraints=join(sort!(unique!(active)),"|"),
            physical_accepted=report.passed,current_capability_accepted=capability,
            local_policy_accepted=validate_local ? policy : nothing))
    end
    (rows=rows,reports=reports,policies=policies,
        accepted=all(r.passed for r in values(reports)) && all(values(policies)))
end

function engineering_local_stages(result)
    Dict(id=>[Dict("stage"=>String(name),"solve"=>solve_status(stage),
        "raw_solver_status"=>stage.operation.raw_solver_status,
        "plant"=>stage.plant,"operation"=>stage.operation,
        "converter_terminal"=>stage.converter_terminal,"grid_phase_current"=>stage.grid_phase_current)
        for (name,stage) in allocation.stages] for (id,allocation) in result.devices)
end
