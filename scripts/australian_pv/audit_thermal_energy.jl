include("thermal_energy.jl")
using Test

function audit_mixed_thermal(output=joinpath(ROOT,"results","thermal_energy"))
    target=joinpath(output,"audit.json")
    writejson(target,Dict("all_planned_accepted"=>false,"reason"=>"AUDIT_INCOMPLETE"))
    count=Ref(0);physical=0
    function require(condition,message)
        count[]+=1;condition || throw(ArgumentError(message))
    end
    try
        env=JSON3.read(read(joinpath(output,"environment.json"),String),Dict{String,Any})
        summary=JSON3.read(read(joinpath(output,"summary.json"),String),Dict{String,Any})
        repo=normpath(joinpath(ROOT,"..",".."))
        for (file,hash) in env["source_sha256"]
            require(bytes2hex(sha256(read(joinpath(repo,file))))==hash,"stale source: $file")
        end
        for file in ("Project","Manifest")
            require(bytes2hex(sha256(read(joinpath(repo,file*".toml"))))==env[file*"_sha256"],"$file mismatch")
        end
        require(bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))))==env["network_sha256"],"network mismatch")
        require(env["controllers"]==collect(MIXED_CONTROLLERS),"controller plan mismatch")
        for key in keys(HARDWARE)
            require(env["hardware"][string(key)]==getproperty(HARDWARE,key),"hardware metadata changed")
        end
        plan=mixed_thermal_cases();statuses=summary["cases"]
        require(Set(keys(statuses)) ⊆ Set(c.case_id for c in plan),"unplanned case")
        require(summary["planned"]==length(plan) && summary["completed"]==length(statuses),"plan coverage mismatch")
        require(summary["accepted"]==Base.count(c->c["accepted"],values(statuses)),"accepted count mismatch")
        require(summary["all_planned_accepted"]==(summary["accepted"]==length(plan)),"completion flag mismatch")
        require(length(env["planned_cases"])==length(plan),"saved plan mismatch")
        for case in plan
            saved=only(c for c in env["planned_cases"] if c["case_id"]==case.case_id)
            for key in keys(case)
                expected=getproperty(case,key);actual=saved[string(key)]
                require(expected isa Tuple ? Tuple(actual)==expected : actual==(expected isa Symbol ? String(expected) : expected),"predeclared plan changed")
            end
        end
        base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
        phasors(v)=ComplexF64[complex(z[1],z[2]) for z in v]
        evidence=Dict{String,Any}()
        for case in plan
            id=case.case_id
            if !haskey(statuses,id) || !statuses[id]["accepted"]
                evidence[id]=Dict("accepted"=>false);continue
            end
            net=voltage_aware_network(base,case);mixed=mixed_thermal_fleet(net,case);spec=mixed.spec
            point=JSON3.read(read(joinpath(output,"points",id*".json"),String),Dict{String,Any})
            active=Set(device for (device,mode) in mixed.modes if mode===:energized)
            require(Set(keys(point["devices"]))==active,"energized coverage mismatch")
            currents=Dict{String,NTuple{3,ComplexF64}}();rows=Dict{String,Any}()
            voltages=Float64[]
            for (device,controlled) in spec.devices
                inv=controlled.device;request=spec.requests[device];ts=point["bus"][inv.bus]
                sensed=[phasor(ts[p])-phasor(ts["n"]) for p in inv.phase_terminals]
                append!(voltages,abs.(sensed))
                if !(device in active)
                    currents[device]=(0im,0im,0im);continue
                end
                record=point["devices"][device]
                ic=phasors(record["converter_current"]);ig=phasors(record["grid_current"]);u=phasors(record["converter_voltage"])
                v=[phasor(ts[p])-(inv.neutral===nothing ? 0im : phasor(ts[inv.neutral])) for p in inv.phase_terminals]
                report=inverter_physical_feasibility(inv;converter_voltage=u,converter_current=ic,
                    grid_voltage=v,grid_current=ig,dc_min=0.,dc_max=request.p_available,
                    pwm_thermal_headroom=inv.pwm_strategy===:NONE ? nothing : PowerOptLab._PWM_ALLOCATION_HEADROOM)
                require(report.checked && report.passed,"physical rejection: $id/$device")
                physical+=length(report.checks)
                power=sum(v.*conj.(ig));pdc=only(c.value for c in report.checks if c.name=="dc_source_power_max")
                if controlled.controller isa PerPhaseVoltVarWatt
                    q,f=PowerOptLab._phase_voltage_requests(controlled.controller,InverterControlMeasurement(sensed),request)
                    checks=PowerOptLab._phase_service_checks(inv,request,v,ig,pdc,q,f,record["phase_reactive_scale"],:rated,request.p_available)
                    require(all(c.passed for c in checks),"phase policy rejected: $id/$device")
                else
                    c=controlled.controller
                    _,u1,_=PowerOptLab._phasor_sequences(Tuple(sensed))
                    f,q,_,_=PowerOptLab._positive_exact(c.positive,u1,abs.(sensed),request)
                    ray=PowerOptLab._inverter_limit_check("Q_ray",:VA,imag(power),min(0.,q),max(0.,q))
                    watt=PowerOptLab._inverter_limit_check("VW",:VA,real(power),-Inf,f*request.p_rated)
                    require(ray.passed && watt.passed,"sequence positive-service policy rejected")
                    require(abs(PowerOptLab._phasor_sequences(Tuple(ic))[3])<=1.01e-5,"unrequested negative sequence")
                end
                require(haskey(point["stages"],device) && length(point["stages"][device])>=2 &&
                    all(s["solve"]["publishable"] for s in point["stages"][device]),"missing/failed allocation stage")
                value(name)=begin;index=findfirst(c->c.name==name,report.checks);index===nothing ? nothing : report.checks[index].value;end
                rows[device]=Dict("controller"=>mixed.labels[device],"p_poc_W"=>real(power),"q_poc_var"=>imag(power),
                    "p_dc_W"=>pdc,"p_available_W"=>request.p_available,"vmax_sensed_V"=>maximum(abs,sensed),
                    "capacitor_thermal_current_A"=>value("capacitor_current"),
                    "capacitor_carrier_audit_current_A"=>value("allocation_pwm_actual_capacitor_current"))
                currents[device]=Tuple(ig)
            end
            proxy=(solve=(publishable=true,),spec=spec,grid_currents=currents,operating_modes=mixed.modes,network=Dict("bus"=>point["bus"]))
            replay=mixed_thermal_replay(net,proxy)
            require(replay.accepted,"mixed fixed-PQ replay failed")
            metrics=Dict("p_poc_W"=>sum(r["p_poc_W"] for r in values(rows)),"q_poc_var"=>sum(r["q_poc_var"] for r in values(rows)),
                "vmax_sensed_V"=>maximum(voltages),"maximum_capacitor_model_A"=>maximum(r["capacitor_thermal_current_A"] for r in values(rows)))
            for (key,value) in metrics
                require(isfinite(value) && abs(value-statuses[id][key])<=1e-5,"case metric mismatch: $id/$key")
            end
            evidence[id]=merge(metrics,Dict("accepted"=>true,"devices"=>rows,"fixed_PQ_replay_difference_V"=>replay.max_difference_V))
        end
        accepted=Base.count(r->r["accepted"],values(evidence))
        report=Dict("all_planned_accepted"=>accepted==length(plan),"accepted"=>accepted,"planned"=>length(plan),
            "conditions"=>count[],"physical_checks"=>physical,"cases"=>evidence,
            "scope"=>"saved-phasor SI energy/PWM/plant/policy reconstruction; same-engine fixed-PQ replay; no global optimality or switching-dynamics certificate")
        writejson(target,report);report
    catch err
        writejson(target,Dict("all_planned_accepted"=>false,"reason"=>"AUDIT_FAILED","error"=>sprint(showerror,err)))
        rethrow()
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @testset "Mixed Australian PWM/energy saved points" begin
        @test audit_mixed_thermal(isempty(ARGS) ? joinpath(ROOT,"results","thermal_energy") : ARGS[1])["all_planned_accepted"]
    end
end
