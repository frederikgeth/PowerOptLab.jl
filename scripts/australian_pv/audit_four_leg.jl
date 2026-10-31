# Fresh SI reconstruction from saved phasors, not optimization-model lifts.
include("four_leg.jl")
using Test

function audit_four_leg(output=joinpath(ROOT,"results","four_leg"))
    path=joinpath(output,"audit.json")
    writejson(path,Dict("all_planned_accepted"=>false,"reason"=>"AUDIT_INCOMPLETE"))
    assertions=Ref(0)
    function require(condition,message)
        assertions[]+=1
        condition || throw(ArgumentError(message))
    end
    try
        env=JSON3.read(read(joinpath(output,"environment.json"),String),Dict{String,Any})
        summary=JSON3.read(read(joinpath(output,"summary.json"),String),Dict{String,Any})
        repo=normpath(joinpath(ROOT,"..",".."))
        for (p,digest) in env["source_sha256"]
            require(bytes2hex(sha256(read(joinpath(repo,p))))==digest,"source fingerprint mismatch: $p")
        end
        for file in ("Project","Manifest")
            require(bytes2hex(sha256(read(joinpath(repo,file*".toml"))))==env[file*"_sha256"],"$file fingerprint mismatch")
        end
        require(bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))))==env["network_sha256"],"network fingerprint mismatch")
        plan=four_leg_cases();ids=Set(c.case_id for c in plan)
        require(length(env["planned_cases"])==length(plan),"plan length mismatch")
        for case in plan
            record=only(c for c in env["planned_cases"] if c["case_id"]==case.case_id)
            for key in keys(case)
                expected=getproperty(case,key);actual=record[string(key)]
                require(expected isa Tuple ? Tuple(actual)==expected : actual==(expected isa Symbol ? String(expected) : expected),
                    "plan mismatch: $(case.case_id)/$key")
            end
        end
        statuses=summary["cases"]
        require(Set(keys(statuses)) ⊆ ids,"unplanned case")
        require(summary["planned"]==length(plan) && summary["completed"]==length(statuses),"summary coverage mismatch")
        require(summary["accepted"]==count(c->c["accepted"],values(statuses)),"accepted count mismatch")
        require(summary["all_planned_accepted"]==(summary["accepted"]==length(plan)),"completion flag mismatch")
        require(Set(summary["unresolved_case_ids"])==Set(id for (id,c) in statuses if !c["accepted"]),"unresolved ids mismatch")
        base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
        evidence=Dict{String,Any}();physical_checks=0;device_count=0
        phasors(v)=ComplexF64[complex(p[1],p[2]) for p in v]
        for case in plan
            id=case.case_id
            if !haskey(statuses,id) || !statuses[id]["accepted"]
                evidence[id]=Dict("accepted"=>false,"reason"=>haskey(statuses,id) ? statuses[id]["reason"] : "MISSING")
                continue
            end
            point=JSON3.read(read(joinpath(output,"points",id*".json"),String),Dict{String,Any})
            net=voltage_aware_network(base,case);spec=four_leg_fleet(net,case.pv_fraction)
            active=case.operating_mode===:energized
            require(Set(keys(point["devices"]))==(active ? Set(keys(spec.devices)) : Set{String}()),"device coverage mismatch: $id")
            currents=Dict{String,NTuple{3,ComplexF64}}();device_evidence=Dict{String,Any}()
            sensed=Float64[];vufs=Float64[]
            for (device,c) in spec.devices
                inv=c.device;ts=point["bus"][inv.bus]
                v=[phasor(ts[p])-phasor(ts["n"]) for p in inv.phase_terminals]
                append!(sensed,abs.(v))
                a=cis(2pi/3);u1=(v[1]+a*v[2]+a^2*v[3])/3;u2=(v[1]+a^2*v[2]+a*v[3])/3
                push!(vufs,100abs(u2)/abs(u1))
                if !active
                    currents[device]=(0im,0im,0im)
                    continue
                end
                record=point["devices"][device]
                ic=phasors(record["converter_current"]);ig=phasors(record["grid_current"])
                report=inverter_physical_feasibility(inv;converter_voltage=phasors(record["converter_voltage"]),
                    converter_current=ic,grid_voltage=v,grid_current=ig,dc_min=0.,dc_max=5000case.pv_fraction)
                require(report.checked && report.passed && report.rail.converged,"physical point rejected: $id/$device")
                physical_checks+=length(report.checks);device_count+=1
                lambda=record["reactive_scale"]
                require(isfinite(lambda) && -.000001<=lambda<=1.000001,"Q coordinate invalid: $id/$device")
                q,f=PowerOptLab._phase_voltage_requests(c.controller,InverterControlMeasurement(v),spec.requests[device])
                pdc=only(ch.value for ch in report.checks if ch.name=="dc_source_power_max")
                policy=PowerOptLab._phase_service_checks(inv,spec.requests[device],v,ig,pdc,q,f,lambda,
                    c.controller.volt_watt_basis,5000case.pv_fraction)
                require(all(ch.passed for ch in policy),"phase policy rejected: $id/$device")
                physical_checks+=length(policy)
                power=sum(v.*conj.(ig));converter=sum(phasors(record["converter_voltage"]).*conj.(ic))
                device_evidence[device]=Dict("p_poc_W"=>real(power),"q_poc_var"=>imag(power),"p_dc_W"=>pdc,
                    "converter_loss_W"=>pdc-real(converter),"filter_loss_W"=>real(converter-power),
                    "neutral_current_A"=>abs(sum(ic)),"maximum_converter_current_A"=>maximum(abs,ic),
                    "vmax_sensed_V"=>maximum(abs,v),"vmin_sensed_V"=>minimum(abs,v),
                    "vuf_percent"=>100abs(u2)/abs(u1),"phase_reactive_service_scale"=>lambda)
                currents[device]=Tuple(ig)
            end
            proxy=(solve=(publishable=true,),spec=spec,grid_currents=currents,
                operating_modes=Dict(device=>case.operating_mode for device in keys(spec.devices)),network=Dict("bus"=>point["bus"]))
            replay=four_leg_pq_replay(net,proxy)
            require(replay.accepted,"fixed-PQ replay rejected: $id")
            p=sum(d["p_poc_W"] for d in values(device_evidence);init=0.)
            q=sum(d["q_poc_var"] for d in values(device_evidence);init=0.)
            metrics=Dict("p_poc_W"=>p,"q_poc_var"=>q,"vmax_sensed_V"=>maximum(sensed),
                "vmin_sensed_V"=>minimum(sensed),"max_vuf_percent"=>maximum(vufs),
                "max_neutral_current_A"=>maximum((d["neutral_current_A"] for d in values(device_evidence));init=0.))
            for (key,value) in metrics
                require(isfinite(value) && abs(value-statuses[id][key])<=1e-5,"summary metric mismatch: $id/$key")
            end
            evidence[id]=merge(metrics,Dict("accepted"=>true,"devices"=>device_evidence,
                "fixed_PQ_replay_difference_V"=>replay.max_difference_V))
        end
        accepted=count(c->c["accepted"],values(evidence))
        report=Dict("planned"=>length(plan),"accepted_saved_points"=>accepted,"all_planned_accepted"=>accepted==length(plan),
            "physical_devices"=>device_count,"physical_checks_recomputed"=>physical_checks,"assertions"=>assertions[],
            "source_hashes_verified"=>length(env["source_sha256"]),"cases"=>evidence,
            "scope"=>"saved-phasor SI hardware and phase-policy audit plus same-engine fixed-PQ replay; not independent EMT or global optimality")
        writejson(path,report);report
    catch err
        writejson(path,Dict("all_planned_accepted"=>false,"reason"=>"AUDIT_FAILED","error"=>sprint(showerror,err)))
        rethrow()
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @testset "Australian four-leg saved points" begin
        report=audit_four_leg(isempty(ARGS) ? joinpath(ROOT,"results","four_leg") : ARGS[1])
        @test report["all_planned_accepted"]
        println("Saved-point audit assertions: ",report["assertions"])
    end
end
