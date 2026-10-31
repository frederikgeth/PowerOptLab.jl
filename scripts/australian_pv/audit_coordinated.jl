# Reconstruct accepted saved points without using the optimizing model's lifts.
include("coordinated.jl")
using Test

function audit_coordinated(output=joinpath(ROOT,"results","coordinated"))
    audit_path=joinpath(output,"audit.json")
    writejson(audit_path,Dict("all_planned_accepted"=>false,"reason"=>"AUDIT_INCOMPLETE"))
    assertions=Ref(0)
    function require(condition,message)
        assertions[]+=1
        condition || throw(ArgumentError(message))
    end
    try
    env=JSON3.read(read(joinpath(output,"environment.json"),String),Dict{String,Any})
    summary=JSON3.read(read(joinpath(output,"summary.json"),String),Dict{String,Any})
    repo=normpath(joinpath(ROOT,"..",".."))
    for (path,digest) in env["source_sha256"]
        require(bytes2hex(sha256(read(joinpath(repo,path))))==digest,"source fingerprint mismatch: $path")
    end
    for file in ("Project","Manifest")
        require(bytes2hex(sha256(read(joinpath(repo,file*".toml"))))==env[file*"_sha256"],"$file fingerprint mismatch")
    end
    require(bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))))==env["network_sha256"],"network fingerprint mismatch")
    plan=env["planned_cases"];statuses=summary["cases"]
    expected=coordinated_cases();ids=Set(c["case_id"] for c in plan)
    require(length(ids)==length(plan)==length(expected),"missing or duplicate planned cases")
    require(Set(keys(statuses)) ⊆ ids,"unplanned completed case")
    for case in expected
        actual=only(c for c in plan if c["case_id"]==case.case_id)
        require(actual["source_pu"]==case.source_pu && actual["pv_fraction"]==case.pv_fraction &&
            actual["current_allocation"]==String(case.current_allocation) && actual["voltage_max_V"]==case.voltage_max_V,
            "predeclared case mismatch: $(case.case_id)")
    end
    require(length(plan)==summary["planned"],"planned count mismatch")
    require(length(statuses)==summary["completed"],"completed count mismatch")
    require(count(c->c["accepted"],values(statuses))==summary["accepted"],"accepted count mismatch")
    require(summary["all_planned_accepted"]==(summary["accepted"]==length(plan)),"completion flag mismatch")
    require(Set(summary["unresolved_case_ids"])==Set(id for (id,s) in statuses if !s["accepted"]),"unresolved ids mismatch")
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    evidence=Dict{String,Any}();physical_checks=0
    phasors(values)=ComplexF64[complex(p[1],p[2]) for p in values]
    for case in plan
        id=case["case_id"]
        if !haskey(statuses,id) || !statuses[id]["accepted"]
            evidence[id]=Dict("accepted"=>false,"reason"=>haskey(statuses,id) ? statuses[id]["reason"] : "MISSING")
            continue
        end
        point=JSON3.read(read(joinpath(output,"points",id*".json"),String),Dict{String,Any})
        net=study_network(base,case["pv_fraction"],case["source_pu"])
        spec=fleet(net,"mean_vv_vw",case["pv_fraction"];settings=TERMINAL_CONTROL_SETTINGS,design=(kappa=0.,phi2=0.))
        require(Set(keys(point["devices"]))==Set(keys(spec.devices)),"device coverage mismatch: $id")
        currents=Dict{String,NTuple{3,ComplexF64}}()
        device_evidence=Dict{String,Any}();vmax=0.
        for (device,record) in point["devices"]
            inv=spec.devices[device].device
            require(!PowerOptLab._lcl_active(inv) && inv.neutral===nothing && iszero(inv.b_filter_shunt),"unsupported saved-point audit topology")
            ts=point["bus"][inv.bus];v=[phasor(ts[p]) for p in inv.phase_terminals]
            ic=phasors(record["converter_current"]);ig=phasors(record["grid_current"])
            report=inverter_physical_feasibility(inv;converter_voltage=phasors(record["converter_voltage"]),
                converter_current=ic,grid_voltage=v,grid_current=ig,dc_min=0.,dc_max=5000case["pv_fraction"])
            require(report.checked && report.passed && report.rail.converged,"physical point rejected: $id/$device")
            physical_checks+=length(report.checks)
            if case["current_allocation"]=="balanced_converter"
                a=cis(2pi/3)
                require(abs(sum(ic)/3)<1e-5,"zero-sequence allocation violated: $id/$device")
                require(abs((ic[1]+a^2*ic[2]+a*ic[3])/3)<1e-5,"negative-sequence allocation violated: $id/$device")
            end
            if case["voltage_max_V"]!==nothing
                check=PowerOptLab._inverter_limit_check("service_voltage",:V,
                    maximum(abs.(v.-phasor(ts["n"]))),0.,case["voltage_max_V"])
                require(check.passed,"customer voltage bound violated: $id/$device")
            end
            power=sum(v.*conj.(ig));converter=sum(phasors(record["converter_voltage"]).*conj.(ic))
            a=cis(2pi/3);current2=abs((ic[1]+a^2*ic[2]+a*ic[3])/3)
            vmax=max(vmax,maximum(abs.(v.-phasor(ts["n"]))))
            device_evidence[device]=Dict("p_poc_W"=>real(power),"q_poc_var"=>imag(power),
                "p_dc_W"=>only(c.value for c in report.checks if c.name=="dc_source_power_max"),
                "converter_loss_W"=>inv.p_loss_fixed+inv.a_loss*sum(abs,ic)+inv.c_loss*sum(abs2,ic),
                "filter_loss_W"=>real(converter-power),"maximum_converter_current_A"=>maximum(abs,ic),
                "converter_negative_sequence_A"=>current2,"continuous_rail_lower_V"=>report.rail.lower_V)
            currents[device]=Tuple(ig)
        end
        proxy=(solve=(publishable=true,),spec=spec,grid_currents=currents,
            operating_modes=Dict(device=>:energized for device in keys(spec.devices)),network=Dict("bus"=>point["bus"]))
        replay=voltage_aware_pq_replay(net,proxy)
        require(replay.accepted,"fixed-PQ replay rejected: $id")
        evidence[id]=Dict("accepted"=>replay.accepted,"physical_devices"=>length(point["devices"]),
            "fixed_PQ_replay_difference_V"=>replay.max_difference_V,"devices"=>device_evidence,
            "p_poc_W"=>sum(d["p_poc_W"] for d in values(device_evidence)),
            "q_poc_var"=>sum(d["q_poc_var"] for d in values(device_evidence)),"vmax_sensed_V"=>vmax)
    end
    accepted=count(c->c["accepted"],values(evidence))
    report=Dict("planned"=>length(plan),"accepted_saved_points"=>accepted,
        "all_planned_accepted"=>accepted==length(plan),"physical_checks_recomputed"=>physical_checks,
        "source_hashes_verified"=>length(env["source_sha256"]),"cases"=>evidence,
        "assertions"=>assertions[],
        "audit_code_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "scope"=>"independent SI inverter/current-allocation/voltage checks from saved phasors; fixed-PQ replay uses same network engine; no global-optimality proof")
    writejson(audit_path,report)
    report
    catch err
        writejson(audit_path,Dict("all_planned_accepted"=>false,"reason"=>"AUDIT_FAILED","error"=>sprint(showerror,err)))
        rethrow()
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @testset "Australian coordinated saved points" begin
        report=audit_coordinated(isempty(ARGS) ? joinpath(ROOT,"results","coordinated") : ARGS[1])
        @test report["all_planned_accepted"]
        println("Saved-point audit assertions: ",report["assertions"])
    end
end
