include("run_engineering_benchmark.jl")

function audit_engineering_benchmark(output;python="python3")
    target=joinpath(output,"audit.json")
    writejson(target,Dict("all_planned_accepted"=>false,"reason"=>"AUDIT_INCOMPLETE"))
    checked=Ref(0);physical=0;network_conditions=0
    function require(condition,message)
        checked[]+=1;condition || throw(ArgumentError(message))
    end
    read_record(file)=JSON3.read(read(file,String),Dict{String,Any})
    raw_accepted(s)=s["termination_status"] in ("LOCALLY_SOLVED","OPTIMAL") &&
        s["primal_status"]=="FEASIBLE_POINT" && all(s[k] for k in ("has_primal","feasible","optimal","publishable"))
    try
        env=read_record(joinpath(output,"environment.json"));summary=read_record(joinpath(output,"summary.json"))
        repo=normpath(joinpath(ROOT,"..",".."))
        for (file,digest) in env["source_sha256"]
            require(bytes2hex(sha256(read(joinpath(repo,file))))==digest,"stale source: $file")
        end
        for name in ("Project","Manifest")
            require(bytes2hex(sha256(read(joinpath(repo,name*".toml"))))==env[name*"_sha256"],"$name mismatch")
        end
        require(bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))))==env["network_sha256"],"network mismatch")
        require(!env["sizing_changed"] && env["s_base_VA"]==VA_BASE,"sizing/base metadata changed")
        require(env["controllers"]==collect(MIXED_CONTROLLERS),"controller coverage changed")
        require(env["planned_modes"]==["local","coordinated"],"mode plan changed")
        require(env["site_availability_fractions"]==collect(MIXED_PV_FRACTIONS),"site availability metadata changed")
        require(env["pwm_thermal_allowance"]==PowerOptLab._PWM_ALLOCATION_HEADROOM,"thermal allowance changed")
        canonical(x)=JSON3.read(JSON3.write(_engineering_json_safe(x)),Dict{String,Any})
        require(env["control_settings"]==canonical(TERMINAL_CONTROL_SETTINGS),"controller settings changed")
        require(env["local_solver_options"]==canonical(ENGINEERING_OPTIONS),"local solver profile changed")
        require(env["coordinated_solver_options"]==canonical(ENGINEERING_COORDINATED_OPTIONS),"coordinated solver profile changed")
        require(env["planning_voltage_max_V"]==253.,"planning ceiling changed")
        for field in keys(HARDWARE)
            require(env["hardware"][string(field)]==getproperty(HARDWARE,field),"hardware changed: $field")
        end
        for field in keys(ENGINEERING_OPERATING_LIMITS)
            require(env["operating_limits"][string(field)]==getproperty(ENGINEERING_OPERATING_LIMITS,field),"engineering ceiling changed")
        end
        plan=engineering_conditions()
        require(length(env["planned_cases"])==length(plan),"case plan changed")
        for condition in plan
            saved=only(c for c in env["planned_cases"] if c["case_id"]==condition.case_id)
            for field in keys(condition)
                expected=getproperty(condition,field);actual=saved[string(field)]
                require(expected isa Tuple ? Tuple(actual)==expected : actual==(expected isa Symbol ? String(expected) : expected),"case condition changed")
            end
        end
        planned=Set(c.case_id*"_"*mode for c in plan for mode in ("local","coordinated"))
        require(Set(keys(summary["cases"]))==planned,"missing or extra result")
        require(summary["planned"]==summary["completed"]==length(planned),"incomplete plan")
        require(summary["accepted"]==count(r->get(r,"accepted",false),values(summary["cases"])),"acceptance count mismatch")
        evidence=Dict{String,Any}();hashes=Dict{String,String}()
        base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
        for condition in plan
            net=voltage_aware_network(base,condition);mixed=engineering_fleet(net,condition)
            baseline_hardware=Dict(id=>engineering_hardware_sha256(c.device) for (id,c) in mixed.spec.devices)
            for mode in ("local","coordinated")
                key=condition.case_id*"_"*mode;status=summary["cases"][key]
                if !get(status,"accepted",false)
                    evidence[key]=Dict("accepted"=>false,"reason"=>get(status,"reason","MODEL_REJECTED"))
                    continue
                end
                path=joinpath(output,"points",key*".json");point=read_record(path)
                details=read_record(joinpath(output,key*"_status.json"))
                digest=bytes2hex(sha256(read(path)));hashes[relpath(path,output)]=digest
                require(digest==status["point_sha256"]==details["point_sha256"],"point digest mismatch")
                require(details["accepted"] && raw_accepted(point["candidate_status"]),"rejected model promoted")
                require(raw_accepted(point["solver_status"]),"raw unsuccessful solve promoted")
                require(point["planning_voltage_max_V"]==env["planning_voltage_max_V"],"point planning ceiling changed")
                for field in keys(condition)
                    expected=getproperty(condition,field);actual=point["case"][string(field)]
                    require(expected isa Tuple ? Tuple(actual)==expected : actual==(expected isa Symbol ? String(expected) : expected),"point condition changed")
                end
                require(point["coordinated"]==(mode=="coordinated"),"wrong control mode")
                require(Set(keys(point["devices"]))==Set(keys(mixed.spec.devices)),"missing or duplicated site")
                require(point["model_hardware_sha256"]==baseline_hardware,"physical models differ between PF and OPF")
                if mode=="local"
                    require(Set(keys(point["allocation_stages"]))==Set(keys(mixed.spec.devices)),"local stage coverage")
                    for stages in values(point["allocation_stages"])
                        require(length(stages)>=2 && all(raw_accepted(s["solve"]) for s in stages),"failed allocation stage")
                    end
                else
                    require(raw_accepted(details["raw"]) && details["replay"]["accepted"],"raw/replay rejection")
                    require(details["temporary_neutral_constraints_remaining"]==0,"initialization rows survived")
                    attempts=details["attempts"]
                    require(last(attempts)["stage"]=="full_capability" && iszero(last(attempts)["weight"]),"wrong final objective/capabilities")
                    require(raw_accepted(last(attempts)["raw"]),"failed final optimization")
                    require(last(attempts)["phasors"]==point["devices"],"initialization point published as final result")
                end
                require(engineering_network_audit(path,joinpath(output,"network_audits",key*".json");python),"independent network rejection")
                network=read_record(joinpath(output,"network_audits",key*".json"))
                network_conditions+=length(network["checks"])
                metrics=engineering_record_metrics(point,mixed;mode)
                require(metrics.accepted,"SI plant/current-capability/local-policy rejection")
                physical+=sum(length(r.checks) for r in values(metrics.reports))
                require(abs(sum(r.p_poc_W for r in metrics.rows)-network["pv_active_power_W"])<=1e-5,"network/device P mismatch")
                currents=Dict(id=>Tuple(ComplexF64[complex(z...) for z in d["grid_current"]]) for (id,d) in point["devices"])
                proxy=(solve=(publishable=true,),spec=mixed.spec,grid_currents=currents,
                    operating_modes=mixed.modes,network=Dict("bus"=>point["bus"]))
                replay=mixed_thermal_replay(net,proxy)
                require(replay.accepted,"saved-point fixed-PQ replay failed")
                evidence[key]=Dict("accepted"=>true,"devices"=>metrics.rows,
                    "p_poc_W"=>sum(r.p_poc_W for r in metrics.rows),"q_poc_var"=>sum(r.q_poc_var for r in metrics.rows),
                    "p_dc_W"=>sum(r.p_dc_W for r in metrics.rows),"curtailment_W"=>sum(r.curtailment_W for r in metrics.rows),
                    "reactive_service_error_var"=>sum(r.reactive_service_error_var for r in metrics.rows),
                    "maximum_sensed_V"=>maximum(r.vmax_sensed_V for r in metrics.rows),
                    "fixed_PQ_replay_difference_V"=>replay.max_difference_V)
            end
            lk=condition.case_id*"_local";ck=condition.case_id*"_coordinated"
            if get(evidence[lk],"accepted",false)&&get(evidence[ck],"accepted",false)
                require(evidence[ck]["p_poc_W"]>=evidence[lk]["p_poc_W"]-.011,"coordinated dispatch below its local-export anchor")
            end
        end
        accepted=count(r->r["accepted"],values(evidence))
        require(summary["all_planned_accepted"]==(accepted==length(planned)),"completion flag mismatch")
        report=Dict("all_planned_accepted"=>accepted==length(planned),"planned"=>length(planned),"accepted"=>accepted,
            "conditions"=>checked[],"physical_checks"=>physical,"network_conditions"=>network_conditions,
            "cases"=>evidence,"record_sha256"=>hashes,
            "environment_sha256"=>bytes2hex(sha256(read(joinpath(output,"environment.json")))),
            "audit_code_sha256"=>bytes2hex(sha256(read(@__FILE__))),
            "scope"=>"saved-phasor SI hardware/current-capability/local-policy reconstruction, independent Python SI feeder oracle and same-engine fixed-PQ replay; free-Q OPF trades away VV/VW obligations; no global optimality or transient guarantee")
        write_engineering_json(target,report);report
    catch err
        writejson(target,Dict("all_planned_accepted"=>false,"reason"=>"AUDIT_FAILED","error"=>sprint(showerror,err)))
        rethrow()
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS) in (1,2) || error("usage: audit_engineering_benchmark.jl output [python]")
    result=audit_engineering_benchmark(ARGS[1];python=length(ARGS)==2 ? ARGS[2] : "python3")
    result["all_planned_accepted"] || exit(1)
end
