include("engineering_records.jl")

function engineering_network_audit(point_path,output;python="python3")
    process=Base.run(ignorestatus(Cmd([python,joinpath(@__DIR__,"network_si_audit.py"),
        joinpath(ROOT,"networks","LV3_55bus.bmopf.json"),point_path,output])))
    report=isfile(output) ? JSON3.read(read(output,String),Dict{String,Any}) : Dict("accepted"=>false)
    success(process) && get(report,"accepted",false)
end

function run_engineering_benchmark(output;python="python3")
    isdir(output) && !isempty(readdir(output)) && error("output must be empty; retain previous runs")
    repo=normpath(joinpath(ROOT,"..",".."));planned=engineering_conditions()
    source=study_source_sha256()
    for file in ("four_leg.jl","thermal_energy.jl","switching_engineering.jl","voltage_aware.jl",
            "coordinated_engineering.jl","engineering_records.jl","run_engineering_benchmark.jl",
            "network_si_audit.py")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    revision=strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String))
    status=read(Cmd(["git","-C",repo,"status","--porcelain"]),String)
    mkpath(joinpath(output,"points"));mkpath(joinpath(output,"network_audits"))
    write_engineering_json(joinpath(output,"environment.json"),Dict("source_sha256"=>source,
        "planned_cases"=>planned,"planned_modes"=>["local","coordinated"],"hardware"=>HARDWARE,
        "operating_limits"=>ENGINEERING_OPERATING_LIMITS,"controllers"=>MIXED_CONTROLLERS,
        "site_availability_fractions"=>MIXED_PV_FRACTIONS,"control_settings"=>TERMINAL_CONTROL_SETTINGS,
        "local_solver_options"=>ENGINEERING_OPTIONS,"coordinated_solver_options"=>ENGINEERING_COORDINATED_OPTIONS,
        "s_base_VA"=>VA_BASE,"planning_voltage_max_V"=>253.,"sizing_changed"=>false,
        "pwm_thermal_allowance"=>PowerOptLab._PWM_ALLOCATION_HEADROOM,
        "git_revision"=>revision,"git_status"=>status,"julia_version"=>string(VERSION),
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Project.toml")))),
        "Manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest.toml")))),
        "network_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json")))),
        "objective"=>"maximum total POC active injection; Q free, original VV/VW commands released; no global optimality bound",
        "current_capabilities"=>"balanced converter current at three-leg sites; free four-leg current with equal per-phase P; all temporary initialization rows removed",
        "scope"=>"matched physical/network limits and hardware; reactive service obligations deliberately differ and their tradeoff must be reported; fixed-point PF and nonconvex OPF, not dynamic stability or regional conformity"))
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    rows=NamedTuple[];cases=NamedTuple[];statuses=Dict{String,Any}()
    function checkpoint()
        write_csv(joinpath(output,"devices.csv"),rows)
        write_csv(joinpath(output,"cases.csv"),cases)
        write_engineering_json(joinpath(output,"summary.json"),Dict("planned"=>2length(planned),
            "completed"=>length(statuses),"accepted"=>count(c->c["accepted"],values(statuses)),
            "all_planned_accepted"=>length(statuses)==2length(planned)&&all(c["accepted"] for c in values(statuses)),
            "cases"=>statuses))
    end
    function save_point!(condition,mixed,point,details,model_accepted,mode)
        key=condition.case_id*"_"*mode
        point["model_hardware_sha256"]=Dict(id=>engineering_hardware_sha256(c.device) for (id,c) in mixed.spec.devices)
        path=joinpath(output,"points",key*".json")
        write_engineering_json(path,point)
        network_ok=engineering_network_audit(path,joinpath(output,"network_audits",key*".json");python)
        metrics=try
            engineering_record_metrics(point,mixed;mode)
        catch err
            details["reconstruction_error"]=sprint(showerror,err)
            nothing
        end
        accepted=model_accepted && network_ok && metrics!==nothing && metrics.accepted
        details["accepted"]=accepted;details["network_accepted"]=network_ok
        details["point_sha256"]=bytes2hex(sha256(read(path)))
        if metrics!==nothing
            details["physical_reports"]=metrics.reports;details["policy_or_capability_checks"]=metrics.policies
            append!(rows,[merge(row,(accepted=accepted,)) for row in metrics.rows])
            push!(cases,(case_id=condition.case_id,mode=mode,accepted=accepted,
                p_poc_W=sum(r.p_poc_W for r in metrics.rows),q_poc_var=sum(r.q_poc_var for r in metrics.rows),
                p_available_W=sum(r.p_available_W for r in metrics.rows),p_dc_W=sum(r.p_dc_W for r in metrics.rows),
                curtailment_W=sum(r.curtailment_W for r in metrics.rows),
                reactive_service_error_var=sum(r.reactive_service_error_var for r in metrics.rows),
                converter_loss_W=sum(r.converter_loss_W for r in metrics.rows),filter_loss_W=sum(r.filter_loss_W for r in metrics.rows),
                maximum_sensed_V=maximum(r.vmax_sensed_V for r in metrics.rows)))
        end
        write_engineering_json(joinpath(output,key*"_status.json"),details)
        statuses[key]=Dict("accepted"=>accepted,"model_accepted"=>model_accepted,"network_accepted"=>network_ok,
            "raw_status"=>details["solve"],"point_sha256"=>details["point_sha256"])
        checkpoint();accepted
    end
    for condition in planned
        net=voltage_aware_network(base,condition);mixed=engineering_fleet(net,condition)
        for mode in ("local","coordinated")
            try
                if mode=="local"
                    local_result=solve_voltage_aware_inverter_fleet(net,mixed.spec;
                        operating_modes=mixed.modes,s_base=VA_BASE,solver_options=ENGINEERING_OPTIONS)
                    point=engineering_point(local_result,condition)
                    point["allocation_stages"]=engineering_local_stages(local_result)
                    details=Dict{String,Any}("solve"=>local_result.solve,"reason"=>local_result.reason,
                        "history"=>local_result.history,"physical_reports"=>local_result.physical_reports)
                    save_point!(condition,mixed,point,details,local_result.solve.publishable,mode)
                else
                    key=condition.case_id*"_local"
                    warm=get(get(statuses,key,Dict()),"accepted",false)
                    local_point=warm ? JSON3.read(read(joinpath(output,"points",key*".json"),String),Dict{String,Any}) : nothing
                    minimum_p=warm ? only(r.p_poc_W for r in cases if r.case_id==condition.case_id&&r.mode=="local") : nothing
                    result=solve_engineering_coordinated(net,mixed;local_point,minimum_p_W=minimum_p,s_base=VA_BASE)
                    point=engineering_point(result,condition;coordinated=true)
                    details=Dict{String,Any}("solve"=>result.solve,"raw"=>result.raw,"raw_solver_status"=>result.raw_solver_status,
                        "replay"=>result.replay,"attempts"=>result.attempts,
                        "temporary_neutral_constraints_remaining"=>result.temporary_neutral_constraints_remaining)
                    save_point!(condition,mixed,point,details,result.model_candidate_accepted,mode)
                end
            catch err
                key=condition.case_id*"_"*mode
                statuses[key]=Dict("accepted"=>false,"reason"=>"STUDY_RUN_ERROR","error"=>sprint(showerror,err))
                checkpoint()
                println(stderr,"Retained failed case ",key,": ",sprint(showerror,err))
            end
        end
    end
    JSON3.read(read(joinpath(output,"summary.json"),String),Dict{String,Any})
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS) in (1,2) || error("usage: run_engineering_benchmark.jl output [python]")
    report=run_engineering_benchmark(ARGS[1];python=length(ARGS)==2 ? ARGS[2] : "python3")
    report["all_planned_accepted"] || exit(1)
end
