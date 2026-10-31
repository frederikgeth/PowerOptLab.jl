include("rated_capability.jl")

function engineering_network_audit(point_path,output;python="python3")
    process=Base.run(ignorestatus(Cmd([python,joinpath(@__DIR__,"network_si_audit.py"),
        joinpath(ROOT,"networks","LV3_55bus.bmopf.json"),point_path,output])))
    report=isfile(output) ? JSON3.read(read(output,String),Dict{String,Any}) : Dict("accepted"=>false)
    success(process) && get(report,"accepted",false)
end

"""Per-limit count of sites where that limit is active (one site can count for several)."""
function engineering_binding_counts(rows)
    counts=Dict{String,Int}()
    for r in rows
        labels=unique(occursin(r"_current_[123]$",n) ? replace(n,r"_[123]$"=>"") : String(n)
            for n in split(r.binding_constraints,"|";keepempty=false))
        isempty(labels) && push!(labels,"none")
        for label in labels
            counts[label]=get(counts,label,0)+1
        end
    end
    counts
end

const ENGINEERING_MODES=("local","opf_equal_service","opf_vv_released")
# Same LV bus set as the OPF planning rows: every non-source bus, phase-neutral.
function engineering_lv_voltage_max(point,net)
    sources=Set(s["bus"] for s in values(net["voltage_source"]))
    maximum(abs(phasor(ts[p])-phasor(ts["n"])) for (b,ts) in point["bus"], p in ("a","b","c")
        if !(b in sources) && haskey(ts,"n") && haskey(ts,p))
end
engineering_planning_ok(point,net)=engineering_lv_voltage_max(point,net)<=253.0+1e-6
# AS/NZS 4777.2 sustained-overvoltage setting (10-minute average). A steady
# snapshot is treated as its own 10-minute average; no time series is implied.
const ENGINEERING_SUSTAINED_OVERVOLTAGE_V=258.

function run_engineering_benchmark(output;python="python3",limits="declared")
    isdir(output) && !isempty(readdir(output)) && error("output must be empty; retain previous runs")
    fixture=engineering_limits(limits)
    repo=normpath(joinpath(ROOT,"..",".."));planned=engineering_conditions()
    source=study_source_sha256()
    for file in ("four_leg.jl","thermal_energy.jl","switching_engineering.jl","voltage_aware.jl",
            "coordinated_engineering.jl","engineering_records.jl","rated_capability.jl",
            "run_engineering_benchmark.jl","network_si_audit.py")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    revision=strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String))
    status=read(Cmd(["git","-C",repo,"status","--porcelain"]),String)
    mkpath(joinpath(output,"points"));mkpath(joinpath(output,"network_audits"))
    cp(joinpath(repo,"Manifest.toml"),joinpath(output,"dependency_manifest.toml"))
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    # R2: the fixture must deliver its nameplate before any feeder comparison.
    gate_net=voltage_aware_network(base,first(planned))
    capability=fleet_rated_capability(engineering_fleet(gate_net,first(planned);limits=fixture))
    write_csv(joinpath(output,"rated_capability.csv"),capability)
    capability_ok=all(r.accepted for r in capability)
    write_engineering_json(joinpath(output,"environment.json"),Dict("source_sha256"=>source,
        "planned_cases"=>planned,"planned_modes"=>collect(ENGINEERING_MODES),"hardware"=>HARDWARE,
        "fixture_limits"=>limits,"operating_limits"=>fixture,"rated_capability_accepted"=>capability_ok,
        "controllers"=>MIXED_CONTROLLERS,
        "site_availability_fractions"=>MIXED_PV_FRACTIONS,"control_settings"=>TERMINAL_CONTROL_SETTINGS,
        "local_solver_options"=>ENGINEERING_OPTIONS,"coordinated_solver_options"=>ENGINEERING_COORDINATED_OPTIONS,
        "s_base_VA"=>VA_BASE,"planning_voltage_max_V"=>253.,"agreement_tolerance"=>ENGINEERING_AGREEMENT_TOL,
        "pwm_thermal_allowance"=>PowerOptLab._PWM_ALLOCATION_HEADROOM,
        "pwm_thermal_model"=>String(PowerOptLab._PWM_ALLOCATION_MODEL),
        "git_revision"=>revision,"git_status"=>status,"julia_version"=>string(VERSION),
        "threads"=>Threads.nthreads(),"cpu"=>Sys.cpu_info()[1].model,
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Project.toml")))),
        "Manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest.toml")))),
        "network_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json")))),
        "objectives"=>Dict("opf_equal_service"=>"maximum total POC active injection with each site's local VV request and VW ceiling imposed at the OPF's own voltages",
            "opf_vv_released"=>"maximum total POC active injection; Q free, local VV/VW commands released"),
        "current_capabilities"=>"balanced converter current at three-leg sites; free four-leg current with equal per-phase P; all temporary initialization rows removed",
        "scope"=>"matched physical/network limits and hardware; fixed-point PF and nonconvex local OPF; no global optimality, dynamic stability or regional conformity claim"))
    capability_ok || limits=="legacy" || error("fixture fails the rated-capability gate; see rated_capability.csv")
    rows=NamedTuple[];cases=NamedTuple[];statuses=Dict{String,Any}();agreement=Dict{String,Any}()
    function checkpoint()
        write_csv(joinpath(output,"devices.csv"),rows)
        write_csv(joinpath(output,"cases.csv"),cases)
        write_engineering_json(joinpath(output,"agreement.json"),agreement)
        n=length(ENGINEERING_MODES)*length(planned)
        write_engineering_json(joinpath(output,"summary.json"),Dict("planned"=>n,
            "completed"=>length(statuses),"accepted"=>count(c->c["accepted"],values(statuses)),
            "all_planned_accepted"=>length(statuses)==n&&all(c["accepted"] for c in values(statuses)),
            "agreement_accepted"=>length(agreement)==length(planned)&&all(a["accepted"] for a in values(agreement)),
            "cases"=>statuses))
    end
    function save_point!(condition,net,mixed,point,details,model_accepted,mode,elapsed)
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
        details["accepted"]=accepted;details["network_accepted"]=network_ok;details["elapsed_s"]=elapsed
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
                maximum_sensed_V=maximum(r.vmax_sensed_V for r in metrics.rows),
                maximum_lv_voltage_V=engineering_lv_voltage_max(point,net),
                planning_253V_exceeded=!engineering_planning_ok(point,net),
                sustained_258V_exceeded=maximum(r.vmax_sensed_V for r in metrics.rows)>ENGINEERING_SUSTAINED_OVERVOLTAGE_V,
                binding_counts=join(sort!(["$(k)=$(v)" for (k,v) in engineering_binding_counts(metrics.rows)]),";"),
                elapsed_s=elapsed,iterations=get(details,"iterations",nothing),
                nlp_variables=get(details,"nlp_variables",nothing),nlp_constraints=get(details,"nlp_constraints",nothing)))
        end
        write_engineering_json(joinpath(output,key*"_status.json"),details)
        statuses[key]=Dict("accepted"=>accepted,"model_accepted"=>model_accepted,"network_accepted"=>network_ok,
            "raw_status"=>details["solve"],"point_sha256"=>details["point_sha256"])
        checkpoint();accepted
    end
    opf_details(result)=Dict{String,Any}("solve"=>result.solve,"raw"=>result.raw,"raw_solver_status"=>result.raw_solver_status,
        "replay"=>result.replay,"attempts"=>result.attempts,"nlp_variables"=>result.model_size.variables,
        "nlp_constraints"=>result.model_size.constraints,
        "temporary_neutral_constraints_remaining"=>result.temporary_neutral_constraints_remaining)
    for condition in planned
        net=voltage_aware_network(base,condition);mixed=engineering_fleet(net,condition;limits=fixture)
        local_point=nothing
        for mode in ENGINEERING_MODES
            key=condition.case_id*"_"*mode
            try
                if mode=="local"
                    elapsed=@elapsed local_result=solve_voltage_aware_inverter_fleet(net,mixed.spec;
                        operating_modes=mixed.modes,s_base=VA_BASE,solver_options=ENGINEERING_OPTIONS)
                    # Local controllers do not enforce the planning limit: an
                    # exceedance is a reported outcome, not a rejected record.
                    point=engineering_point(local_result,condition;voltage_max=nothing)
                    point["allocation_stages"]=engineering_local_stages(local_result)
                    details=Dict{String,Any}("solve"=>local_result.solve,"reason"=>local_result.reason,
                        "iterations"=>local_result.iterations,
                        "history"=>local_result.history,"physical_reports"=>local_result.physical_reports)
                    save_point!(condition,net,mixed,point,details,local_result.solve.publishable,mode,elapsed) &&
                        (local_point=point)
                elseif local_point===nothing
                    statuses[key]=Dict("accepted"=>false,"reason"=>"LOCAL_POINT_UNAVAILABLE");checkpoint()
                elseif mode=="opf_equal_service"
                    elapsed=@elapsed result=solve_engineering_equal_service(net,mixed;local_point,s_base=VA_BASE)
                    point=engineering_point(result,condition;coordinated=true,mode)
                    point["service_rounds"]=result.service_rounds
                    details=merge(opf_details(result),Dict{String,Any}("service_rounds"=>result.service_rounds,
                        "service_converged"=>result.service_converged,"iterations"=>length(result.service_rounds)))
                    save_point!(condition,net,mixed,point,details,result.model_candidate_accepted,mode,elapsed)
                    # Agreement: the same OPF without the planning row versus the local fixed point.
                    agree=solve_engineering_equal_service(net,mixed;local_point,voltage_max=nothing,s_base=VA_BASE)
                    agreement[condition.case_id]=agree.model_candidate_accepted ?
                        engineering_agreement(local_point,engineering_point(agree,condition;coordinated=true),mixed) :
                        Dict("accepted"=>false,"reason"=>"AGREEMENT_SOLVE_REJECTED")
                    checkpoint()
                else
                    # Anchor at the local export only where the local point meets the planning limit.
                    anchor=engineering_planning_ok(local_point,net) ?
                        only(r.p_poc_W for r in cases if r.case_id==condition.case_id&&r.mode=="local") : nothing
                    elapsed=@elapsed result=solve_engineering_coordinated(net,mixed;local_point,minimum_p_W=anchor,s_base=VA_BASE)
                    point=engineering_point(result,condition;coordinated=true,mode)
                    details=merge(opf_details(result),Dict{String,Any}("local_export_anchor_W"=>anchor))
                    save_point!(condition,net,mixed,point,details,result.model_candidate_accepted,mode,elapsed)
                end
            catch err
                statuses[key]=Dict("accepted"=>false,"reason"=>"STUDY_RUN_ERROR","error"=>sprint(showerror,err))
                checkpoint()
                println(stderr,"Retained failed case ",key,": ",sprint(showerror,err))
            end
        end
    end
    JSON3.read(read(joinpath(output,"summary.json"),String),Dict{String,Any})
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS) in (1,2,3) || error("usage: run_engineering_benchmark.jl output [python] [declared|legacy]")
    report=run_engineering_benchmark(ARGS[1];python=length(ARGS)>=2 ? ARGS[2] : "python3",
        limits=length(ARGS)==3 ? ARGS[3] : "declared")
    report["all_planned_accepted"] && report["agreement_accepted"] || exit(1)
end
