# Reproducible quasi-static experiments for the staged physical allocator.
# Usage: julia --project=. scripts/australian_pv/voltage_aware.jl [snapshots|chronology] [output]
include("run.jl")

const VA_OPTIONS=("tol"=>1e-9,"max_iter"=>1000,"bound_relax_factor"=>0.,
    "mu_strategy"=>"adaptive","nlp_scaling_method"=>"none")
const VA_BASE=1e5
const VA_DAY_GRID=TimeGrid([6.,3.,3.,3.,3.,6.])
const VA_DAY_PV=[0.,.15,.65,1.,.45,0.]
const VA_DAY_LOAD=[.6,.8,.9,1.,1.1,1.2]

function voltage_aware_cases(suite::Symbol)
    cases=NamedTuple[]
    function add(id,group,mode;pv=1.,source=1.03,vuf=.01,angle=0.,
            shares=(.60,.25,.15),operation=:energized,step=0,start=0.,duration=0.,load_scale=1.)
        push!(cases,(case_id=id,group=group,controller=mode,pv_fraction=pv,source_pu=source,
            source_vuf=vuf,source_negative_angle_rad=angle,phase_shares=shares,
            operating_mode=operation,step=step,start_h=start,duration_h=duration,load_scale=load_scale))
    end
    if suite===:snapshots
        for (si,source) in enumerate((1.,1.03)),(pvi,pv) in enumerate((0.,.5,1.)),mode in TERMINAL_MODES
            add("matched_s$(si)_p$(pvi)_$mode","matched",mode;source,pv)
        end
        # Magnitude, relative negative-sequence angle and load phase allocation
        # vary separately. No hardware or controller gain is changed.
        for vuf in (0.,.02),mode in TERMINAL_MODES
            add("vuf$(round(Int,100vuf))_$mode","source_unbalance",mode;vuf)
        end
        for (i,angle) in enumerate((pi/2,pi,3pi/2)),mode in TERMINAL_MODES
            add("angle$(i)_$mode","source_angle",mode;angle)
        end
        for (i,shares) in enumerate(((1/3,1/3,1/3),(.15,.60,.25))),mode in TERMINAL_MODES
            add("loads$(i)_$mode","load_unbalance",mode;shares)
        end
        for (i,source) in enumerate((1.,1.03))
            add("disconnected_s$i","disconnected","mean_vv_vw";pv=0.,source,operation=:disconnected)
        end
    elseif suite===:chronology
        for policy in (:energized_night,:disconnected_night),mode in TERMINAL_MODES
            start=0.
            for step in eachindex(VA_DAY_PV)
                pv=VA_DAY_PV[step];duration=VA_DAY_GRID[step]
                operation=policy===:disconnected_night && iszero(pv) ? :disconnected : :energized
                add("$(policy)_$(mode)_t$step",String(policy),mode;pv,operation,step,start,duration,
                    load_scale=VA_DAY_LOAD[step])
                start+=duration
            end
        end
    else
        throw(ArgumentError("suite must be :snapshots or :chronology"))
    end
    length(unique(c.case_id for c in cases))==length(cases) || error("duplicate case ids")
    cases
end

function voltage_aware_network(base,case)
    net=study_network(base,case.pv_fraction,case.source_pu;
        source_vuf=case.source_vuf,phase_shares=case.phase_shares)
    # Rotate U2 relative to U1 without altering |U1| or adding U0.
    a=cis(2pi/3)
    for source in values(net["voltage_source"])
        v=source["v_magnitude"].*cis.(source["v_angle"])
        u1=(v[1]+a*v[2]+a^2*v[3])/3
        u2=case.source_vuf*u1*cis(case.source_negative_angle_rad)
        v=u1.*[1,a^2,a].+u2.*[1,a,a^2]
        source["v_magnitude"]=abs.(v);source["v_angle"]=angle.(v)
    end
    for load in values(net["load"])
        load["p_nom"].*=case.load_scale;load["q_nom"].*=case.load_scale
    end
    net
end

# A different injection formulation in the same network engine. This is a
# fixed-PQ cross-check, not an independent circuit or external simulator.
function voltage_aware_pq_replay(net,result)
    result.solve.publishable || return (accepted=false,termination="NOT_RUN",max_difference_V=NaN)
    fixed=deepcopy(net)
    for id in keys(result.spec.devices)
        delete!(fixed["ibr"],id)
        result.operating_modes[id]===:disconnected && continue
        inv=result.spec.devices[id].device
        # This campaign's fixed three-wire L filter has no POC shunt.
        iszero(inv.b_filter_shunt) && inv.neutral===nothing || error("replay scope requires three-wire/no-shunt study plant")
        current=result.grid_currents[id]
        abs(sum(current))<1e-5 || error("nonzero three-wire current sum")
        ts=result.network["bus"][inv.bus];vn=phasor(ts["n"])
        for (i,phase) in enumerate(inv.phase_terminals)
            s=(phasor(ts[phase])-vn)*conj(current[i])
            fixed["load"]["replay_$(id)_$phase"]=Dict("bus"=>inv.bus,
                "terminal_map"=>[phase,"n"],"configuration"=>"SINGLE_PHASE",
                "model"=>"constant_power","v_nom"=>[HARDWARE.v_nom],
                "p_nom"=>[-real(s)],"q_nom"=>[-imag(s)])
        end
    end
    pf=solve_pf(fixed;per_unit=true,s_base=VA_BASE,solver_options=VA_OPTIONS)
    status=String(pf["termination_status"])
    if !(status in ("LOCALLY_SOLVED","OPTIMAL") && get(pf,"feasible",false))
        return (accepted=false,termination=status,max_difference_V=NaN)
    end
    delta=maximum(abs(phasor(t)-phasor(pf["bus"][bus][terminal]))
        for (bus,ts) in result.network["bus"] for (terminal,t) in ts)
    (accepted=isfinite(delta) && delta<=.001,termination=status,max_difference_V=delta)
end

function voltage_aware_energy(cases,devices,planned)
    output=NamedTuple[]
    for group in ("energized_night","disconnected_night"),mode in TERMINAL_MODES
        plan=[c for c in planned if c.group==group && c.controller==mode]
        isempty(plan) && continue
        actual=[c for c in cases if c.group==group && c.controller==mode]
        accepted=length(actual)==length(plan) && all(c.accepted for c in actual)
        ids=Set(c.case_id for c in plan)
        selected=[d for d in devices if d.case_id in ids]
        if !accepted
            push!(output,(group=group,controller=mode,accepted=false,periods_completed=length(actual),
                available_kWh=nothing,export_kWh=nothing,import_kWh=nothing,net_poc_kWh=nothing,
                unused_dc_kWh=nothing,converter_loss_kWh=nothing,filter_loss_kWh=nothing,
                reactive_support_kvarh=nothing,balance_residual_Wh=nothing))
            continue
        end
        hours=Dict(c.case_id=>c.duration_h for c in plan)
        energy(f)=sum(hours[d.case_id]*f(d) for d in selected)/1000
        available=energy(d->d.p_available_W)
        exported=energy(d->max(d.p_poc_W,0.));imported=energy(d->max(-d.p_poc_W,0.))
        # DC standby is undefined when disconnected. These intervals have zero
        # PV availability and no AC model; do not assign an invented standby loss.
        unused=energy(d->d.model_present ? d.p_available_W-d.p_dc_W : d.p_available_W)
        converter=energy(d->d.model_present ? d.converter_loss_W : 0.)
        filter=energy(d->d.model_present ? d.filter_loss_W : 0.)
        residual=1000*(available-(exported-imported)-unused-converter-filter)
        abs(residual)<=.001 || error("chronological energy balance failed: $residual Wh")
        push!(output,(group=group,controller=mode,accepted=true,periods_completed=length(actual),
            available_kWh=available,export_kWh=exported,import_kWh=imported,net_poc_kWh=exported-imported,
            unused_dc_kWh=unused,converter_loss_kWh=converter,filter_loss_kWh=filter,
            reactive_support_kvarh=energy(d->d.q_poc_var),balance_residual_Wh=residual))
    end
    output
end

function run_voltage_aware_study(;suite=:snapshots,
        output=joinpath(ROOT,"results","voltage_aware_$(suite)"))
    planned=voltage_aware_cases(suite)
    isdir(output) && !isempty(readdir(output)) && throw(ArgumentError("output must be empty: $output"))
    mkpath(joinpath(output,"network"))
    repo=normpath(joinpath(ROOT,"..",".."))
    source=study_source_sha256()
    source[relpath(@__FILE__,repo)]=bytes2hex(sha256(read(@__FILE__)))
    setup_tests=joinpath(@__DIR__,"test_voltage_aware.jl")
    source[relpath(setup_tests,repo)]=bytes2hex(sha256(read(setup_tests)))
    writejson(joinpath(output,"environment.json"),Dict(
        "suite"=>String(suite),"source_sha256"=>source,
        "git_revision"=>strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String)),
        "git_status"=>read(Cmd(["git","-C",repo,"status","--porcelain"]),String),
        "julia"=>string(VERSION),"PowerOptLab"=>string(Base.pkgversion(PowerOptLab)),
        "BMOPFTools"=>string(Base.pkgversion(BMOPFTools)),"Ipopt"=>string(Base.pkgversion(Ipopt)),
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Project.toml")))),
        "Manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest.toml")))),
        "network_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json")))),
        "planned_cases"=>planned,"hardware"=>HARDWARE,"curves"=>AUS_A,
        "control_settings"=>TERMINAL_CONTROL_SETTINGS,"solver_options"=>Dict(VA_OPTIONS),
        "s_base_VA"=>VA_BASE,"max_iterations"=>30,"relaxation"=>1.,
        "voltage_tolerance_V"=>1e-7,"current_tolerance_A"=>1e-6,
        "reactive_mode"=>"within_request","pq_replay_tolerance_V"=>.001,
        "dc_source_model"=>"non-charging PV, open at line frequency; fundamental/2omega averaged model, no PWM",
        "disconnection_model"=>"open grid-side contactor; complete AC model absent; DC standby excluded",
        "chronology"=>"synthetic piecewise-constant quasi-static intervals; no measured yield, DC state transitions or dynamic response",
        "network_audit_scope"=>"native SI equality residuals plus fixed-PQ cross-formulation replay in the same engine"))
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    cases=NamedTuple[];devices=NamedTuple[];history=NamedTuple[];checks=NamedTuple[]
    warm=Dict{Tuple{String,String},Any}()
    function checkpoint()
        for (file,rows) in (("cases.csv",cases),("devices.csv",devices),("iterations.csv",history),
                ("physical_checks.csv",checks),("energy.csv",voltage_aware_energy(cases,devices,planned)))
            write_csv(joinpath(output,file),rows)
        end
        writejson(joinpath(output,"summary.json"),Dict("planned"=>length(planned),
            "completed"=>length(cases),"accepted"=>count(c->c.accepted,cases),
            "all_planned_accepted"=>length(cases)==length(planned) && all(c.accepted for c in cases),
            "unresolved_case_ids"=>[c.case_id for c in cases if !c.accepted]))
    end
    try
        for case in planned
            net=voltage_aware_network(base,case)
            spec=fleet(net,case.controller,case.pv_fraction;settings=TERMINAL_CONTROL_SETTINGS,
                design=(kappa=0.,phi2=0.))
            key=(case.group,case.controller)
            initial=case.step>1 ? get(warm,key,nothing) : nothing
            result=solve_voltage_aware_inverter_fleet(net,spec;
                operating_modes=case.operating_mode,initial_currents=initial,
                s_base=VA_BASE,solver_options=VA_OPTIONS)
            # Warm starts affect only numerics. Never retry a failed case.
            result.solve.publishable && case.step>0 && (warm[key]=result.grid_currents)
            replay=voltage_aware_pq_replay(net,result)
            accepted=result.solve.publishable && replay.accepted
            rows=voltage_aware_inverter_rows(result)
            prefix=(case_id=case.case_id,group=case.group,controller=case.controller)
            append!(devices,[merge(prefix,row,(accepted=accepted,)) for row in rows])
            append!(history,[merge(prefix,row) for row in result.history])
            for (id,report) in sort!(collect(result.physical_reports);by=first),check in report.checks
                push!(checks,merge(prefix,(device_id=id,name=check.name,unit=check.unit,
                    value=check.value,lower=check.lower,upper=check.upper,tolerance=check.tolerance,
                    margin=check.margin,passed=check.passed)))
            end
            sumfield(name)=sum(getproperty(row,name) for row in rows)
            push!(cases,merge(case,(accepted=accepted,coupling_accepted=result.solve.publishable,
                termination=result.termination_status,reason=result.reason,iterations=result.iterations,
                warm_start_used=initial!==nothing,pq_replay_accepted=replay.accepted,
                pq_replay_termination=replay.termination,pq_replay_difference_V=replay.max_difference_V,
                network_check_count=length(result.network_equation_checks),
                network_checks_passed=all(c.passed for c in result.network_equation_checks),
                network_max_tolerance_ratio=isempty(result.network_equation_checks) ? NaN :
                    maximum(c.value/c.tolerance for c in result.network_equation_checks),
                physical_device_count=length(result.physical_reports),
                failed_physical_devices=count(r->!r.passed,values(result.physical_reports)),
                p_poc_W=sumfield(:p_poc_W),q_poc_var=sumfield(:q_poc_var),
                vmin_sensed_V=minimum(r.vmin_sensed_V for r in rows),
                vmax_sensed_V=maximum(r.vmax_sensed_V for r in rows),
                max_vuf_percent=maximum(r.vuf_percent for r in rows))))
            writejson(joinpath(output,"network",case.case_id*".json"),result.network["bus"])
            checkpoint()
            println(case.case_id," ",result.reason," accepted=",accepted,
                " iterations=",result.iterations," P=",last(cases).p_poc_W," Q=",last(cases).q_poc_var)
            flush(stdout)
        end
    finally
        checkpoint()
    end
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    suite=isempty(ARGS) ? :snapshots : Symbol(ARGS[1])
    output=length(ARGS)>=2 ? ARGS[2] : joinpath(ROOT,"results","voltage_aware_$(suite)")
    run_voltage_aware_study(;suite,output)
end
