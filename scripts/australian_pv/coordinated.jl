# Feasible coordinated capability references; not global upper bounds.
include("voltage_aware.jl")

function solve_australian_coordinated(net,spec;current_allocation=:balanced_converter,
        voltage_max=253.,s_base=VA_BASE,solver_options=VA_OPTIONS)
    current_allocation in (:balanced_converter,:free) || throw(ArgumentError("unsupported benchmark current allocation"))
    voltage_max===nothing || isfinite(voltage_max) && voltage_max>0 || throw(ArgumentError("voltage_max must be positive or nothing"))
    isfinite(s_base) && s_base>0 || throw(ArgumentError("s_base must be positive"))
    PowerOptLab._validate_controlled_fleet(net,spec)
    devices=Dict(id=>ConstrainedInverter(c.device,InverterOperatingConstraints();
        current_allocation,dc_source=DCSourceBounds(p_min=0.)) for (id,c) in spec.devices)
    handles=Dict{String,Any}()
    builder=BMOPFTools.OpfDeviceBuilder(:PowerOptLab,function(ctx,ids)
        for id in ids
            PowerOptLab._disable_native_ibr_port!(ctx,id)
            handles[id]=stamp_device!(ctx,devices[id];request=spec.requests[id])
        end
        ctx
    end)
    build_spec=BMOPFTools.OpfBuildSpec(component_builders=Dict((:ibr,id)=>builder for id in keys(devices)))
    hook! = ctx -> begin
        m=PowerOptLab._opf_model(ctx)
        if voltage_max!==nothing
            for c in values(spec.devices),phase in c.device.phase_terminals
                bus=c.device.bus;vb=BMOPFTools.opf_ac_coordinate_bases(ctx,bus).voltage
                re,im=PowerOptLab._dv(ctx,bus,phase,"n")
                @constraint(m,(re^2+im^2)/(voltage_max/vb)^2<=1.)
            end
        end
        rating=sum(d.device.s_max for d in values(devices))
        @objective(m,Max,sum(h.plant.p_poc*h.plant.sb for h in values(handles))/rating)
    end
    ctx=PowerOptLab.build_opf_model(net;per_unit=true,s_base=Float64(s_base),add_objective=false,
        build_spec,model_hook! = hook!,optimizer=Ipopt.Optimizer)
    PowerOptLab._set_inverter_solver_options!(PowerOptLab._opf_model(ctx),solver_options)
    PowerOptLab._optimize_inverter_rails!(ctx,[(devices[id].device,h.plant) for (id,h) in handles])
    raw=solve_status(PowerOptLab._opf_model(ctx))
    extracted=Dict(id=>extract_device(devices[id],h,raw) for (id,h) in handles)
    native,checked=PowerOptLab._voltage_aware_network_checks(ctx,raw)
    outcome=PowerOptLab._solve_outcome(PowerOptLab._opf_model(ctx))
    network=PowerOptLab._network_without_replaced_ibrs(PowerOptLab._extract_result(ctx,outcome),keys(devices))
    voltage_checks=InverterLimitCheck[]
    if raw.has_primal && voltage_max!==nothing
        for (id,c) in spec.devices,phase in c.device.phase_terminals
            ts=network["bus"][c.device.bus]
            value=abs(phasor(ts[phase])-phasor(ts["n"]))
            push!(voltage_checks,PowerOptLab._inverter_limit_check("service_voltage/$id/$phase",:V,value,0.,voltage_max))
        end
    end
    plant_pass=all(d.plant.physical_feasibility.passed for d in values(extracted))
    native_pass=checked && all(c.passed for c in native)
    voltage_pass=(voltage_max===nothing || length(voltage_checks)==3length(devices)) && all(c.passed for c in voltage_checks)
    accepted=raw.publishable && plant_pass && native_pass && voltage_pass
    reason=!raw.publishable ? :SOLVER_FAILED : !plant_pass ? :INVERTER_PHYSICAL_REJECTED :
        !native_pass ? :NETWORK_EQUATION_REJECTED : !voltage_pass ? :SERVICE_VOLTAGE_REJECTED : :ACCEPTED_FEASIBLE
    currents=Dict(id=>d.grid_phase_current for (id,d) in extracted)
    result=(solve=SolveStatus(raw.termination_status,raw.primal_status,raw.has_primal,raw.feasible,raw.optimal,accepted),
        raw=raw,reason=reason,devices=extracted,current_allocation=current_allocation,
        voltage_max=voltage_max,voltage_checks=voltage_checks,network_checks=native,
        spec=spec,operating_modes=Dict(id=>:energized for id in keys(devices)),grid_currents=currents,network=network)
    replay=voltage_aware_pq_replay(net,result)
    final=accepted && replay.accepted
    merge(result,(replay=replay,accepted=final,
        reason=accepted && !replay.accepted ? :PQ_REPLAY_REJECTED : reason,
        solve=SolveStatus(raw.termination_status,raw.primal_status,raw.has_primal,raw.feasible,raw.optimal,final)))
end

function coordinated_cases()
    [(case_id="s$(si)_p$(li)_$(allocation)_$(cap===nothing ? "unrestricted" : "253V")",
        source_pu=source,pv_fraction=level,current_allocation=allocation,voltage_max_V=cap)
        for (si,source) in enumerate((1.,1.03)) for (li,level) in enumerate((0.,.5,1.))
        for allocation in (:balanced_converter,:free) for cap in (nothing,253.)]
end

function run_coordinated(;output=joinpath(ROOT,"results","coordinated"))
    planned=coordinated_cases()
    isdir(output) && !isempty(readdir(output)) && throw(ArgumentError("output must be empty"))
    mkpath(joinpath(output,"points"))
    repo=normpath(joinpath(ROOT,"..",".."));source=study_source_sha256()
    for file in ("voltage_aware.jl","test_voltage_aware.jl","coordinated.jl","test_coordinated.jl","audit_coordinated.jl")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    writejson(joinpath(output,"environment.json"),Dict("source_sha256"=>source,"planned_cases"=>planned,
        "git_revision"=>strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String)),
        "git_status"=>read(Cmd(["git","-C",repo,"status","--porcelain"]),String),
        "Manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest.toml")))),
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Project.toml")))),
        "network_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json")))),
        "hardware"=>HARDWARE,"control_settings"=>TERMINAL_CONTROL_SETTINGS,
        "solver_options"=>Dict(VA_OPTIONS),"s_base_VA"=>VA_BASE,
        "objective"=>"maximum total POC P; free Q; feasible local nonconvex benchmark, no global bound",
        "voltage_bound"=>"illustrative 253 V upper phase-neutral planning bound at all 12 customers; not regional conformity",
        "scope"=>"same three-leg hardware and DC-link PV bounds; controllers released; balanced/free sequence capabilities labeled separately",
        "source_vuf"=>.01,"source_negative_angle_rad"=>0.,"phase_shares"=>(.60,.25,.15)))
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    cases=NamedTuple[];rows=NamedTuple[];checks=NamedTuple[];network_checks=NamedTuple[]
    function checkpoint()
        for (file,data) in (("cases.csv",cases),("devices.csv",rows),("physical_checks.csv",checks),("network_checks.csv",network_checks))
            write_csv(joinpath(output,file),data)
        end
        writejson(joinpath(output,"summary.json"),Dict("planned"=>length(planned),"completed"=>length(cases),
            "accepted"=>count(c->c.accepted,cases),"unresolved_case_ids"=>[c.case_id for c in cases if !c.accepted],
            "cases"=>Dict(c.case_id=>Dict("accepted"=>c.accepted,"reason"=>String(c.reason),"termination"=>c.termination) for c in cases),
            "all_planned_accepted"=>length(cases)==length(planned) && all(c.accepted for c in cases)))
    end
    try
        for case in planned
            net=study_network(base,case.pv_fraction,case.source_pu)
            spec=fleet(net,"mean_vv_vw",case.pv_fraction;settings=TERMINAL_CONTROL_SETTINGS,design=(kappa=0.,phi2=0.))
            result=solve_australian_coordinated(net,spec;current_allocation=case.current_allocation,voltage_max=case.voltage_max_V)
            prefix=(case_id=case.case_id,)
            for (id,d) in sort!(collect(result.devices);by=first)
                plant=d.plant
                sensed=try
                    PowerOptLab._network_phase_state(result.network["bus"],spec.devices[id])
                catch err
                    err isa ArgumentError || rethrow()
                    ntuple(_->ComplexF64(NaN,NaN),3)
                end
                report=plant.physical_feasibility
                push!(rows,merge(prefix,(device_id=id,accepted=result.accepted,p_available_W=spec.requests[id].p_available,
                    p_poc_W=plant.p_poc,q_poc_var=plant.q_poc,p_dc_W=plant.p_dc,
                    converter_loss_W=plant.p_loss,capacitor_loss_W=plant.p_cap_loss,filter_loss_W=plant.p_filter_loss,
                    physical_checked=report.checked,physical_feasible=report.passed,
                    maximum_converter_current_A=maximum(plant.i_mag),converter_negative_sequence_A=plant.i_negative,
                    dc_ripple_voltage_V=plant.dv2,capacitor_thermal_current_A=plant.i_cap_thermal,
                    continuous_rail_lower_V=report.rail.lower_V,vmin_sensed_V=minimum(abs,sensed),vmax_sensed_V=maximum(abs,sensed))))
                for check in report.checks
                    push!(checks,merge(prefix,(device_id=id,name=check.name,unit=check.unit,value=check.value,
                        lower=check.lower,upper=check.upper,tolerance=check.tolerance,margin=check.margin,passed=check.passed)))
                end
            end
            for check in vcat(result.network_checks,result.voltage_checks)
                push!(network_checks,merge(prefix,(name=check.name,unit=check.unit,value=check.value,
                    lower=check.lower,upper=check.upper,tolerance=check.tolerance,margin=check.margin,passed=check.passed)))
            end
            selected=[r for r in rows if r.case_id==case.case_id]
            push!(cases,merge(case,(accepted=result.accepted,reason=result.reason,termination=result.raw.termination_status,
                raw_feasible=result.raw.feasible,raw_optimal=result.raw.optimal,pq_replay_accepted=result.replay.accepted,
                pq_replay_difference_V=result.replay.max_difference_V,
                p_poc_W=sum(r.p_poc_W for r in selected),q_poc_var=sum(r.q_poc_var for r in selected),
                vmax_sensed_V=maximum(r.vmax_sensed_V for r in selected),
                maximum_negative_sequence_A=maximum(r.converter_negative_sequence_A for r in selected))))
            finite=all(all(isfinite,d.converter_terminal.phase_voltage) &&
                all(isfinite,d.converter_terminal.phase_current) && all(isfinite,d.grid_phase_current)
                for d in values(result.devices))
            points=finite ? Dict("bus"=>result.network["bus"],"devices"=>Dict(id=>Dict(
                "converter_voltage"=>[collect(reim(v)) for v in d.converter_terminal.phase_voltage],
                "converter_current"=>[collect(reim(v)) for v in d.converter_terminal.phase_current],
                "grid_current"=>[collect(reim(v)) for v in d.grid_phase_current]) for (id,d) in result.devices)) :
                Dict("candidate_finite"=>false,"termination"=>result.raw.termination_status)
            writejson(joinpath(output,"points",case.case_id*".json"),points)
            checkpoint()
            println(case.case_id," ",result.reason," accepted=",result.accepted," P=",last(cases).p_poc_W," Q=",last(cases).q_poc_var)
            flush(stdout)
        end
    finally
        checkpoint()
    end
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_coordinated(;output=isempty(ARGS) ? joinpath(ROOT,"results","coordinated") : ARGS[1])
end
