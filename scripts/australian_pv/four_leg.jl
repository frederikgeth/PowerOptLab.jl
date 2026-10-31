# Explicit four-wire counterfactual, separate from the three-wire campaign.
include("voltage_aware.jl")

function four_leg_fleet(net,availability)
    devices=Dict{String,ControlledDevice}();requests=Dict{String,InverterControlRequest}()
    vv=PiecewiseLinearLaw(AUS_A.vv_v,AUS_A.vv_q;smoothing_epsilon=.05)
    vw=PiecewiseLinearLaw(AUS_A.vw_v,AUS_A.vw_p;smoothing_epsilon=.05)
    for (id,record) in net["ibr"]
        # A declared topology change, not a retrofit or a sizing optimization.
        record["topology"]="FOUR_LEG";record["terminal_map"]=["a","b","c","n"]
        hw=HARDWARE
        inv=AdvancedInverter(id=id,bus=record["bus"],topology=:FOUR_LEG,
            phase_terminals=["a","b","c"],neutral="n",s_max=hw.s_max,
            power_limit_location=:poc,i_max=hw.i_max,In_max=hw.i_max,
            v_dc=hw.v_dc,c_dc=hw.c_dc,dv2_max=hw.dv2_max,i_cap_max=hw.i_cap_max,
            m_max=hw.m_max,r_filter=filter_r(hw),x_filter=filter_x(hw),
            r_filter_neutral=filter_r(hw),x_filter_neutral=filter_x(hw),
            p_loss_fixed=hw.p_loss_fixed,a_loss=hw.a_loss,c_loss=hw.c_loss)
        devices[id]=ControlledDevice(inv,PerPhaseVoltVarWatt(volt_var=vv,volt_watt=vw,
            volt_watt_basis=AUS_A.vw_basis,current_target=:converter))
        requests[id]=InverterControlRequest(p_available=5000availability,p_rated=hw.s_max,q_scale=hw.s_max)
    end
    ControlledInverterFleetSpec(devices,requests)
end

function four_leg_cases()
    cases=NamedTuple[]
    function add(id;pv=1.,source=1.03,vuf=.01,shares=(.60,.25,.15),operation=:energized)
        push!(cases,(case_id=id,pv_fraction=pv,source_pu=source,source_vuf=vuf,
            source_negative_angle_rad=0.,phase_shares=shares,load_scale=1.,operating_mode=operation))
    end
    for (si,source) in enumerate((1.,1.03)),(pi,pv) in enumerate((0.,.5,1.))
        add("matched_s$(si)_p$(pi)";source,pv)
    end
    for vuf in (0.,.02);add("vuf$(round(Int,100vuf))";vuf);end
    add("balanced_loads";shares=(1/3,1/3,1/3))
    for (i,source) in enumerate((1.,1.03))
        add("disconnected_s$i";source,pv=0.,operation=:disconnected)
    end
    cases
end

# Fixed phase-neutral PQ loads give an independent injection formulation in
# the same engine. Four-wire return current is carried by their neutral ports.
function four_leg_pq_replay(net,result)
    result.solve.publishable || return (accepted=false,termination="NOT_RUN",max_difference_V=NaN)
    fixed=deepcopy(net)
    for (id,c) in result.spec.devices
        delete!(fixed["ibr"],id)
        result.operating_modes[id]===:disconnected && continue
        inv=c.device
        inv.topology===:FOUR_LEG && inv.neutral=="n" && iszero(inv.b_filter_shunt) ||
            throw(ArgumentError("four-leg replay requires explicit n and no POC shunt"))
        ts=result.network["bus"][inv.bus];vn=phasor(ts["n"])
        for (k,phase) in enumerate(inv.phase_terminals)
            s=(phasor(ts[phase])-vn)*conj(result.grid_currents[id][k])
            fixed["load"]["replay_$(id)_$phase"]=Dict("bus"=>inv.bus,
                "terminal_map"=>[phase,"n"],"configuration"=>"SINGLE_PHASE",
                "model"=>"constant_power","v_nom"=>[HARDWARE.v_nom],
                "p_nom"=>[-real(s)],"q_nom"=>[-imag(s)])
        end
    end
    pf=solve_pf(fixed;per_unit=true,s_base=VA_BASE,solver_options=VA_OPTIONS)
    status=String(pf["termination_status"])
    good=status in ("LOCALLY_SOLVED","OPTIMAL") && get(pf,"feasible",false)
    delta=good ? maximum(abs(phasor(t)-phasor(pf["bus"][b][p]))
        for (b,ts) in result.network["bus"] for (p,t) in ts) : NaN
    (accepted=good && isfinite(delta) && delta<=.001,termination=status,max_difference_V=delta)
end

function run_four_leg(;output=joinpath(ROOT,"results","four_leg"))
    isdir(output) && !isempty(readdir(output)) && throw(ArgumentError("output must be empty"))
    mkpath(joinpath(output,"points"));repo=normpath(joinpath(ROOT,"..",".."))
    source=study_source_sha256()
    for file in ("voltage_aware.jl","four_leg.jl","test_four_leg.jl","audit_four_leg.jl")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    planned=four_leg_cases()
    writejson(joinpath(output,"environment.json"),Dict("planned_cases"=>planned,"source_sha256"=>source,
        "git_revision"=>strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String)),
        "git_status"=>read(Cmd(["git","-C",repo,"status","--porcelain"]),String),
        "Manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest.toml")))),
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Project.toml")))),
        "network_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json")))),
        "hardware"=>HARDWARE,"curves"=>AUS_A,"solver_options"=>Dict(VA_OPTIONS),"s_base_VA"=>VA_BASE,
        "topology_assumptions"=>"FOUR_LEG with explicit neutral, 10 A neutral RMS rating, neutral R/X equal to phase arm; phase/DC sizing unchanged",
        "policy"=>"exact per-phase PN VV/VW, shared Q-vector service fraction first, equal phase P second, shared DC source",
        "scope"=>"counterfactual four-wire research control; not Australian connection conformity, switching or dynamics validation"))
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    cases=NamedTuple[];devices=NamedTuple[];checks=NamedTuple[];history=NamedTuple[]
    statuses=Dict{String,Any}()
    function checkpoint()
        for (file,rows) in (("cases.csv",cases),("devices.csv",devices),("physical_checks.csv",checks),("iterations.csv",history))
            write_csv(joinpath(output,file),rows)
        end
        writejson(joinpath(output,"summary.json"),Dict("planned"=>length(planned),"completed"=>length(cases),
            "accepted"=>count(c->c.accepted,cases),"all_planned_accepted"=>length(cases)==length(planned) && all(c.accepted for c in cases),
            "cases"=>statuses,"unresolved_case_ids"=>[c.case_id for c in cases if !c.accepted]))
    end
    try
        for case in planned
            net=voltage_aware_network(base,case);spec=four_leg_fleet(net,case.pv_fraction)
            result=solve_voltage_aware_inverter_fleet(net,spec;operating_modes=case.operating_mode,
                s_base=VA_BASE,solver_options=VA_OPTIONS)
            replay=four_leg_pq_replay(net,result);accepted=result.solve.publishable && replay.accepted
            rows=voltage_aware_inverter_rows(result);prefix=(case_id=case.case_id,)
            append!(devices,[merge(prefix,row,(accepted=accepted,)) for row in rows])
            append!(history,[merge(prefix,row) for row in result.history])
            for (id,report) in result.physical_reports,c in report.checks
                push!(checks,merge(prefix,(device_id=id,name=c.name,unit=c.unit,value=c.value,
                    lower=c.lower,upper=c.upper,tolerance=c.tolerance,margin=c.margin,passed=c.passed)))
            end
            push!(cases,merge(case,(accepted=accepted,coupling_accepted=result.solve.publishable,
                termination=result.termination_status,reason=result.reason,iterations=result.iterations,
                pq_replay_accepted=replay.accepted,pq_replay_difference_V=replay.max_difference_V,
                p_poc_W=sum(r.p_poc_W for r in rows),q_poc_var=sum(r.q_poc_var for r in rows),
                vmax_sensed_V=maximum(r.vmax_sensed_V for r in rows),vmin_sensed_V=minimum(r.vmin_sensed_V for r in rows),
                max_vuf_percent=maximum(r.vuf_percent for r in rows),
                max_neutral_current_A=maximum(r.neutral_total_rms_current_A for r in rows),
                min_reactive_service_scale=isempty(result.devices) ? nothing : minimum(r.reactive_scale for r in values(result.devices)))))
            statuses[case.case_id]=last(cases)
            pairs(v)=[[real(z),imag(z)] for z in v]
            saved=Dict(id=>Dict("converter_voltage"=>pairs(p.state.converter_voltage),
                "converter_current"=>pairs(p.state.converter_current),"grid_current"=>pairs(p.state.grid_current),
                "reactive_scale"=>haskey(result.devices,id) ? result.devices[id].reactive_scale : nothing) for (id,p) in result.replayed)
            writejson(joinpath(output,"points",case.case_id*".json"),Dict("bus"=>result.network["bus"],"devices"=>saved))
            checkpoint()
            println(case.case_id," ",result.reason," accepted=",accepted," P=",last(cases).p_poc_W,
                " neutral=",last(cases).max_neutral_current_A);flush(stdout)
        end
    finally
        checkpoint()
    end
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_four_leg(;output=isempty(ARGS) ? joinpath(ROOT,"results","four_leg") : ARGS[1])
end
