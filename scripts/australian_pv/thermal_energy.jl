# Matched fixed-sizing heterogeneous controllers in one coupled feeder case.
include("four_leg.jl")
const MIXED_CONTROLLERS=("mean_vv_vw","worst_vv_vw","sequence_vv_vw","phase_vv_vw")
const MIXED_PV_FRACTIONS=(.25,.6,1.)
_mixed_json_safe(x::Real)=isfinite(x) ? x : nothing
_mixed_json_safe(x::NamedTuple)=Dict(string(k)=>_mixed_json_safe(v) for (k,v) in pairs(x))
_mixed_json_safe(x::AbstractDict)=Dict(string(k)=>_mixed_json_safe(v) for (k,v) in x)
_mixed_json_safe(x::Union{AbstractArray,Tuple})=[_mixed_json_safe(v) for v in x]
_mixed_json_safe(x::Complex)=Dict("re"=>_mixed_json_safe(real(x)),"im"=>_mixed_json_safe(imag(x)))
_mixed_json_safe(x)=x
write_mixed_json(path,x)=writejson(path,_mixed_json_safe(x))

function mixed_thermal_cases()
    cases=NamedTuple[]
    for model in (:fundamental,:pwm,:energy),source in (1.,1.03)
        push!(cases,(case_id="$(model)_s$(source==1. ? 1 : 2)",model=model,pv_fraction=1.,
            source_pu=source,source_vuf=.01,source_negative_angle_rad=0.,phase_shares=(.60,.25,.15),
            load_scale=1.,operation=:energized))
    end
    for operation in (:night,:mixed_disconnection)
        push!(cases,(case_id="energy_$(operation)",model=:energy,pv_fraction=operation===:night ? 0. : 1.,
            source_pu=1.03,source_vuf=.01,source_negative_angle_rad=0.,phase_shares=(.60,.25,.15),
            load_scale=1.,operation=operation))
    end
    cases
end

function mixed_thermal_fleet(net,case)
    four=four_leg_fleet(net,case.pv_fraction)
    devices=Dict{String,ControlledDevice}();requests=Dict{String,InverterControlRequest}()
    modes=Dict{String,Symbol}();labels=Dict{String,String}()
    for (index,id) in enumerate(sort!(collect(keys(net["ibr"]))))
        label=MIXED_CONTROLLERS[mod1(index,4)];labels[id]=label
        base=label=="phase_vv_vw" ? four.devices[id].device : plant(id,net["ibr"][id]["bus"];settings=TERMINAL_CONTROL_SETTINGS)
        fields=NamedTuple{fieldnames(AdvancedInverter)}(Tuple(getfield(base,k) for k in fieldnames(AdvancedInverter)))
        inv=AdvancedInverter(;merge(fields,(pwm_strategy=case.model===:fundamental ? :NONE : :CENTERED,
            f_sw=case.model===:fundamental ? nothing : 10000.,
            dc_energy_model=case.model===:energy ? :nonlinear : :linearized))...)
        law=label=="phase_vv_vw" ? four.devices[id].controller :
            controller(label;design=(kappa=0.,phi2=0.),settings=TERMINAL_CONTROL_SETTINGS)
        devices[id]=ControlledDevice(inv,law)
        net["ibr"][id]["topology"]=label=="phase_vv_vw" ? "FOUR_LEG" : "THREE_LEG"
        net["ibr"][id]["terminal_map"]=label=="phase_vv_vw" ? ["a","b","c","n"] : ["a","b","c"]
        requests[id]=InverterControlRequest(p_available=5000case.pv_fraction*MIXED_PV_FRACTIONS[mod1(index,3)],
            p_rated=HARDWARE.s_max,q_scale=HARDWARE.s_max)
        modes[id]=case.operation===:mixed_disconnection && index%3==0 ? :disconnected : :energized
    end
    (spec=ControlledInverterFleetSpec(devices,requests),modes=modes,labels=labels)
end

function mixed_thermal_replay(net,result)
    result.solve.publishable || return (accepted=false,termination="NOT_RUN",max_difference_V=NaN)
    fixed=deepcopy(net)
    for (id,controlled) in result.spec.devices
        delete!(fixed["ibr"],id)
        result.operating_modes[id]===:disconnected && continue
        inv=controlled.device
        iszero(inv.b_filter_shunt) || error("mixed replay requires zero POC shunt")
        current=result.grid_currents[id]
        inv.neutral===nothing && abs(sum(current))>=1e-5 && error("three-wire KCL failed")
        ts=result.network["bus"][inv.bus]
        for (k,phase) in enumerate(inv.phase_terminals)
            s=(phasor(ts[phase])-phasor(ts["n"]))*conj(current[k])
            fixed["load"]["mixed_replay_$(id)_$phase"]=Dict("bus"=>inv.bus,
                "terminal_map"=>[phase,"n"],"configuration"=>"SINGLE_PHASE","model"=>"constant_power",
                "v_nom"=>[230.],"p_nom"=>[-real(s)],"q_nom"=>[-imag(s)])
        end
    end
    pf=solve_pf(fixed;per_unit=true,s_base=VA_BASE,solver_options=VA_OPTIONS)
    status=String(pf["termination_status"])
    good=status in ("LOCALLY_SOLVED","OPTIMAL") && get(pf,"feasible",false)
    delta=good ? maximum(abs(phasor(t)-phasor(pf["bus"][b][p]))
        for (b,ts) in result.network["bus"] for (p,t) in ts) : NaN
    (accepted=good && isfinite(delta) && delta<=.001,termination=status,max_difference_V=delta)
end

function run_mixed_thermal(;output=joinpath(ROOT,"results","thermal_energy"))
    isdir(output) && !isempty(readdir(output)) && error("output must be empty")
    mkpath(joinpath(output,"points"));planned=mixed_thermal_cases()
    repo=normpath(joinpath(ROOT,"..",".."));source=study_source_sha256()
    for file in ("four_leg.jl","voltage_aware.jl","thermal_energy.jl","audit_thermal_energy.jl","test_thermal_energy.jl")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    writejson(joinpath(output,"environment.json"),Dict("source_sha256"=>source,"planned_cases"=>planned,
        "git_revision"=>strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String)),
        "git_status"=>read(Cmd(["git","-C",repo,"status","--porcelain"]),String),
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Project.toml")))),
        "Manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest.toml")))),
        "network_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"networks","LV3_55bus.bmopf.json")))),
        "hardware"=>HARDWARE,"controllers"=>MIXED_CONTROLLERS,"site_PV_fractions"=>MIXED_PV_FRACTIONS,
        "solver_options"=>Dict(VA_OPTIONS),"s_base_VA"=>VA_BASE,"switching_Hz"=>10000.,
        "pwm_thermal_headroom"=>PowerOptLab._PWM_ALLOCATION_HEADROOM,
        "scope"=>"coupled mixed-controller three/four-leg feeder equilibria; phase/DC sizing unchanged; no added 253 V constraint, regional conformity or dynamic-stability claim",
        "source_assumptions"=>"open carrier-frequency source, constant-power input for nonlinear energy; mean converter loss; representative EMT required"))
    cases=NamedTuple[];devices=NamedTuple[];checks=NamedTuple[];history=NamedTuple[]
    statuses=Dict{String,Any}();base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    function checkpoint()
        for (file,rows) in (("cases.csv",cases),("devices.csv",devices),("physical_checks.csv",checks),("iterations.csv",history))
            write_csv(joinpath(output,file),rows)
        end
        write_mixed_json(joinpath(output,"summary.json"),Dict("planned"=>length(planned),"completed"=>length(cases),
            "accepted"=>count(c->c.accepted,cases),"all_planned_accepted"=>length(cases)==length(planned)&&all(c.accepted for c in cases),
            "unresolved_case_ids"=>[c.case_id for c in cases if !c.accepted],"cases"=>statuses))
    end
    try
        for case in planned
            net=voltage_aware_network(base,case);mixed=mixed_thermal_fleet(net,case)
            r=solve_voltage_aware_inverter_fleet(net,mixed.spec;operating_modes=mixed.modes,
                s_base=VA_BASE,solver_options=VA_OPTIONS)
            replay=mixed_thermal_replay(net,r);accepted=r.solve.publishable && replay.accepted
            rows=voltage_aware_inverter_rows(r)
            append!(devices,[merge((case_id=case.case_id,controller=mixed.labels[row.device_id]),row,(accepted=accepted,)) for row in rows])
            append!(history,[merge((case_id=case.case_id,),row) for row in r.history])
            for (id,report) in r.physical_reports,c in report.checks
                push!(checks,(case_id=case.case_id,device_id=id,name=c.name,unit=c.unit,value=c.value,
                    lower=c.lower,upper=c.upper,tolerance=c.tolerance,margin=c.margin,passed=c.passed))
            end
            push!(cases,merge(case,(accepted=accepted,reason=r.reason,termination=r.termination_status,
                iterations=r.iterations,pq_replay_accepted=replay.accepted,pq_replay_difference_V=replay.max_difference_V,
                p_poc_W=sum(row.p_poc_W for row in rows),q_poc_var=sum(row.q_poc_var for row in rows),
                vmax_sensed_V=maximum(row.vmax_sensed_V for row in rows),
                maximum_capacitor_model_A=maximum((row.capacitor_thermal_current_A for row in rows if row.capacitor_thermal_current_A!==nothing);init=0.))))
            statuses[case.case_id]=last(cases)
            pairs(v)=[[real(z),imag(z)] for z in v]
            records=Dict(id=>Dict("converter_voltage"=>pairs(p.state.converter_voltage),
                "converter_current"=>pairs(p.state.converter_current),"grid_current"=>pairs(p.state.grid_current),
                "phase_reactive_scale"=>haskey(r.devices,id)&&r.devices[id] isa PhaseInverterServiceAllocationResult ? r.devices[id].reactive_scale : nothing)
                for (id,p) in r.replayed)
            stages=Dict(id=>[(name=n,solve=solve_status(s),operation=s.operation) for (n,s) in a.stages] for (id,a) in r.devices)
            write_mixed_json(joinpath(output,"points",case.case_id*".json"),Dict("bus"=>r.network["bus"],"devices"=>records,"stages"=>stages))
            checkpoint()
            println(case.case_id," ",r.reason," accepted=",accepted," P=",last(cases).p_poc_W," Q=",last(cases).q_poc_var);flush(stdout)
        end
    finally
        checkpoint()
    end
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_mixed_thermal(;output=isempty(ARGS) ? joinpath(ROOT,"results","thermal_energy") : ARGS[1])
end
