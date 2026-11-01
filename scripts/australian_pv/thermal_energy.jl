# Mixed-controller fleet factory: mean, worst-phase, sequence and four-leg sites.
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
