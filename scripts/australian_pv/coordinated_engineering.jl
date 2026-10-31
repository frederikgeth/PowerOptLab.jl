# Matched capabilities: retain balanced converter current at three-leg sites
# and equal per-phase active P at four-leg sites. Coordination releases the
# local VV/VW requests, not any sizing or additional sequence/phase-P freedom.
include("switching_engineering.jl")

_engineering_json_safe(x::Real)=isfinite(x) ? x : nothing
_engineering_json_safe(x::Complex)=[_engineering_json_safe(real(x)),_engineering_json_safe(imag(x))]
_engineering_json_safe(x::Union{AbstractArray,Tuple})=[_engineering_json_safe(v) for v in x]
_engineering_json_safe(x::NamedTuple)=Dict(string(k)=>_engineering_json_safe(v) for (k,v) in pairs(x))
_engineering_json_safe(x::AbstractDict)=Dict(string(k)=>_engineering_json_safe(v) for (k,v) in x)
_engineering_json_safe(x::Union{SolveStatus,InverterLimitCheck,InverterFeasibilityReport,ConverterTerminalResult})=
    Dict(string(k)=>_engineering_json_safe(getfield(x,k)) for k in fieldnames(typeof(x)))
_engineering_json_safe(x)=x
write_engineering_json(path,x)=writejson(path,_engineering_json_safe(x))
const ENGINEERING_COORDINATED_OPTIONS=merge(ENGINEERING_OPTIONS,
    Dict("mu_strategy"=>"adaptive","nlp_scaling_method"=>"gradient-based","max_iter"=>2000))

function engineering_conditions()
    [(case_id=id,model=:energy,pv_fraction=pv,source_pu=source,source_vuf=vuf,
        source_negative_angle_rad=angle,phase_shares=(.60,.25,.15),load_scale=1.,
        operation=iszero(pv) ? :night : :energized)
        for (id,source,pv,vuf,angle) in (("nominal",1.,1.,.01,0.),("high_source",1.03,1.,.01,0.),
            ("rotated_unbalance",1.,1.,.02,pi/2),("night",1.03,0.,.01,0.))]
end

function engineering_bus_voltages(ctx)
    vr,vi=PowerOptLab._opf_voltage_maps(ctx)
    buses=Dict{String,Any}()
    for ((bus,terminal),re) in vr
        vb=BMOPFTools.opf_ac_coordinate_bases(ctx,bus).voltage
        v=complex(JuMP.value(re),JuMP.value(vi[(bus,terminal)]))*vb
        get!(buses,bus,Dict{String,Any}())[terminal]=Dict(
            "vr"=>real(v),"vi"=>imag(v),"vm"=>abs(v),"va"=>angle(v))
    end
    buses
end

function engineering_device_phasors(handles)
    pairs(x)=[[real(z),imag(z)] for z in x]
    phasors(re,im,scale)=complex.(JuMP.value.(re),JuMP.value.(im)).*scale
    Dict(id=>Dict(
        "converter_voltage"=>pairs(phasors(h.plant.vrint,h.plant.viint,h.plant.vb)),
        "converter_current"=>pairs(phasors(h.plant.cri,h.plant.cii,h.plant.ib)),
        "grid_current"=>pairs(phasors(h.plant.gri,h.plant.gii,h.plant.ib))) for (id,h) in handles)
end

function engineering_night_seed(inv,voltage)
    if inv.topology===:THREE_LEG
        a=cis(2pi/3);positive=ComplexF64[1,a^2,a]
        u1=sum(voltage.*conj.(positive))/3
        current=conj(-inv.p_loss_fixed/(3u1)).*positive
    else
        # Equal phase P and zero neutral current, with the least-norm Q seed.
        # Only the fixed mean loss is estimated; the actual NLP still enforces
        # the complete nonlinear source/loss balance and all engineering rows.
        A=zeros(5,6)
        for k in 1:3
            A[k,2k-1]=real(voltage[k]);A[k,2k]=imag(voltage[k])
            A[4,2k-1]=1.;A[5,2k]=1.
        end
        x=A\[-fill(inv.p_loss_fixed/3,3);0.;0.]
        current=complex.(x[1:2:6],x[2:2:6])
    end
    R,X=PowerOptLab._series_filter_matrices(inv,3,1.;side=:converter)
    primitive=complex.(R,X)*(inv.neutral===nothing ? current : [current;-sum(current)])
    drop=inv.neutral===nothing ? primitive : primitive[1:3].-primitive[4]
    pairs(x)=[[real(z),imag(z)] for z in x]
    Dict("converter_voltage"=>pairs(voltage+drop),"converter_current"=>pairs(current),"grid_current"=>pairs(current))
end

function engineering_point(result,case;coordinated=false,voltage_max=253.)
    pairs(x)=[[real(z),imag(z)] for z in x]
    devices=if coordinated
        result.raw.has_primal ? engineering_device_phasors(result.handles) : Dict()
    else
        Dict(id=>Dict("converter_voltage"=>pairs(d.state.converter_voltage),
            "converter_current"=>pairs(d.state.converter_current),"grid_current"=>pairs(d.state.grid_current))
            for (id,d) in result.replayed)
    end
    buses=coordinated && result.raw.has_primal ? engineering_bus_voltages(result.model_context) : result.network["bus"]
    Dict("case"=>case,"planning_voltage_max_V"=>voltage_max,"bus"=>buses,"devices"=>devices,
        "solver_status"=>coordinated ? result.raw : result.solve,"candidate_status"=>result.solve,
        "coordinated"=>coordinated)
end

function solve_engineering_coordinated(net,mixed;local_point=nothing,minimum_p_W=nothing,
        voltage_max=253.,s_base=VA_BASE,solver_options=ENGINEERING_COORDINATED_OPTIONS,
        continuation_weights=(0.,),initialize_neutral=true,objective_scale=1.)
    !isempty(continuation_weights) && last(continuation_weights)==0. ||
        error("the final objective must be unregularized total active export")
    all(w->isfinite(w)&&w>=0.,continuation_weights) || error("invalid continuation weights")
    isfinite(objective_scale) && objective_scale>0 || error("objective scale must be finite and positive")
    any(string(k)=="obj_scaling_factor" for (k,_) in
        (solver_options isa NamedTuple ? pairs(solver_options) : solver_options)) &&
        error("use objective_scale: a raw Ipopt obj_scaling_factor can override the MAX objective sign")
    PowerOptLab._validate_controlled_fleet(net,mixed.spec)
    all(==(:energized),values(mixed.modes)) || error("coordinated benchmark requires energized sites")
    devices=Dict(id=>ConstrainedInverter(c.device,InverterOperatingConstraints();
        current_allocation=c.device.topology===:THREE_LEG ? :balanced_converter : :free,
        dc_source=DCSourceBounds(p_min=0.)) for (id,c) in mixed.spec.devices)
    handles=Dict{String,Any}()
    objectives=Dict{Symbol,Any}()
    builder=BMOPFTools.OpfDeviceBuilder(:PowerOptLab,function(ctx,ids)
        for id in ids
            PowerOptLab._disable_native_ibr_port!(ctx,id)
            handles[id]=stamp_device!(ctx,devices[id];request=mixed.spec.requests[id])
        end
        ctx
    end)
    build_spec=BMOPFTools.OpfBuildSpec(component_builders=Dict((:ibr,id)=>builder for id in keys(devices)))
    source_buses=Set(s["bus"] for s in values(net["voltage_source"]))
    lv_buses=sort!([b for b in keys(net["bus"]) if !(b in source_buses)])
    hook! = ctx -> begin
        m=PowerOptLab._opf_model(ctx)
        for bus in lv_buses,phase in ("a","b","c")
            vb=BMOPFTools.opf_ac_coordinate_bases(ctx,bus).voltage
            re,im=PowerOptLab._dv(ctx,bus,phase,"n")
            @constraint(m,(re^2+im^2)/(voltage_max/vb)^2<=1.)
        end
        for (id,device) in devices
            inv=device.device;h=handles[id].plant
            if inv.topology===:FOUR_LEG
                iszero(inv.b_filter_shunt) || error("matched phase-P scope excludes POC shunts")
                p=[@expression(m,h.vrpoc[k]*h.gri[k]+h.vipoc[k]*h.gii[k]) for k in 1:3]
                for k in 2:3
                    @constraint(m,(p[k]-p[1])/(inv.s_max/h.sb)==0.)
                end
            end
        end
        rating=sum(d.device.s_max for d in values(devices))
        total=@expression(m,sum(h.plant.p_poc*h.plant.sb for h in values(handles)))
        penalty=@expression(m,sum(((h.plant.vipoc[k]*h.plant.gri[k]-
            h.plant.vrpoc[k]*h.plant.gii[k])*h.plant.sb/(devices[id].device.s_max/3))^2
            for (id,h) in handles for k in 1:3)/(3length(handles)))
        objectives[:export]=total/rating;objectives[:regularizer]=penalty
        if minimum_p_W!==nothing
            @constraint(m,total/rating>=(minimum_p_W-.01)/rating)
        end
        @objective(m,Max,total/rating)
    end
    ctx=PowerOptLab.build_opf_model(net;per_unit=true,s_base=Float64(s_base),add_objective=false,
        build_spec,model_hook! = hook!,optimizer=Ipopt.Optimizer)
    m=PowerOptLab._opf_model(ctx)
    PowerOptLab._set_inverter_solver_options!(m,solver_options)
    if local_point!==nothing
        vr,vi=PowerOptLab._opf_voltage_maps(ctx)
        seed(variable,value)=variable isa JuMP.VariableRef ? JuMP.set_start_value(variable,value) : nothing
        for (bus,ts) in local_point["bus"],(phase,record) in ts
            vb=BMOPFTools.opf_ac_coordinate_bases(ctx,bus).voltage
            seed(vr[(bus,phase)],record["vr"]/vb)
            seed(vi[(bus,phase)],record["vi"]/vb)
        end
        for (id,record) in local_point["devices"]
            h=handles[id].plant;inv=devices[id].device
            if iszero(mixed.spec.requests[id].p_available)
                ts=local_point["bus"][inv.bus]
                neutral=inv.neutral===nothing ? 0im : phasor(ts[inv.neutral])
                voltage=ComplexF64[phasor(ts[phase])-neutral for phase in inv.phase_terminals]
                record=engineering_night_seed(inv,voltage)
            end
            for (key,re,im,scale) in (("converter_voltage",h.vrint,h.viint,h.vb),
                    ("converter_current",h.cri,h.cii,h.ib),("grid_current",h.gri,h.gii,h.ib))
                for k in 1:3
                    seed(re[k],record[key][k][1]/scale)
                    seed(im[k],record[key][k][2]/scale)
                end
            end
            u=ComplexF64[complex(x...) for x in record["converter_voltage"]]
            i=ComplexF64[complex(x...) for x in record["converter_current"]]
            e=PowerOptLab._dc_energy_state(inv,sum(u.*i))
            JuMP.set_start_value(h.dc_energy.center,e.center_V2/inv.v_dc^2)
            JuMP.set_start_value(h.dc_energy.radius,sqrt(e.center_V2^2-(2inv.v_dc*abs(e.energy_coordinate))^2)/inv.v_dc^2)
            JuMP.set_start_value(h.dc_energy.er,real(e.energy_coordinate))
            JuMP.set_start_value(h.dc_energy.ei,imag(e.energy_coordinate))
            JuMP.set_start_value(h.dre,real(e.voltage_2omega));JuMP.set_start_value(h.dim,imag(e.voltage_2omega))
        end
    end
    function run_stage!(name,weight)
        @objective(m,Max,objective_scale*(objectives[:export]-weight*objectives[:regularizer]))
        PowerOptLab._optimize_inverter_rails!(ctx,[(d.device,handles[id].plant) for (id,d) in devices])
        stage=solve_status(m)
        stage_devices=Dict(id=>extract_device(d,handles[id],stage) for (id,d) in devices)
        result=(stage=name,weight=weight,raw=stage,raw_solver_status=JuMP.raw_status(m),
            devices=stage_devices,bus=stage.has_primal ? engineering_bus_voltages(ctx) : Dict(),
            phasors=stage.has_primal ? engineering_device_phasors(handles) : Dict(),
            total_p_W=stage.has_primal ? sum(JuMP.value(h.plant.p_poc)*h.plant.sb for h in values(handles)) : NaN)
        println("OPF stage=",name," weight=",weight," raw=",stage.termination_status,
            " P=",result.total_p_W);flush(stdout)
        result
    end
    # The fitted fourth-leg norm is very sharp at cancellation. Start on its
    # zero-current subspace, then remove every temporary row before solving
    # the original problem. This changes neither the final feasible set nor
    # its objective; initialization records cannot be accepted as the result.
    neutral_groups=Dict{String,Vector{JuMP.ConstraintRef}}()
    if initialize_neutral
        for (id,d) in devices
            d.device.topology===:FOUR_LEG || continue
            h=handles[id].plant
            group=JuMP.ConstraintRef[@constraint(m,sum(h.cri)==0.),@constraint(m,sum(h.cii)==0.)]
            neutral_groups[id]=group
        end
    end
    order=sort!(collect(keys(neutral_groups));by=id->(-mixed.spec.requests[id].p_available,id))
    continuation=PowerOptLab._inverter_solver_continuation!(m,
        [id=>neutral_groups[id] for id in order],continuation_weights;
        solve_stage! = run_stage!,accepted=r->r.raw.publishable &&
            all(d.plant.physical_feasibility.passed for d in values(r.devices)))
    attempts=continuation.attempts
    raw=solve_status(m)
    extracted=Dict(id=>extract_device(d,handles[id],raw) for (id,d) in devices)
    outcome=PowerOptLab._solve_outcome(m)
    network=PowerOptLab._network_without_replaced_ibrs(PowerOptLab._extract_result(ctx,outcome),keys(devices))
    equations,checked=PowerOptLab._voltage_aware_network_checks(ctx,raw)
    voltage_checks=InverterLimitCheck[]
    if raw.has_primal
        candidate_buses=engineering_bus_voltages(ctx)
        for bus in lv_buses,phase in ("a","b","c")
            ts=candidate_buses[bus]
            push!(voltage_checks,PowerOptLab._inverter_limit_check("planning_voltage/$bus/$phase",:V,
                abs(phasor(ts[phase])-phasor(ts["n"])),0.,voltage_max))
        end
    end
    physical=all(d.plant.physical_feasibility.passed for d in values(extracted))
    candidate=raw.publishable && physical && checked && all(c.passed for c in [equations;voltage_checks])
    result=(raw=raw,solve=SolveStatus(raw.termination_status,raw.primal_status,raw.has_primal,raw.feasible,raw.optimal,candidate),
        devices=extracted,spec=mixed.spec,operating_modes=mixed.modes,
        grid_currents=Dict(id=>d.grid_phase_current for (id,d) in extracted),network=network,
        network_checks=equations,voltage_checks=voltage_checks,
        raw_solver_status=JuMP.raw_status(m),model_candidate_accepted=candidate,
        model_context=ctx,handles=handles,attempts=attempts,
        temporary_neutral_constraints_remaining=continuation.temporary_constraints_remaining,
        scope="matched local current capabilities; released VV/VW commands; independent network audit still required")
    replay=mixed_thermal_replay(net,result)
    accepted=candidate && replay.accepted
    merge(result,(replay=replay,model_candidate_accepted=accepted,
        solve=SolveStatus(raw.termination_status,raw.primal_status,raw.has_primal,raw.feasible,raw.optimal,accepted)))
end
