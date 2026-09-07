include("prepare.jl")
include(joinpath(@__DIR__,"..","enwl_advanced_ibr","common.jl"))
using .ENWLAdvancedIBR: write_csv
const MODES = ["unity", "mean_vv", "mean_vv_vw", "worst_vv_vw", "sequence_vv_vw", "sequence_droop"]
function controller(mode)
    # Illustrative curves, NOT a declaration of AS/NZS 4777 compliance.
    vv = PiecewiseLinearLaw([207.,220.,240.,258.],[0.44,0.,0.,-0.60]; smoothing_epsilon=0.05)
    vw = PiecewiseLinearLaw([250.,260.],[1.,0.]; smoothing_epsilon=0.05)
    positive = if mode == "unity"
        AverageVoltageVoltVarWatt()
    elseif mode == "mean_vv"
        AverageVoltageVoltVarWatt(volt_var=vv)
    elseif mode == "mean_vv_vw"
        AverageVoltageVoltVarWatt(volt_var=vv,volt_watt=vw)
    elseif mode == "worst_vv_vw"
        WorstPhaseVoltVarWatt(volt_var=vv,volt_watt=vw)
    elseif mode in ("sequence_vv_vw","sequence_droop")
        PositiveSequenceVoltVarWatt(volt_var=vv,volt_watt=vw)
    else
        error("Unknown mode $mode")
    end
    unbalance = mode == "sequence_droop" ? NegativeSequenceAdmittanceDroop(
        PiecewiseLinearLaw([0.,0.0005,0.003,0.01],[0.,0.,1.5,3.0]; smoothing_epsilon=0.000025);
        impedance_angle=atan(0.3), ripple_blend=0.0) : NoUnbalanceControl()
    SequenceController(positive; unbalance, current_target=:grid)
end
function study_network(base, availability, source_pu)
    net = deepcopy(base)
    for source in values(net["voltage_source"])
        # Controlled 1% upstream negative-sequence voltage; no zero sequence.
        # This is an explicit study perturbation, not a measured site condition.
        v=Float64.(source["v_magnitude"]) .* cis.(Float64.(source["v_angle"]))
        a=cis(2pi/3); u1=(v[1]+a*v[2]+a^2*v[3])/3*source_pu
        v[1:3] = u1 .* [1,a^2,a] .+ 0.01u1 .* [1,a,a^2]
        source["v_magnitude"]=abs.(v); source["v_angle"]=angle.(v)
    end
    # Upgrade the supply connection at each of the 12 customers. Three
    # explicit constant-PQ phase loads retain deliberately unbalanced demand.
    # Total demand is 1.5 kW/customer, PF=.95; A/B/C shares are 60/25/15%.
    net["load"] = Dict{String,Any}()
    net["ibr"] = Dict{String,Any}()
    for (customer,ld) in base["load"]
        bus=ld["bus"]
        @assert all(t in net["bus"][bus]["terminal_names"] for t in ["a","b","c","n"])
        for (phase,share) in zip(["a","b","c"],[0.60,0.25,0.15])
            record=deepcopy(ld)
            record["terminal_map"]=[phase,"n"]
            record["p_nom"]=[1500.0*share]
            record["q_nom"]=[1500.0*share*tan(acos(0.95))]
            net["load"][customer*"_"*phase]=record
        end
        net["ibr"]["pv_"*customer] = Dict("bus"=>bus,"topology"=>"THREE_LEG",
            "terminal_map"=>["a","b","c"],"s_max"=>fill(6000/3,3),
            "p_min"=>zeros(3),"p_max"=>fill(6000/3,3),
            "q_min"=>fill(-6000/3,3),"q_max"=>fill(6000/3,3))
    end
    net
end
function fleet(net, mode, availability)
    devices=Dict{String,ControlledDevice}(); requests=Dict{String,InverterControlRequest}()
    for (id,record) in net["ibr"]
        plant=AdvancedInverter(id=id,bus=record["bus"],topology=:THREE_LEG,
            phase_terminals=["a","b","c"],neutral=nothing,s_max=6000.,i_max=10.,
            v_dc=750.,c_dc=0.002,m_max=0.96,r_filter=0.05,x_filter=0.15)
        devices[id]=ControlledDevice(plant,controller(mode))
        requests[id]=InverterControlRequest(p_available=5000.0*availability,p_rated=6000.,q_scale=6000.)
    end
    ControlledInverterFleetSpec(devices,requests)
end
function customer_rows(base, buses, mode, availability, source_pu)
    a=cis(2pi/3)
    [begin
        ts=buses[ld["bus"]]; v=[phasor(ts[p]) for p in ("a","b","c")]
        vn=phasor(ts["n"]); phase=first(ld["terminal_map"])
        v1=(v[1]+a*v[2]+a^2*v[3])/3; v2=(v[1]+a^2*v[2]+a*v[3])/3
        (mode=mode,availability=availability,source_pu=source_pu,source_vuf_percent=1.0,customer=id,bus=ld["bus"],phase=phase,
         original_phase_voltage_V=abs(phasor(ts[phase])-vn),
         customer_voltage_V=maximum(abs.(v .- vn)),
         customer_vmin_V=minimum(abs.(v .- vn)),
         va_pn_V=abs(v[1]-vn),vb_pn_V=abs(v[2]-vn),vc_pn_V=abs(v[3]-vn),
         neutral_V=abs(vn),vuf_percent=100abs(v2)/abs(v1),
         vmax_pg_V=maximum(abs,v),vmax_pn_V=maximum(abs.(v .- vn)))
    end for (id,ld) in sort!(collect(base["load"]);by=first)]
end
function run(; modes=MODES, levels=[0.0,0.5,1.0], sources=[1.0,1.03])
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    rows=NamedTuple[]; devices=NamedTuple[]; cases=NamedTuple[]
    output=joinpath(ROOT,"results"); mkpath(output)
    writejson(joinpath(output,"environment.json"), Dict(
        "julia"=>string(VERSION),"PowerOptLab"=>string(Base.pkgversion(PowerOptLab)),
        "BMOPFTools"=>string(Base.pkgversion(BMOPFTools)),
        "Ipopt"=>string(Base.pkgversion(Ipopt)),
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"..","..","Project.toml")))),
        "source_vuf_percent"=>1.0,"modes"=>modes,"availability"=>levels,
        "source_positive_sequence_pu"=>sources))
    for source in sources, level in levels, mode in modes
        println("Solving ", (source,level,mode)); flush(stdout)
        net=study_network(base,level,source)
        result=solve_controlled_inverter_fleet(net,fleet(net,mode,level);
            s_base=100000.,solver_options=OPTIONS)
        dr=controlled_inverter_rows(result)
        cr=customer_rows(base,result.network["bus"],mode,level,source)
        # Main predates PR49: add an independent finite/firmware check here.
        # PWM is NONE: these are fundamental-frequency steady states only.
        finite=all(isfinite(getproperty(r,k)) for r in cr for k in
            (:customer_voltage_V,:neutral_V,:vuf_percent,:vmax_pn_V))
        finite &= all(isfinite(getproperty(r,k)) for r in dr for k in
            (:p_poc_W,:q_poc_var,:converter_loss_W,:maximum_converter_current_A,
             :maximum_grid_current_A,:ripple_power_VA,:dc_ripple_voltage_V,
             :capacitor_current_A,:power_scale,:current_scale))
        residual=maximum(r.exact_smooth_current_residual_A for r in dr)
        accepted=solve_status(result).publishable && finite && isfinite(residual) && residual<0.02 &&
            all(r.maximum_converter_current_A <= 10.0+1e-5 for r in dr)
        push!(cases,(mode=mode,availability=level,source_pu=source,
            source_vuf_percent=1.0,termination=result.termination_status,accepted=accepted,exact_residual_A=residual))
        if accepted
            append!(rows,cr)
            append!(devices,[(;mode,availability=level,source_pu=source,r...) for r in dr])
        end
        write_csv(joinpath(output,"cases.csv"),cases)
        write_csv(joinpath(output,"customers.csv"),rows)
        write_csv(joinpath(output,"devices.csv"),devices)
    end
    @assert all(r.accepted for r in cases) "Rejected cases recorded; omitted from scientific plots"
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run()
end
