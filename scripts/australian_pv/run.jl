include("prepare.jl")
include(joinpath(@__DIR__,"..","enwl_advanced_ibr","common.jl"))
using .ENWLAdvancedIBR: write_csv
using JuMP

const MODES = ["unity", "mean_vv", "mean_vv_vw", "worst_vv_vw", "sequence_vv_vw",
               "sequence_droop", "sequence_droop_power", "gfm_virtual_impedance",
               "single_phase_native"]

# Hardware and control parameters. Every value is derived, sourced or assumed as
# documented in docs/src/ibr/parameter_selection.md and the study README.
const HARDWARE = (
    s_max=6000.0,        # assumed: residential three-phase inverter rating (VA)
    v_nom=230.0,         # sourced: AS 60038 nominal LV phase voltage (V)
    i_max=10.0,          # derived: >= s_max/(3*207 V) = 9.66 A, ~15% above I_rated
    v_dc=750.0,          # derived: >= sqrt(6)*266 V/0.96 = 679 V (three-leg hull)
    m_max=0.96,          # assumed: modulation utilisation
    c_dc=300e-6,         # assumed: 50 µF/kVA film link; swept in sensitivity
    dv2_max=7.5,         # assumed: 1% of v_dc peak 2ω ripple allowance (V)
    x_pu=0.05,           # assumed: 5% series filter reactance (Liserre et al. ≤ 10%)
    x_over_r=10.0,       # assumed: inductor X/R at 50 Hz
    p_loss_fixed=15.0,   # assumed loss curve: ~2.2% at rated, ~3.5% at 10% load
    a_loss=2.0,
    c_loss=0.3,
    single_phase_i_max=28.7,   # derived: 1.1*6 kVA/230 V for the status-quo inverter
)
i_rated(hw=HARDWARE) = hw.s_max/(3hw.v_nom)
z_base(hw=HARDWARE) = hw.v_nom/i_rated(hw)
filter_x(hw=HARDWARE) = hw.x_pu*z_base(hw)
filter_r(hw=HARDWARE) = filter_x(hw)/hw.x_over_r

# AS/NZS 4777.2:2020 Australia A (BMOPFTools Aus_A preset; to be verified against
# the licensed text). Volt-watt is a fraction of the rating (p_ref = S_MAX).
const AUS_A = (vv_v=[207.0, 220.0, 240.0, 258.0], vv_q=[0.44, 0.0, 0.0, -0.60],
               vw_v=[253.0, 260.0], vw_p=[1.0, 0.2], vw_basis=:rated)
# The earlier pedagogical curve, kept as a named sensitivity case only.
const PEDAGOGICAL = (vv_v=AUS_A.vv_v, vv_q=AUS_A.vv_q,
                     vw_v=[250.0, 260.0], vw_p=[1.0, 0.0], vw_basis=:available)
const SENSING_REFERENCE = "n"   # assumed: three-wire power stage, neutral sensing lead
const TARGET_ATTENUATION = 0.8  # design: |U2|/|E2| for the aligned fleet response

function controller(mode; curves=AUS_A, design)
    vv = PiecewiseLinearLaw(curves.vv_v, curves.vv_q; smoothing_epsilon=0.05)
    vw = PiecewiseLinearLaw(curves.vw_v, curves.vw_p; smoothing_epsilon=0.05)
    basis = curves.vw_basis
    positive = if mode == "unity"
        AverageVoltageVoltVarWatt(volt_watt_basis=basis)
    elseif mode == "mean_vv"
        AverageVoltageVoltVarWatt(volt_var=vv, volt_watt_basis=basis)
    elseif mode == "mean_vv_vw"
        AverageVoltageVoltVarWatt(volt_var=vv, volt_watt=vw, volt_watt_basis=basis)
    elseif mode == "worst_vv_vw"
        WorstPhaseVoltVarWatt(volt_var=vv, volt_watt=vw, volt_watt_basis=basis)
    elseif mode in ("sequence_vv_vw", "sequence_droop", "sequence_droop_power",
                    "gfm_virtual_impedance")
        PositiveSequenceVoltVarWatt(volt_var=vv, volt_watt=vw, volt_watt_basis=basis)
    else
        error("Unknown controlled mode $mode")
    end
    unbalance = if mode in ("sequence_droop", "sequence_droop_power")
        NegativeSequenceAdmittanceDroop(
            PiecewiseLinearLaw([0.0, 0.0005, 0.002], [0.0, 0.0, design.kappa];
                               smoothing_epsilon=2.5e-5);
            impedance_angle=design.phi2, ripple_blend=0.0)
    elseif mode == "gfm_virtual_impedance"
        negative_sequence_virtual_impedance(filter_r(), filter_x())
    else
        NoUnbalanceControl()
    end
    limiter = mode in ("sequence_droop_power", "gfm_virtual_impedance") ?
        SequencePriorityLimiter(priority=:power) : CommonScaleLimiter()
    SequenceController(positive, unbalance, limiter;
                       current_target=:grid, sensing_reference=SENSING_REFERENCE)
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
            "terminal_map"=>["a","b","c"],"s_max"=>fill(HARDWARE.s_max/3,3),
            "p_min"=>zeros(3),"p_max"=>fill(HARDWARE.s_max/3,3),
            "q_min"=>fill(-HARDWARE.s_max/3,3),"q_max"=>fill(HARDWARE.s_max/3,3))
    end
    net
end

function plant(id, bus; hw=HARDWARE)
    AdvancedInverter(id=id, bus=bus, topology=:THREE_LEG,
        phase_terminals=["a","b","c"], neutral=nothing, s_max=hw.s_max,
        i_max=hw.i_max, v_dc=hw.v_dc, c_dc=hw.c_dc, dv2_max=hw.dv2_max,
        m_max=hw.m_max, r_filter=filter_r(hw), x_filter=filter_x(hw),
        p_loss_fixed=hw.p_loss_fixed, a_loss=hw.a_loss, c_loss=hw.c_loss)
end

function fleet(net, mode, availability; curves=AUS_A, hw=HARDWARE,
               design=negative_sequence_design())
    devices=Dict{String,ControlledDevice}(); requests=Dict{String,InverterControlRequest}()
    for (id,record) in net["ibr"]
        devices[id]=ControlledDevice(plant(id, record["bus"]; hw=hw),
                                     controller(mode; curves=curves, design=design))
        requests[id]=InverterControlRequest(p_available=5000.0*availability,
            p_rated=hw.s_max, q_scale=hw.s_max)
    end
    ControlledInverterFleetSpec(devices,requests)
end

# Negative-sequence design (docs: "Negative-sequence gain and angle"). The fleet
# voltage sensitivity in sequence coordinates gives each device's self and
# fleet-aggregate negative-sequence impedance. φ2 is the median aggregate angle
# and κ meets TARGET_ATTENUATION for the median aggregate magnitude.
const _DESIGN = Ref{Any}(nothing)
function negative_sequence_design(; recompute::Bool=false)
    !recompute && _DESIGN[] !== nothing && return _DESIGN[]
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    net=study_network(base, 0.0, 1.0)
    spec=fleet(net, "unity", 0.0; design=(kappa=0.0, phi2=0.0))
    ids=sort!(collect(keys(spec.devices)))
    zero_current=Dict(id=>ntuple(_->0.0im,3) for id in ids)
    sens=controlled_inverter_fleet_network_voltage_sensitivity(
        net, spec, zero_current; step=0.5, s_base=100000., solver_options=OPTIONS)
    @assert sens.base_status.publishable && all(s->s.publishable, sens.plus_status)
    Z=sens.sensitivity
    block(i,j)=complex(Z[4(i-1)+3,4(j-1)+3], Z[4(i-1)+4,4(j-1)+3])
    n=length(ids)
    self=[block(i,i) for i in 1:n]
    aggregate=[sum(block(i,j) for j in 1:n) for i in 1:n]
    median(x)=sort(x)[cld(length(x),2)]
    phi2=median(angle.(aggregate))
    g=1/TARGET_ATTENUATION-1
    kappa=g/median(abs.(aggregate))/n
    design=(kappa=kappa, phi2=phi2, ids=ids,
        self_impedance=self, aggregate_impedance=aggregate,
        kappa_per_unit=kappa*z_base(), target_attenuation=TARGET_ATTENUATION)
    _DESIGN[]=design
    design
end

function design_json(design)
    Dict("kappa_A_per_V"=>design.kappa, "kappa_per_unit"=>design.kappa_per_unit,
         "phi2_rad"=>design.phi2, "phi2_deg"=>rad2deg(design.phi2),
         "target_attenuation"=>design.target_attenuation,
         "devices"=>[Dict("id"=>id,
             "self_Z2_ohm"=>[real(z),imag(z)], "self_Z2_angle_deg"=>rad2deg(angle(z)),
             "aggregate_Z2_ohm"=>[real(za),imag(za)],
             "aggregate_Z2_angle_deg"=>rad2deg(angle(za)),
             "aligned_delta_deg"=>rad2deg(angle(za)-design.phi2))
             for (id,z,za) in zip(design.ids,design.self_impedance,design.aggregate_impedance)])
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

# Status quo (C3): each customer keeps a single-phase 6 kVA PV inverter on its
# original phase with native BMOPFTools AS/NZS 4777.2 droops (phase-to-neutral
# sensing). BMOPFTools' solve_pf strips IBR s_max and thereby zeroes S_MAX-based
# droop bases (dev/notes/ibr_review_resolution.md, A10), so this uses the OPF
# builder with bus voltage bounds removed and maximum PV export: P is the
# tighter of availability and the Volt-watt cap, Q follows Volt-var.
function single_phase_network(base, availability, source_pu; curves=AUS_A)
    net=study_network(base, availability, source_pu)
    net["ibr"]=Dict{String,Any}()
    for b in values(net["bus"]); delete!(b,"v_min"); delete!(b,"v_max"); end
    net["control_profile"]=Dict{String,Any}("as4777_a"=>Dict{String,Any}(
        "volt_var"=>Dict{String,Any}("voltage_reference"=>"PN_PER_PHASE",
            "breakpoints"=>curves.vv_v, "q_limits"=>[minimum(curves.vv_q), maximum(curves.vv_q)],
            "q_unit"=>"VA_FRACTION", "q_ref"=>"VAR_MAX"),
        "volt_watt"=>Dict{String,Any}("voltage_reference"=>"PN_PER_PHASE",
            "breakpoints"=>curves.vw_v, "p_limits"=>sort(curves.vw_p),
            "p_unit"=>"VA_FRACTION",
            "p_ref"=>curves.vw_basis === :rated ? "S_MAX" : "P_AVAILABLE")))
    for (customer,ld) in base["load"]
        net["ibr"]["pv_"*customer]=Dict{String,Any}("bus"=>ld["bus"],
            "topology"=>"SINGLE_PHASE", "terminal_map"=>[first(ld["terminal_map"]),"n"],
            "s_max"=>[HARDWARE.s_max], "i_max"=>[HARDWARE.single_phase_i_max],
            "p_min"=>[0.0], "p_max"=>[5000.0*availability],
            "p_avail"=>5000.0*availability, "control_profile"=>"as4777_a",
            "prime_mover"=>"PV")
    end
    net
end

function solve_single_phase(net)
    ids=sort!(collect(keys(net["ibr"])))
    hook!(ctx)=begin
        obj=zero(QuadExpr)
        for id in ids
            rec=net["ibr"][id]
            vr,vi=PowerOptLab._dv(ctx, rec["bus"], rec["terminal_map"][1], "n")
            cr=PowerOptLab._opf_ibr_current(ctx,id,1)
            ci=PowerOptLab._opf_ibr_current(ctx,id,1;component=:imag)
            add_to_expression!(obj, -(vr*cr+vi*ci))
        end
        @objective(PowerOptLab._opf_model(ctx), Min, obj)
    end
    ctx=BMOPFTools.build_opf_model(net; per_unit=true, s_base=100000., add_objective=false,
        model_hook! = hook!, optimizer=Ipopt.Optimizer)
    PowerOptLab._set_solver_options!(PowerOptLab._opf_model(ctx), OPTIONS)
    BMOPFTools.optimize_opf!(ctx)
    outcome=PowerOptLab._solve_outcome(PowerOptLab._opf_model(ctx))
    (status=string(outcome.termination_status), publishable=SolveStatus(outcome).publishable,
     result=PowerOptLab._extract_result(ctx,outcome))
end

function single_phase_device_rows(net, result, mode, availability, source_pu)
    [begin
        rec=net["ibr"][id]; d=only(values(result["ibr"][id]))
        p=d["pg"]; q=d["qg"]; i=hypot(d["cri"],d["cii"])
        merge((mode=mode, availability=availability, source_pu=source_pu,
               device_id=id),
              NamedTuple{DEVICE_COLUMNS}((rec["p_avail"], p, q, rec["p_avail"]-p,
                  NaN, i, i, NaN, NaN, NaN, NaN, NaN, 0.0)))
    end for id in sort!(collect(keys(net["ibr"])))]
end

const DEVICE_COLUMNS = (:p_available_W, :p_poc_W, :q_poc_var,
    :poc_active_power_shortfall_W, :converter_loss_W, :maximum_converter_current_A,
    :maximum_grid_current_A, :negative_sequence_current_A, :negative_sequence_scale,
    :dc_availability_margin_W, :dc_ripple_voltage_V, :capacitor_current_A,
    :exact_smooth_current_residual_A)

function run(; modes=MODES, levels=[0.0,0.5,1.0], sources=[1.0,1.03],
             curves=AUS_A, hw=HARDWARE, output=joinpath(ROOT,"results"),
             label="main")
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    rows=NamedTuple[]; devices=NamedTuple[]; cases=NamedTuple[]
    mkpath(output)
    design=negative_sequence_design()
    writejson(joinpath(output,"negative_sequence_design.json"), design_json(design))
    writejson(joinpath(output,"environment.json"), Dict(
        "label"=>label,
        "julia"=>string(VERSION),"PowerOptLab"=>string(Base.pkgversion(PowerOptLab)),
        "BMOPFTools"=>string(Base.pkgversion(BMOPFTools)),
        "Ipopt"=>string(Base.pkgversion(Ipopt)),
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"..","..","Project.toml")))),
        "source_vuf_percent"=>1.0,"modes"=>modes,"availability"=>levels,
        "source_positive_sequence_pu"=>sources,
        "hardware"=>Dict(String(k)=>v for (k,v) in pairs(hw)),
        "filter_r_ohm"=>filter_r(hw), "filter_x_ohm"=>filter_x(hw),
        "curves"=>Dict(String(k)=>v for (k,v) in pairs(curves)),
        "sensing_reference"=>SENSING_REFERENCE))
    for source in sources, level in levels, mode in modes
        println("Solving ", (label,source,level,mode)); flush(stdout)
        if mode == "single_phase_native"
            net=single_phase_network(base,level,source; curves=curves)
            solved=solve_single_phase(net)
            dr=single_phase_device_rows(net,solved.result,mode,level,source)
            cr=customer_rows(base,solved.result["bus"],mode,level,source)
            accepted=solved.publishable &&
                all(isfinite(r.vmax_pn_V) && isfinite(r.vuf_percent) for r in cr)
            push!(cases,(mode=mode,availability=level,source_pu=source,
                source_vuf_percent=1.0,termination=solved.status,accepted=accepted,
                exact_residual_A=0.0,min_dc_availability_margin_W=NaN))
        else
            net=study_network(base,level,source)
            result=solve_controlled_inverter_fleet(net,
                fleet(net,mode,level; curves=curves, hw=hw, design=design);
                s_base=100000.,solver_options=OPTIONS)
            full=controlled_inverter_rows(result)
            dr=[merge((mode=mode, availability=level, source_pu=source,
                       device_id=r.device_id),
                      NamedTuple{DEVICE_COLUMNS}(Tuple(getproperty(r,k)
                          for k in DEVICE_COLUMNS))) for r in full]
            cr=customer_rows(base,result.network["bus"],mode,level,source)
            finite=all(isfinite(getproperty(r,k)) for r in cr for k in
                (:customer_voltage_V,:neutral_V,:vuf_percent,:vmax_pn_V))
            finite &= all(isfinite(getproperty(r,k)) for r in full for k in
                (:p_poc_W,:q_poc_var,:converter_loss_W,:maximum_converter_current_A,
                 :maximum_grid_current_A,:ripple_power_VA,:dc_ripple_voltage_V,
                 :capacitor_current_A,:power_scale,:current_scale))
            residual=maximum(r.exact_smooth_current_residual_A for r in full)
            margin=minimum(r.dc_availability_margin_W for r in full)
            accepted=solve_status(result).publishable && finite && isfinite(residual) &&
                residual<0.02 && margin > -1e-2 &&
                all(r.maximum_converter_current_A <= hw.i_max+1e-5 for r in full) &&
                all(r.dc_ripple_voltage_V <= hw.dv2_max+1e-4 for r in full)
            push!(cases,(mode=mode,availability=level,source_pu=source,
                source_vuf_percent=1.0,termination=result.termination_status,
                accepted=accepted,exact_residual_A=residual,
                min_dc_availability_margin_W=margin))
        end
        if last(cases).accepted
            append!(rows,cr); append!(devices,dr)
        end
        write_csv(joinpath(output,"cases.csv"),cases)
        write_csv(joinpath(output,"customers.csv"),rows)
        write_csv(joinpath(output,"devices.csv"),devices)
    end
    @assert all(r.accepted for r in cases) "Rejected cases recorded; omitted from scientific plots"
    cases
end

# Sensitivity campaigns at the stressed operating point (full PV, 1.03 pu
# source). They change one assumption each and write to results/sensitivity_*.
const SENSITIVITIES = [
    ("pedagogical_volt_watt", (curves=PEDAGOGICAL,)),
    ("large_dc_link_2mF", (hw=merge(HARDWARE, (c_dc=2e-3,)),)),
    ("small_dc_link_150uF", (hw=merge(HARDWARE, (c_dc=150e-6,)),)),
]
function run_sensitivities()
    for (name, override) in SENSITIVITIES
        run(; levels=[1.0], sources=[1.03], label=name,
            output=joinpath(ROOT,"results","sensitivity_"*name), override...)
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run()
    run_sensitivities()
end
