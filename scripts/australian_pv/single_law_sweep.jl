# Single-law sweep: the same LV3_55bus feeder, every inverter running the same
# volt-var/volt-watt law, swept over the upstream source voltage.
#
# Study changes relative to the matched benchmark (declared, not measured):
#   - every one of the 12 customers has a 6 kVA three-phase inverter with 5 kW PV;
#   - each customer bus also carries a fixed single-phase PV export on phase C
#     (unity PF, no volt-var), standing in for neighbouring single-phase rooftops;
#   - source voltage swept; 2 % negative-sequence source voltage rotated -120 deg so
#     that phase C is the high phase (the benchmark uses 1 % at 0 deg).
# Local fixed point only (no OPF). Illustrative, one feeder, one snapshot per point.
#
# usage: julia --project=. scripts/australian_pv/single_law_sweep.jl output [phase_c_pv_W] [source_pu,...]
include("run_engineering_benchmark.jl")

const SWEEP_LAWS=("mean_vv_vw","worst_vv_vw","phase_vv_vw")

function single_law_fleet(net,case,label;limits=ENGINEERING_OPERATING_LIMITS,p_available=5000.)
    four=four_leg_fleet(net,case.pv_fraction)
    devices=Dict{String,ControlledDevice}();requests=Dict{String,InverterControlRequest}()
    modes=Dict{String,Symbol}();labels=Dict{String,String}()
    for id in sort!(collect(keys(net["ibr"])))
        labels[id]=label
        base=label=="phase_vv_vw" ? four.devices[id].device : plant(id,net["ibr"][id]["bus"];settings=TERMINAL_CONTROL_SETTINGS)
        fields=NamedTuple{fieldnames(AdvancedInverter)}(Tuple(getfield(base,k) for k in fieldnames(AdvancedInverter)))
        inv=AdvancedInverter(;merge(fields,(pwm_strategy=:CENTERED,f_sw=10000.,dc_energy_model=:nonlinear))...)
        law=label=="phase_vv_vw" ? four.devices[id].controller :
            controller(label;design=(kappa=0.,phi2=0.),settings=TERMINAL_CONTROL_SETTINGS)
        devices[id]=ControlledDevice(engineering_plant(inv;limits),law)
        net["ibr"][id]["topology"]=label=="phase_vv_vw" ? "FOUR_LEG" : "THREE_LEG"
        net["ibr"][id]["terminal_map"]=label=="phase_vv_vw" ? ["a","b","c","n"] : ["a","b","c"]
        requests[id]=InverterControlRequest(p_available=p_available,p_rated=HARDWARE.s_max,q_scale=HARDWARE.s_max)
        modes[id]=:energized
    end
    (spec=ControlledInverterFleetSpec(devices,requests),modes=modes,labels=labels)
end

function add_phase_c_pv!(net,watts)
    for id in collect(keys(net["load"]))
        endswith(id,"_c") || continue
        rec=deepcopy(net["load"][id]);rec["p_nom"]=[-watts];rec["q_nom"]=[0.]
        net["load"][replace(id,r"_c$"=>"_pv1ph_c")]=rec
    end
    net
end

function run_single_law_sweep(output;phase_c_pv=4000.,sources=(1.00,1.01,1.02,1.03,1.04,1.05),source_vuf=.02,source_angle=-2pi/3)
    mkpath(joinpath(output,"points"))
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    rows=NamedTuple[];cases=NamedTuple[]
    for src in sources, label in SWEEP_LAWS
        cond=(case_id="sweep_$(label)_$(round(Int,1000src))",model=:energy,pv_fraction=1.,source_pu=src,
            source_vuf=source_vuf,source_negative_angle_rad=source_angle,phase_shares=(.60,.25,.15),load_scale=1.,operation=:energized)
        net=add_phase_c_pv!(voltage_aware_network(base,cond),phase_c_pv)
        fleet=single_law_fleet(net,cond,label)
        t=@elapsed res=solve_voltage_aware_inverter_fleet(net,fleet.spec;operating_modes=fleet.modes,
            s_base=VA_BASE,solver_options=ENGINEERING_OPTIONS)
        ok=res.solve.publishable
        println("source=$(src) law=$(label) publishable=$(ok) reason=$(res.reason) t=$(round(t,digits=1))s")
        ok || (push!(cases,(source_pu=src,law=label,accepted=false,elapsed_s=t));continue)
        point=engineering_point(res,cond;voltage_max=nothing)
        write_engineering_json(joinpath(output,"points",cond.case_id*".json"),point)
        m=engineering_record_metrics(point,fleet;mode="local")
        for r in m.rows
            push!(rows,merge((source_pu=src,law=label),r))
        end
        push!(cases,(source_pu=src,law=label,accepted=m.accepted,elapsed_s=t,
            p_poc_W=sum(r.p_poc_W for r in m.rows),q_poc_var=sum(r.q_poc_var for r in m.rows),
            curtailment_W=sum(r.curtailment_W for r in m.rows),
            maximum_lv_voltage_V=engineering_lv_voltage_max(point,net),
            maximum_sensed_V=maximum(r.vmax_sensed_V for r in m.rows),
            maximum_vuf_percent=maximum(r.vuf_percent for r in m.rows)))
        write_csv(joinpath(output,"devices.csv"),rows);write_csv(joinpath(output,"cases.csv"),cases)
    end
    write_engineering_json(joinpath(output,"environment.json"),Dict("laws"=>SWEEP_LAWS,"sources_pu"=>sources,
        "phase_c_single_phase_pv_W_per_bus"=>phase_c_pv,"p_available_W_per_site"=>5000.,
        "load_phase_shares"=>(.60,.25,.15),"source_vuf"=>source_vuf,"source_negative_angle_rad"=>source_angle,"control_settings"=>TERMINAL_CONTROL_SETTINGS,
        "git_revision"=>strip(read(`git -C $(normpath(joinpath(ROOT,"..",".."))) rev-parse HEAD`,String))))
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    pv=length(ARGS)>=2 ? parse(Float64,ARGS[2]) : 4000.
    src=length(ARGS)>=3 ? Tuple(parse.(Float64,split(ARGS[3],","))) : (1.00,1.01,1.02,1.03,1.04,1.05)
    run_single_law_sweep(ARGS[1];phase_c_pv=pv,sources=src)
end
