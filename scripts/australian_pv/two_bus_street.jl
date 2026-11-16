# Two-bus teaching case: the LV3_55bus source and 500 kVA transformer feed one
# 4-wire overhead line to a single "street" bus. At the street bus sit `homes`
# identical customers, each with a 6 kVA three-phase inverter (5 kW PV) and the
# study's 1.5 kW 60/25/15 % A/B/C load, plus neighbours' single-phase PV on
# phase C (unity PF, no volt-var). Balanced source: all unbalance is local.
#
# Declared, illustrative construction: one bus, one snapshot per point, local
# fixed point only. Not a generalised feeder result.
#
# usage: julia --project=. scripts/australian_pv/two_bus_street.jl output
include("single_law_sweep.jl")

const STREET_LAWS=("mean_vv_vw","worst_vv_vw","phase_vv_vw")

function two_bus_base(base;homes=8,length_m=250.,linecode="abc4x95_lv_oh_4w_bundled")
    net=deepcopy(base)
    template=first(values(base["line"]))
    net["bus"]=Dict(b=>base["bus"][b] for b in ("b804","b641"))
    net["bus"]["street"]=Dict("perfectly_grounded_terminals"=>String[],"terminal_names"=>["a","b","c","n"])
    line=deepcopy(template);line["bus_from"]="b641";line["bus_to"]="street";line["length"]=length_m
    line["linecode"]=linecode
    net["line"]=Dict("street_line"=>line);net["switch"]=Dict{String,Any}()
    net["shunt"]=Dict("street_grounding"=>Dict("B_1_1"=>0,"G_1_1"=>0.1,"bus"=>"street","terminal_map"=>["n"]))
    ld=first(values(base["load"]))
    net["load"]=Dict{String,Any}()
    for k in 1:homes
        rec=deepcopy(ld);rec["bus"]="street";net["load"]["home$(lpad(k,2,'0'))"]=rec
    end
    net
end

street_fleet(net,cond,label)=single_law_fleet(net,cond,label)

function run_two_bus_street(output;homes=8,length_m=250.,source_pu=1.0,
        neighbour_pv_W=(0.,250.,500.,750.,1000.,1250.,1500.),laws=STREET_LAWS)
    mkpath(joinpath(output,"points"))
    base=two_bus_base(parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"));homes,length_m)
    rows=NamedTuple[];cases=NamedTuple[]
    vmag(ts,p)=abs(phasor(ts[p])-phasor(ts["n"]))
    for pv in neighbour_pv_W, label in laws
        cond=(case_id="street_$(label)_$(round(Int,pv))",model=:energy,pv_fraction=1.,source_pu=source_pu,
            source_vuf=0.,source_negative_angle_rad=0.,phase_shares=(.60,.25,.15),load_scale=1.,operation=:energized)
        net=add_phase_c_pv!(voltage_aware_network(base,cond),pv)
        fleet=street_fleet(net,cond,label)
        # Damped retries for the fixed point near the volt-watt knee (253 V).
        res=nothing;t=0.;relax=1.
        for (relax,iters) in ((1.,30),(.5,120),(.25,300))
            t+=@elapsed res=solve_voltage_aware_inverter_fleet(net,fleet.spec;operating_modes=fleet.modes,
                s_base=VA_BASE,solver_options=ENGINEERING_OPTIONS,relaxation=relax,max_iterations=iters)
            res.solve.publishable && break
        end
        ok=res.solve.publishable
        println("neighbour_pv=$(pv) law=$(label) publishable=$(ok) reason=$(res.reason) iterations=$(res.iterations) t=$(round(t,digits=1))s")
        ok || continue
        point=engineering_point(res,cond;voltage_max=nothing)
        write_engineering_json(joinpath(output,"points",cond.case_id*".json"),point)
        m=engineering_record_metrics(point,fleet;mode="local")
        street=point["bus"]["street"];head=point["bus"]["b641"]
        id=first(sort!(collect(keys(fleet.spec.devices))))
        d=point["devices"][id]
        gv=[complex(z...) for z in d["grid_current"]]
        vv=[phasor(street[p])-phasor(street["n"]) for p in ("a","b","c")]
        sp=vv.*conj.(gv)
        for r in m.rows
            push!(rows,merge((neighbour_pv_W=pv,law=label),r))
        end
        push!(cases,(neighbour_pv_W=pv,law=label,accepted=m.accepted,
            va=vmag(street,"a"),vb=vmag(street,"b"),vc=vmag(street,"c"),
            head_va=vmag(head,"a"),head_vb=vmag(head,"b"),head_vc=vmag(head,"c"),
            p_home_a=real(sp[1]),p_home_b=real(sp[2]),p_home_c=real(sp[3]),
            q_home_a=imag(sp[1]),q_home_b=imag(sp[2]),q_home_c=imag(sp[3]),
            p_poc_W=sum(r.p_poc_W for r in m.rows),q_poc_var=sum(r.q_poc_var for r in m.rows),
            curtailment_W=sum(r.curtailment_W for r in m.rows),
            maximum_vuf_percent=maximum(r.vuf_percent for r in m.rows)))
        write_csv(joinpath(output,"devices.csv"),rows);write_csv(joinpath(output,"cases.csv"),cases)
    end
    write_engineering_json(joinpath(output,"environment.json"),Dict("laws"=>laws,"homes"=>homes,
        "line_length_m"=>length_m,"linecode"=>"abc4x95_lv_oh_4w_bundled","source_pu"=>source_pu,
        "neighbour_single_phase_pv_W_per_home_on_phase_c"=>neighbour_pv_W,"fixed_point_retries"=>"relaxation 1/0.5/0.25, max iterations 30/120/300","p_available_W_per_home"=>5000.,
        "load_per_home"=>"1.5 kW PF 0.95, 60/25/15 % A/B/C","source_vuf"=>0.,
        "git_revision"=>strip(read(`git -C $(normpath(joinpath(ROOT,"..",".."))) rev-parse HEAD`,String))))
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_two_bus_street(ARGS[1])
end
