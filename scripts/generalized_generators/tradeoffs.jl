# Reproducible steady-state trade-offs on a loaded, unbalanced four-wire feeder.
# Run: julia --project=. scripts/generalized_generators/tradeoffs.jl
using PowerOptLab, BMOPFTools, JuMP, Ipopt, LinearAlgebra, Printf

function generator_tradeoff_network()
    parse_bmopf("""
    {"bus":{
      "grid":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"]},
      "poc":{"terminal_names":["a","b","c","n"],"v_min":[180.0,180.0,180.0],"v_max":[280.0,280.0,280.0]}},
     "voltage_source":{"grid":{"bus":"grid","terminal_map":["a","b","c"],
       "v_magnitude":[240.0,225.0,232.0],"v_angle":[0.03,-2.12,2.10]}},
     "linecode":{"feeder":{"R_series_1_1":0.1,"R_series_2_2":0.1,"R_series_3_3":0.1,"R_series_4_4":0.2,
       "X_series_1_1":0.3,"X_series_2_2":0.3,"X_series_3_3":0.3,"X_series_4_4":0.1}},
     "line":{"feeder":{"bus_from":"grid","bus_to":"poc","terminal_map_from":["a","b","c","n"],
       "terminal_map_to":["a","b","c","n"],"linecode":"feeder","length":1.0}},
     "load":{
       "a":{"bus":"poc","terminal_map":["a","n"],"configuration":"SINGLE_PHASE","p_nom":[1500.0],"q_nom":[300.0]},
       "b":{"bus":"poc","terminal_map":["b","n"],"configuration":"SINGLE_PHASE","p_nom":[3500.0],"q_nom":[600.0]},
       "c":{"bus":"poc","terminal_map":["c","n"],"configuration":"SINGLE_PHASE","p_nom":[800.0],"q_nom":[100.0]}}}
    """;from_string=true)
end

function generator_tradeoff_study()
    ports=[("a","n"),("b","n"),("c","n")]
    emf=230 .* cis.([0,-2pi/3,2pi/3])
    Z=Matrix(Diagonal(fill(0.3+0.6im,3)))+(0.1+0.1im)*ones(3,3)
    cap=GeneratorCapability(i_max=60.0,terminal_i_max=[60.0,60.0,60.0,40.0])
    common=GeneratorVoltageLaw(:common_magnitude;magnitude_min=180.0,magnitude_max=280.0)
    cases=[
        ("a: phase PQ",GeneratorVoltageLaw(),GeneratorControl(p=[600.0,300.0,900.0],q=[100.0,100.0,100.0]),cap),
        ("b: fixed phasor",GeneratorVoltageLaw(:fixed_phasor;phasor=emf),GeneratorControl(),cap),
        ("c: fixed excitation",GeneratorVoltageLaw(:fixed_magnitudes;phasor=emf),GeneratorControl(p=1800.0),cap),
        ("d: common magnitude",common,GeneratorControl(p=1800.0,q=300.0),cap),
        ("e: phase magnitudes",GeneratorVoltageLaw(:phase_magnitudes;magnitude_min=180.0,magnitude_max=280.0),
            GeneratorControl(p=1800.0,q=[100.0,100.0,100.0]),cap),
        ("terminal PV",common,GeneratorControl(p=1800.0,voltage_target=230.0),cap),
        ("f: sequence dispatch",GeneratorVoltageLaw(),GeneratorControl(p=1800.0,q=300.0),
            GeneratorCapability(i_max=60.0,terminal_i_max=[60.0,60.0,60.0,40.0],
                sequence=GeneratorSequenceLimits(current_max=[2.0,40.0,4.0])))
    ]
    devices=Pair{String,AbstractDevice}[]
    for (label,law,control,capability) in cases
        push!(devices,label=>GeneralizedGenerator(id="device",bus="poc",connections=ports,
            impedance=Z,voltage=law,control=control,capability=capability))
    end
    for (label,ground) in (("source: open earth",nothing),("source: finite earth",0.8+0.1im),("source: ideal earth",0.0))
        push!(devices,label=>SourceGenerator(id="device",bus="poc",grounding=ground,
            impedance=ComplexF64[0.3+0.6im,0.3+0.6im,0.3+0.6im,0.2],
            voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=emf),capability=cap))
    end
    # Change only the electrical primitive in this pair of PQ cases: internal
    # excitation is unconstrained, so terminal dispatch/voltage should agree.
    push!(devices,"a: ideal series"=>GeneralizedGenerator(id="device",bus="poc",connections=ports,
        control=GeneratorControl(p=[600.0,300.0,900.0],q=[100.0,100.0,100.0]),capability=cap))
    options=Pair{String,Any}["tol"=>1e-9,"constr_viol_tol"=>1e-9,"bound_relax_factor"=>0.0,"max_iter"=>1500]
    rows=NamedTuple[]
    for (label,d) in devices
        result=solve_generator_opf(generator_tradeoff_network(),[d];per_unit=true,s_base=1e4,
            objective=:loss,solver_options=options)
        result.solve.publishable || error("tutorial case '$label' failed: $(result.solve)")
        g=result.devices["device"]
        push!(rows,(case=label,p_kw=g.p/1000,q_kvar=g.q/1000,
            v1=abs(g.voltage_sequence[2]),vuf_percent=100abs(g.voltage_sequence[3])/abs(g.voltage_sequence[2]),
            i_phase=maximum(abs.(g.port_current)),i_neutral=abs(g.terminal_current[end]),
            i_earth=abs(g.earth_current),loss_w=g.series_loss+g.ground_loss,
            power_error_w=abs(g.power_balance_error),result=g))
    end
    rows
end

function print_generator_tradeoffs(rows)
    @printf("%-23s %7s %7s %7s %7s %7s %7s %7s %8s\n",
        "Case","P kW","Q kvar","V1 V","VUF %","Iph A","In A","Ig A","Loss W")
    for r in rows
        @printf("%-23s %7.3f %7.3f %7.2f %7.3f %7.2f %7.2f %7.2f %8.2f\n",
            r.case,r.p_kw,r.q_kvar,r.v1,r.vuf_percent,r.i_phase,r.i_neutral,r.i_earth,r.loss_w)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    print_generator_tradeoffs(generator_tradeoff_study())
end
