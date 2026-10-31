_coupled_test_plant(id="pv";kwargs...)=AdvancedInverter(;id,bus="poc",phase_terminals=["a","b","c"],
    neutral=nothing,topology=:THREE_LEG,s_max=6000.,i_max=10.,v_dc=750.,c_dc=300e-6,m_max=.96,
    r_filter=.13225,x_filter=1.3225,p_loss_fixed=15.,a_loss=2.,c_loss=.3,
    dv2_max=7.5,i_cap_max=1.,power_limit_location=:poc,kwargs...)

function _coupled_test_case(;plants=[_coupled_test_plant()],available=5000.,net=inv_grid3_bal(),floor=1.)
    vv=PiecewiseLinearLaw([180.,280.],[-.1,-.1];formulation=LocalC2Formulation(.01))
    controller=SequenceController(AverageVoltageVoltVarWatt(volt_var=vv);
        limiter=CommonScaleLimiter(pq_priority=:var),power_voltage_floor=floor)
    net=deepcopy(net)
    net["ibr"]=Dict(inv.id=>Dict{String,Any}("bus"=>inv.bus,
        "terminal_map"=>inv.neutral===nothing ? inv.phase_terminals : [inv.phase_terminals;inv.neutral],
        "topology"=>inv.neutral===nothing ? "THREE_LEG" : "FOUR_LEG",
        "s_max"=>fill(inv.s_max/3,3),"p_min"=>zeros(3),"p_max"=>fill(inv.s_max/3,3),
        "q_min"=>fill(-inv.s_max/3,3),"q_max"=>fill(inv.s_max/3,3)) for inv in plants)
    devices=Dict(inv.id=>ControlledDevice(inv,controller) for inv in plants)
    requests=Dict(inv.id=>InverterControlRequest(p_available=available,p_rated=6000.,q_scale=6000.) for inv in plants)
    net,ControlledInverterFleetSpec(devices,requests)
end

@testset "Voltage-aware feeder coupling: physical closure and numerical bases" begin
    net,spec=_coupled_test_case()
    original=deepcopy(net)
    reference=nothing
    for base in (1e4,1e5,1e6)
        r=solve_voltage_aware_inverter_fleet(net,spec;s_base=base)
        @test solve_status(r).publishable
        @test r.converged && r.reason===:CONVERGED
        @test 1<=r.iterations<15
        @test !isempty(r.network_equation_checks) && all(c.passed for c in r.network_equation_checks)
        @test all(p.passed for p in values(physical_feasibility(r)))
        @test r.history[end].voltage_residual_inf_V<=1e-7
        @test r.history[end].current_residual_inf_A<=1e-6
        @test !haskey(r.network["ibr"],"pv")
        p=r.replayed["pv"];inv=spec.devices["pv"].device
        @test p.p_dc ≈ 5000. atol=.01
        @test p.q_poc ≈ -600. atol=.01
        @test p.p_conv-p.p_poc ≈ inv.r_filter*sum(abs2,p.state.converter_current) atol=1e-6
        @test p.p_dc-p.p_conv ≈ inv.p_loss_fixed+inv.a_loss*sum(abs,p.state.converter_current)+
            inv.c_loss*sum(abs2,p.state.converter_current) atol=1e-6
        @test maximum(abs,p.state.grid_current.-collect(r.grid_currents["pv"]))<1e-12
        # Independent Ohm/KCL calculation for this diagonal 0.05-ohm feeder.
        phasor(b,t)=complex(r.network["bus"][b][t]["vr"],r.network["bus"][b][t]["vi"])
        @test maximum(abs((phasor("poc",t)-phasor("grid",t))/.05-r.grid_currents["pv"][k])
            for (k,t) in enumerate(["a","b","c"]))<1e-5
        @test any(startswith(c.name,"equilibrium_current") for c in physical_feasibility(r)["pv"].checks)
        rows=voltage_aware_inverter_rows(r)
        @test length(rows)==1 && rows[1].accepted && rows[1].model_present
        @test rows[1].p_poc_W==p.p_poc
        if reference!==nothing
            @test abs(p.p_poc-reference.replayed["pv"].p_poc)<.01
            @test maximum(abs,p.state.grid_voltage-reference.replayed["pv"].state.grid_voltage)<1e-5
        end
        reference=r
    end
    @test net==original
    warm=solve_voltage_aware_inverter_fleet(net,spec;initial_currents=reference.grid_currents)
    @test warm.solve.publishable && warm.iterations==1
    @test abs(warm.replayed["pv"].p_poc-reference.replayed["pv"].p_poc)<.01
    damped=solve_voltage_aware_inverter_fleet(net,spec;relaxation=.5)
    @test damped.solve.publishable
    @test abs(damped.replayed["pv"].p_poc-reference.replayed["pv"].p_poc)<.01
end

@testset "Voltage-aware feeder coupling: unselected native device stays active" begin
    net,spec=_coupled_test_case()
    net["ibr"]["native"]=deepcopy(net["ibr"]["pv"])
    native=net["ibr"]["native"]
    native["p_min"]=native["p_max"]=fill(100.,3)
    native["q_min"]=native["q_max"]=zeros(3)
    original=deepcopy(net)
    r=solve_voltage_aware_inverter_fleet(net,spec)
    @test r.solve.publishable
    @test Set(keys(r.network["ibr"]))==Set(["native"])
    @test sum(t["pg"] for t in values(r.network["ibr"]["native"])) ≈ 300. atol=.001
    phasor(b,t)=complex(r.network["bus"][b][t]["vr"],r.network["bus"][b][t]["vi"])
    @test maximum(abs((phasor("poc",t)-phasor("grid",t))/.05-r.grid_currents["pv"][k]-
        conj(100/phasor("poc",t))) for (k,t) in enumerate(["a","b","c"]))<1e-5
    @test net==original
end

@testset "Voltage-aware feeder coupling: energized zero PV and open contactor" begin
    net,spec=_coupled_test_case(;available=0.)
    on=solve_voltage_aware_inverter_fleet(net,spec)
    @test on.solve.publishable
    @test on.replayed["pv"].p_poc < -15.
    @test on.replayed["pv"].q_poc ≈ -600. atol=.01
    @test on.replayed["pv"].p_dc ≈ 0. atol=.01
    off=solve_voltage_aware_inverter_fleet(net,spec;operating_modes=:disconnected)
    @test off.solve.publishable && off.converged
    @test off.reason===:DISCONNECTED_POWER_FLOW && off.iterations==0
    @test isempty(off.devices) && isempty(off.replayed) && isempty(physical_feasibility(off))
    row=only(voltage_aware_inverter_rows(off))
    @test !row.model_present && row.p_poc_W==0. && row.q_poc_var==0.
    @test row.p_dc_W===nothing && row.physical_feasible===nothing
    no_device=deepcopy(net);empty!(no_device["ibr"])
    baseline=PowerOptLab.BMOPFTools.solve_pf(no_device;per_unit=true,
        solver_options=("tol"=>1e-9,"bound_relax_factor"=>0.,))
    @test maximum(abs(complex(v["vr"],v["vi"])-complex(baseline["bus"][b][t]["vr"],baseline["bus"][b][t]["vi"]))
        for (b,ts) in off.network["bus"] for (t,v) in ts)<1e-6
    mixed_net,mixed_spec=_coupled_test_case(plants=[_coupled_test_plant(),_coupled_test_plant("off")])
    mixed=solve_voltage_aware_inverter_fleet(mixed_net,mixed_spec;
        operating_modes=Dict("pv"=>:energized,"off"=>:disconnected))
    @test mixed.solve.publishable
    @test Set(keys(mixed.devices))==Set(["pv"])
    @test Set(keys(mixed.physical_reports))==Set(["pv"])
    @test mixed.grid_currents["off"]==ntuple(_->0.0im,3)
end

@testset "Voltage-aware feeder coupling: LCL arm and neutral sensing" begin
    inv=_coupled_test_plant(neutral="n",r_filter=.05,x_filter=.15,r_filter_grid=.03,x_filter_grid=.1,
        c_filter_mid=20e-6,r_filter_damping=.5,b_filter_shunt=1e-4)
    net,spec=_coupled_test_case(plants=[inv],net=inv_grid3_floating_neutral())
    r=solve_voltage_aware_inverter_fleet(net,spec)
    @test r.solve.publishable
    @test r.replayed["pv"].physical_feasibility.passed
    @test abs(sum(r.replayed["pv"].state.converter_current))<1e-5
    @test abs(sum(r.replayed["pv"].state.grid_current))>.001
    @test abs(complex(r.network["bus"]["poc"]["n"]["vr"],r.network["bus"]["poc"]["n"]["vi"]))>.01
    @test r.replayed["pv"].q_poc ≈ -600. atol=.01
    # Disconnection removes the POC shunt too, rather than merely fixing Ig=0.
    off=solve_voltage_aware_inverter_fleet(net,spec;operating_modes=:disconnected)
    baseline=deepcopy(net);empty!(baseline["ibr"])
    pf=PowerOptLab.BMOPFTools.solve_pf(baseline;per_unit=true,solver_options=("tol"=>1e-9,))
    @test off.solve.publishable
    @test maximum(abs(complex(v["vr"],v["vi"])-complex(pf["bus"][b][t]["vr"],pf["bus"][b][t]["vi"]))
        for (b,ts) in off.network["bus"] for (t,v) in ts)<1e-6
end

@testset "Voltage-aware feeder coupling: failures stay explicit" begin
    net,spec=_coupled_test_case()
    limited=solve_voltage_aware_inverter_fleet(net,spec;max_iterations=1)
    @test !limited.solve.publishable && !limited.converged
    @test limited.reason===:ITERATION_LIMIT
    @test limited.solve.optimal # last virtual network solve, not controller acceptance
    @test !only(voltage_aware_inverter_rows(limited)).accepted
    floor_net,floor_spec=_coupled_test_case(floor=300.)
    domain=solve_voltage_aware_inverter_fleet(floor_net,floor_spec)
    @test domain.reason===:CONTROLLER_VOLTAGE_DOMAIN
    @test !domain.solve.publishable
    @test domain.physical_reports["pv"].passed # independent plant witness at the initial point
    @test domain.history[end].failed_allocation_ids==["pv"]
    bad_net,bad_spec=_coupled_test_case(plants=[_coupled_test_plant(v_dc=550.,r_filter=0.,x_filter=0.)])
    bad=solve_voltage_aware_inverter_fleet(bad_net,bad_spec)
    @test bad.reason===:ALLOCATION_FAILED && !bad.solve.publishable
    @test bad.history[end].failed_allocation_ids==["pv"]
    @test bad.devices["pv"].policy_feasibility===:unresolved
    @test_throws ArgumentError solve_voltage_aware_inverter_fleet(net,spec;relaxation=0.)
    @test_throws ArgumentError solve_voltage_aware_inverter_fleet(net,spec;operating_modes=:idle)
    @test_throws ArgumentError solve_voltage_aware_inverter_fleet(net,spec;operating_modes=Dict("other"=>:energized))
    @test_throws ArgumentError solve_voltage_aware_inverter_fleet(net,spec;initial_currents=Dict("pv"=>(1.0+0im,0.0im,0.0im)))
end
