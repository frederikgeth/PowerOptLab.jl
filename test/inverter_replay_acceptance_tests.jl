_replay_test_plant(id="pv";kwargs...)=AdvancedInverter(;
    id,bus="poc",phase_terminals=["a","b","c"],neutral=nothing,topology=:THREE_LEG,
    s_max=20000.,v_dc=750.,c_dc=.05,m_max=.96,n_samples=6,c_loss=.01,kwargs...)

function _replay_test_fleet(plants;available=0.)
    net=inv_grid3_bal()
    net["ibr"]=Dict(inv.id=>Dict{String,Any}("bus"=>"poc","terminal_map"=>["a","b","c"],
        "topology"=>"THREE_LEG","s_max"=>fill(inv.s_max/3,3),
        "p_min"=>zeros(3),"p_max"=>fill(inv.s_max/3,3),
        "q_min"=>fill(-inv.s_max/3,3),"q_max"=>fill(inv.s_max/3,3)) for inv in plants)
    devices=Dict(inv.id=>ControlledDevice(inv,SequenceController(AverageVoltageVoltVarWatt())) for inv in plants)
    requests=Dict(inv.id=>InverterControlRequest(p_available=available,p_rated=inv.s_max,q_scale=inv.s_max) for inv in plants)
    net,ControlledInverterFleetSpec(devices,requests)
end

@testset "Current replays: virtual iterations cannot publish an overdrawn final point" begin
    inv=_replay_test_plant(i_max=40.,r_filter=.05,x_filter=.15)
    controlled=ControlledDevice(inv,SequenceController(AverageVoltageVoltVarWatt()))
    request=InverterControlRequest(p_available=5000.,p_rated=5000.,q_scale=0.)
    low=ComplexF64[180.,180cis(-2pi/3),180cis(2pi/3)]
    r=solve_inverter_control_network_fixed_point(inv_grid3_bal(),controlled,request;
        initial_voltage=low,max_iterations=1,solver_options=("tol"=>1e-9,))
    @test r.solve.optimal
    @test !r.solve.publishable
    @test !r.converged
    @test !physical_feasibility(r).passed
    @test any(c.name=="replay_dc_availability" && !c.passed for c in physical_feasibility(r).checks)
    net,spec=_replay_test_fleet([inv];available=5000.)
    fleet=solve_controlled_inverter_fleet_network_fixed_point(net,spec;
        initial_voltage=Dict("pv"=>low),max_iterations=1,solver_options=("tol"=>1e-9,))
    @test fleet.solve.optimal
    @test !fleet.solve.publishable
    @test !fleet.converged
    @test !physical_feasibility(fleet)["pv"].passed
end

@testset "Current replays: sampled rail false success cannot publish" begin
    bad=_replay_test_plant(v_dc=550.)
    zero=ntuple(_->0.0im,3)
    sampled=PowerOptLab._solve_exact_current_plant(inv_grid3_bal(),bad,zero,
        ConverterCurrentTarget();rail_refinement=false)
    @test sampled.status.optimal && sampled.status.feasible
    @test sampled.status.termination_status in ("OPTIMAL","LOCALLY_SOLVED")
    @test !sampled.status.publishable
    @test !sampled.plant.physical_feasibility.passed
    @test sampled.plant.physical_feasibility.rail.upper_V < -35.
    @test !isempty(sampled.network["raw_solver_status"])
    refined=PowerOptLab._solve_exact_current_plant(inv_grid3_bal(),bad,zero,ConverterCurrentTarget())
    @test !refined.status.publishable
    # One bad device rejects the whole fleet while preserving each report.
    net,spec=_replay_test_fleet([bad,_replay_test_plant("good")])
    currents=Dict(id=>zero for id in keys(spec.devices))
    fleet=PowerOptLab._solve_exact_current_fleet_plant(net,spec,currents;rail_refinement=false)
    @test fleet.status.optimal && fleet.status.feasible
    @test !fleet.status.publishable
    @test !fleet.plants["pv"].physical_feasibility.passed
    @test fleet.plants["good"].physical_feasibility.passed
    @test !isempty(fleet.network["raw_solver_status"])
    @test !haskey(fleet.network["ibr"],"pv")
    refined_fleet=PowerOptLab._solve_exact_current_fleet_plant(net,spec,currents)
    @test !refined_fleet.status.publishable
end

@testset "Current replays: physical commands and source boundary" begin
    a=cis(2pi/3);current=Tuple(4 .* [1,a^2,a])
    inv=_replay_test_plant(r_filter=.05,x_filter=.15)
    for base in (1e4,1e6,1e8),target in (ConverterCurrentTarget(),GridCurrentTarget())
        r=PowerOptLab._solve_exact_current_plant(inv_grid3_bal(),inv,current,target;s_base=base,
            dc_source=DCSourceBounds(p_min=0.,p_max=5000.),solver_options=("tol"=>1e-9,))
        @test r.status.publishable
        @test r.plant.physical_feasibility.passed
        actual=target isa ConverterCurrentTarget ? r.converter_terminal.phase_current : r.grid_phase_current
        @test maximum(abs,actual.-current)<1e-5
        @test count(c->startswith(c.name,"replayed_current_"),r.plant.physical_feasibility.checks)==3
    end
    insufficient=PowerOptLab._solve_exact_current_plant(inv_grid3_bal(),inv,current,ConverterCurrentTarget();
        dc_source=DCSourceBounds(p_max=100.))
    @test !insufficient.status.publishable
    net,spec=_replay_test_fleet([inv];available=0.)
    currents=Dict("pv"=>current)
    # A feeder sensitivity probe is an imposed-current experiment, not a claim
    # that zero available PV can supply this active injection.
    probe=PowerOptLab._solve_exact_current_fleet_plant(net,spec,currents)
    @test probe.status.publishable
    @test probe.plants["pv"].p_dc>2000.
    controlled=PowerOptLab._solve_exact_current_fleet_plant(net,spec,currents;respect_availability=true)
    @test !controlled.status.publishable
    no_charge=PowerOptLab._solve_exact_current_plant(inv_grid3_bal(),inv,Tuple(-collect(current)),
        ConverterCurrentTarget();dc_source=DCSourceBounds(p_min=0.))
    @test !no_charge.status.publishable
end

@testset "Current replays: rejected sensitivity base remains unresolved" begin
    inv=_replay_test_plant(v_dc=550.)
    controlled=ControlledDevice(inv,SequenceController(AverageVoltageVoltVarWatt()))
    zero=ntuple(_->0.0im,3)
    result=inverter_control_network_voltage_sensitivity(inv_grid3_bal(),controlled,zero)
    @test !result.base_status.publishable
    @test all(isnan,result.sensitivity)
    net,spec=_replay_test_fleet([inv])
    fleet=controlled_inverter_fleet_network_voltage_sensitivity(net,spec,Dict("pv"=>zero))
    @test !fleet.base_status.publishable
    @test all(isnan,fleet.sensitivity)
end
