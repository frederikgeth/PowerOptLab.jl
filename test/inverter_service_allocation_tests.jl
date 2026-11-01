_service_plant(;kwargs...)=AdvancedInverter(;id="services",bus="poc",topology=:THREE_LEG,
    phase_terminals=["a","b","c"],neutral=nothing,s_max=6000.,i_max=10.,v_dc=900.,c_dc=.001,kwargs...)
_service_balanced(v)=v.*[1,cis(-2pi/3),cis(2pi/3)]
const _SERVICE_REQUEST=InverterControlRequest(p_available=6000.,p_rated=6000.,q_scale=6000.)

@testset "Service allocation: Q then P, analytic capability and numerical bases" begin
    for base in (1e4,1e5,1e6),q in (-3000.,3000.)
        r=allocate_inverter_services(_service_plant(),_SERVICE_REQUEST;grid_voltage=_service_balanced(230.),
            reactive_request=q,s_base=base)
        @test solve_status(r).publishable
        @test r.policy_feasibility===:feasible
        @test r.plant_feasibility===:feasible
        @test physical_feasibility(r).passed
        @test r.current_command!==nothing
        @test first.(r.stages)==[:reactive_tracking,:active_export]
        @test r.allocation.plant.q_poc ≈ q atol=.01
        @test r.allocation.plant.p_poc ≈ sqrt(6000^2-q^2) atol=.02
        @test abs(r.allocation.plant.q_poc-r.stages[1].second.plant.q_poc)<.01
    end
    saturated=allocate_inverter_services(_service_plant(),_SERVICE_REQUEST;
        grid_voltage=_service_balanced(230.),reactive_request=8000.)
    @test solve_status(saturated).publishable
    @test saturated.allocation.plant.q_poc ≈ 6000. atol=.02
    @test abs(saturated.allocation.plant.p_poc)<1.
end

@testset "Service allocation: preserve an active Australian reactive request without duplicate bounds" begin
    inv=_service_plant(v_dc=750.,c_dc=300e-6,m_max=.96,r_filter=.13225,x_filter=1.3225,
        p_loss_fixed=15.,a_loss=2.,c_loss=.3,dv2_max=7.5,i_cap_max=1.,power_limit_location=:poc)
    v=_service_balanced(253.) .+ 2.53cis(pi/3).*[1,cis(2pi/3),cis(-2pi/3)]
    r=allocate_inverter_services(inv,InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.);
        grid_voltage=v,reactive_request=-2854.888024653047)
    @test solve_status(r).publishable
    @test physical_feasibility(r).passed
    @test r.allocation.plant.q_poc ≈ -2854.888024653047 atol=.01
    @test r.allocation.plant.p_poc > 4800.
end

@testset "Service allocation: voltage rails, explicit override and unresolved strict policy" begin
    inv=_service_plant(v_dc=620.,c_dc=300e-6,m_max=.96,r_filter=.13225,x_filter=1.3225,
        p_loss_fixed=20.,c_loss=.01,dv2_max=7.5,i_cap_max=1.)
    request=InverterControlRequest(p_available=5000.,p_rated=4800.,q_scale=6000.)
    v=_service_balanced(253.)
    preserved=allocate_inverter_services(inv,request;grid_voltage=v,reactive_request=-6000.)
    @test solve_status(preserved).publishable
    @test preserved.allocation.plant.q_poc ≈ -6000. atol=.01
    @test 1700. < preserved.allocation.plant.p_poc < 1800.
    @test physical_feasibility(preserved).rail.lower_V>=-1e-5
    strict=allocate_inverter_services(inv,request;grid_voltage=v,reactive_request=-3600.)
    @test !solve_status(strict).publishable
    @test strict.current_command===nothing
    @test strict.policy_feasibility===:unresolved
    @test strict.plant_feasibility===:feasible
    @test strict.plant_diagnostic!==nothing
    @test solve_status(strict.plant_diagnostic).publishable
    override=allocate_inverter_services(inv,request;grid_voltage=v,reactive_request=-3600.,
        reactive_mode=:protective_override)
    @test solve_status(override).publishable
    @test override.allocation.plant.q_poc < -5000.
    @test override.reactive_mode===:protective_override
    @test all(physical_feasibility(s).passed for (_,s) in override.stages)
end

@testset "Service allocation: available/rated ceiling and energized zero PV" begin
    inv=_service_plant(p_loss_fixed=100.,r_filter=.05,x_filter=.15,b_filter_shunt=1e-4)
    request=InverterControlRequest(p_available=3000.,p_rated=6000.,q_scale=6000.)
    rated=allocate_inverter_services(inv,request;grid_voltage=_service_balanced(230.),
        reactive_request=-600.,active_fraction=.5)
    available=allocate_inverter_services(inv,request;grid_voltage=_service_balanced(230.),
        reactive_request=-600.,active_fraction=.5,volt_watt_basis=:available)
    @test solve_status(rated).publishable && solve_status(available).publishable
    @test rated.allocation.plant.p_poc > 2800.
    @test 1400. < available.allocation.plant.p_poc < 1500.
    @test rated.allocation.plant.q_poc ≈ -600. atol=.01
    night=allocate_inverter_services(inv,InverterControlRequest(p_available=0.,p_rated=6000.,q_scale=6000.);
        grid_voltage=_service_balanced(230.),reactive_request=-600.)
    @test solve_status(night).publishable
    @test night.allocation.plant.p_dc ≈ 0. atol=.01
    @test night.allocation.plant.p_poc < -100.
    @test night.allocation.plant.q_poc ≈ -600. atol=.01
end

@testset "Service allocation: negative support retains both higher-priority services" begin
    r=allocate_inverter_services(_service_plant(),_SERVICE_REQUEST;grid_voltage=_service_balanced(230.),
        reactive_request=-1000.,negative_sequence_request=.5im)
    @test solve_status(r).publishable
    @test first.(r.stages)==[:reactive_tracking,:active_export,:negative_sequence_support]
    @test r.allocation.plant.q_poc ≈ r.stages[2].second.plant.q_poc atol=.01
    @test r.allocation.plant.p_poc ≈ r.stages[2].second.plant.p_poc atol=.01
    @test 0. <= r.negative_scale <= 1.
    @test abs(r.allocation.converter_terminal.current_sequence[3]-.5im*r.negative_scale)<1e-5
    @test abs(r.allocation.plant.q_poc+1000.)<.01
end

@testset "Voltage-aware law adapter: matched sensing and explicit boundaries" begin
    vv=PiecewiseLinearLaw([207.,220.,240.,258.],[.44,0.,0.,-.60];formulation=LocalC2Formulation(.01))
    vw=PiecewiseLinearLaw([253.,260.],[1.,.2];formulation=LocalC2Formulation(.01))
    measurement=InverterControlMeasurement(_service_balanced(250.))
    for policy in (AverageVoltageVoltVarWatt(volt_var=vv,volt_watt=vw),
            WorstPhaseVoltVarWatt(volt_var=vv,volt_watt=vw),
            PositiveSequenceVoltVarWatt(volt_var=vv,volt_watt=vw))
        c=ControlledDevice(_service_plant(),SequenceController(policy;
            limiter=CommonScaleLimiter(pq_priority=:var),current_target=:converter))
        r=evaluate_voltage_aware(c,measurement,_SERVICE_REQUEST;grid_voltage=_service_balanced(230.))
        @test solve_status(r).publishable
        @test r.reactive_request_var ≈ -2000. atol=1e-6
        @test r.allocation.plant.q_poc ≈ -2000. atol=.01
        @test collect(r.allocation.plant.poc_phase_voltage) ≈ _service_balanced(230.) atol=1e-8
    end
    bad=ControlledDevice(_service_plant(),SequenceController(AverageVoltageVoltVarWatt()))
    @test_throws ArgumentError evaluate_voltage_aware(bad,measurement,_SERVICE_REQUEST;
        grid_voltage=_service_balanced(230.))
    @test_throws ArgumentError allocate_inverter_services(_service_plant(),_SERVICE_REQUEST;
        grid_voltage=_service_balanced(230.),reactive_request=NaN)
    failed=allocate_inverter_services(_service_plant(),_SERVICE_REQUEST;
        grid_voltage=_service_balanced(230.),reactive_request=3000.,diagnose_plant=false,
        solver_options=("max_iter"=>0,))
    @test !solve_status(failed).publishable
    @test failed.current_command===nothing
    @test failed.plant_feasibility===:unresolved
end

@testset "Voltage-aware allocation: LCL grid command and measured admittance support" begin
    inv=_service_plant(r_filter=.05,x_filter=.15,r_filter_grid=.03,x_filter_grid=.1,
        c_filter_mid=20e-6,r_filter_damping=.5,b_filter_shunt=1e-4)
    v=_service_balanced(230.) .+ (3+1im)
    r=allocate_inverter_services(inv,_SERVICE_REQUEST;grid_voltage=v,reactive_request=-600.,target=:grid)
    @test solve_status(r).publishable
    @test physical_feasibility(r).passed
    @test abs(r.allocation.converter_terminal.current_sequence[1])<1e-5
    grid0=abs(sum(r.allocation.grid_phase_current)/3)
    @test grid0 > .001
    @test r.current_command==r.allocation.grid_phase_current
    @test r.allocation.plant.q_poc ≈ -600. atol=.01
    gain=PiecewiseLinearLaw([0.,.1],[.1,.1];formulation=LocalC2Formulation(.001))
    positive=AverageVoltageVoltVarWatt()
    controller=SequenceController(positive;unbalance=NegativeSequenceAdmittanceDroop(gain;
        ripple_blend=0.,impedance_angle=.3),limiter=SequencePriorityLimiter(pq_priority=:var),current_target=:grid)
    u=_service_balanced(230.) .+ 2 .* [1,cis(2pi/3),cis(-2pi/3)]
    support=evaluate_voltage_aware(ControlledDevice(_service_plant(),controller),
        InverterControlMeasurement(u),_SERVICE_REQUEST;grid_voltage=u)
    @test solve_status(support).publishable
    @test support.negative_sequence_request_A ≈ -.2cis(-.3) atol=1e-10
    @test support.target===:grid
    unsupported=SequenceController(positive;unbalance=NegativeSequenceAdmittanceDroop(gain;ripple_blend=.5),
        limiter=CommonScaleLimiter(pq_priority=:var))
    @test_throws ArgumentError evaluate_voltage_aware(ControlledDevice(_service_plant(),unsupported),
        InverterControlMeasurement(u),_SERVICE_REQUEST;grid_voltage=u)
end
