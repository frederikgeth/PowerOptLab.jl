_phase_test_balanced(x)=ComplexF64.(x.*[1,cis(-2pi/3),cis(2pi/3)])
_phase_test_plant(;kwargs...)=AdvancedInverter(;id="pv",bus="poc",topology=:FOUR_LEG,
    phase_terminals=["a","b","c"],neutral="n",s_max=6000.,i_max=10.,In_max=10.,
    v_dc=750.,c_dc=300e-6,m_max=.96,kwargs...)
const _PHASE_TEST_REQUEST=InverterControlRequest(p_available=5000.,p_rated=6000.,q_scale=6000.)

@testset "Composite magnitude: cancelling current sums stay finite" begin
    for staged in (false,true)
        target=staged ? PowerOptLab.build_opf_model(inv_grid3_bal();add_objective=false) : Model()
        m=staged ? PowerOptLab._opf_model(target) : target
        @variable(m,x[1:3]);@variable(m,y[1:3])
        for ib in (1.,4000.)
            expression=magnitude_expression(target,(sum(x),sum(y)),
                MagnitudeApproximation(1e-8;unit=:A,scale=10.);component_scale=ib)
            for radius in (1.,5.1,11.3),theta in (0.,.1,1.2,2.5)
                current=_phase_test_balanced(radius*cis(theta))
                current[3]=-(current[1]+current[2])
                point=Dict(v=>value for k in 1:3 for (v,value) in
                    ((x[k],real(current[k])/ib),(y[k],imag(current[k])/ib)))
                result=value(v->point[v],expression)
                @test isfinite(result)
                @test abs(result)<=1e-8
            end
        end
    end
end

function _phase_test_fleet(;available=5000.,mixed=false,plant_kwargs=(;))
    net=inv_grid3_floating_neutral()
    inv=_phase_test_plant(;power_limit_location=:poc,r_filter=.13225,x_filter=1.3225,
        p_loss_fixed=15.,a_loss=2.,c_loss=.3,dv2_max=7.5,i_cap_max=1.,plant_kwargs...)
    law=PiecewiseLinearLaw([180.,280.],[-.1,-.1];smoothing_epsilon=.01)
    devices=Dict{String,ControlledDevice}("pv"=>ControlledDevice(inv,PerPhaseVoltVarWatt(volt_var=law)))
    if mixed
        three=_phase_test_plant(id="three",topology=:THREE_LEG,neutral=nothing,
            power_limit_location=:poc,p_loss_fixed=15.,a_loss=2.,c_loss=.3)
        devices["three"]=ControlledDevice(three,SequenceController(AverageVoltageVoltVarWatt(volt_var=law);
            limiter=CommonScaleLimiter(pq_priority=:var)))
    end
    net["ibr"]=Dict(id=>Dict{String,Any}("bus"=>"poc","topology"=>c.device.topology===:FOUR_LEG ? "FOUR_LEG" : "THREE_LEG",
        "terminal_map"=>c.device.neutral===nothing ? ["a","b","c"] : ["a","b","c","n"],
        "s_max"=>fill(2000.,3),"p_min"=>zeros(3),"p_max"=>fill(2000.,3),
        "q_min"=>fill(-2000.,3),"q_max"=>fill(2000.,3)) for (id,c) in devices)
    requests=Dict(id=>InverterControlRequest(p_available=available,p_rated=6000.,q_scale=6000.) for id in keys(devices))
    net,ControlledInverterFleetSpec(devices,requests)
end

@testset "Four-leg feeder coupling: full LCL and neutral loss accounting" begin
    net,spec=_phase_test_fleet(plant_kwargs=(r_filter_neutral=.07,x_filter_neutral=.2,
        r_filter_grid=.03,x_filter_grid=.1,r_filter_grid_neutral=.04,x_filter_grid_neutral=.1,
        c_filter_mid=20e-6,r_filter_damping=.5,b_filter_shunt=1e-4))
    r=solve_voltage_aware_inverter_fleet(net,spec)
    @test r.solve.publishable
    point=r.replayed["pv"];inv=spec.devices["pv"].device
    ic,ig=point.state.converter_current,point.state.grid_current
    @test abs(sum(ic)-sum(ig))>.01
    losses=inv.r_filter*sum(abs2,ic)+inv.r_filter_neutral*abs2(sum(ic))+
        inv.r_filter_grid*sum(abs2,ig)+inv.r_filter_grid_neutral*abs2(sum(ig))+
        inv.r_filter_damping*sum(abs2,ic.-ig)
    @test point.p_conv-point.p_poc ≈ losses atol=.0001
    power=point.state.grid_voltage.*conj.(ig).+im*inv.b_filter_shunt.*abs2.(point.state.grid_voltage)
    @test imag.(power) ≈ fill(-200.,3) atol=.01
    @test maximum(real.(power))-minimum(real.(power))<.01
    @test all(c.passed for c in physical_feasibility(r)["pv"].checks)
    @test point.p_dc<=5000.01
end

@testset "Four-leg feeder coupling: neutral KCL and mixed controller families" begin
    for mixed in (false,true)
        net,spec=_phase_test_fleet(;mixed);original=deepcopy(net)
        r=solve_voltage_aware_inverter_fleet(net,spec)
        @test r.solve.publishable && r.converged
        @test all(report.passed for report in values(r.physical_reports))
        @test r.devices["pv"] isa PhaseInverterServiceAllocationResult
        @test isempty(r.network["ibr"])
        ph(bus,t)=complex(r.network["bus"][bus][t]["vr"],r.network["bus"][bus][t]["vi"])
        voltage=[ph("poc",t)-ph("poc","n") for t in ("a","b","c")]
        current=r.grid_currents["pv"]
        power=voltage.*conj.(collect(current))
        @test imag.(power) ≈ fill(-200.,3) atol=.01
        @test maximum(real.(power))-minimum(real.(power))<.01
        @test abs(sum(current))>.01
        @test abs(ph("poc","n"))>.1
        total=[sum(i[k] for i in values(r.grid_currents)) for k in 1:3]
        load=conj(6000/voltage[1])
        for (k,t) in enumerate(("a","b","c"))
            @test abs((ph("poc",t)-ph("grid",t))/.05-total[k]+(k==1 ? load : 0))<1e-5
        end
        @test abs(ph("poc","n")/.2+sum(total)-load)<1e-5
        rows=Dict(row.device_id=>row for row in voltage_aware_inverter_rows(r))
        @test rows["pv"].negative_sequence_scale===nothing
        @test rows["pv"].phase_reactive_service_scale ≈ 1. atol=1e-6
        @test rows["pv"].neutral_total_rms_current_A<=10.00001
        @test net==original
    end
end

@testset "Four-leg feeder coupling: energized, disconnected and domain failures" begin
    net,spec=_phase_test_fleet(available=0.)
    on=solve_voltage_aware_inverter_fleet(net,spec)
    @test on.solve.publishable
    @test on.replayed["pv"].p_poc < -15.
    @test on.replayed["pv"].p_dc<=.01
    off=solve_voltage_aware_inverter_fleet(net,spec;operating_modes=:disconnected)
    @test off.solve.publishable && off.reason===:DISCONNECTED_POWER_FLOW
    @test isempty(off.devices) && isempty(off.physical_reports)
    @test only(voltage_aware_inverter_rows(off)).neutral_total_rms_current_A==0.
    inv=spec.devices["pv"].device
    floor=ControlledInverterFleetSpec(Dict("pv"=>ControlledDevice(inv,PerPhaseVoltVarWatt(power_voltage_floor=300.))),spec.requests)
    failed=solve_voltage_aware_inverter_fleet(net,floor)
    @test failed.reason===:CONTROLLER_VOLTAGE_DOMAIN && !failed.solve.publishable
    @test_throws ArgumentError solve_voltage_aware_inverter_fleet(net,spec;reactive_mode=:protective_override)
end

@testset "Four-leg phase allocation: analytic neutral limit and shared source" begin
    inv=_phase_test_plant(In_max=2.)
    q=[0.,-600.,600.];v=_phase_test_balanced(230.)
    reference=nothing
    for base in (1e4,1e5,1e6)
        r=allocate_phase_inverter_services(inv,_PHASE_TEST_REQUEST;grid_voltage=v,reactive_request=q,s_base=base)
        @test r isa AbstractInverterServiceAllocationResult
        @test solve_status(r).publishable && physical_feasibility(r).passed
        @test r.plant_feasibility===r.policy_feasibility===:feasible
        @test r.reactive_scale ≈ 2*230/(sqrt(3)*600) atol=1e-6
        @test r.allocation.plant.p_poc ≈ 5000. atol=.01
        @test r.allocation.plant.p_dc<=5000.01
        current=collect(r.allocation.grid_phase_current)
        phase=v.*conj.(current)
        @test real.(phase) ≈ fill(5000/3,3) atol=.01
        @test imag.(phase) ≈ r.reactive_scale.*q atol=.01
        @test abs(sum(current)) ≈ 2. atol=1e-5
        @test maximum(abs,current)<10.
        @test r.current_command!==nothing
        @test [n for (n,_) in r.stages]==[:phase_reactive_tracking,:phase_active_export]
        @test all(solve_status(s).publishable for (_,s) in r.stages)
        reference===nothing || @test collect(r.current_command) ≈ collect(reference.current_command) atol=1e-5
        reference=r
    end
    # A phase VA rating is not inferred from total VA or from a missing current
    # rating: opposing Q can cancel in the declared total POC apparent power.
    free=allocate_phase_inverter_services(_phase_test_plant(i_max=nothing,In_max=nothing),_PHASE_TEST_REQUEST;
        grid_voltage=v,reactive_request=[0.,-2000.,2000.])
    @test solve_status(free).publishable
    @test free.allocation.plant.p_poc ≈ 5000. atol=.01
    @test maximum(free.allocation.plant.i_mag)>10.
end

@testset "Four-leg phase allocation: losses, watt bases and sensed requests" begin
    v=_phase_test_balanced(230.)
    inv=_phase_test_plant(r_filter=.13225,x_filter=1.3225,p_loss_fixed=15.,a_loss=2.,c_loss=.3,
        dv2_max=7.5,i_cap_max=1.)
    for available in (0.,5000.)
        request=InverterControlRequest(p_available=available,p_rated=6000.,q_scale=6000.)
        r=allocate_phase_inverter_services(inv,request;grid_voltage=v,reactive_request=fill(-600.,3))
        @test solve_status(r).publishable
        @test r.reactive_scale ≈ 1. atol=1e-6
        phase=v.*conj.(collect(r.allocation.grid_phase_current))
        @test imag.(phase) ≈ fill(-600.,3) atol=.01
        @test maximum(real.(phase))-minimum(real.(phase))<.01
        current=collect(r.allocation.converter_terminal.phase_current)
        exact_loss=15+2*(sum(abs,current)+abs(sum(current)))+.3*(sum(abs2,current)+abs2(sum(current)))
        @test r.allocation.plant.p_loss ≈ exact_loss atol=.00011
        @test r.allocation.plant.p_dc<=available+.01
        iszero(available) ? (@test r.allocation.plant.p_poc < -15.) : (@test r.allocation.plant.p_poc>4800.)
    end
    for (basis,expected) in ((:rated,3000.),(:available,2500.))
        r=allocate_phase_inverter_services(_phase_test_plant(),_PHASE_TEST_REQUEST;grid_voltage=v,
            reactive_request=zeros(3),active_fraction=(1.,.5,1.),volt_watt_basis=basis)
        @test solve_status(r).publishable
        @test r.allocation.plant.p_poc ≈ expected atol=.01
    end
    law=PiecewiseLinearLaw([220.,240.],[0.,-.5];smoothing_epsilon=.01)
    controlled=ControlledDevice(_phase_test_plant(),PerPhaseVoltVarWatt(volt_var=law))
    sensed=ComplexF64[241,230cis(-2pi/3),219cis(2pi/3)]
    r=evaluate_voltage_aware(controlled,InverterControlMeasurement(sensed),_PHASE_TEST_REQUEST;grid_voltage=v)
    @test solve_status(r).publishable
    @test r.phase_reactive_request_var==(-1000.,-500.,0.)
    @test imag.(v.*conj.(collect(r.allocation.grid_phase_current))) ≈ r.reactive_scale.*[-1000.,-500.,0.] atol=.01
    @test_throws ArgumentError evaluate_voltage_aware(controlled,InverterControlMeasurement(
        ComplexF64[.5,230cis(-2pi/3),230cis(2pi/3)]),_PHASE_TEST_REQUEST;grid_voltage=v)
    # A four-wire phase policy does not divide by positive-sequence voltage.
    # In this extreme synthetic U0-only boundary, the 10 A neutral, rather than
    # a missing U1 denominator, limits lossless active export to 230 V * 10 A.
    zero_sequence=fill(230.0+0im,3)
    r=evaluate_voltage_aware(ControlledDevice(_phase_test_plant(),PerPhaseVoltVarWatt()),
        InverterControlMeasurement(zero_sequence),_PHASE_TEST_REQUEST;grid_voltage=zero_sequence)
    @test solve_status(r).publishable
    @test r.allocation.plant.p_poc ≈ 2300. atol=.01
    @test r.allocation.plant.i_neutral ≈ 10. atol=1e-5
end

@testset "Four-leg phase allocation: explicit failures and domains" begin
    inv=_phase_test_plant();v=_phase_test_balanced(230.)
    r=allocate_phase_inverter_services(inv,_PHASE_TEST_REQUEST;grid_voltage=v,reactive_request=zeros(3),
        solver_options=("max_iter"=>0,),diagnose_plant=false)
    @test !solve_status(r).publishable && r.current_command===nothing
    @test r.policy_feasibility===:unresolved
    @test_throws ArgumentError allocate_phase_inverter_services(_phase_test_plant(topology=:THREE_LEG),
        _PHASE_TEST_REQUEST;grid_voltage=v,reactive_request=zeros(3))
    @test_throws DimensionMismatch allocate_phase_inverter_services(inv,_PHASE_TEST_REQUEST;
        grid_voltage=v,reactive_request=[0.])
    @test_throws ArgumentError allocate_phase_inverter_services(inv,_PHASE_TEST_REQUEST;
        grid_voltage=v,reactive_request=zeros(3),reactive_mode=:protective_override)
    @test_throws ArgumentError validate_device(ControlledDevice(_phase_test_plant(neutral=nothing),PerPhaseVoltVarWatt()))
    @test_throws ArgumentError PerPhaseVoltVarWatt(power_voltage_floor=0.)
    @test_throws ArgumentError PerPhaseVoltVarWatt(volt_watt_basis=:wrong)
    @test_throws ArgumentError solve_controlled_inverter(inv_grid3_bal(),
        ControlledDevice(inv,PerPhaseVoltVarWatt()),_PHASE_TEST_REQUEST)
end
