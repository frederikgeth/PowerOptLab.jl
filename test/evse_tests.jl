@testset "EV electrical isolation and conductor capability" begin
    ev = EVDevice(id="off", bus="poc", phase_terminals=["a","b","c"],
        i_max=16.0, p_charge_max=11e3, energy_max=40e3, energy_init=20e3,
        available=[false], departure_energy=20e3)
    for pu in (true, false)
        r = solve_multiperiod_opf([inv_grid3()], [ev]; per_unit=pu)
        @test solve_status(r).publishable
        d = r.dispatch["off"]
        @test maximum(d.current_magnitude_a[1]) < 1e-6
        @test d.energy_wh ≈ [20e3,20e3] atol=0.1
        @test maximum(d.diagnostics.disconnected_current_a) < 1e-6
    end
    # The former counterexample imposed a nonzero conductor current while the
    # disconnected aggregate P,Q remained zero. It must now be infeasible.
    h = Ref{Any}()
    multi = build_multi_context([inv_grid3()]; per_unit=false,
        hook_factory=t->ctx->(h[]=stamp_device!(ctx,ev;period=t)))
    @constraint(multi.model, h[].currents[1][1] == 10.0)
    foreach(PowerOptLab.enforce_kcl!,multi.contexts)
    optimize!(multi.model)
    @test termination_status(multi.model) ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)

    # On an unbalanced feeder, equal power is not equal phase current.
    on = EVDevice(id="on", bus="poc", phase_terminals=["a","b","c"],
        i_max=16.0, neutral_current_max=10.0, s_max=11e3,
        p_charge_max=11e3, energy_max=40e3, energy_init=0.0,
        available=[true], departure_energy=6000.0, energy_final=6000.0)
    for pu in (true,false)
        r = solve_multiperiod_opf([inv_grid3_unbal()], [on]; per_unit=pu)
        @test solve_status(r).publishable
        d = r.dispatch["on"]
        @test all(d.phase_power_w[1] .≈ -2000.0)
        @test maximum(d.current_magnitude_a[1]) < 16.01
        @test d.neutral_current_a[1] < 10.01
        @test maximum(d.diagnostics.phase_sharing_residual_w) < 0.01
    end
end

@testset "EV modes: negative-price artificial dissipation" begin
    function full_ev(operation; tolerance=1e-8)
        EVDevice(id="full",bus="bus1",p_charge_max=10e3,p_discharge_max=10e3,
            energy_max=40e3,energy_init=40e3,eff_charge=0.9,eff_discharge=0.9,
            available=[true],departure_energy=40e3,energy_final=40e3,
            operation=operation,complementarity_tolerance=tolerance)
    end
    nets = [single_bus_net(src_cost=-0.1,pload=10e3)]
    loose = solve_multiperiod_opf(nets,[full_ev(:independent)])
    @test solve_status(loose).publishable
    @test loose.dispatch["full"].p_charge[1] ≈ 10e3 atol=0.1
    @test loose.dispatch["full"].p_discharge[1] ≈ 8.1e3 atol=0.1
    # With a*b <= tau and b=.81a, the analytic maximum a=sqrt(tau/.81).
    for pu in (true,false)
        tau = 1e-6
        r = solve_multiperiod_opf(nets,[full_ev(:relaxed;tolerance=tau)];
            per_unit=pu,solver_options=(tol=1e-10,bound_relax_factor=0.0))
        @test solve_status(r).publishable
        d=r.dispatch["full"]
        @test d.p_charge[1] ≈ 10e3*sqrt(tau/0.81) atol=0.2
        @test d.diagnostics.complementarity_product[1] <= tau+1e-9
        @test maximum(abs.(d.diagnostics.energy_balance_wh)) < 0.01
    end
end

@testset "EVSE sessions: intersection, permissions, assignment and outage" begin
    v=EV(id="car",energy_max=40e3,onboard_charge_max=7e3,
        onboard_discharge_max=7e3,eff_charge=0.9)
    e=EVSE(id="outlet",bus="poc",s_max=22e3,i_max=32.0)
    s=ChargingSession(id="visit",ev=v,evse=e,available=[false,true,true,false],
        energy_init=10e3,departure_energy=16.3e3)
    @test validate_device(s,fill(inv_grid(),4)) === nothing
    r=solve_multiperiod_opf(fill(inv_grid(),4),[s];time_grid=TimeGrid([1.0,0.5,0.5,1.0]))
    @test solve_status(r).publishable
    d=r.dispatch["visit"]
    @test d.vehicle_id == "car" && d.evse_id == "outlet"
    @test d.energy_wh[4] ≈ 16.3e3 atol=0.1
    @test sum(d.p_charge[2:3])*0.5 ≈ 7e3 atol=0.1
    @test maximum(abs.(d.p_discharge)) < 1e-8
    @test maximum(d.current_magnitude_a[2]) <= 32.01
    @test maximum(abs.(d.diagnostics.energy_balance_wh)) < 0.01
    @test solve_diagnostics(r).device_diagnostics["visit"] == d.diagnostics
    @test r.snapshots[2]["custom_injection"]["p"] ≈ d.p_net[2]

    # Two different vehicles can use one outlet successively, but not overlap;
    # an outage still occupies the physical outlet.
    v2=EV(id="car2",energy_max=40e3,onboard_charge_max=7e3)
    first=ChargingSession(id="first",ev=v,evse=e,available=[true,false],
        energy_init=10e3,departure_energy=10e3)
    second=ChargingSession(id="second",ev=v2,evse=e,available=[false,true],
        energy_init=20e3,departure_energy=20e3)
    @test solve_status(solve_multiperiod_opf(fill(inv_grid(),2),[first,second])).publishable
    overlap=ChargingSession(id="overlap",ev=v2,evse=e,available=[true,false],
        energy_init=20e3,departure_energy=20e3)
    @test_throws ArgumentError solve_multiperiod_opf(fill(inv_grid(),2),[first,overlap])
    reset=ChargingSession(id="reset",ev=v,evse=e,available=[false,true],
        energy_init=0.0,departure_energy=0.0)
    @test_throws ArgumentError solve_multiperiod_opf(fill(inv_grid(),2),[first,reset])
    eoff=EVSE(id="outage",bus="poc",s_max=7e3,i_max=32.0,available=[true,false])
    outage=ChargingSession(id="outage",ev=v,evse=eoff,available=[true,true],
        energy_init=0.0,departure_energy=1000.0)
    out=solve_multiperiod_opf(fill(inv_grid(),2),[outage])
    @test solve_status(out).publishable
    @test out.dispatch["outage"].p_charge[2] ≈ 0.0 atol=1e-5
    @test out.dispatch["outage"].occupied[2]
    @test out.dispatch["outage"].energy_wh[2] ≈ out.dispatch["outage"].energy_wh[3] atol=0.01

    # Raw invalid inputs may not be hidden by disabled permissions or min().
    bad=EV(id="bad",energy_max=40e3,onboard_charge_max=7e3,onboard_discharge_max=-1.0)
    @test_throws ArgumentError validate_device(ChargingSession(id="bad",ev=bad,evse=e,
        available=[true],energy_init=0.0,departure_energy=0.0),[inv_grid()])
    @test_throws ArgumentError validate_device(ChargingSession(id="bad",ev=v,evse=e,
        available=[true,false,true],energy_init=0.0,departure_energy=0.0),fill(inv_grid(),3))
    @test_throws ArgumentError validate_device(ChargingSession(id="bad",ev=v,evse=e,
        available=[true,true],departure_period=1,energy_init=0.0,departure_energy=0.0),fill(inv_grid(),2))
    @test_throws ArgumentError validate_device(EVDevice(id="bad",bus="poc",
        phase_terminals=["a","b","c"],p_charge_max=7e3,energy_max=40e3,
        energy_init=0.0,available=[true],departure_energy=0.0),[inv_grid3()])
end

include(joinpath(@__DIR__,"..","examples","ev_charging.jl"))
@testset "EV pedagogical examples and taper error budget" begin
    ratings=EVChargingExamples.ratings()
    @test solve_status(ratings).publishable
    @test ratings.objective ≈ 1.4 atol=1e-4  # 7 kWh at 0.20 currency/kWh
    @test ratings.dispatch["visit"].energy_wh[end] ≈ 16.3e3 atol=0.1
    cycle=EVChargingExamples.arbitrage()
    @test solve_status(cycle).publishable
    @test cycle.dispatch["cycle"].p_charge[1] ≈ 10e3 atol=0.5
    @test cycle.dispatch["cycle"].p_discharge[2] ≈ 8.1e3 atol=0.5
    @test cycle.dispatch["cycle"].energy_wh[end] ≈ 20e3 atol=0.1
    depleted=EVChargingExamples.arbitrage(replenish=false)
    @test solve_status(depleted).publishable
    @test depleted.dispatch["cycle"].energy_wh[end] < 1.0
    @test solve_status(EVChargingExamples.successive_visits()).publishable
    for f in (:auto,LocalC2Formulation(0.01),SoftplusFormulation(0.01))
        taper=EVChargingExamples.taper(f)
        @test solve_status(taper).publishable
        d=taper.dispatch["taper"]
        @test maximum(d.diagnostics.acceptance_violation_w) < 0.1
        @test d.p_charge[1] <= 7000/1.35 + 0.1
        if f === :auto
            @test d.p_charge[1] ≈ 7000/1.35 atol=0.1
        end
    end
    fine=EVChargingExamples.taper(dt_h=0.25,periods=4)
    @test solve_status(fine).publishable
    # Endpoint-conservative discretization approaches continuous max charging
    # z(t)=1-0.5exp(-0.35t) from below as the time grid is refined.
    coarse=EVChargingExamples.taper()
    @test fine.dispatch["taper"].energy_wh[end] > coarse.dispatch["taper"].energy_wh[end]
    @test fine.dispatch["taper"].energy_wh[end] < 40e3*(1-0.5exp(-0.35))
end

@testset "EV equipment ratings bind and permissions intersect" begin
    # A 10 A unity-PF load on a stiff 230 V source is limited to 2300 W.
    for (smax,imax,expected) in ((7000.0,10.0,2300.0),(2000.0,32.0,2000.0))
        e=EVSE(id="limited",bus="poc",s_max=smax,i_max=imax)
        v=EV(id="car",energy_max=40e3,onboard_charge_max=7e3)
        s=ChargingSession(id="limited",ev=v,evse=e,available=[true],
            energy_init=0.0,departure_energy=0.0)
        r=solve_multiperiod_opf([EVChargingExamples.grid(-0.1)],[s])
        @test solve_status(r).publishable
        @test r.dispatch["limited"].p_charge[1] ≈ expected atol=0.1
        @test maximum(r.dispatch["limited"].diagnostics.current_limit_violation_a) < 0.01
        @test maximum(r.dispatch["limited"].diagnostics.apparent_power_violation_va) < 0.1
    end
    # Session permission, the vehicle, and the EVSE independently gate export.
    for (vp,ep,permit,expected) in ((3000.,2000.,false,0.),(0.,2000.,true,0.),
            (3000.,0.,true,0.),(3000.,2000.,true,2000.))
        v=EV(id="car",energy_max=40e3,onboard_charge_max=7e3,onboard_discharge_max=vp)
        e=EVSE(id="outlet",bus="poc",s_max=7e3,i_max=32.0,p_discharge_max=ep)
        s=ChargingSession(id="v2g",ev=v,evse=e,available=[true],energy_init=20e3,
            departure_energy=0.0,allow_v2g=permit)
        r=solve_multiperiod_opf([EVChargingExamples.grid(0.2)],[s])
        @test solve_status(r).publishable
        @test r.dispatch["v2g"].p_net[1] ≈ expected atol=0.1
    end
    # P/Q must share a kVA budget; the box alone would allow 2 kW + 1.2 kvar.
    ev=EVDevice(id="pq",bus="poc",s_max=2000.0,i_max=32.0,
        p_charge_max=7000.0,q_min=1200.0,q_max=1200.0,
        energy_max=40e3,energy_init=0.0,available=[true],departure_energy=0.0)
    r=solve_multiperiod_opf([EVChargingExamples.grid(-0.1)],[ev])
    @test solve_status(r).publishable
    @test r.dispatch["pq"].p_charge[1] ≈ 1600.0 atol=0.1
end

@testset "EV acceptance domains and explicit formulation contract" begin
    vcurve(xs,ys;unit=:W) = EV(id="car",energy_max=40e3,onboard_charge_max=7e3,
        charge_acceptance=PWLFunction(xs,ys;input_unit=:unitless,output_unit=unit))
    outlet=EVSE(id="outlet",bus="poc",s_max=7e3,i_max=32.0)
    session(v;form=:auto)=ChargingSession(id="visit",ev=v,evse=outlet,
        available=[true],energy_init=10e3,departure_energy=10e3,acceptance_formulation=form)
    @test_throws ArgumentError validate_device(session(vcurve([0.,1.],[0.,7e3])),[inv_grid()])
    @test_throws ArgumentError validate_device(session(vcurve([0.,0.9],[7e3,0.])),[inv_grid()])
    @test_throws ArgumentError validate_device(session(vcurve([0.,1.],[7e3,0.];unit=:kW)),[inv_grid()])
    nonconcave=vcurve([0.,0.5,1.],[7e3,2e3,0.])
    @test_throws PowerOptLab.UnsupportedFormulation validate_device(session(nonconcave),[inv_grid()])
    @test validate_device(session(nonconcave;form=LocalC2Formulation(0.01)),[inv_grid()]) === nothing
    @test_throws ArgumentError validate_device(session(nonconcave;form=ExactPWLGraph()),[inv_grid()])
    bare=build_multi_context([inv_grid()];optimizer=nothing)
    @test length(bare.contexts) == 1
end
