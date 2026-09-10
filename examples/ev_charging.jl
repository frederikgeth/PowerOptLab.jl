"""Small, analytic AC charging examples. Prices below are currency/kWh."""
module EVChargingExamples
using PowerOptLab
using BMOPFTools: parse_bmopf

# Ideal 230 V point of supply; no line losses, so analytical expectations are
# exact for the selected device model. Real feeders require a separate study.
grid(price=0.0) = parse_bmopf("""
{"bus":{"poc":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"source":{"bus":"poc","terminal_map":["1"],
 "v_magnitude":[230.0],"v_angle":[0.0],"cost":[$price]}}}
""";from_string=true)

function ratings()
    car=EV(id="car",energy_max=40e3,onboard_charge_max=7e3,eff_charge=0.9)
    outlet=EVSE(id="outlet",bus="poc",s_max=22e3,i_max=32.0)
    visit=ChargingSession(id="visit",ev=car,evse=outlet,
        available=[false,true,true,false],energy_init=10e3,departure_energy=16.3e3)
    solve_multiperiod_opf([grid(0.2) for _ in 1:4],[visit];
        time_grid=TimeGrid([1.0,0.5,0.5,1.0]))
end

function negative_price(operation=:relaxed; tolerance=1e-6)
    ev=EVDevice(id="full",bus="poc",p_charge_max=10e3,p_discharge_max=10e3,
        energy_max=40e3,energy_init=40e3,eff_charge=0.9,eff_discharge=0.9,
        available=[true],departure_energy=40e3,energy_final=40e3,
        operation=operation,complementarity_tolerance=tolerance)
    solve_multiperiod_opf([grid(-0.1)],[ev];
        solver_options=(tol=1e-10,bound_relax_factor=0.0))
end

function arbitrage(; replenish=true)
    ev=EVDevice(id="cycle",bus="poc",p_charge_max=10e3,p_discharge_max=10e3,
        energy_max=40e3,energy_init=20e3,eff_charge=0.9,eff_discharge=0.9,
        available=[true,true],departure_energy=replenish ? 20e3 : 0.0,
        energy_final=replenish ? 20e3 : nothing)
    solve_multiperiod_opf([grid(p) for p in (0.05,0.25)],[ev];
        solver_options=(tol=1e-10,bound_relax_factor=0.0))
end

function taper(formulation=:auto; dt_h=1.0, periods=1)
    # Illustrative envelope, not fitted vehicle data. Domain is PE energy
    # fraction; the falling branch is 14000*(1-z) W.
    curve=PWLFunction([0.0,0.5,1.0],[7000.0,7000.0,0.0];
        input_unit=:unitless,output_unit=:W)
    car=EV(id="car",energy_max=40e3,onboard_charge_max=7e3,charge_acceptance=curve)
    outlet=EVSE(id="outlet",bus="poc",s_max=7e3,i_max=32.0)
    visit=ChargingSession(id="taper",ev=car,evse=outlet,available=fill(true,periods),
        energy_init=20e3,departure_energy=20e3,acceptance_formulation=formulation)
    solve_multiperiod_opf([grid(-0.1) for _ in 1:periods],[visit];dt_h=dt_h,
        solver_options=(tol=1e-10,bound_relax_factor=0.0))
end

function successive_visits()
    outlet=EVSE(id="shared",bus="poc",s_max=7e3,i_max=32.0)
    first=ChargingSession(id="morning",ev=EV(id="a",energy_max=40e3,onboard_charge_max=7e3),
        evse=outlet,available=[true,false],energy_init=0.0,departure_energy=3e3)
    second=ChargingSession(id="afternoon",ev=EV(id="b",energy_max=40e3,onboard_charge_max=3e3),
        evse=outlet,available=[false,true],energy_init=0.0,departure_energy=3e3)
    solve_multiperiod_opf([grid(0.1),grid(0.2)],[first,second])
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    for (name,run) in (("ratings",EVChargingExamples.ratings),
            ("arbitrage",EVChargingExamples.arbitrage),
            ("taper",EVChargingExamples.taper),
            ("successive visits",EVChargingExamples.successive_visits))
        r=run()
        println((example=name,status=r.termination_status,objective=r.objective))
        for (id,d) in r.dispatch
            println((id=id,energy_wh=d.energy_wh,p_charge=d.p_charge,p_discharge=d.p_discharge))
        end
    end
end
