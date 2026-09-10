"""Illustrative workplace study: fixed assignments, unequal phases, real line losses.
The neutral is ideal/equipotential; this example does not study neutral displacement.
"""
module EVWorkplaceExample
using PowerOptLab
using BMOPFTools: parse_bmopf

function feeder(price)
    parse_bmopf("""
    {"bus":{
      "grid":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"]},
      "work":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"],
              "v_min":[200.0,200.0,200.0],"v_max":[250.0,250.0,250.0]}},
     "voltage_source":{"supply":{"bus":"grid","terminal_map":["a","b","c"],
        "v_magnitude":[245.0,215.0,230.0],"v_angle":[0.0,-2.0943951023931953,2.0943951023931953],
        "cost":[$price,$price,$price]}},
     "linecode":{"cable":{"R_series_1_1":0.2,"R_series_2_2":0.2,"R_series_3_3":0.2}},
     "line":{"service":{"bus_from":"grid","bus_to":"work","terminal_map_from":["a","b","c"],
        "terminal_map_to":["a","b","c"],"linecode":"cable","length":1.0}},
     "load":{
       "office_a":{"bus":"work","terminal_map":["a","n"],"configuration":"SINGLE_PHASE","p_nom":[500.0],"q_nom":[0.0]},
       "office_b":{"bus":"work","terminal_map":["b","n"],"configuration":"SINGLE_PHASE","p_nom":[1000.0],"q_nom":[0.0]},
       "office_c":{"bus":"work","terminal_map":["c","n"],"configuration":"SINGLE_PHASE","p_nom":[1500.0],"q_nom":[0.0]}}}
    """; from_string=true)
end

"""Four one-hour intervals. Targets are battery-energy gains (Wh).
`weak_current` changes equipment, keeping every vehicle request fixed.
`phase` changes the fixed installation connection, not an optimization decision.
"""
function study(; weak_current=16.0, visitor_gain=5500.0, phase="b", per_unit=true)
    outlet = EVSE(id="visitors",bus="work",phase_terminals=[phase],
        s_max=7e3,i_max=weak_current,p_charge_max=7e3)
    staff_outlet = EVSE(id="staff",bus="work",phase_terminals=["a"],s_max=7e3,i_max=32.0)
    car(id) = EV(id=id,energy_max=40e3,onboard_charge_max=7e3,eff_charge=0.9)
    sessions = [
        ChargingSession(id="morning",ev=car("visitor_a"),evse=outlet,
            available=[true,true,false,false],energy_init=10e3,departure_energy=10e3+visitor_gain),
        ChargingSession(id="afternoon",ev=car("visitor_b"),evse=outlet,
            available=[false,false,true,true],energy_init=10e3,departure_energy=10e3+visitor_gain),
        ChargingSession(id="staff",ev=car("staff_car"),evse=staff_outlet,
            available=fill(true,4),energy_init=10e3,departure_energy=20e3)]
    nets = [feeder(p) for p in (0.30,0.10,0.10,0.30)]
    result = solve_multiperiod_opf(nets,sessions;per_unit,
        solver_options=(tol=1e-9,bound_relax_factor=0.0))
    (; result, sessions)
end

# The bus phase-to-neutral voltage is <=250 V in this ideal-neutral fixture.
# Thus this is a necessary bound, regardless of the local NLP solver's status.
visitor_energy_bound(current=16.0) = 2.0 * 0.9 * min(7000.0,250.0*current)

function summary(case)
    r=case.result
    (status=r.termination_status, total_source_cost=r.objective,
     charging_w=Dict(id=>d.p_charge for (id,d) in r.dispatch),
     stored_energy_wh=Dict(id=>d.energy_wh for (id,d) in r.dispatch))
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    for (label,case) in (("16 A visitors",EVWorkplaceExample.study()),
            ("32 A upgrade",EVWorkplaceExample.study(weak_current=32.0)),
            ("8 kWh request at 16 A",EVWorkplaceExample.study(visitor_gain=8e3)),
            ("8 kWh request at 32 A",EVWorkplaceExample.study(weak_current=32.0,visitor_gain=8e3)))
        println((; label, EVWorkplaceExample.summary(case)...))
    end
end
