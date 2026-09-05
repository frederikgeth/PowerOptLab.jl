# Small runnable case for scripts/validate_doe_from_dsse.jl.
# Replace this file with an adapter for a real DSSE snapshot and meter export.

using PowerOptLab
using BMOPFTools: parse_bmopf, solve_pf, add_statcom!, augment_case

function _demo_net(; include_loads=true, include_ibrs=true)
    bounds = (include_loads || include_ibrs) ? ",\"v_min\":[216.0],\"v_max\":[245.0]" : ""
    loads = include_loads ? """
    ,"load":{
      "d1":{"bus":"bus1","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[1200.0],"q_nom":[150.0]},
      "d2":{"bus":"bus2","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[800.0],"q_nom":[100.0]}}""" : ""
    ibrs = include_ibrs ? """
    ,"ibr":{
      "pv1":{"bus":"bus1","terminal_map":["1","n"],"topology":"SINGLE_PHASE","prime_mover":"PV","s_max":[10000.0],"p_min":[0.0],"p_max":[10000.0],"q_min":[0.0],"q_max":[0.0]},
      "pv2":{"bus":"bus2","terminal_map":["1","n"],"topology":"SINGLE_PHASE","prime_mover":"PV","s_max":[10000.0],"p_min":[0.0],"p_max":[10000.0],"q_min":[0.0],"q_max":[0.0]}}""" : ""
    return parse_bmopf("""
    {"bus":{
      "src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
      "bus1":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]$bounds},
      "bus2":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]$bounds}},
     "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.30,"X_series_1_1":0.30}},
     "line":{
       "l1":{"bus_from":"src","bus_to":"bus1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
       "l2":{"bus_from":"bus1","bus_to":"bus2","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}
     $loads $ibrs}
    """; from_string=true)
end

function _demo_materialize_estimate(estimate, template)
    # This fixture has one single-phase load per bus and zero pre-DOE DER P/Q.
    # Consequently the estimated net injection identifies load demand directly.
    # A real feeder with active DER needs a separate injection-disaggregation model.
    for load in values(template["load"])
        bus, terminal = load["bus"], first(load["terminal_map"])
        for (kind, field) in ((:pinj, "p_nom"), (:qinj, "q_nom"))
            row = only(r for r in estimate.residuals if
                r.kind == kind && r.bus == bus && r.terminal == terminal)
            isfinite(row.estimated) || error("missing estimated load injection")
            load[field] = [-row.estimated]
        end
    end
    net, _ = augment_case(template)
    return net
end

function doe_validation_case()
    physics = _demo_net(include_loads=false, include_ibrs=false)
    operational = _demo_net()
    truth = deepcopy(operational)
    for inv in values(truth["ibr"])
        inv["p_max"] = [0.0] # measurement-generating pre-DOE PV dispatch
    end
    pf = solve_pf(truth; per_unit=false)
    measurements = Measurement[]
    for (bus, p, q) in (("bus1", -1200.0, -150.0), ("bus2", -800.0, -100.0))
        vm = pf["bus"][bus]["1"]["vm"]
        append!(measurements, [Measurement(kind=:vmag, bus=bus, value=vm, sigma=0.5),
                               Measurement(kind=:pinj, bus=bus, value=p, sigma=20.0),
                               Measurement(kind=:qinj, bus=bus, value=q, sigma=20.0)])
    end
    statcom = deepcopy(operational)
    add_statcom!(statcom, "bus2"; s_max=5000.0)
    stat = statcom["ibr"]["statcom_bus2"]
    stat["p_min"] = [0.0]; stat["p_max"] = [0.0]
    # The runner must freeze the Q selected by DOE issuance in the separate PF.
    stat["q_min"] = [-5000.0]; stat["q_max"] = [5000.0]
    statcom, _ = augment_case(statcom)
    return (physics_net=physics, operational_net=operational, truth_net=truth,
            with_statcom_net=statcom, measurements=measurements,
            materialize_estimate=_demo_materialize_estimate,
            connection_points=[ConnectionPoint(id="pv1", bus="bus1", ibr_id="pv1", export_max=10e3),
                               ConnectionPoint(id="pv2", bus="bus2", ibr_id="pv2", export_max=10e3)],
            doe_keywords=(security=:bound_point, control_policy=IssuePlusLocalLaws()))
end
