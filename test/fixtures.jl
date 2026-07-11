# Shared test fixtures: tiny single-phase resistive feeders with a grounded
# neutral, whose solutions are easy to reason about.

using BMOPFTools: parse_bmopf

# source ──line── bus1 , unity-PF load at bus1. `src_cost` sets the slack import
# price used by the OPF objective; `pload` the load (W).
single_bus_net(; src_cost=0.0, pload=100000.0) = parse_bmopf("""
{"bus":{
    "sourcebus":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus1":     {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],
                 "v_min":[900.0],"v_max":[1100.0]}},
 "voltage_source":{"vs":{"bus":"sourcebus","terminal_map":["1"],
     "v_magnitude":[1000.0],"v_angle":[0.0],"cost":[$src_cost]}},
 "linecode":{"lc":{"R_series_1_1":0.1}},
 "line":{"l1":{"bus_from":"sourcebus","bus_to":"bus1",
     "terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}},
 "load":{"ld1":{"bus":"bus1","terminal_map":["1","n"],
     "configuration":"SINGLE_PHASE","p_nom":[$pload],"q_nom":[0.0]}}}
"""; from_string=true)

# source ──l1── bus1 ──l2── bus2 , unity-PF loads at bus1 and bus2. Physics net
# for state estimation (loads optional; pass load=false for a load-free net).
two_bus_net(; load=true) = parse_bmopf("""
{"bus":{
    "src": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus1":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus2":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],
     "v_magnitude":[1000.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.5}},
 "line":{
    "l1":{"bus_from":"src","bus_to":"bus1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
    "l2":{"bus_from":"bus1","bus_to":"bus2","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}
 $(load ? """,
 "load":{
    "d1":{"bus":"bus1","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[20000.0],"q_nom":[0.0]},
    "d2":{"bus":"bus2","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[20000.0],"q_nom":[0.0]}}""" : "")
}
"""; from_string=true)

# Stiff single-phase grid for advanced-inverter tests: slack at "grid", short line
# to the inverter POC bus, so the inverter's own limits (not the network) bind.
inv_grid(; vmax=250.0) = parse_bmopf("""
{"bus":{
    "grid":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "poc": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],"v_min":[200.0],"v_max":[$vmax]}},
 "voltage_source":{"vs":{"bus":"grid","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.05}},
 "line":{"l1":{"bus_from":"grid","bus_to":"poc","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

# Balanced three-phase stiff grid (for grid-forming tests).
inv_grid3() = parse_bmopf("""
{"bus":{
    "grid":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"]},
    "poc": {"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"],"v_min":[200.0,200.0,200.0],"v_max":[250.0,250.0,250.0]}},
 "voltage_source":{"vs":{"bus":"grid","terminal_map":["a","b","c"],"v_magnitude":[230.0,230.0,230.0],"v_angle":[0.0,-2.0944,2.0944]}},
 "linecode":{"lc":{"R_series_1_1":0.05,"R_series_2_2":0.05,"R_series_3_3":0.05,"R_series_4_4":0.05}},
 "line":{"l1":{"bus_from":"grid","bus_to":"poc","terminal_map_from":["a","b","c","n"],"terminal_map_to":["a","b","c","n"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

# Parametrised three-phase stiff grid for advanced-inverter topology tests:
# slack at "grid", short line to the POC. `mags`/`angs` are the per-phase source
# magnitude (V) / angle (rad); wide v_max so the inverter's own limits bind.
inv_grid3_src(; mags, angs, vmax=270.0) = parse_bmopf("""
{"bus":{
    "grid":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"]},
    "poc": {"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"],"v_min":[180.0,180.0,180.0],"v_max":[$vmax,$vmax,$vmax]}},
 "voltage_source":{"vs":{"bus":"grid","terminal_map":["a","b","c"],"v_magnitude":$(mags),"v_angle":$(angs)}},
 "linecode":{"lc":{"R_series_1_1":0.05,"R_series_2_2":0.05,"R_series_3_3":0.05,"R_series_4_4":0.05}},
 "line":{"l1":{"bus_from":"grid","bus_to":"poc","terminal_map_from":["a","b","c","n"],"terminal_map_to":["a","b","c","n"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

# Balanced (230 V) and magnitude+angle-unbalanced three-phase grids.
inv_grid3_bal(v=230.0) = inv_grid3_src(mags=[v, v, v], angs=[0.0, -2.0944, 2.0944])
inv_grid3_unbal() = inv_grid3_src(mags=[245.0, 215.0, 230.0], angs=[0.05, -2.15, 2.0])

# LV radial feeder for operating-envelope tests: source ──l1── bus1 ──l2── bus2,
# with v_max on the DER buses so simultaneous export is voltage-limited. `p1`,`p2`
# set the baseline loads (W) that define the interval's headroom.
doe_feeder(; p1, p2, vmax=245.0) = parse_bmopf("""
{"bus":{
    "src": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus1":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],"v_min":[216.0],"v_max":[$vmax]},
    "bus2":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],"v_min":[216.0],"v_max":[$vmax]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.4}},
 "line":{
    "l1":{"bus_from":"src","bus_to":"bus1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
    "l2":{"bus_from":"bus1","bus_to":"bus2","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}},
 "load":{
    "d1":{"bus":"bus1","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$p1],"q_nom":[0.0]},
    "d2":{"bus":"bus2","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$p2],"q_nom":[0.0]}}}
"""; from_string=true)
