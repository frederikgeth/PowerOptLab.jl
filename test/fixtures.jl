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
