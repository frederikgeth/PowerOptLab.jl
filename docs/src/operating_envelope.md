# Dynamic operating envelopes

A **dynamic operating envelope** (DOE) is the time-varying active-power export
limit a network operator allocates to each connection point so that, when every
point exports at its allocated limit at once, the network's operational limits
(voltage bounds, thermal ratings) still hold. [`solve_operating_envelope`](@ref)
computes one.

Unlike [state estimation](state_estimation.md), a DOE **keeps** the engine's
operational bounds — they are exactly what limit the envelope. Each interval is
an independent constrained allocation solved on `solve_opf` with a `model_hook!`
that stamps a free export at each connection point and replaces the
generation-cost objective with a fairness objective; a `solution_hook!` reads the
allocation back out. "Dynamic" refers to recomputing per interval as the baseline
load (and thus the available headroom) changes; the intervals are otherwise
independent.

## Fairness rules

Two allocation rules trade equity against efficiency:

- `:equal` maximises a single export level assigned to **every** point. The
  result is equitable but capped by the weakest point and the tightest
  constraint.
- `:sum` maximises the **total** allocated export. The result is efficient but
  can be uneven — an electrically stronger point (nearer the source, more
  headroom) may receive most of the allocation.

## Worked example

An LV feeder with DER at two buses. As baseline load rises the local consumption
absorbs export and frees voltage headroom, so the envelope grows — the dynamic
response a static export limit misses.

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

feeder(p1, p2) = parse_bmopf("""
{"bus":{
    "src": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus1":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],"v_min":[216.0],"v_max":[245.0]},
    "bus2":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],"v_min":[216.0],"v_max":[245.0]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.4}},
 "line":{
    "l1":{"bus_from":"src","bus_to":"bus1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
    "l2":{"bus_from":"bus1","bus_to":"bus2","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}},
 "load":{
    "d1":{"bus":"bus1","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$p1],"q_nom":[0.0]},
    "d2":{"bus":"bus2","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$p2],"q_nom":[0.0]}}}
"""; from_string=true)

cps  = [ConnectionPoint(id="der1", bus="bus1", export_max=10e3),
        ConnectionPoint(id="der2", bus="bus2", export_max=10e3)]
nets = [feeder(200.0, 200.0),      # low load  → tight envelope
        feeder(5000.0, 5000.0)]    # high load → more headroom

env = solve_operating_envelope(nets, cps; fairness=:equal)

env.envelope["der1"]   # e.g. [≈3.2e3, ≈8.0e3]  — grows with load (W)
env.total_export       # sum across connection points, per interval (W)
```

With `:equal` both points receive the same limit and the far bus sits at its
`v_max`; switching to `fairness=:sum` raises the total but skews it toward the
stronger point — the efficiency/equity tradeoff a DOE policy must choose between.

See the API reference for [`ConnectionPoint`](@ref),
[`solve_operating_envelope`](@ref), and [`OperatingEnvelopeResult`](@ref).

Reactive power is held at zero (unity-PF export); import envelopes and
power-factor-flexible allocation are natural extensions.
