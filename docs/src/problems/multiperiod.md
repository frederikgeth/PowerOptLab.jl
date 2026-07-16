# Multi-period OPF

> **Kind:** Problem specification · **Maturity:** promotion candidate · **Direction:** forward · **Temporal:** inter-temporal

[`solve_multiperiod_opf`](@ref) co-optimises a sequence of network snapshots in
one JuMP model, with storage/EV devices whose state of charge links each period
to the next. This is the formulation the single-snapshot `solve_opf` cannot
express, because state of charge couples time steps.

Under the hood it uses the BMOPFTools staged API: every snapshot is built into
one shared model with `build_opf_model(add_objective=false)`, each device's
state of charge is linked across the snapshots, the per-snapshot generation-cost
rates are multiplied by the corresponding [`TimeGrid`](@ref) duration and summed
into one interval-cost objective, KCL is enforced per snapshot, and the model is
solved once. The `dt_h` keyword is the uniform-duration shorthand.

## Worked example: battery arbitrage

Two periods with a time-varying slack import price (set via each snapshot's
`voltage_source` `cost`). A cyclic battery discharges into the expensive period
and recharges in the cheap one:

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

net(price) = parse_bmopf("""
{"bus":{
    "sourcebus":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus1":     {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],
                 "v_min":[900.0],"v_max":[1100.0]}},
 "voltage_source":{"vs":{"bus":"sourcebus","terminal_map":["1"],
     "v_magnitude":[1000.0],"v_angle":[0.0],"cost":[$price]}},
 "linecode":{"lc":{"R_series_1_1":0.1}},
 "line":{"l1":{"bus_from":"sourcebus","bus_to":"bus1",
     "terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}},
 "load":{"ld1":{"bus":"bus1","terminal_map":["1","n"],
     "configuration":"SINGLE_PHASE","p_nom":[100000.0],"q_nom":[0.0]}}}
"""; from_string=true)

nets = [net(0.20), net(0.05)]      # expensive, then cheap
bat  = StorageDevice(id="bat", bus="bus1",
                     p_charge_max=40e3, p_discharge_max=40e3,
                     energy_max=100e3, energy_init=40e3, cyclic=true)

res = solve_multiperiod_opf(nets, [bat]; dt_h=1.0)

res.dispatch["bat"].p_net   # ≈ [ +40e3, −40e3 ]  discharge then charge (W)
res.dispatch["bat"].soc     # ≈ [ 40e3, 0.0, 40e3 ]  (Wh, cyclic)
```

For nonuniform telemetry or tariff intervals, supply one positive duration per
snapshot. The same duration weights both the objective rate and state update:

```julia
grid = TimeGrid([0.25, 0.75])
res = solve_multiperiod_opf(nets, [bat]; time_grid=grid)
```

## Multiple devices

Pass any mix of [`StorageDevice`](@ref) and [`EVDevice`](@ref); each gets its own
state-of-charge trajectory and constraints, all co-optimised in the same model:

```julia
res = solve_multiperiod_opf(nets, [bat, ev1, ev2]; dt_h=0.5)
```

Per-period economics come from the snapshots themselves — a time-varying import
price, differing loads, or any other per-snapshot data the BMOPFTools objective
sees.

## Result

[`solve_multiperiod_opf`](@ref) returns a [`MultiperiodResult`](@ref): the
per-period BMOPFTools result dicts (`res.snapshots`) plus each device's SI
`p_charge`, `p_discharge`, `p_net`, `q`, and `soc` trajectories in
`res.dispatch`. See the API reference for [`solve_multiperiod_opf`](@ref) and
[`MultiperiodResult`](@ref).

Numerical values are published only when the shared solve reaches `OPTIMAL` or
`LOCALLY_SOLVED` with a feasible primal point. Iteration-limited, infeasible, or
relaxed-tolerance outcomes retain their termination status and return `NaN`
objective/trajectory values rather than presenting a candidate iterate as an
optimum.
