# State estimation

[`solve_state_estimation`](@ref) is a *different problem specification* over the
same network physics: given noisy measurements of an energised network, find the
bus voltage state that best fits them in a weighted-least-squares (WLS) sense.

It reuses the BMOPFTools device model but with three changes, all expressed
through the public seams:

1. **No operational bounds.** Because `build_opf_model` adds limits only where the
   net declares them, a bounds-free net yields a pure physics model with free bus
   voltages.
2. **Free injections instead of fixed loads.** A `model_hook!` adds a free
   injection current at each measured bus (so KCL closes with the voltages free to
   fit the data). Buses without an injection measurement are treated as
   zero-injection — the classic zero-injection pseudo-measurement.
3. **A residual objective instead of generation cost.** The hook sets the WLS
   objective ``\sum_i (z_i - h_i(\text{state}))^2 / \sigma_i^2`` over the
   measurements.

## Measurements

Each [`Measurement`](@ref) is an SI scalar with a standard deviation `sigma`
(WLS weight ``1/\sigma^2``):

- `:vmag` — voltage magnitude at `(bus, terminal)` in volts.
- `:pinj` — active power injected into the network at `(bus, terminal)` in watts.
- `:qinj` — reactive power injection in vars.

## Worked example

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf, solve_pf

# A physics-only net: buses, lines, source; no operational limits.
net = parse_bmopf("""
{"bus":{
    "src": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus1":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus2":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],
     "v_magnitude":[1000.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.5}},
 "line":{
    "l1":{"bus_from":"src","bus_to":"bus1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
    "l2":{"bus_from":"bus1","bus_to":"bus2","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

meas = [
    Measurement(kind=:vmag, bus="bus1", value=979.5,     sigma=2.0),
    Measurement(kind=:pinj, bus="bus1", value=-20_000.0, sigma=400.0),
    Measurement(kind=:qinj, bus="bus1", value=0.0,       sigma=400.0),
    Measurement(kind=:vmag, bus="bus2", value=969.2,     sigma=2.0),
    Measurement(kind=:pinj, bus="bus2", value=-20_000.0, sigma=400.0),
    Measurement(kind=:qinj, bus="bus2", value=0.0,       sigma=400.0),
]

se = solve_state_estimation(net, meas)

se.bus["bus1"]["1"]["vm"]   # estimated |V| at bus1 (V)
se.residuals                # per-measurement measured/estimated/residual/normalized
se.objective                # optimal weighted-residual sum
```

With more measurements than unknowns (here six measurements for four voltage
components), the fused estimate filters measurement noise: it is typically closer
to the true state than the raw voltage readings, because it reconciles every
measurement through the network physics.

See the API reference for [`Measurement`](@ref), [`solve_state_estimation`](@ref),
and [`StateEstimationResult`](@ref).
