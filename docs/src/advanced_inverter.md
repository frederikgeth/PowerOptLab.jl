# Advanced inverter (prototype)

[`AdvancedInverter`](@ref) is a more detailed inverter-based-resource (IBR) than
the BMOPFTools engine's built-in current-injection IBR. The engine models an IBR
as a bounded current source at the point of connection (POC); this prototype adds
the core structural idea from the BMOPFTools
[IBR model extensions design doc](https://github.com/frederikgeth/BMOPFTools.jl/blob/main/docs/ibr_model_extensions.md)
— an explicit **internal AC node** behind the converter — and the five feature
phases layered on it.

```
 POC bus ──[filter r+jx (+ grid shunt b)]── internal node ──[converter]── DC
   network sets V here                       EMF lives here          losses/ripple here
```

It is built entirely on the BMOPFTools staged API through a `model_hook!`; it
does **not** modify the engine. The prototype runs in SI so device parameters
(ohms, volts, amps, watts, siemens) enter directly — the per-unit base coupling
the design doc flags as highest-risk (AC `v_base` vs DC `v_dc_base`) is an
engine-integration concern the prototype sidesteps.

## The five phases

| Phase | Feature | Model |
|---|---|---|
| 0 | Output filter (L/LC) | series `r+jx` from POC to the internal node, optional grid-side shunt `b` |
| 1 | Internal EMF bounds | `|V_int|` box, or DC-link modulation `|V_int| ≤ modulation_max·v_dc/√3` |
| 2 | Grid-forming | balanced 120° internal EMF with a bounded magnitude decision variable |
| 3 | Converter losses | non-branching `P_dc = P_ac + P_loss`, `P_loss = p_loss_fixed + a_loss·|I| + c_loss·|I|²` |
| 4 | Double-frequency ripple | `|Σ_k V_int_k·I_k|² ≤ p_ripple_max²` (2ω pulsation) |

Every feature is opt-in: with only `id`, `bus`, and `s_max` the device is a plain
grid-following converter, and the internal node collapses onto the POC when the
filter is zero.

## Key modelling choices (from the design doc)

- **Limits on the converter side.** The apparent-power circle `s_max` and the
  current limit `i_max` are applied on the converter quantities (internal-node
  voltage × current), matching real nameplate — so an output filter reduces the
  power actually delivered to the grid below the converter rating.
- **Non-branching losses.** With AC power positive = injected to grid and DC power
  positive = drawn from the DC source, the single equation `P_dc = P_ac + P_loss`
  (`P_loss ≥ 0`) holds for both discharge and charge — no direction `if`-branch.
- **Grid-forming ≠ slack.** A grid-forming inverter holds a balanced, bounded
  internal EMF *behind the filter*, but does not replace the network's reference;
  the surrounding grid still needs a slack source.

## Worked example

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

# A stiff grid: slack at "grid", short line to the inverter POC.
net = parse_bmopf("""
{"bus":{
    "grid":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "poc": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],"v_min":[200.0],"v_max":[250.0]}},
 "voltage_source":{"vs":{"bus":"grid","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.05}},
 "line":{"l1":{"bus_from":"grid","bus_to":"poc","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

# A converter with an output filter and a three-term loss curve; minimise loss
# while delivering 3 kW to the grid.
inv = AdvancedInverter(id="inv", bus="poc", s_max=5000.0,
                       r_filter=0.2, x_filter=0.5,
                       p_loss_fixed=20.0, a_loss=0.3, c_loss=0.02)

r = solve_advanced_inverter(net, inv; objective=:min_loss, p_set=3000.0)

r.p_poc     # ≈ 3000 W delivered at the POC
r.p_conv    # converter-side active power (> p_poc: filter losses)
r.p_loss    # 20 + 0.3·|I| + 0.02·|I|²
r.p_dc      # = p_conv + p_loss  (the non-branching DC-link balance)
r.v_int_mag # internal EMF magnitude per phase (V)
```

Switch `objective=:max_export` to maximise POC active power and watch the
converter rating, filter, EMF/modulation, or ripple limits bind. For a
three-phase `grid_forming=true` inverter the solved internal EMF magnitudes are
equal across phases (balanced 120°) and the 2ω ripple is ≈ 0.

See the API reference for [`AdvancedInverter`](@ref),
[`solve_advanced_inverter`](@ref), and [`InverterResult`](@ref).

## Scope

This is a **prototype** for experimentation, not a validated engine feature. It
implements the design doc's Phases 0–4 as a hook-stamped device; the reactive/
active priority modes, sequence-current limits, and grid-forming-as-reference
capabilities listed in the doc's backlog are not included. If a piece of this
matures, it can be folded back into the engine per the design doc's plan.
