# Parameter estimation (calibration)

[`solve_parameter_estimation`](@ref) calibrates uncertain **network parameters**
— line lengths and transformer/regulator tap ratios — from smart-meter data. It
is the *shared-parameter dual* of [state estimation](state_estimation.md):

| | fixes | estimates |
|---|---|---|
| **State estimation** | the parameters | the per-snapshot **state** |
| **Parameter estimation** | nothing about the uncertain elements | the shared, time-invariant **parameters** |

Distribution GIS records are frequently wrong — line lengths are approximate and
tap positions drift or are mislogged — which quietly biases every downstream
power-flow, hosting-capacity, and OPF study. Given enough metered snapshots, the
true parameters are recoverable by fitting the physics to the data.

## Why multiple time steps

A single snapshot cannot separate a *long line* from a *heavy load*: both pull the
downstream voltage down by the same amount, so infinitely many
(length, load) pairs explain one reading. With `T` snapshots of **diverse loads**
the shared parameters become over-determined — the one length (or tap) that
explains *every* snapshot at once is pinned. Calibration therefore needs a time
series; this is the essential difference from a one-shot state estimate.

## How it works

Every snapshot is built into **one shared JuMP model** (the multi-period build
pattern). The uncertain elements are the unknowns, so they are *not* part of the
per-snapshot physics nets — those carry the known source, the known lines, and
each snapshot's metered loads. The unknowns are stamped into the shared model by a
`model_hook!`, each with a single free variable reused across all snapshots:

- a **variable-length line** — ``V_f - V_t = (r_0 + j x_0)\,\ell\,I`` with the
  length ``\ell`` free (`CalibLine`);
- a **variable-ratio ideal transformer** — ``V_f = a\,V_t``, ``I_s = a\,I_p`` with
  ``a = a_0\,\tau`` and the tap multiplier ``\tau`` free (`CalibTap`).

Both introduce bilinear terms (``\ell\cdot I``, ``\tau\cdot V``), so the calibration
is a smooth nonlinear program solved by Ipopt. The objective is the combined
weighted-least-squares voltage residual over every snapshot,
``\sum_{t}\sum_i (z_{t,i} - |V|_{t,i})^2 / \sigma_i^2``.

Smart meters supply both halves of the data: the **loads** (baked into each
snapshot net as fixed injections) and the **voltage magnitudes** (the
[`Measurement`](@ref)s of kind `:vmag` to fit).

## Worked example

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf, solve_pf

# Physics net for one snapshot: known source + metered loads, but NOT the
# uncertain line l1 (src→b1) — that length is the unknown we estimate.
r0, x0 = 0.4, 0.25
calnet(p) = parse_bmopf("""
{"bus":{
   "src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
   "b1": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "load":{"d1":{"bus":"b1","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$p],"q_nom":[0.0]}}}
"""; from_string=true)

# Diverse-load snapshots + the voltage the meter at b1 reported in each.
loads = [1000.0, 4000.0, 7000.0, 2500.0]
nets  = [calnet(p) for p in loads]
meas  = [[Measurement(kind=:vmag, bus="b1", value=v, sigma=0.5)]
         for v in (229.3, 226.8, 224.1, 228.0)]   # smart-meter readings (V)

r = solve_parameter_estimation(nets, meas;
        lines=[CalibLine(id="l1", bus_from="src", bus_to="b1",
                         r_per_length=r0, x_per_length=x0)])

r.line_length["l1"]   # estimated length of l1
r.residual_rms        # RMS voltage misfit (V) — near the meter noise floor if the fit is good
r.snapshots           # per-snapshot fitted state (SI)
```

Add `CalibTap`s to the `taps` keyword to calibrate transformer/regulator taps in
the same solve; lengths and taps are estimated jointly. A `residual_rms` well
above the known meter noise is the diagnostic that something structural is still
wrong (a misidentified tap, a missing line, an unmetered load).

## Scope

This is a **prototype** capability. It handles single-phase elements (one terminal
per element end — add one `CalibLine`/`CalibTap` per phase for a multi-phase
feeder), models taps as lossless ideal transformers, and runs in SI units so the
residual is directly in volts. Estimating series impedance *type* (not just
length), regularising toward GIS priors, and per-unit conditioning of the
transformer base referral are natural extensions left for later.

See the API reference for [`CalibLine`](@ref), [`CalibTap`](@ref),
[`solve_parameter_estimation`](@ref), and [`ParameterEstimationResult`](@ref).
