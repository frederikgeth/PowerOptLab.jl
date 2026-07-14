# Parameter estimation (calibration)

> **Kind:** Problem specification · **Maturity:** prototype · **Direction:** inverse · **Temporal:** inter-temporal (shared parameters)

[`solve_parameter_estimation`](@ref) calibrates uncertain **network parameters**
— line lengths and transformer/regulator tap ratios — from smart-meter data. It
is the *shared-parameter dual* of [state estimation](state_estimation.md), and a
multi-time-step extension of it:

| | fixes | estimates |
|---|---|---|
| **State estimation** | the parameters | the per-snapshot **state** |
| **Parameter estimation** | nothing about the uncertain elements | the shared, time-invariant **parameters** |

The formulation follows Vanin, Geth, Heidari & Van Hertem, *Distribution System
State and Impedance Estimation Augmented with Carson's Equations*
([arXiv:2506.04949](https://arxiv.org/abs/2506.04949)): a weighted
measurement-residual objective over a time series, the multiconductor IVR
(rectangular current–voltage) power flow as the physics, time-invariant impedances
written as `nominal · length`, and the smart-meter data as noisy `(P, Q, |V|)`
triples per user.

## Why multiple time steps

A single snapshot cannot separate a *long line* from a *heavy load*: both pull the
downstream voltage down by the same amount. With `T` snapshots of **diverse
loads** the shared parameters become over-determined — the one length (or tap)
that explains *every* snapshot at once is pinned. Calibration therefore needs a
time series; this is the essential difference from a one-shot state estimate.

!!! note "Identifiability caveat"
    As Vanin et al. discuss, in realistic feeders (meters only at load buses,
    smart-meter noise) *individual* branch parameters generally cannot be
    recovered exactly — multiple equivalent solutions fit the data equally well.
    What is recoverable, and what matters for downstream power-flow/OPF, is the
    *cumulative* impedance along each path and a model that reproduces the measured
    voltages. The small radial examples here are deliberately well-conditioned so
    the estimates land close to the truth.

## Measurements

Each metered user contributes noisy `(P, Q, |V|)` readings as
[`Measurement`](@ref)s — the same struct as state estimation:

- `:vmag` — phase-to-neutral voltage magnitude [V].
- `:pinj`, `:qinj` — active/reactive power **injected into the network** [W]/[var]
  (negative for a load). A free injection current is added at each measured bus and
  fit to these; buses with no injection measurement are zero-injection.

`sigma` is the standard deviation (WLS weight `1/σ²`). Choose the objective with
`objective=:wls` (smooth least squares) or `objective=:wlav` (weighted least
*absolute* value — the robust, bad-data-rejecting choice used in Vanin et al.).

## How the uncertain elements are stamped

The two element classes are handled differently, reflecting how the BMOPFTools
engine exposes them:

- **Transformer taps** use the engine's **native free-tap variable**. A transformer
  carrying `tap_min`/`tap_max` gets a `:tap` decision variable that the engine
  threads through its per-unit-correct, base-referred winding constraints.
  [`CalibTap`](@ref) keeps the transformer in every snapshot net, sets those
  bounds, and the solver couples the tap equal across snapshots — the engine's
  transformer physics (leakage, losses, referral) is reused unchanged.
- **Line lengths** have no native free variable, so [`CalibLine`](@ref)s are
  **omitted** from the physics nets and re-stamped here with a shared free length
  `ℓ`: `V_f − V_t = (r₀+jx₀)·ℓ·I`. This is the one genuine "replace Ohm's law"
  step; per-unit scaling by `z_base` is applied so SI and per-unit solves agree.

## Worked example

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf, solve_pf

# Physics net for one snapshot: known source; the uncertain line l1 (src→b1) is
# OMITTED (its length is the unknown); NO loads (injections come from the meters).
calnet() = parse_bmopf("""
{"bus":{
   "src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
   "b1": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}}}
"""; from_string=true)

# Per-snapshot meter readings (P, Q, |V|) at b1 — here for four diverse loads.
meas = [
  [Measurement(kind=:vmag, bus="b1", value=v,  sigma=0.3),
   Measurement(kind=:pinj, bus="b1", value=-p, sigma=50.0),   # negative = load
   Measurement(kind=:qinj, bus="b1", value=-q, sigma=50.0)]
  for (v, p, q) in [(229.3, 1000.0, 200.0), (226.8, 4000.0, 800.0),
                    (224.1, 7000.0, 1400.0), (228.0, 2500.0, 500.0)]
]
nets = [calnet() for _ in meas]

r = solve_parameter_estimation(nets, meas;
        lines=[CalibLine(id="l1", bus_from="src", bus_to="b1",
                         r_per_length=0.4, x_per_length=0.25)],
        objective=:wls)      # or :wlav

r.line_length["l1"]   # estimated length of l1
r.residual_rms        # RMS voltage misfit (V) — near the meter noise floor if good
r.snapshots           # per-snapshot fitted state (SI)
```

Add [`CalibTap`](@ref)s (naming transformers present in the nets) to the `taps`
keyword to calibrate taps jointly. Both `per_unit=true` and `per_unit=false` give
identical estimates. A `residual_rms` well above the known meter noise is the
diagnostic that something structural is still wrong (a misidentified tap, a missing
line, an unmetered load).

## Scope

This is a **prototype**. It handles single-phase elements (one terminal per element
end — add one [`CalibLine`](@ref) per phase for a multi-phase feeder). Estimating
the series-impedance *type* via Carson's equations (conductor geometry/material) as
in Vanin et al., regularising toward GIS priors, and four-wire mutual coupling are
natural extensions left for later.

See the API reference for [`CalibLine`](@ref), [`CalibTap`](@ref),
[`solve_parameter_estimation`](@ref), and [`ParameterEstimationResult`](@ref).
