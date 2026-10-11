# Parameter estimation (calibration)

> **Kind:** Problem specification · **Maturity:** prototype · **Direction:** inverse · **Temporal:** inter-temporal (shared parameters)

[`solve_parameter_estimation`](@ref) calibrates uncertain **network parameters**
— line lengths and transformer/regulator tap ratios — from smart-meter data. It
shares those parameters across multiple snapshots while estimating the network
state in each snapshot. The current uncertain-line primitive connects one
terminal per end with a fixed nominal R/X ratio. It does not estimate a full
four-wire impedance matrix, mutual coupling, topology or conductor identities.

| | fixes | estimates |
|---|---|---|
| **State estimation** | topology and network parameters | the per-snapshot **state** |
| **Parameter estimation** | topology, terminal maps, nominal per-length R/X and all unselected parameters | selected shared **lengths/taps** and per-snapshot **states** |

The shared-parameter formulation is motivated by Vanin, Geth, Heidari & Van Hertem, *Distribution System
State and Impedance Estimation Augmented with Carson's Equations*
([arXiv:2506.04949](https://arxiv.org/abs/2506.04949)). This prototype uses a
weighted measurement-residual objective, rectangular current–voltage equations,
impedances written as `nominal · length`, and noisy `(P, Q, |V|)` readings. It
does not implement the paper's full multiconductor/Carson estimation problem.

## Why multiple time steps

Load variation supplies additional equations for the same lengths and taps and
can improve identifiability. It does not guarantee that the parameters are
uniquely determined. Identifiability depends on which quantities are measured,
meter placement, network structure, excitation and noise. A single snapshot may
identify a parameter in a sufficiently instrumented small circuit; multiple
snapshots may still leave equivalent parameter combinations.

!!! note "Identifiability caveat"
    As Vanin et al. discuss, in realistic feeders (meters only at load buses,
    smart-meter noise) *individual* branch parameters generally cannot be
    recovered exactly — multiple equivalent solutions fit the data equally well.
    The returned local fit and a small residual do not establish unique parameter
    recovery, identification of cumulative impedance, or accuracy under new
    operating conditions. The small radial examples are illustrative. This API
    does not return a parameter-identifiability or global-uniqueness certificate.

## Measurements

Each metered user contributes noisy `(P, Q, |V|)` readings as
[`Measurement`](@ref)s — the same struct as state estimation:

- `:vmag` — phase-to-neutral voltage magnitude [V].
- `:pinj`, `:qinj` — active/reactive power **injected into the network** [W]/[var]
  (negative for a load). A free injection current is added at each measured bus and
  fit to these; buses with no injection measurement are zero-injection.

`sigma` is the standard deviation (WLS weight `1/σ²`). Choose the objective with
`objective=:wls` (smooth least squares) or `objective=:wlav` (weighted least
*absolute* value, which is less sensitive to isolated large residuals and is
used in Vanin et al.). WLAV does not identify or remove bad measurements; inspect
the fitted residuals before attributing a discrepancy to a particular meter.

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
keyword to calibrate taps jointly. Both `per_unit=true` and `per_unit=false`
use SI inputs and return SI results; scaling and local convergence can affect
the numerical fit. A `residual_rms` well above the expected voltage noise is a
reason to investigate data, modeling assumptions and convergence, rather than a
diagnosis of a particular structural error.

## Scope

This is a **prototype**. It handles single-phase elements (one terminal per element
end). Separate [`CalibLine`](@ref)s model independent conductor branches; they do
not reproduce a coupled multiconductor line or infer a shared physical length
across those branches. Estimating
the series-impedance *type* via Carson's equations (conductor geometry/material) as
in Vanin et al., regularising toward GIS priors, and four-wire mutual coupling are
natural extensions left for later.

See the API reference for [`CalibLine`](@ref), [`CalibTap`](@ref),
[`solve_parameter_estimation`](@ref), and [`ParameterEstimationResult`](@ref).

## API reference

The [calibration API](../api/parameter_estimation.md) lists the public types, functions and result fields
for this guide.
