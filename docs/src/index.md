# PowerOptLab.jl

A staging ground for **component models** and **problem specifications** built on
top of the [BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl) reference
optimal-power-flow engine.

BMOPFTools ships one small, correct four-wire rectangular current–voltage OPF and
deliberately refuses to become a "model zoo". PowerOptLab is where the zoo lives:
experimental devices and alternative problem formulations that reuse the engine's
device physics, per-unit handling, and result extraction **through its public
extension seams** — `model_hook!` / `solution_hook!` and the staged
`build_opf_model` / `enforce_kcl!` / `generation_cost` / `extract_result` API —
without forking the engine. Anything here that matures into accepted practice can
later be folded back into the BMOPF spec.

## Capabilities

| Capability | Entry point | Reuses |
|---|---|---|
| Storage / battery devices with state of charge | [`StorageDevice`](@ref) | `model_hook!` current injection + KCL |
| EV charging (V1G / V2G) with availability & departure energy | [`EVDevice`](@ref) | storage device + per-period masking |
| Multi-period OPF co-optimising many snapshots | [`solve_multiperiod_opf`](@ref) | staged API (one shared model) |
| State estimation (weighted least squares) | [`solve_state_estimation`](@ref) | bounds-free physics + custom objective |
| Parameter estimation (calibrate line lengths / taps) | [`solve_parameter_estimation`](@ref) | shared parameters across snapshots + WLS objective |
| Dynamic operating envelopes (DER export limits) | [`solve_operating_envelope`](@ref) | operational bounds + fairness objective |
| Advanced inverter (internal-node IBR prototype) | [`AdvancedInverter`](@ref) | `model_hook!` internal node + filter/EMF/loss/ripple |

## Installation

BMOPFTools is not yet registered, so develop both from local checkouts (or Git
URLs):

```julia
using Pkg
Pkg.activate(".")
Pkg.develop(path="../BMOPFTools.jl")   # or Pkg.develop(url="https://github.com/frederikgeth/BMOPFTools.jl")
Pkg.instantiate()
```

Everything is SI at the interface (watts, vars, watt-hours, volts); per-unit
conditioning inside each solve is handled by the engine's `ctx.bases`.

## How it fits together

Each capability is a thin layer over the BMOPFTools staged API:

- A **device** ([`StorageDevice`](@ref), [`EVDevice`](@ref)) is stamped into a
  network snapshot by a `model_hook!` as a per-phase current injection added to
  the engine's Kirchhoff-current-law accumulators, plus a charge/discharge power
  split that a state-of-charge variable links across snapshots.
- **Multi-period OPF** ([`solve_multiperiod_opf`](@ref)) builds several snapshots
  into one JuMP model with `build_opf_model(add_objective=false)`, links each
  device's state of charge, sums `generation_cost` for one objective, enforces
  KCL per snapshot, and solves once.
- **State estimation** ([`solve_state_estimation`](@ref)) is a *different problem
  specification* on the same physics: a bounds-free net gives free bus voltages,
  a `model_hook!` adds free injections at measured buses and a weighted
  least-squares residual objective, and the solve returns the fitted state.
- **Parameter estimation** ([`solve_parameter_estimation`](@ref)) inverts that:
  many snapshots share one model, uncertain lines/taps are stamped with free
  parameters common to every snapshot, and a WLS voltage objective calibrates the
  shared line lengths and tap ratios from smart-meter time series.

See the per-topic pages for worked examples.
