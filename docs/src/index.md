# PowerOptLab.jl

A staging ground for **component models**, **problem specifications**, and
**bespoke algorithms** built on top of the
[BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl) reference
optimal-power-flow engine.

BMOPFTools ships one small, correct four-wire rectangular current–voltage OPF and
deliberately refuses to become a "model zoo". PowerOptLab is where the zoo lives:
experimental devices and alternative problem formulations that reuse the engine's
device physics, per-unit handling, and result extraction **through its public
extension seams** — without forking the engine. Anything here that matures into
accepted practice can later be folded back into the BMOPF spec.

New here? Start with [Concepts](concepts.md) for the three kinds of contribution
and the engine seams they reuse, then browse the capabilities below.

## Capabilities

### Component models

New network elements, stamped via `model_hook!` / `solution_hook!`.

| Capability | Entry point | Maturity |
|---|---|---|
| [Storage / battery](components/devices.md) with state of charge | [`StorageDevice`](@ref) | promotion candidate |
| [EV charging](components/devices.md) (V1G / V2G) with availability & departure energy | [`EVDevice`](@ref) | promotion candidate |
| [Advanced inverter](components/advanced_inverter.md) (internal-node IBR) | [`AdvancedInverter`](@ref) | prototype |

### Problem specifications

New formulations over the same physics, via the staged `build_opf_model` /
`enforce_kcl!` / `generation_cost` / `extract_result` API.

| Capability | Entry point | Direction | Maturity |
|---|---|---|---|
| [Multi-period OPF](problems/multiperiod.md) co-optimising many snapshots | [`solve_multiperiod_opf`](@ref) | forward | promotion candidate |
| [Legacy WLS state estimation](problems/state_estimation.md) | [`solve_state_estimation`](@ref) | inverse | promotion candidate |
| [Constrained NLLS state estimation](problems/constrained_state_estimation.md) | [`solve_sparse_state_estimator`](@ref) | inverse | prototype |
| [Parameter estimation](problems/parameter_estimation.md) (line lengths / taps) | [`solve_parameter_estimation`](@ref) | inverse | prototype |
| [Dynamic operating envelopes](problems/operating_envelope.md) (active import/export capacity) | [`solve_operating_envelope`](@ref) | forward | research prototype |

### Bespoke algorithms

New solution methods (custom solve loops). None yet — see
[Bespoke algorithms](algorithms/index.md) for the reserved slot.

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

See [Concepts](concepts.md) for how each capability layers over the BMOPFTools
staged API, and [Contributing](contributing.md) for how to add your own and the
promotion path back to the BMOPF spec. Each capability page carries a worked
example.
