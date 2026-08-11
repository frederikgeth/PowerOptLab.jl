# PowerOptLab.jl

PowerOptLab is a research laboratory for **four-wire distribution-network
decisions when the network state or model is uncertain**. Its organizing
questions are:

1. What states and models are compatible with telemetry and metadata?
2. What safe measurement or intervention would distinguish material
   alternatives?
3. What operating decision is justified by the remaining alternatives?

The package builds on the
[BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl) reference
current–voltage OPF engine and reuses its neutral-explicit device physics,
per-unit handling, and result extraction through public extension seams.

The complete loop is not implemented yet. Current capabilities provide its
foundations as separate, explicit research prototypes; the documentation does
not treat proposed model-forensics, active-probing, or global-certification APIs
as existing functionality. Start with [Concepts](concepts.md) for the modelling
contracts, then read the [Research program](research_program.md) for the current
boundary and priorities.

## Capabilities

### Component models

New network elements, stamped via `model_hook!` / `solution_hook!`.

| Capability | Entry point | Maturity |
|---|---|---|
| [Storage / battery](components/devices.md) with state of charge | [`StorageDevice`](@ref) | promotion candidate |
| [EV charging](components/devices.md) (V1G / V2G) with availability & departure energy | [`EVDevice`](@ref) | promotion candidate |
| [Advanced inverter](ibr/index.md) (circuit-aware IBR) | [`AdvancedInverter`](@ref) | experimental |
| [IVQ battery](components/ivq_battery.md) (current–voltage–charge model) | [`IVQBattery`](@ref) | prototype |

### Problem specifications

New formulations over the same physics, via the staged `build_opf_model` /
`enforce_kcl!` / `generation_cost` / `extract_result` API.

| Capability | Entry point | Direction | Maturity |
|---|---|---|---|
| [Multi-period OPF](problems/multiperiod.md) co-optimising many snapshots | [`solve_multiperiod_opf`](@ref) | forward | promotion candidate |
| [Legacy WLS state estimation](problems/state_estimation.md) | [`solve_state_estimation`](@ref) | inverse | prototype |
| [Constrained NLLS state estimation](problems/constrained_state_estimation.md) | [`solve_sparse_state_estimator`](@ref) | inverse | prototype |
| [Parameter estimation](problems/parameter_estimation.md) (line lengths / taps) | [`solve_parameter_estimation`](@ref) | inverse | prototype |
| [Inverse Carson reconstruction](problems/inverse_carson.md) (compatible overhead constructions) | [`solve_inverse_carson`](@ref) | inverse | prototype |
| [Dynamic operating envelopes](problems/operating_envelope.md) (active import/export capacity) | [`solve_operating_envelope`](@ref) | forward | research prototype |
| [Bilevel PV/tap POC](problems/bilevel.md) (DiffOpt lower-level response) | [`solve_bilevel_pv_tap`](@ref) | hierarchical | proof of concept |

### Bespoke algorithms

Custom solution methods, currently including [HELM power
flow](algorithms/helm.md). See [Bespoke algorithms](algorithms/index.md) for the
method roadmap.

## Installation

BMOPFTools is not yet registered. Develop both from local checkouts; automated
builds pin BMOPFTools commit
`b7aa9a1bb48bcc8b790d3bcf5417d6a32036352a` (the PowerIO 0.7-compatible
snapshot with the semantic IBR power and monitored-voltage keys):

```julia
using Pkg
Pkg.activate(".")
Pkg.develop(path="../BMOPFTools.jl")
Pkg.instantiate()
```

Everything is SI at the interface (watts, vars, watt-hours, volts); per-unit
conditioning inside each solve is handled through the engine's
`opf_bases(ctx)` accessor.

## How it fits together

See [Concepts](concepts.md) for how each capability layers over the BMOPFTools
staged API, and [Contributing](contributing.md) for how to add your own and the
promotion path back to the BMOPF spec. Each capability page carries a worked
example.
