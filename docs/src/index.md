# PowerOptLab.jl

PowerOptLab is a research laboratory for **four-wire distribution-network
decisions when the network state or model is uncertain**.

## Start here

Follow [Your first PowerOptLab study](getting_started.md) to install from a
checkout, build a complete network, schedule a battery, and check status and
physical residuals. The guide uses synthetic data and the standard solver, and
its Julia examples execute during documentation builds.

For an existing checkout, install the pinned dependencies from its root:

```sh
julia --project=. scripts/instantiate_pinned.jl
```

Then choose a topic from the guide's next-study table or the capabilities below.

## Research context

The organizing questions are:

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
not treat proposed joint model-forensics, active-probing, or general
uncertainty-set certification APIs as existing functionality.
[Concepts](concepts.md) explains the shared modeling
contracts; the [Research program](research_program.md) distinguishes current
capabilities from the intended integrated workflow.

## Capabilities

### Component models

Start with the [component model families](components/index.md) for generators,
battery storage and electric vehicles. Converter physics and controls have a
dedicated [inverter-based resources guide](ibr/index.md).

| Capability | Entry point | Maturity |
|---|---|---|
| [Generalized generators](components/generators.md) | [`GeneralizedGenerator`](@ref), [`SourceGenerator`](@ref) | research prototype |
| [Battery storage](components/battery_storage.md) with state of charge | [`StorageDevice`](@ref) | promotion candidate |
| [EV charging](ev/index.md) (V1G / V2G) with availability & departure energy | [`EVDevice`](@ref), [`ChargingSession`](@ref) | promotion candidate |
| [Advanced inverter](ibr/index.md) (circuit-aware IBR) | [`AdvancedInverter`](@ref) | experimental |
| [Phase-aware inverter controls](components/inverter_controls.md) (local three-leg laws) | [`SequenceController`](@ref) | research prototype |
| [IVQ battery](components/ivq_battery.md) (current–voltage–charge model) | [`IVQBattery`](@ref) | prototype |

### Problem specifications

Network-level questions built on the engine physics. Optimization wrappers use
the staged OPF API; the compiled constrained estimator and static operability
checker use their own residual and diagnostic workflows.

| Capability | Entry point | Direction | Maturity |
|---|---|---|---|
| [Multi-period OPF](problems/multiperiod.md) co-optimising many snapshots | [`solve_multiperiod_opf`](@ref) | forward | promotion candidate |
| [Legacy WLS state estimation](problems/state_estimation.md) | [`solve_state_estimation`](@ref) | inverse | prototype |
| [Constrained NLLS state estimation](problems/constrained_state_estimation.md) | [`solve_sparse_state_estimator`](@ref) | inverse | prototype |
| [Post-OPF operability](problems/operability.md) (native static snapshot diagnostics and sufficient certificates) | [`check_opf_operability`](@ref) | forward | research prototype |
| [Parameter estimation](problems/parameter_estimation.md) (line lengths / taps) | [`solve_parameter_estimation`](@ref) | inverse | prototype |
| [Inverse Carson reconstruction](problems/inverse_carson.md) (compatible overhead constructions) | [`solve_inverse_carson`](@ref) | inverse | prototype |
| [Dynamic operating envelopes](tutorials/doe_getting_started.md) (active import/export capacity) | [`solve_operating_envelope`](@ref) | forward | research prototype |
| [Bilevel PV/tap POC](problems/bilevel.md) (DiffOpt lower-level response) | [`solve_bilevel_pv_tap`](@ref) | hierarchical | proof of concept |

### Bespoke algorithms

Custom solution methods include [HELM power flow](algorithms/helm.md) and
natural-parameter or pseudo-arclength static operability continuation. See
[Bespoke algorithms](algorithms/index.md) for implemented methods and their
scope.

## How it fits together

See [Concepts](concepts.md) for how each capability layers over the BMOPFTools
staged API, and [Contributing](contributing.md) for how to add your own and the
promotion path back to the BMOPF spec. The topic guides link equations, worked
examples, API references and validation evidence.
