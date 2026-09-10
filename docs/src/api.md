# API reference

All exported names, grouped by capability.

## Module

```@docs
PowerOptLab
```

## Shared contracts

```@docs
SolveOutcome
SolveStatus
AbstractSolveResult
solve_status
solve_diagnostics
TimeGrid
MultiContext
build_multi_context
```

## Extension interfaces

```@docs
AbstractDevice
device_id
validate_device
stamp_device!
link_device!
extract_device
AbstractMeasurement
measurement_kind
measurement_value
measurement_sigma
measurement_prediction
```

## Devices

Storage and charging types are documented in the [EV and storage API](ev/api.md):
[`StorageDevice`](@ref), [`EVDevice`](@ref), [`EV`](@ref), [`EVSE`](@ref), and
[`ChargingSession`](@ref).

## Multi-period OPF

```@docs
solve_multiperiod_opf
MultiperiodResult
```

## Legacy WLS state estimation

```@docs
Measurement
solve_state_estimation
StateEstimationResult
```

## Constrained NLLS state estimation

```@docs
TerminalConnection
ExactInjectionSpecification
ConstantPowerDevice
ConstantCurrentDevice
ZIPDevice
ExactDeviceEquation
BranchMeasurement
SEPreflightFinding
SEPreflightReport
SEUnsupportedNetwork
state_estimator_preflight
initial_state_estimator
SEStructure
SEParameters
SEEvaluation
compile_state_estimator
evaluate_state_estimator
residual_jacobian
constraint_jacobian
ConstrainedStateEstimationResult
solve_compiled_state_estimator
SparseConstrainedStateEstimationResult
solve_sparse_state_estimator
ContinuationStateEstimationResult
solve_with_continuation
SEObservability
observability_diagnostics
unobservable_directions
selected_state_covariance
derived_covariance
StatePrior
set_state_prior!
TimeSeriesStateEstimationResult
solve_time_series_state_estimator
```

## Parameter estimation

```@docs
CalibLine
CalibTap
solve_parameter_estimation
ParameterEstimationResult
```

## Inverse Carson reconstruction

```@docs
SequenceLineObservation
OverheadCarsonCandidate
solve_inverse_carson
profile_inverse_carson
InverseCarsonFit
InverseCarsonResult
InverseCarsonProfileInterval
materialize_inverse_carson
```

## Bilevel distribution-network proof of concept

```@docs
BilevelPVResult
BilevelPVResponse
SingleLevelPVResult
solve_bilevel_pv_tap
solve_bilevel_pv_response
solve_single_level_pv_tap
bilevel_demo_network
```

## Advanced inverter

```@docs
AdvancedInverter
solve_advanced_inverter
InverterResult
```

## Phase-aware inverter controls

```@docs
AbstractInverterControlLaw
AbstractPositiveSequencePolicy
AbstractUnbalancePolicy
AbstractLimiterPolicy
AbstractCurrentTarget
ConverterCurrentTarget
GridCurrentTarget
PiecewiseLinearLaw
WorstPhaseVoltVarWatt
AverageVoltageVoltVarWatt
PositiveSequenceVoltVarWatt
NoUnbalanceControl
NegativeSequenceAdmittanceDroop
CommonScaleLimiter
SequenceController
InverterControlMeasurement
InverterControlRequest
InverterControlRatings
InverterControlResult
ConverterTerminalResult
ControlledDevice
ControlledInverterResult
evaluate_exact
evaluate_smooth
stamp_smooth_control!
solve_controlled_inverter
inverter_spec
inverter_handles
```

## Current–voltage (IVQ) battery

```@docs
IVQBattery
solve_ivq_battery
IVQBatteryResult
solve_multiperiod_ivq
MultiperiodIVQResult
BatteryChemistry
thevenin_chemistry
linear_chemistry
tabulated_chemistry
illustrative_lfp
illustrative_nmc
illustrative_nca
illustrative_lead_acid
illustrative_leaf
```

## HELM power flow

```@docs
solve_pf_helm
helm_series
HelmResult
```

## Explicit-neutral Kron reduction

```@docs
kron_reduce_bmopf
```

`kron_reduce_bmopf(net; neutral_terminals=nothing, as_json=false)` returns a
new BMOPF dictionary (or schema-valid JSON when `as_json=true`) and never
mutates `net`. A JSON string is accepted as the first argument. Neutral labels
come from `terminal_conventions.neutral` when that block is present; otherwise
the conservative `n`/`N` naming convention is used. A
`Dict("bus_id" => "label")` `neutral_terminals` override is available for
datasets whose labels are not declared.

For each affected line, the series impedance is reduced exactly with
`Zᵣ = Zₚₚ − Zₚₙ Zₙₙ⁻¹ Zₙₚ`, with the explicit neutral voltage constrained to
zero. Inline matrices and referenced linecodes are supported, including
line lengths, line π-shunt projections, shared-linecode cloning when retained
conductor orders differ, and conductor-indexed current ratings. The output
removes neutral bus names and AC terminal-map entries across lines, sources,
loads, generators, switches, shunts, capacitors, IBRs, and transformer
windings. Grounded WYE/SINGLE_PHASE capacitors become equivalent phase-to-ground
shunts. Phase-to-neutral voltage bounds are intersected with phase-to-ground
bounds; `vn_max` and neutral-conductor current limits are dropped with an audit
entry because the neutral is now ideal ground.

| Field/family | Reduction rule | Audit or refusal |
|---|---|---|
| `bus.terminal_names`, `perfectly_grounded_terminals` | Remove the resolved neutral; filter grounds by value | Classifies each bus as exact or forced-ground projection |
| `bus.v_min/v_max`, `vpn_min/vpn_max` | Select terminal-indexed phase entries and intersect bounds | Contradictory intersections error; `vn_max` is dropped |
| Line `R/X_series_*` | Exact `Zpp − Zpn Znn⁻¹ Znp` Schur complement | Stores `K = −Znn⁻¹Znp`, conductor maps, and neutral π-row coefficients |
| Line `G/B_*`, shunts | Remove grounded-neutral rows/columns; neutral-only ground shunts are removed | Neutral current ratings are explicitly reported as lost |
| WYE/SINGLE_PHASE capacitors | Materialize equivalent phase-to-ground `shunt` | Time-varying capacitor conversion is refused |
| Loads, sources, generators, switches | Remove neutral map entries and slice only conductor-indexed vectors | Ambiguous role/order errors; phase arrays are preserved |
| IBR and transformers | Active neutral-leg IBRs and neutral-bearing transformers are refused | No invalid shortened BMOPF subtype is emitted |

DC buses, DC terminal maps, and DC grounding are intentionally outside the
transformation. Structural, matrix, or bound time-series targets are refused;
ordinary snapshot data remains untouched. The reduced output is BMOPF JSON
schema-valid when serialized with `as_json=true`; upstream 0.1.0 may still
report the known `WYE` arity warning for a grounded WYE whose implicit midpoint
is no longer listed.

The transformation is refused for singular or near-singular neutral blocks,
ambiguous maps, DELTA elements containing a neutral, active IBR
`FOUR_LEG` neutral physics, unsupported center-tap/single-phase transformer
reductions, and time-varying structural/impedance data whose commutation is not
proven. Matrix keys follow BMOPFTools' sparse convention: absent entries are
zero, a missing reciprocal is mirrored, and explicitly stored full entries take
precedence. Every reduction records original/kept conductor ordering, the
neutral-current recovery vector `K = −Zₙₙ⁻¹ Zₙₚ`, bound rewrites, dropped
constraints, and whether a bus was already perfectly grounded or was forced to
ideal ground in `net["_meta"]["kron_reduction"]` and
`net["extras"]["kron_reduction"]`.

The package test target requires OpenDSSDirect.jl and runs an independent
paired oracle: an explicit four-wire node-0-ground deck and its manually
projected three-wire deck. It asserts convergence and compares retained phase
voltages and line currents at both ends, source P/Q, total losses, and the
neutral-current recovery relation `Iₙ = K Iₚ` against the reduced BMOPF solve.
The older bundled finite-0.001-ohm grounding fixture remains useful as a
near-ground comparison, but is not the exact-equivalence gate.

```julia
reduced = kron_reduce_bmopf(net)
json = kron_reduce_bmopf(net; as_json=true)
```
