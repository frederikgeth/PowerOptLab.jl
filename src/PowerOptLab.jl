"""
    PowerOptLab

A staging ground for component models and problem specifications built on top of
the BMOPFTools reference OPF engine, using its public extension seams
(`model_hook!` / `solution_hook!` and the staged `build_opf_model` /
`enforce_kcl!` / `generation_cost` / `extract_result` API) rather than forking
the engine.

The BMOPFTools engine deliberately stays small and correct — it ships one
four-wire rectangular current–voltage OPF and refuses to be a "model zoo".
PowerOptLab is where the zoo lives: experimental devices (storage, EVs) and
alternative problem specifications (multi-period co-optimisation, weighted
least-squares state estimation) that reuse the engine's device physics,
per-unit handling, and result extraction unchanged. Anything here that matures
into accepted practice can later be folded back into the spec and the engine.

Contributions are organised by *what layer of the engine they extend*:

### Component models — new network elements (`src/components/`)

Stamped into a solve through `model_hook!` / `solution_hook!`.

- **Storage / EV** ([`StorageDevice`](@ref), [`EVDevice`](@ref)) — battery/EV
  inverter ports stamped as current injections with an inter-temporal
  state-of-charge state (an energy/power "PE" model with fixed efficiency).
- **IVQ battery** ([`IVQBattery`](@ref)) — the current–voltage counterpart: cells
  of a [`BatteryChemistry`](@ref) modelled in the voltage–current–charge space
  (`v = OCV(soc) − i·R`), so voltage and current limits bind individually and a
  current-dependent cell efficiency emerges from the physics (a Rint proxy, not a
  full round-trip energy efficiency). Reuses the [`AdvancedInverter`](@ref) for
  the AC↔DC converter.
- **Advanced inverter** ([`AdvancedInverter`](@ref)) — a prototype internal-AC-node
  IBR with an output filter, internal-EMF/DC-modulation bounds, grid-forming
  operation, converter losses, and double-frequency ripple limits.

### Problem specifications — new formulations over the staged API (`src/problems/`)

A different objective/variable/constraint structure on the same physics.

- **Multi-period OPF** ([`solve_multiperiod_opf`](@ref)) — several network
  snapshots co-optimised in one model with storage/EV state linking each step to
  the next.
- **State estimation** ([`solve_state_estimation`](@ref)) — weighted
  least-squares estimation of the network state from noisy measurements (an
  *inverse* problem on the same physics).
- **Parameter estimation** ([`solve_parameter_estimation`](@ref)) — calibration
  of uncertain line lengths and transformer tap ratios from smart-meter data
  across multiple time steps (the shared-parameter dual of state estimation).
- **Inverse Carson reconstruction** ([`solve_inverse_carson`](@ref)) — screens
  discrete overhead construction candidates against diagonal sequence data,
  retaining ambiguity and local identifiability diagnostics while reconstructing
  a full primitive, neutral-explicit line model.
- **Dynamic operating envelopes** ([`solve_operating_envelope`](@ref)) —
  per-connection-point active-power import/export capacity with parameterized
  fairness, forecast/model scenarios, explicit corner-security semantics, and
  prescribed IBR Q-V controls retained from the network model.

### Bespoke algorithms — new solution methods (`src/algorithms/`)

Custom solve loops (decomposition, sequential linearization, warm-start
schemes) and alternative solution methods.

- **HELM** ([`solve_pf_helm`](@ref)) — the Holomorphic Embedding Load-flow
  Method: a non-iterative power flow that expands each voltage as a power series
  in a load-scaling parameter and evaluates it by Padé analytic continuation.
  Physical mismatch, Padé spread, and coefficient-tail diagnostics distinguish
  convergence from finite-order numerical divergence without treating the
  latter as a non-existence certificate.

Everything is SI at the interface; per-unit conditioning inside the solve is
handled via the engine's `ctx.bases`.
"""
module PowerOptLab

using BMOPFTools
using Dates
using ForwardDiff
using JuMP
using Ipopt
using LinearAlgebra
using SparseArrays

# Shared validation and solver-result contracts.
include("contracts.jl")
include("interfaces.jl")

# Component models — new network elements stamped via model_hook! / solution_hook!
include("components/devices.jl")
include("components/advanced_inverter.jl")
include("components/battery_chemistry.jl")
include("components/ivq_battery.jl")

# Problem specifications — new objective/constraint structures over the staged API
include("problems/multiperiod.jl")
include("problems/state_estimation.jl")
include("problems/constrained_state_estimation.jl")
include("problems/parameter_estimation.jl")
include("problems/inverse_carson.jl")
include("problems/operating_envelope.jl")

# Bespoke algorithms — new solution methods (custom solve loops)
include("algorithms/pade.jl")
include("algorithms/helm.jl")

# Shared extension interfaces
export AbstractDevice, AbstractMeasurement, AbstractSolveResult
export TimeGrid, MultiContext, build_multi_context
export SolveStatus, solve_status, solve_diagnostics
export device_id, validate_device, stamp_device!, link_device!, extract_device
export measurement_kind, measurement_value, measurement_sigma, measurement_prediction

# Devices
export StorageDevice, EVDevice

# Shared solve-result contract
export SolveOutcome

# Multi-period OPF
export solve_multiperiod_opf, MultiperiodResult

# State estimation
export Measurement, BranchMeasurement, solve_state_estimation, StateEstimationResult
export TerminalID, ExactInjectionSpecification, NoExactInjection,
       ExactZeroInjection, ExactDeviceEquation, TerminalConnection,
       ConstantPowerDevice, ConstantCurrentDevice, ZIPDevice, SEStructure, SEParameters,
       SEEvaluation, compile_state_estimator, evaluate_state_estimator,
       residual_jacobian, constraint_jacobian,
       ConstrainedStateEstimationResult, solve_compiled_state_estimator,
       ContinuationStateEstimationResult, solve_with_continuation,
       SparseConstrainedStateEstimationResult, solve_sparse_state_estimator,
       SEObservability, observability_diagnostics, unobservable_directions,
       selected_state_covariance, derived_covariance,
       StatePrior, set_state_prior!, TimeSeriesStateEstimationResult,
       solve_time_series_state_estimator

# Parameter estimation (calibration of line lengths / transformer taps)
export CalibLine, CalibTap, solve_parameter_estimation, ParameterEstimationResult

# Inverse Carson reconstruction from diagonal sequence data
export SequenceLineObservation, OverheadCarsonCandidate,
       InverseCarsonFit, InverseCarsonResult, InverseCarsonProfileInterval,
       solve_inverse_carson, profile_inverse_carson,
       materialize_inverse_carson

# Dynamic operating envelopes
export ConnectionPoint, FairnessPolicy, solve_operating_envelope,
       verify_operating_envelope, compare_operating_envelope_policies,
       OperatingEnvelopeResult, OperatingEnvelopeVerification

# Advanced inverter (prototype internal-node IBR)
export AdvancedInverter, solve_advanced_inverter, InverterResult

# Current–voltage (IVQ) battery storage + chemistry library
export BatteryChemistry, thevenin_chemistry, linear_chemistry, tabulated_chemistry,
       illustrative_lfp, illustrative_nmc, illustrative_nca,
       illustrative_lead_acid, illustrative_leaf
export IVQBattery, solve_ivq_battery, IVQBatteryResult
export solve_multiperiod_ivq, MultiperiodIVQResult

# HELM power flow (holomorphic embedding load-flow, a bespoke solution method)
export helm_series, HelmResult, solve_pf_helm

end # module PowerOptLab
