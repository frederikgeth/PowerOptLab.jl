"""
    PowerOptLab

A research laboratory for four-wire distribution-network decisions when the
network state or model is uncertain. It connects model evidence and forensics,
informative interventions, and verification of operating decisions. The complete
workflow is a research direction; current APIs provide separate foundations and
state their local, scenario, or prototype limits explicitly.

PowerOptLab builds on the BMOPFTools reference current–voltage OPF engine, using
its public extension seams (`model_hook!` / `solution_hook!` and the staged
`build_opf_model` / `enforce_kcl!` / `generation_cost` / `extract_result` API)
rather than forking the engine. Experimental devices and formulations reuse the
engine's neutral-explicit physics, per-unit handling, and result extraction.
Stable foundational work can later be proposed for the BMOPF spec.

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
- **Advanced inverter** ([`AdvancedInverter`](@ref)) — an experimental internal-AC-node
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
- **Constrained NLLS state estimation** ([`solve_sparse_state_estimator`](@ref))
  — a compiled neutral-explicit residual/constraint model with tangent-space
  observability, selected local covariance, and sparse or dense reference solves.
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
- **Bilevel PV/tap POC** ([`solve_bilevel_pv_tap`](@ref)) — differentiated
  aggregate-export and local-controller lower levels with native Volt-var/
  Volt-watt controls, compared with a centralized single-level tap/export solve.
- **Network-scale inverter-control studies**
  ([`solve_controlled_inverter_fleet`](@ref)) — selected dataset IBRs are
  replaced by phase-aware advanced inverters in one simultaneous snapshot,
  with explicit ownership provenance and table-ready hardware-stress results.

### Bespoke algorithms — new solution methods (`src/algorithms/`)

Question-driven custom solve loops and alternative solution methods.

- **HELM** ([`solve_pf_helm`](@ref)) — the Holomorphic Embedding Load-flow
  Method: a non-iterative power flow that expands each voltage as a power series
  in a load-scaling parameter and evaluates it by Padé analytic continuation.
  Physical mismatch, Padé spread, and coefficient-tail diagnostics distinguish
  convergence from finite-order numerical divergence without treating the
  latter as a non-existence certificate.

Everything is SI at the interface; per-unit conditioning inside the solve is
handled via the engine's `opf_bases(ctx)` accessor.
"""
module PowerOptLab

using BMOPFTools
using Dates
using DiffOpt
using ForwardDiff
using JuMP
using Ipopt
using LinearAlgebra
using Random
using SHA
using SparseArrays

# Shared validation and solver-result contracts.
include("contracts.jl")
include("interfaces.jl")

# One isolated compatibility adapter for the load decomposition that
# BMOPFTools 0.1.0 does not yet expose publicly.
include("upstream.jl")

# Component models — new network elements stamped via model_hook! / solution_hook!
include("components/devices.jl")
include("components/advanced_inverter.jl")
include("components/inverter_controls.jl")
include("components/battery_chemistry.jl")
include("components/ivq_battery.jl")

# Problem specifications — new objective/constraint structures over the staged API
include("problems/multiperiod.jl")
include("problems/state_estimation.jl")
include("problems/constrained_state_estimation.jl")
include("problems/parameter_estimation.jl")
include("problems/inverse_carson.jl")
include("problems/operating_envelope.jl")
include("problems/bilevel.jl")
include("problems/inverter_control_study.jl")
include("problems/inverter_control_experiments.jl")
include("problems/inverter_control_sizing.jl")
include("problems/closed_loop_evidence.jl")
include("problems/operability.jl")

# Bespoke algorithms — new solution methods (custom solve loops)
include("algorithms/pade.jl")
include("algorithms/helm.jl")
include("algorithms/operability_continuation.jl")

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
export ConnectionPoint, FairnessPolicy, DOEScenario, DOEScenarioSet,
       DOEUncertaintySample, DOEUncertaintySampleSet,
       sample_doe_gaussian_uncertainty, doe_uncertainty_manifest,
       materialize_doe_scenarios,
       select_doe_scenarios, DOEScenarioTimeSplit,
       split_doe_scenarios_by_time, DOEScenarioCalibrationAudit,
       audit_doe_scenario_calibration, DOECovariateShiftResult,
       test_doe_covariate_shift, test_doe_time_series_covariate_shift,
       DOEProbabilityObservation,
       DOEProbabilityCalibrationResult, doe_probability_observations,
       evaluate_doe_probability_calibration, DOEControlRegistration,
       DOEControlRule, DOEControlPolicy,
       PerfectRecourse, IssueFixedControls, IssuePlusLocalLaws,
       DOEStudySpec, doe_study_manifest, doe_benchmark_rows,
       doe_context_benchmark_rows,
       solve_operating_envelope, solve_operating_envelope_multistart,
       verify_operating_envelope, search_operating_envelope_utilizations,
       search_operating_envelope_adversarial,
       confirm_operating_envelope_counterexample,
       evaluate_operating_envelope_coverage,
       evaluate_operating_envelope_coverage_curve,
       compare_doe_coverage_shift,
       solve_adversarial_search_stable_operating_envelope,
       solve_search_stable_operating_envelope,
       compare_operating_envelope_policies,
       OperatingEnvelopeResult, OperatingEnvelopeVerification,
       OperatingEnvelopeContextResult, OperatingEnvelopeSearchResult,
       DOEAdversarialSearchResult, DOECounterexampleConfirmationResult,
       DOECoverageResult, DOECoverageCurveResult, DOECoverageShiftResult,
       AdversarialSearchStableOperatingEnvelopeResult,
       SearchStableOperatingEnvelopeResult,
       OperatingEnvelopeMultistartResult

# Post-OPF static voltage operability (first slice: native ybus load scope)
export OperabilitySpec, OperabilityCheck, OperabilityModelError, OperabilityResult,
       check_opf_operability, OperabilityStressDirection,
       operability_scope_audit, operability_upstream_audit,
       operability_stress_network, operability_stress_rows,
       operability_stress_summary, operability_stress_ensemble_rows,
       operability_snapshot_row, operability_snapshot_rows

# Bilevel distribution-network proof of concept
export BilevelPVResult, BilevelPVResponse, SingleLevelPVResult,
       solve_bilevel_pv_tap, solve_bilevel_pv_response,
       solve_single_level_pv_tap, bilevel_demo_network

# Advanced inverter (prototype internal-node IBR)
export AdvancedInverter, solve_advanced_inverter, InverterResult

# Local phase-aware inverter controls
export AbstractInverterControlLaw, AbstractPositiveSequencePolicy
export AbstractUnbalancePolicy, AbstractLimiterPolicy, AbstractCurrentTarget
export ConverterCurrentTarget, GridCurrentTarget
export PiecewiseLinearLaw, WorstPhaseVoltVarWatt, AverageVoltageVoltVarWatt
export PositiveSequenceVoltVarWatt
export NoUnbalanceControl, NegativeSequenceAdmittanceDroop, CommonScaleLimiter
export SequenceController, InverterControlMeasurement, InverterControlRequest
export InverterControlRatings, InverterControlResult, ConverterTerminalResult
export ControlledDevice, ControlledInverterResult
export evaluate_exact, evaluate_smooth, stamp_smooth_control!, solve_controlled_inverter
export inverter_spec, inverter_handles
export ControlledInverterFleetSpec, ControlledInverterFleetResult
export solve_controlled_inverter_fleet, controlled_inverter_rows
export controlled_inverter_phase_rows
export InverterControlStudyCase, InverterControlStudyCaseResult
export InverterControlStudyResult, run_inverter_control_study
export inverter_control_study_case_rows, inverter_control_study_device_rows
export inverter_control_study_phase_rows, inverter_control_study_summary_rows
export inverter_control_paired_rows, inverter_control_paired_summary_rows
export InverterHardwareSweepPoint, resize_controlled_inverter_fleet
export validate_inverter_control_campaign, expand_inverter_hardware_cases
export inverter_control_hardware_requirement_rows
export FixedPointIterationResult, FixedPointGainScreen
export fixed_point_oracle, finite_difference_jacobian, screen_fixed_point_gain
export InverterControlScalingAudit, inverter_control_scaling_audit
export InverterControlFixedPointResult, inverter_control_fixed_point_oracle
export InverterControlNetworkFixedPointResult
export solve_inverter_control_network_fixed_point
export ControlledInverterFleetNetworkFixedPointResult
export solve_controlled_inverter_fleet_network_fixed_point
export controlled_inverter_network_fixed_point_rows
export ControlledInverterFleetMultiStartResult
export solve_controlled_inverter_fleet_multistart
export controlled_inverter_fleet_multistart_rows
export inverter_control_current_jacobian, inverter_control_loop_gain
export InverterControlNetworkSensitivityResult
export inverter_control_network_voltage_sensitivity
export ControlledInverterFleetNetworkSensitivityResult
export controlled_inverter_fleet_network_voltage_sensitivity
export controlled_inverter_fleet_network_sensitivity_rows
export controlled_inverter_fleet_loop_gain

# Current–voltage (IVQ) battery storage + chemistry library
export BatteryChemistry, thevenin_chemistry, linear_chemistry, tabulated_chemistry,
       illustrative_lfp, illustrative_nmc, illustrative_nca,
       illustrative_lead_acid, illustrative_leaf
export IVQBattery, solve_ivq_battery, IVQBatteryResult
export solve_multiperiod_ivq, MultiperiodIVQResult

# HELM power flow (holomorphic embedding load-flow, a bespoke solution method)
export helm_series, HelmResult, solve_pf_helm

# Static operability continuation and local fold evidence
export OperabilityContinuationSpec, OperabilityContinuationResult,
       continue_opf_operability, OperabilityPseudoArclengthSpec,
       continue_opf_operability_pseudo_arclength, OperabilityFoldResult,
       locate_opf_operability_fold, operability_continuation_margin,
       operability_continuation_rows

end # module PowerOptLab
