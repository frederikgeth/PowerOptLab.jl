# API reference

All exported names, grouped by capability.

## Module

```@docs
PowerOptLab
```

## Devices

```@docs
StorageDevice
EVDevice
```

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

## Dynamic operating envelopes

```@docs
ConnectionPoint
solve_operating_envelope
OperatingEnvelopeResult
```

## Advanced inverter

```@docs
AdvancedInverter
solve_advanced_inverter
InverterResult
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

## Index

```@index
```
