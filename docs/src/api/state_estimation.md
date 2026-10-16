# State estimation API

Select an estimator using the [state estimation overview](../estimation/index.md)
and [comparison](../estimation/comparison.md). The
[legacy WLS](../problems/state_estimation.md) and
[constrained NLLS](../problems/constrained_state_estimation.md) guides describe
the network assumptions, measurement equations and convergence criteria.

## Legacy WLS estimator

```@docs
Measurement
solve_state_estimation
StateEstimationResult
```

## Compiled constrained estimator

```@docs
TerminalConnection
ExactInjectionSpecification
ConstantPowerDevice
ConstantCurrentDevice
ZIPDevice
ExactDeviceEquation
BranchMeasurement
NetInjectionMeasurement
TransformerMeasurement
ShuntMeasurement
FunctionalMeasurement
SEPreflightFinding
SEPreflightReport
SEUnsupportedNetwork
state_estimator_preflight
initial_state_estimator
injection_state_indices
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
IpoptStateEstimationResult
solve_ipopt_state_estimator
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
