# Shared interfaces and API index

This reference describes solution status, time intervals and the interfaces
used to add devices and measurements. The links below lead to the functions
and input parameters for each study. The
[symbol index](symbol_index.md) provides alphabetical lookup across all chapters.

## Module

```@docs
PowerOptLab
```

## Results and time intervals

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

[Generator models](components/generalized_generator_api.md),
[battery storage](components/battery_storage_api.md) and
[electric vehicles](ev/api.md) each have their own API reference.

## Multi-period OPF

See the [multi-period scheduling API](api/multiperiod.md).

## Legacy WLS state estimation

See the [state estimation API](api/state_estimation.md#Legacy-WLS-estimator).

## Constrained NLLS state estimation

See the [state estimation API](api/state_estimation.md#Compiled-constrained-estimator).

## Parameter estimation

See the [parameter calibration API](api/parameter_estimation.md).

## Inverse Carson reconstruction

See the [line reconstruction API](api/inverse_carson.md).

## Bilevel distribution-network proof of concept

See the [bilevel PV and tap control API](api/bilevel.md).

## Advanced inverter

See the [advanced inverter API](api/advanced_inverter.md).

## Phase-aware inverter controls

See the [phase-aware control API](api/inverter_controls.md),
[fleet study reference](ibr/network_control_studies.md) and
[closed-loop evidence reference](ibr/closed_loop_evidence.md).

## Current–voltage (IVQ) battery

See the [battery storage API](components/battery_storage_api.md).

## HELM power flow

See the [HELM API](api/helm.md).

## Explicit-neutral Kron reduction

See the [neutral reduction reference and API](api/network_reduction.md).

## Post-OPF operability

See the [operability API](api/operability.md).

## Numerical formulations

See the [function formulation API](formulations/api.md).
