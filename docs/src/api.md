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

## State estimation

```@docs
Measurement
solve_state_estimation
StateEstimationResult
```

## Parameter estimation

```@docs
CalibLine
CalibTap
solve_parameter_estimation
ParameterEstimationResult
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
