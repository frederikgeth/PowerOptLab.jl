# Battery storage API

Start with the [battery storage overview](battery_storage.md) for model selection,
units and the distinction between energy and charge-based SoC.

## Power-energy storage

See the [PE model](devices.md) for equations and usage, and
[Multi-period OPF](../problems/multiperiod.md) for scheduling and solve results.
The [scheduling API](../api/multiperiod.md) documents the shared solver and result
type used by PE storage and EV sessions.

```@docs
StorageDevice
```

## Current-voltage-charge storage

See the [IVQ model](ivq_battery.md) for the cell/converter boundary, chemistry
assumptions and numerical limits.

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
