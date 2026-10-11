# Battery storage

PowerOptLab provides a power–energy (PE) storage model and a
current–voltage–charge (IVQ) model. PE is suitable for energy scheduling when
fixed efficiencies and power ratings are adequate. IVQ represents cell voltage
and current explicitly and requires the corresponding cell data.

| Model | Stored state | Electrical representation | Entry point |
|:--|:--|:--|:--|
| Power-energy (PE), [`StorageDevice`](@ref) | Energy in Wh | AC port with fixed one-way efficiencies and supplied power/current ratings | [`solve_multiperiod_opf`](@ref) |
| Current-voltage-charge (IVQ), [`IVQBattery`](@ref) | Charge-based, dimensionless SoC | Cell voltage/current model coupled to an [`AdvancedInverter`](@ref) | [`solve_ivq_battery`](@ref), [`solve_multiperiod_ivq`](@ref) |

PE storage is a promotion candidate; IVQ is a prototype. PE's `soc` result is a
compatibility alias for energy in Wh. It is not the dimensionless SoC used by IVQ.

## Models and examples

| Your task | Guide |
|:--|:--|
| Select a model and the data needed for comparison | [PE versus IVQ](../tutorials/battery_storage_models.md) |
| Schedule storage with energy balance and AC port limits | [Power-energy model](devices.md) |
| Study cell-limited power and charge-linked dispatch | [Current-voltage-charge model](ivq_battery.md) |
| Look up constructors, chemistry helpers and solve results | [Battery storage API](battery_storage_api.md) |

The PE–IVQ comparison explains the assumptions, data requirements and losses
in each model. The individual guides give the equations and examples.
[Multi-period OPF](../problems/multiperiod.md) describes the PE scheduling problem.

## EV charging and converter models

[Electric vehicles](../ev/index.md) reuse the PE energy account but add vehicle,
equipment and departure requirements. The [operating-mode guide](../ev/modes.md)
also applies to bidirectional PE storage: positive efficiencies alone do not
exclude simultaneous charging and discharging. The
[compatibility guide](../ev/migration.md) records shared PE port changes.

IVQ uses the converter physics documented in
[Inverter-based resources](../ibr/index.md). Its current implementation couples
battery and converter through DC power. The IVQ guide explains the resulting
limitations when converter ratings depend on the battery terminal voltage.
