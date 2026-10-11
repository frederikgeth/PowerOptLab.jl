# Battery storage

PowerOptLab provides two storage abstractions. Choose the least detailed one that
retains the mechanism relevant to your study and that your data can support.

| Model | Stored state | Physical boundary | Entry point |
|:--|:--|:--|:--|
| Power-energy (PE), [`StorageDevice`](@ref) | Energy in Wh | AC port with fixed one-way efficiencies and supplied power/current ratings | [`solve_multiperiod_opf`](@ref) |
| Current-voltage-charge (IVQ), [`IVQBattery`](@ref) | Charge-based, dimensionless SoC | Cell voltage/current model coupled to an [`AdvancedInverter`](@ref) | [`solve_ivq_battery`](@ref), [`solve_multiperiod_ivq`](@ref) |

PE storage is a promotion candidate; IVQ is a prototype. PE's `soc` result is a
compatibility alias for energy in Wh. It is not the dimensionless SoC used by IVQ.

## Choose a starting point

| Your task | Guide |
|:--|:--|
| Choose model fidelity, data and a defensible comparison | [PE versus IVQ](../tutorials/battery_storage_models.md) |
| Schedule storage with energy balance and AC port limits | [Power-energy model](devices.md) |
| Study cell-limited power and charge-linked dispatch | [Current-voltage-charge model](ivq_battery.md) |
| Look up constructors, chemistry helpers and solve results | [Battery storage API](battery_storage_api.md) |

The model comparison is the recommended first read. The model pages provide
equations, usage and numerical limits; the API lists the callable interface.
For PE scheduling context, see [Multi-period OPF](../problems/multiperiod.md).

## Connections to other families

[Electric vehicles](../ev/index.md) reuse the PE energy account but add vehicle,
equipment and departure requirements. The [operating-mode guide](../ev/modes.md)
also applies to bidirectional PE storage: positive efficiencies alone do not
exclude simultaneous charging and discharging. The
[compatibility guide](../ev/migration.md) records shared PE port changes.

IVQ uses the converter physics documented in
[Inverter-based resources](../ibr/index.md). Its current implementation couples
battery and converter at DC power; read the IVQ model's coupling limits before
interpreting voltage-dependent converter capability.
