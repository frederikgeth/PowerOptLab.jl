# Device models and scheduling

The device model determines which electrical limits and stored states appear in
the optimisation. The guides below cover generators, stationary batteries and
EV charging. The [multi-period OPF guide](../problems/multiperiod.md) shows how
storage and EV energy balances link successive operating points.

| Family | Study question | Start here |
|:--|:--|:--|
| Generators | Study the effects of source impedance, voltage regulation, winding connection and grounding. | [Generator overview](generators.md) |
| Inverters | Compare forward current control with voltage-dependent operating constraints. | [Constraint-based inverter operation](inverter_operating_constraints.md) |
| Battery storage | Compare power-energy scheduling with a model of cell voltage and current limits. | [Battery storage overview](battery_storage.md) |
| Electric vehicles | Calculate charging schedules subject to outlet ratings, availability and departure energy. | [Electric vehicle overview](../ev/index.md) |

The [scheduling API](../api/multiperiod.md) describes the multi-period solver.
The generator and IVQ models have separate solvers, listed in their API
references.

The [inverter chapter](../ibr/index.md) covers converter circuits, losses,
ratings and local control laws for studies in which these determine the feasible
operating region.

## Model assembly

Storage and EV devices contribute currents to the KCL equations at each operating
point and energy-balance constraints between intervals. They implement
[`AbstractDevice`](@ref) using [`validate_device`](@ref),
[`stamp_device!`](@ref), [`link_device!`](@ref) and [`extract_device`](@ref).
[`solve_multiperiod_opf`](@ref) assembles and solves the complete schedule.

See [Concepts](../concepts.md) for units, sign conventions and interpretation of
solver status, and the [shared API](../api.md) for the common programming
interfaces. The [neutral reduction reference](../api/network_reduction.md)
describes the assumptions and changes involved in eliminating a grounded neutral.
