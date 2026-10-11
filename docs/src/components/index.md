# Device modeling and scheduling

Choose a component family from the physical device and the question your study
asks. Each family brings its models, worked examples, validation or literature
evidence, and API together. For storage and EV schedules, the shared
[multi-period OPF guide](../problems/multiperiod.md) explains the chronological
solve, duration weighting and result contract.

| Family | Study question | Start here |
|:--|:--|:--|
| Generators | How do source impedance, voltage regulation, connection and grounding affect a network decision? | [Generator overview](generators.md) |
| Battery storage | Is a power-energy schedule sufficient, or do cell voltage and current limits determine dispatch? | [Battery storage overview](battery_storage.md) |
| Electric vehicles | Can a vehicle meet its departure requirement through the available charging equipment? | [Electric vehicle overview](../ev/index.md) |

Use the [scheduling API](../api/multiperiod.md) for the shared solver and result
types, and the family APIs for the devices passed to that solver. Generalized
generators and IVQ batteries use dedicated solve entry points.

Converter physics and local controls have a dedicated
[Inverter-based resources](../ibr/index.md) chapter. A generator's voltage law,
a battery's stored state, and an inverter's conversion constraints describe
different parts of a device; use the linked guides when those parts interact.

## Shared interfaces

The storage and EV scheduling devices implement [`AbstractDevice`](@ref) through
[`validate_device`](@ref), [`stamp_device!`](@ref), [`link_device!`](@ref), and
[`extract_device`](@ref). [`solve_multiperiod_opf`](@ref) uses this lifecycle to
stamp each snapshot and link stored energy across intervals. Generator and IVQ
studies use their dedicated solve entry points; do not assume every component
can be passed to the same builder.

See [Concepts](../concepts.md) for units, engine extension seams and solve-status
contracts, and the [shared API](../api.md) for result, time-grid and extension
interfaces. Neutral-grounding model transformations are documented in the
[neutral reduction reference](../api/network_reduction.md).
