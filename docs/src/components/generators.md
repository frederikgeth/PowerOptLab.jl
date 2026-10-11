# Generators

PowerOptLab's generalized generator models represent source voltage behind
impedance, connection topology, operating controls and capability limits.
[`GeneralizedGenerator`](@ref) describes oriented electrical ports;
[`SourceGenerator`](@ref) adds an explicit source neutral and grounding path.
These are research prototypes for fundamental-frequency network studies.

## Choose a starting point

| Your task | Guide |
|:--|:--|
| Understand what each modeling choice changes | [Worked examples](../tutorials/generalized_generator_models.md) |
| Derive connections, source laws, grounding and power balances | [Scientific model](generalized_generator.md) |
| Compare source assumptions on the same unbalanced feeder | [Trade-offs under unbalance](../tutorials/generalized_generator_tradeoffs.md) |
| Check independent circuit calculations and equality rank | [Circuit and rank validation](../tutorials/generalized_generator_validation.md) |
| Construct a device and interpret its results | [Implementation and API](generalized_generator_api.md) |
| Assess precedent and the limits of the evidence | [Literature evidence](generalized_generator_evidence.md) |
| Review proposed upstream fields and serialization | [Specification drafts](../specifications/generators/index.md) |

Start with the worked examples, then use the scientific model and API for your
chosen connection and source law. The specification drafts describe proposed
upstream conventions; consult their status before treating them as an existing
data interface.

For switching capability, converter losses, filter circuits or DC-link limits,
see [Inverter-based resources](../ibr/index.md). A generalized voltage-source
model does not supply those converter mechanisms.
