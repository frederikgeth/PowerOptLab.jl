# Generators

The generator models describe an internal voltage source behind impedance,
winding connections, voltage and power controls, and operating limits.
[`GeneralizedGenerator`](@ref) describes oriented electrical ports;
[`SourceGenerator`](@ref) adds an explicit source neutral and grounding path.
These are research prototypes for fundamental-frequency network studies.

## Models and examples

| Your task | Guide |
|:--|:--|
| Compare generator modelling assumptions | [Worked examples](../tutorials/generalized_generator_models.md) |
| Derive connections, source laws, grounding and power balances | [Scientific model](generalized_generator.md) |
| Compare source assumptions on the same unbalanced feeder | [Trade-offs under unbalance](../tutorials/generalized_generator_tradeoffs.md) |
| Check independent circuit calculations and equality rank | [Circuit and rank validation](../tutorials/generalized_generator_validation.md) |
| Construct a device and interpret its results | [Implementation and API](generalized_generator_api.md) |
| Review the supporting literature and validation | [Literature evidence](generalized_generator_evidence.md) |
| Review the proposed data format and JSON representation | [Specification drafts](../specifications/generators/index.md) |

Start with the worked examples, then use the scientific model and API for your
chosen winding connection and voltage-source assumptions. The specification drafts propose
a data format for upstream adoption. PowerOptLab already implements the experimental
[local schema and JSON reader/writer](../specifications/generators/serialization.md);
these records are not an adopted upstream data interface.

For a converter whose switching limits, losses, filter or DC link affect the
study, use the models in [Inverter-based resources](../ibr/index.md).
