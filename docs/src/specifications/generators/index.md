# Generator specification drafts

These pages stage mathematical and JSON data-model specifications for the
`GeneralizedGenerator` and `SourceGenerator` components. They are proposals for
later contributions to the IEEE PES Task Force specification, **not an adopted
schema or a JSON reader implemented by PowerOptLab**. The executable interface
remains the [Julia API](../../components/generalized_generator_api.md).

The drafts follow the upstream component structure: **1. Data model; 2. Input
symbols; 3. Variables; 4. Equality constraints; 5. Inequality constraints**, followed
by implementation notes and examples. The electrical primitive is specified once
for each component. Voltage laws and capabilities are composable physical
assumptions; they are not competing numerical formulations of the same device.

- [Generalized generators](generator.md): independently connected source ports,
  series impedance, terminal injection, and power conservation.
- [Source generators](source.md): explicit internal star, conductor impedance,
  independent neutral and earth returns, and complete power accounting.
- [Shared voltage laws, controls, and capabilities](operating-model.md): every
  field and equation shared by the two components, including modes a–f.
- [Implementation and upstream migration](integration.md): rectangular equations,
  exact limits, API mapping, upstream inconsistencies, and proposed PR boundaries.

## Reference revision and attribution

Reviewed on 6 September 2026 against upstream commit
[`73fae2b6bae2663d9a2e901c41a4c062457bf834`](https://github.com/distribution-system-opt/math-and-data-model-specifications/tree/73fae2b6bae2663d9a2e901c41a4c062457bf834).
The source repository identifies its current content as work in progress without
a validated release. This commit is an audit reference, not a released standard.

The organization, field-table conventions, notation, and overlapping baseline
model definitions are adapted from the **IEEE PES Task Force on Benchmarking
Multiconductor OPF for Distribution Systems and contributors**, under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). This draft adds source
impedance, voltage-law selection, grounding, controls, capabilities, and migration
analysis; the Task Force has not endorsed these changes. Preserve this attribution
when moving the adapted material upstream or elsewhere.

The review included these revision-pinned pages:

| Upstream source | Convention or finding used here |
|:--|:--|
| [Generators](https://github.com/distribution-system-opt/math-and-data-model-specifications/blob/73fae2b6bae2663d9a2e901c41a4c062457bf834/docs/src/spec/generator.md) | Wye dispatch, phase P/Q/S bounds, conductor current bounds, per-phase cost, implicit neutral return |
| [Voltage sources](https://github.com/distribution-system-opt/math-and-data-model-specifications/blob/73fae2b6bae2663d9a2e901c41a4c062457bf834/docs/src/spec/source.md) | Fixed terminal voltage, grounded source-bus neutral, unrestricted slack current, exactly-one-source scope |
| [Notation](https://github.com/distribution-system-opt/math-and-data-model-specifications/blob/73fae2b6bae2663d9a2e901c41a4c062457bf834/docs/src/spec/notation.md) | Colored symbol definitions, complex foundational equations, rectangular implementation, 1/3 Fortescue transform |
| [Data formatting](https://github.com/distribution-system-opt/math-and-data-model-specifications/blob/73fae2b6bae2663d9a2e901c41a4c062457bf834/docs/src/spec/data-format.md) | SI, radians, string IDs, ordered vectors, real-field complex encoding, row-first matrix entries, omission semantics |
| [Grounding](https://github.com/distribution-system-opt/math-and-data-model-specifications/blob/73fae2b6bae2663d9a2e901c41a4c062457bf834/docs/src/spec/grounding.md) and [buses](https://github.com/distribution-system-opt/math-and-data-model-specifications/blob/73fae2b6bae2663d9a2e901c41a4c062457bf834/docs/src/spec/bus.md) | Common zero-potential earth, explicit neutral, bus grounding, source-fixed terminal assumptions |
| [Objective](https://github.com/distribution-system-opt/math-and-data-model-specifications/blob/73fae2b6bae2663d9a2e901c41a4c062457bf834/docs/src/spec/objective.md) | Injection-positive power, currency/kWh, W-to-kW conversion, duration weighting |
| [Contributing](https://github.com/distribution-system-opt/math-and-data-model-specifications/blob/73fae2b6bae2663d9a2e901c41a4c062457bf834/docs/src/contributing.md) | Separate foundational physics from solver choices; normative review and coordinated schema changes |

## What studying the generator specification changes

The existing `generator` is a dispatchable wye injection. Its `p_min`, `q_max`, and
`s_max` arrays refer to phases; `i_max` refers to conductors, potentially including
the neutral. In the current PowerOptLab API, by contrast, `GeneratorCapability.i_max`
refers to **ports**, while `terminal_i_max` refers to conductors. The draft retains
the upstream meaning of JSON `i_max` and introduces `i_port_max` for port currents.

The existing `cost` is a required **per-phase vector**. PowerOptLab currently offers
a scalar price on **complete PCC power**. The draft names this `cost_total` and
requires it explicitly, including a literal zero for an unpriced source. It does
not reinterpret the upstream `cost` field. Equal phase prices reproduce a total
price for a closed-return generator; this equivalence needs an additional neutral/
earth power condition for a grounded source.

The existing source is an ideal terminal-voltage reference. A fixed **internal**
phasor behind finite impedance does not fix the PCC voltages. Neither that
assumption nor the dedicated earth return can be introduced as an editorial
change to the source page. They require changes to grounding, bus-limit application,
source/reference scope, notation, objective semantics, and the schema.

## Draft envelope and field conventions

For staging, records are keyed by string IDs in two new, provisional top-level
objects: `generalized_generator` and `source_generator`. Existing `generator` and
`voltage_source` objects retain their meaning. Whether upstream instead chooses
explicit subtypes of those objects is a review decision; an importer must never
silently reinterpret an old record as a new one.

All numeric data are finite real JSON numbers. Phasors are RMS fundamental-frequency
quantities. Units are fixed by the field: V, A, Ω, W, var, VA, rad, and currency/kWh.
No per-unit bases or solver initialization fields belong to this proposed physical
data model. Numerical frequency-dependent parameters describe one declared study
frequency; there is no frequency state or automatic frequency rescaling.

A missing optional bound adds no constraint. A missing optional impedance entry
is zero. Required enum fields have no implicit default. `null`, NaN, infinity,
complex literals, and strings containing units are not numeric encodings here.
Vectors have their declared full length; there is no implicit scalar broadcast.
An omitted whole bound vector is distinct from a vector containing zeros.

Matrix entry ``(p,q)`` is encoded by `R_series_p_q` and `X_series_p_q`, **1-indexed,
row first**, in ohms. Dimension comes from the component's ports or conductors,
not the largest supplied index. Entries outside that dimension are invalid.
There is no implicit symmetry: a nonzero reciprocal pair must be given in both
positions. A zero matrix, an exactly singular nonzero matrix, and a small nonzero
matrix have their literal meanings.

Arrays of distinct terminal strings select terminals of the host bus in their
listed order. No phase or neutral is inferred from a label. The common ground
reference is not an entry of a port or conductor matrix; the reserved upstream
name `"g"` must not be used as a component terminal. Grounding is declared by the
source's grounding model or by the network's bus/shunt data.

## Status of the examples

The JSON blocks are complete **component-collection fragments**, not complete
network cases and not documents validated against the published upstream schema.
Their hosts and terminals must exist in an assembled case. They have an explicit
mapping to the current Julia constructors described in [Integration](integration.md).
The example values are illustrative model assumptions, not equipment ratings.
