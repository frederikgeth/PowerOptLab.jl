# Reading and writing generator data

PowerOptLab implements the flat records in these specification drafts with a
local JSON Schema and semantic validation. The format remains experimental;
it is not a released IEEE PES Task Force or `dsopt-schema` format. Native
BMOPFTools network records keep their existing meaning. The generator extension
is read separately, then its components enter the public staged engine.

## 1. An executable delta example

This envelope describes three unequal winding P/Q injections. Voltages are
line-line RMS values, `i_port_max` limits winding current, and `i_max` limits
external line current. There is no neutral conductor on this component.

```@example generator_data
using PowerOptLab
text = """
{
  "generalized_generator": {
    "roof": {
      "bus": "poc",
      "terminal_map": ["a", "b", "c"],
      "configuration": "DELTA",
      "voltage_model": "NONE",
      "power_setpoint_location": "POC",
      "p_set": [600, 300, 900],
      "q_set": [50, 100, 150],
      "i_port_max": [10, 10, 10],
      "i_max": [20, 20, 20],
      "cost_total": 0.15
    }
  }
}
"""
data = read_generator_data(text; from_string=true)
d = only(data.devices)
@assert d.connections == [("a", "b"), ("b", "c"), ("c", "a")]
@assert d.control.p == [600.0, 300.0, 900.0]
(data.identifiers, generator_data(d)["configuration"])
```

No impedance entries were supplied, so ``Z=0`` exactly. The source law is
unconstrained (`NONE`), and the explicit winding P/Q setpoints make this a PQ
model. Incidence automatically enforces terminal voltage closure and zero total
external current. The three winding currents can still have a nonzero mean.

For a network `net` parsed by BMOPFTools, the solve sequence is:

```julia
data = read_generator_data("generators.json"; net)
result = solve_generator_opf(net, data.devices; per_unit=true, s_base=1e4)
write_generator_data("generators-canonical.json", data)
```

Supply `net` to check bus and terminal membership immediately. Without it, the
importer validates the component against a temporary terminal inventory and
leaves network membership for model construction. It does not infer bus grounds,
source replacement ownership, or voltage references. Those remain explicit
network/build choices.

## 2. Two validation layers

[`generator_data_schema`](@ref) returns the bundled
`data/schema/generators.schema.json` as a fresh dictionary. This Draft 7 schema
validates object structure, allowed fields, types, enumerations, required and
forbidden mode fields, mutually exclusive aggregate/per-port representations,
fixed sequence-vector length, and nonnegative engineering magnitudes. All
references are local; reading the format does not fetch remote schemas.

The semantic translator additionally checks:

- Vector lengths derived from the actual connection and conductor inventory.
- Full impedance dimensions, finite entries, and positive-semidefinite Hermitian
  part; omitted entries are exactly zero, not reflected across the diagonal.
- Port orientation, conductor participation, supported forests/delta cycles,
  ideal-loop template consistency, and positive angle-law magnitude domains.
- Bound ordering, common-magnitude interval intersection, sequence applicability,
  and hard-PV control compatibility.
- Grounding semantics and, when `net` is supplied, network terminal membership
  and duplicate ideal source/engine grounding.

A schema cannot establish passivity, solve a magnitude-space nullspace, or prove
network feasibility. Passing validation means the data describes a supported
component; it does not guarantee an operating point or a globally optimal solve.
JSON numbers are converted to `Float64` before entering the typed API. Infinities,
NaNs, null bounds, partial bound vectors, and unknown fields are rejected.

## 3. Ordering, identities and canonicalization

Collection IDs are strings and unique within their family. Envelope imports use
`generalized_generator:<id>` or `source_generator:<id>` as runtime IDs for every
record. Prefixing all records is injective even when an original ID itself
contains a colon. The `GeneratorDataSet.identifiers` dictionary retains the
original `(family, id)`, and dataset export restores those collection keys.
Single-record [`generator_from_data`](@ref) keeps the explicit ID supplied by its
caller. No automatic migration of native generators or sources is performed.

`PORTS` can list conductors in an order different from the runtime's order of
outgoing terminals followed by unseen returns. Import permutes the conductor
`i_max` vector into runtime order. Port controls, port ratings and port impedance
keep their declared winding order. Canonical export writes the runtime conductor
order and corresponding `port_map`, preserving the physical incidence relation.
Result phasors use `GeneratorResult.terminals`; they must be matched by labels
when comparing to the original data order. This codec exports component data,
not a separate result-file schema.

Canonical export expands per-entry scalar ratings to full vectors, writes angle
defaults explicitly, emits full nonzero matrix entries, and normalizes zero
source-ground impedance to `IDEAL`. Common-magnitude lower bounds collapse to
the greatest input lower bound and upper bounds to the least input upper bound;
export repeats that equivalent interval for each port. Fixed phasors are written
as magnitude and principal angle. Fixed-template magnitude consistency permits
only a floating-point roundoff allowance (eight machine epsilons times the
greater of 1 V and template magnitude), so an exact boundary survives conversion. Unused location enums may become explicit
`POC`. An empty sequence overlay can disappear because it imposes no bounds.
These are physical round trips, not promises to preserve JSON text or spelling.
`voltage_scale` is a numerical option and is deliberately excluded from the data.

`write_generator_data` validates the full envelope before opening its destination,
so a validation failure leaves an existing file intact. Passing a vector of
components preserves their IDs within their respective families; passing a parsed
dataset restores its original IDs.

## 4. Verification and upstream boundary

Tests cover all five source laws for both component families under unbalance,
P/Q measurement locations, every capability field's translation, delta and
single/split-phase connections, conductor permutations, original-ID preservation,
file round trips and malformed records. The solve comparisons check physical
voltages, power and loss accounting; where controls leave multiple valid
solutions, tests do not require an arbitrary optimizer tie-break to be preserved.

The mathematical pages remain suitable as staged upstream proposals. The local
schema must be reviewed and adapted alongside the upstream schema version and
collection/discriminator decisions before it is claimed as conforming there.

## API

```@docs
GeneratorDataSet
generator_data_schema
generator_from_data
read_generator_data
generator_data
write_generator_data
```
