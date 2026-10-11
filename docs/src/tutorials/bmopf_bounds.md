# Deriving and filling BMOPF bounds

`propagate_bounds` provides the first solver-free slice of the
[bound propagation design](../problems/bmopf_bound_propagation.md). It derives
outer enclosures from native BMOPFTools physics and declared operational limits,
records how each endpoint was obtained, and reports equations it could not use.
It does not run optimization-based bound tightening yet.

```julia
using PowerOptLab, BMOPFTools, JSON3

net = parse_bmopf("network.bmopf.json")
report = propagate_bounds(net; max_steps=100_000, time_limit=30.0)

# Inspect a phase-to-earth voltage, in volts. Replace the ID/terminal with yours.
key = QuantityKey(:bus, "bus1", :voltage, ("a",))
domain = bound(report, key)
domain.lower, domain.upper
explain_bound(report, key; endpoint=:upper)
coverage(report)

# Both operations leave the input and report unchanged.
filled, patch = materialize_bounds(net, report)
write_bmopf(filled, "network.with_bounds.bmopf.json")
open("network.bounds.json", "w") do io
    JSON3.pretty(io, bounds_sidecar(report))
end
```

For a runnable small feeder example, see `examples/bmopf_bounds.jl`. A report with
`status == :conflict` cannot be materialized. A non-conflicting report is not a
proof of AC feasibility, and its box is not an operating region in which every
point is feasible.

## Current coverage

| Family | Implemented | Boundary |
|---|---|---|
| Bus | Terminal real/imaginary/magnitude domains, neutral and phase-pair voltage bounds, perfect grounding, KCL | Sequence and angle constraints are relaxed. Native source/ground exemptions are preserved. |
| Line | Full mutually coupled matrix KVL, inline or linecode coefficients, length scaling, both π-shunts, both terminal current/power ratings | Unresolved impedance models retain unknown injections. Series current and terminal current are different quantities. |
| Switch | Fixed open/closed state, mapped voltage equality, current and apparent-power limits | No topology optimization. Ideal cycles are retained and may leave unbounded current domains. |
| Load | Wye, delta and single-phase connection currents and power identities; constant P/Q | Voltage-dependent laws are relaxed, with nominal power left unconstrained. |
| Generator | Connection P/Q boxes, current and apparent-power limits, supported neutral-current sums | Native generator model only; no PowerOptLab extension-device hooks. |
| Voltage source | Magnitude reference, exact zero-angle components, source grounding and P/Q boxes | Nonzero angles are relaxed. No approximate sine/cosine value is used as an exact bound. |
| Shunt | Full terminal admittance and current/power domains | Uses terminal-matrix semantics. |
| Transformer, IBR, capacitor, other AC devices | Object inventory and free AC injections at referenced buses | Device physics and internal quantities are not yet supported. Unknown attachment layouts retain free injections at all AC buses. |
| DC-only objects | Unsupported-object inventory | No direct AC port; converter effects are covered by the unsupported IBR's free AC injections. |
| Linecode, geometry, wire data, control profile, time series | Dependency inventory | These are input data, not independent operating variables. |

`coverage(report)` lists `:supported`, `:relaxed`, `:unsupported`, `:data` or
`:metadata` for each object. Support describes compiled relations; it does not
mean every endpoint is finite. In particular, unresolved bounds are not proofs
of unboundedness. An unsupported device is never replaced by zero injection.

## Quantity keys and units

`QuantityKey(family, id, quantity, index; component=:magnitude)` uses string
tuples for ports and terminal labels. Real and imaginary components use
`:real`/`:imag`. Power components use `:active`/`:reactive`; their units are W and
var, while power magnitude is VA. Voltages are V and currents are A.

| Key family/quantity | Index convention | Meaning |
|---|---|---|
| `:bus`, `:voltage` | `(terminal,)` | Terminal-to-earth voltage. |
| `:bus`, `:connection_voltage` | `(positive, negative)` | Difference of the two bus terminal voltages; pairs follow declaration order. |
| `:bus`, `:ground_current` | `(terminal,)` | Current entering the terminal from earth. |
| `:line`, `:current` | `("series", conductor_position)` | Series current directed from the from-side to the to-side. Position is a string such as `"1"`. |
| `:line`, `:current` or `:power` | `("from", terminal)` / `("to", terminal)` | Current/power entering the line at that terminal, including its shunt. Power uses terminal-to-earth voltage. |
| `:line`, `:shunt_current` | `(side, terminal)` | Current leaving the bus into that end's shunt. |
| `:switch`, `:current` or `:power` | `("from", terminal)` | From-side current/power entering the switch. |
| `:load`, `:generator`, `:voltage_source` | `("coil", coil_position)` | `:connection_voltage`, `:current`, or `:power`; current is load consumption or generation/source injection. Coil order follows the native connection map. |
| `:generator`, `:current` | `("neutral", terminal)` | Return-current sum when a separate neutral rating is supplied. |
| `:shunt`, `:current` or `:power` | `(terminal,)` | Current/power entering the shunt. |
| `:unsupported`, `:current` | `(bus, terminal)` | Unconstrained injection retained for an unsupported object; its ID is the object's JSON pointer. |

Use `keys(report.domains)` to discover the registered quantities. Grounding or
an explicit source can establish finite phase-to-earth domains. Phase-to-neutral
bounds alone do not establish finite neutral-to-earth bounds in a floating system.

## Numerical and model contract

Input is copied and normalized through BMOPFTools, including terminal aliases and
migrations. Schema validation runs by default when the input declares a
recognized schema. `report.validation` records whether it ran; normalization is
not a full physical-model validation. Malformed supported terminal references,
array sizes and numerical coefficients are rejected. Conflicting domains are
returned with endpoint evidence rather than silently corrected.

The kernel interprets normalized binary64 numbers as exact real inputs, and
forms its supported algebraic relations from them. Addition, multiplication,
division and coefficient arithmetic use exact rational comparisons to round
outwards. Square-root endpoints are checked by exact squaring. This does not
certify unknown physical parameters, the original decimal-to-binary conversion,
or an independently rounded solver model. There is no artificial voltage floor,
large default current cap, or local-solver bound.

`max_steps` and `time_limit` cap propagation **after compilation**. `rtol` and
`atol` skip small improvements; they never move an endpoint inward. Default
`atol=0` avoids one absolute tolerance across different SI quantities. Status is
`:fixed_point`, `:tolerance_stopped`, `:budget_exhausted`, or `:conflict`.
All accepted domains remain valid after a budget stop. These contractors can be
weak, especially without angle information or in an ungrounded network.

Evidence is recorded separately for each endpoint and includes parent evidence
IDs, a rule and source pointers into the normalized input. Some relation source
pointers identify a whole device or a derived grounding/KCL relation rather than
an original scalar field. Natural nonnegative-magnitude domains have their own
evidence. `explain_bound` returns the endpoint's derivation history in dependency
order. The sidecar schema is `data/schema/bounds.schema.json`; JSON `null` with
`status="unresolved"` represents an infinite endpoint.

## Scope and export

Use `BoundScope(snapshot=t)` to resolve native time-series references with
`BMOPFTools.get_snapshot`. Bounds apply only to that snapshot. Control-profile
time-series semantics remain unsupported along with IBR controls.
`materialize_bounds` returns the static snapshot, never a network with those
bounds attached unconditionally to every time step. Source/profile changes
invalidate an earlier report, even if it happened to produce the same bounds.

`materialize_bounds(net, report; mode=:fill_missing)` fills absent bus voltage
and line current/power limits. `mode=:intersect` can also tighten existing fields.
Line bounds cover both ends and respect inherited linecode ratings. Existing
input limits are never widened. A vector field is written only if all its entries
are finite and its phase entries are active in the native model. Fields involving
source-fixed or grounded phases are left unchanged. Partially bounded vectors and all other quantities stay in the
sidecar. Finite computed limits are redundant constraints for this snapshot,
not inferred equipment nameplates or permissible ratings for another study.

The returned dictionary is normalized and can contain BMOPFTools' internal
metadata and derived terminal fields. Use `BMOPFTools.write_bmopf` for an
exchange JSON file; it strips internal fields and retains supported provenance.
The accompanying patch records changed fields, original and inherited values,
quantity keys and evidence IDs. Keep the sidecar with the exported network.

The next increments are dedicated transformer/IBR/DC and voltage-dependent-load
contractors, a public quantity-to-JuMP adapter, and selective OBBT over a valid
outer relaxation. Those APIs from the design proposal are not implemented yet.

## API

```@docs
BoundScope
QuantityKey
BoundDomain
BoundEvidence
BoundReport
propagate_bounds
bound
coverage
explain_bound
bounds_sidecar
materialize_bounds
```
