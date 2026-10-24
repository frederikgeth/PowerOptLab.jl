# Neutral reduction reference and API

Use this transformation when the assumed neutral grounding is appropriate to
the study. The [power-flow and operability overview](../algorithms/index.md)
places model reduction alongside the available solution and diagnostic methods.

```@docs
kron_reduce_bmopf
```

`kron_reduce_bmopf(net; neutral_terminals=nothing, as_json=false, include_coordinates=false)` returns a
new BMOPF dictionary (or schema-valid JSON when `as_json=true`) and never
mutates `net`. A JSON string is accepted as the first argument. Neutral labels
come from `terminal_conventions.neutral` when that block is present; otherwise
the conservative `n`/`N` naming convention is used. A
`Dict("bus_id" => "label")` `neutral_terminals` override is available for
datasets whose labels are not declared.

Dictionary output preserves bus coordinates. JSON output omits the non-spec
`longitude`/`latitude` fields by default to retain the strict-schema contract.
Use `as_json=true, include_coordinates=true` for an explicit coordinate-preserving
extension; that output is not strict-schema JSON when coordinates are present.
The serialization choice applies both when neutrals are reduced and when the
input already has no explicit neutral; it does not change the electrical model.

For each affected line, the series impedance is reduced exactly with
`Zᵣ = Zₚₚ − Zₚₙ Zₙₙ⁻¹ Zₙₚ`, with the explicit neutral voltage constrained to
zero. Inline matrices and referenced linecodes are supported, including
line lengths, line π-shunt projections, shared-linecode cloning when retained
conductor orders differ, and conductor-indexed current ratings. The output
removes neutral bus names and AC terminal-map entries across lines, sources,
loads, generators, switches, shunts, capacitors, IBRs, and transformer
windings. Grounded WYE/SINGLE_PHASE capacitors become equivalent phase-to-ground
shunts. Phase-to-neutral voltage bounds are intersected with phase-to-ground
bounds; `vn_max` and neutral-conductor current limits are dropped with an audit
entry because the neutral is now ideal ground.

| Field/family | Reduction rule | Audit or refusal |
|---|---|---|
| `bus.terminal_names`, `perfectly_grounded_terminals` | Remove the resolved neutral; filter grounds by value | Classifies each bus as exact or forced-ground projection |
| `bus.v_min/v_max`, `vpn_min/vpn_max` | Select terminal-indexed phase entries and intersect bounds | Contradictory intersections error; `vn_max` is dropped |
| Line `R/X_series_*` | Exact `Zpp − Zpn Znn⁻¹ Znp` Schur complement | Stores `K = −Znn⁻¹Znp`, conductor maps, and neutral π-row coefficients |
| Line `G/B_*`, shunts | Remove grounded-neutral rows/columns; neutral-only ground shunts are removed | Neutral current ratings are explicitly reported as lost |
| WYE/SINGLE_PHASE capacitors | Materialize equivalent phase-to-ground `shunt` | Time-varying capacitor conversion is refused |
| Loads, sources, generators, switches | Remove neutral map entries and slice only conductor-indexed vectors | Ambiguous role/order errors; phase arrays are preserved |
| IBR and transformers | Active neutral-leg IBRs and neutral-bearing transformers are refused | No invalid shortened BMOPF subtype is emitted |

DC buses, DC terminal maps, and DC grounding are intentionally outside the
transformation. Structural, matrix, or bound time-series targets are refused;
ordinary snapshot data remains untouched. The reduced output is BMOPF JSON
schema-valid when serialized with `as_json=true, include_coordinates=false`; upstream 0.1.0 may still
report the known `WYE` arity warning for a grounded WYE whose implicit midpoint
is no longer listed.

The transformation is refused for singular or near-singular neutral blocks,
ambiguous maps, DELTA elements containing a neutral, active IBR
`FOUR_LEG` neutral physics, unsupported center-tap/single-phase transformer
reductions, and time-varying structural/impedance data whose commutation is not
proven. Matrix keys follow BMOPFTools' sparse convention: absent entries are
zero, a missing reciprocal is mirrored, and explicitly stored full entries take
precedence. Every reduction records original/kept conductor ordering, the
neutral-current recovery vector `K = −Zₙₙ⁻¹ Zₙₚ`, bound rewrites, dropped
constraints, and whether a bus was already perfectly grounded or was forced to
ideal ground in `net["_meta"]["kron_reduction"]` and
`net["extras"]["kron_reduction"]`.

The package test target requires OpenDSSDirect.jl and runs an independent
comparison: an explicit four-wire node-0-ground deck and its manually
projected three-wire deck. It asserts convergence and compares retained phase
voltages and line currents at both ends, source P/Q, total losses, and the
neutral-current recovery relation `Iₙ = K Iₚ` against the reduced BMOPF solve.
The older bundled finite-0.001-ohm grounding fixture remains useful as a
near-ground comparison, but is not the test of exact equivalence.

```julia
reduced = kron_reduce_bmopf(net)
json = kron_reduce_bmopf(net; as_json=true)
```
