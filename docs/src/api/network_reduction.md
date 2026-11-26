# Neutral reduction reference and API

The explicit-neutral (conductor) Kron reduction now lives in BMOPFTools as
`kron_reduce_neutrals`, and PowerOptLab re-exports it. It eliminates each
bus's neutral terminal from a four-wire BMOPF case, holding the neutral at
earth potential: every line's series impedance becomes the Schur complement
`Zᵣ = Zₚₚ − Zₚₙ Zₙₙ⁻¹ Zₙₚ`, and every change is recorded in
`net["_meta"]["kron_reduction"]`.

Use it only when the assumed neutral grounding suits the study. The reduction
is exact only where the neutral is already perfectly earthed at both ends of a
line; elsewhere it forces the neutral to ideal earth and reports where. The
[power-flow and operability overview](../algorithms/index.md) places model
reduction alongside the available solution and diagnostic methods.

The arguments are unchanged:
`kron_reduce_neutrals(net; neutral_terminals=nothing, as_json=false, include_coordinates=false)`.

```julia
reduced = kron_reduce_neutrals(net)
json = kron_reduce_neutrals(net; as_json=true)
```

The full reference covers the reduction rules, refusals, provenance, and the
OpenDSS cross-check. BMOPFTools also compares this reduction with network
simplification and error-bounded bus reduction (Opti-KRON). See:

- [`kron_reduce_neutrals` in the BMOPFTools API reference](https://frederikgeth.github.io/BMOPFTools.jl/dev/api/#BMOPFTools.kron_reduce_neutrals)
- [Network reduction: which tool, and what it preserves](https://frederikgeth.github.io/BMOPFTools.jl/dev/network_reduction/)

## Deprecated names

`kron_reduce_bmopf` and `reduce_bmopf_neutrals` remain as deprecated aliases of
`kron_reduce_neutrals` and will be removed in a future release.

```@docs
kron_reduce_bmopf
reduce_bmopf_neutrals
```
