# Contributing

PowerOptLab is a staging ground: new work lands here first, and anything that
matures folds back into the [BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl)
spec. Before adding something, decide which of the three
[kinds of contribution](concepts.md) it is — that determines where the code and
docs go and which engine seam you reuse.

## Adding a component model

A new network *element*. Put it in `src/components/`.

1. Define a struct describing the element (ports, ratings, physics parameters).
2. Stamp it into a solve with a `model_hook!`: add its per-phase current injection
   to the engine's KCL accumulators, plus any extra variables/constraints (e.g. a
   state-of-charge state).
3. Read results back with a `solution_hook!`.
4. Add a page under `docs/src/components/` with the badge line and a worked
   example; register it in `docs/make.jl` and export public names from
   `PowerOptLab.jl`.

See [Storage & EVs](components/devices.md) as the reference example.

## Adding a problem specification

A new *formulation* over the same physics. Put it in `src/problems/`.

1. Build the network snapshot(s) via the staged API — `build_opf_model(add_objective=false)`,
   then compose your own objective and constraints, `enforce_kcl!` per snapshot,
   and `extract_result`.
2. Return a bespoke result struct.
3. Add a page under `docs/src/problems/` (badge line + worked example), register in
   `docs/make.jl`, export from `PowerOptLab.jl`.

See [State estimation](problems/state_estimation.md) for the *inverse*-problem
pattern and [Multi-period OPF](problems/multiperiod.md) for inter-temporal linking.

## Adding a bespoke algorithm

A new *solution method*. Put it in `src/algorithms/` and add a page under
`docs/src/algorithms/`. See [Bespoke algorithms](algorithms/index.md).

## Badge line convention

Each capability page opens with a one-line metadata block so readers can place it
at a glance:

```markdown
> **Kind:** Problem specification · **Maturity:** promotion candidate ·
> **Direction:** inverse · **Temporal:** single-snapshot
```

- **Kind** — Component model / Problem specification / Bespoke algorithm.
- **Maturity** — *prototype* (expect churn) or *promotion candidate* (ready to
  fold into the BMOPF spec).
- **Direction** — *forward* (dispatch/OPF) or *inverse* (estimation).
- **Temporal** — *single-snapshot* or *inter-temporal*.

## The promotion path

When a contribution stabilises — a settled interface, tests, and a worked example
— it becomes a *promotion candidate*. Promotion means proposing the model or
formulation for the BMOPF spec upstream, at which point the engine may absorb it
and the PowerOptLab version becomes a thin re-export or is retired.

## Conventions

- SI at the interface; per-unit only inside a solve (via `ctx.bases`).
- `Manifest.toml` is intentionally not committed (library convention).
- Every feature is opt-in and covered by a test under `test/`.
