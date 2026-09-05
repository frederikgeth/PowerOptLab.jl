# Contributing

PowerOptLab develops experimental methods for four-wire model inference,
informative interventions, and verified network decisions. Stable foundational
physics may later fold back into the
[BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl) spec.

Before adding a capability, answer two separate questions:

1. **Research fit:** which inference, experiment-design, decision-verification,
   or solution-validity question does it enable?
2. **Implementation kind:** is it a component model, problem specification, or
   bespoke algorithm?

A standalone device or generic algorithm is not a sufficient contribution by
itself. Its documentation should name the research question, the evidence or
decision it connects to, and the validity limits of the result. Once that fit is
clear, the [kind of contribution](concepts.md) determines where the code and docs
go and which engine seam it reuses.

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
- **Maturity** — *proof of concept*, *prototype*, *research prototype*, or
  *promotion candidate*. Use the definitions in [Concepts](concepts.md).
- **Direction** — *forward* (dispatch/OPF), *inverse* (estimation), or
  *hierarchical* (multiple decision/response levels).
- **Temporal** — *single-snapshot* or *inter-temporal*.

## The promotion path

When a contribution stabilises — a settled interface, tests, and a worked example
— it becomes a *promotion candidate*. Promotion means proposing the model or
formulation for the BMOPF spec upstream, at which point the engine may absorb it
and the PowerOptLab version becomes a thin re-export or is retired.

## Conventions

- SI at the interface; per-unit only inside a solve (via `opf_bases(ctx)`).
- `Manifest.toml` is intentionally not committed (library convention).
- Every feature is opt-in and covered by a test under `test/`.

## BMOPFTools compatibility contract

BMOPFTools is unregistered. The root and docs `Project.toml` source tables pin
commit `5b51d2f361dab91bd7c16711019584407da79ed8`; package compatibility remains
0.1.0. Update both source pins together, then run the complete tests and
documentation build on the Julia compatibility floor and current stable.

Run `julia --project=. scripts/instantiate_pinned.jl` for local tests and
`julia --project=docs scripts/instantiate_pinned.jl` for documentation. CI uses
these commands too. The script checks agreement between the two declarations
and explicitly installs the commit, including on Julia 1.10 where `[sources]`
is not consumed automatically.

Both manifest files remain local. A `Pkg.develop` path override means tests
use that checkout's working tree. Re-run the setup script to restore the tested
commit before diagnosing a discrepancy with CI. Record `versioninfo()`,
`Pkg.status()`, the Ipopt version, and explicit solver tolerances with published
results.

Most integrations use BMOPFTools' public staged-model and admittance APIs. HELM
also needs the engine's parsed constant-power and constant-impedance sub-loads,
which BMOPFTools 0.1.0 does not expose publicly. Those imports are isolated in
`src/upstream.jl`; they are the only supported private compatibility adapter.
The upstream API needed to remove it is a public, read-only load-decomposition
function returning each sub-load's terminal pair, complex power, and shunt
admittance. Until BMOPFTools provides that seam, an upstream upgrade must treat
changes to the adapter's imported names as breaking.
