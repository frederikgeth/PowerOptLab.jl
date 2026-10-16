# Contributing

PowerOptLab develops experimental methods for four-wire model inference,
informative interventions, and verified network decisions. Stable foundational
physics may later fold back into the
[BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl) spec.

Before adding a capability, answer two separate questions:

1. **Research fit:** which inference, experiment-design, decision-verification,
   or solution-validity question does it enable?
2. **Implementation kind:** is it a component model, problem specification, or
   solution method?

A standalone device or generic algorithm is not a sufficient contribution by
itself. Its documentation should name the research question, the evidence or
decision it connects to, and the validity limits of the result. Once that fit is
clear, the [kind of contribution](concepts.md) determines where the code and docs
go and which model-building functions it uses.

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

Register the new page once under the study topic and device family it supports.
The [device modeling and scheduling overview](components/index.md) links the
generator, storage and EV families. See
[PE battery storage](components/devices.md) as the reference example.

## Adding a problem specification

A new *formulation* over the same physics. Put it in `src/problems/`.

1. Build the network snapshot(s) via the staged API — `build_opf_model(add_objective=false)`,
   then compose your own objective and constraints, `optimize_opf!(ctxs...)` once,
   and `extract_result`.
2. Return a bespoke result struct.
3. Add a page under `docs/src/problems/` (badge line + worked example), register in
   `docs/make.jl`, export from `PowerOptLab.jl`.

See [State estimation](problems/state_estimation.md) for the *inverse*-problem
pattern and [Multi-period OPF](problems/multiperiod.md) for inter-temporal linking.

## Adding a solution method

A new *solution method*. Put it in `src/algorithms/` and add a page under
`docs/src/algorithms/`. Register it with the study question it answers; the
[power-flow and operability chapter](algorithms/index.md) groups the current
HELM and static-continuation methods.

## Documentation navigation and API ownership

The published chapters follow study topics, while source folders describe
implementation layers. Add a guide or tutorial to its study chapter in
`docs/make.jl`; avoid separate global lists of components, problems, algorithms
or tutorials. Register each published page once.

Place canonical `@docs` entries in the relevant topic API reference. The
[shared API](api.md) describes result status, time intervals and extension interfaces, and links to
the topic references. Use `@ref` to link a symbol from guides and examples,
rather than reproducing its docstring. The documentation build checks exported
symbols and cross-references, and the [symbol index](symbol_index.md) spans all
chapters.

## Executable study tutorials

Run `julia --project=. scripts/run_study_tutorials.jl` after changing the
battery comparison, constrained state estimation, three focused inverter
studies, HELM comparison or post-OPF operability tutorial. The runner executes
each page's plain `julia` blocks in order in a fresh module, with the repository
root as its working directory. CI runs these checks on Julia 1.10 and current
stable; missing pages, missing executable markers, runtime errors and failed
assertions fail the check. The DOE tutorials have their own runner,
`scripts/run_doe_tutorials.jl`.

Keep setup data, expected results and scientific assertions together in the
published page. Check physical balances, limits and inference claims as well as
solver status. Label synthetic data, explain numerical acceptance budgets, and
update the expected-result tables when a justified model change alters them.
Do not maintain a separate copy of tutorial code for testing. Mark schematic
API guides clearly and link to a complete worked example.

### Regenerate the battery and state-estimation figures

The plots use solved data exported from the tutorial modules. Install the
plotting dependency in a Python 3.12 environment, then run from the repository
root:

```sh
python -m pip install -r scripts/requirements-study-figures.txt
julia --project=. scripts/run_study_tutorials.jl --figure-data docs/build/study-figure-data
python scripts/render_study_figures.py
python scripts/render_study_figures.py --check
```

Commit the SVGs and plotting-data snapshot in `docs/src/assets/studies/` together
with any renderer changes. Current-stable CI compares fresh solved quantities
with that snapshot (absolute tolerances of 0.1 W, 0.0001 V or normalized
residual, and 0.00001 SoC fraction for the corresponding arrays),
checks renderer/input fingerprints and exercises Matplotlib. This avoids
pixel comparisons that depend on platform rendering. The SVGs have no generation
timestamp and can be used directly in a paper or presentation. Use
`--preview-dir docs/build/study-figure-preview` to produce PNGs for visual review.

## Writing for power-engineering readers

Start with the network problem, the assumptions and the quantities being
calculated. Explain the equations, units and sign conventions before describing
the software implementation. Use standard power-system terminology and short,
direct sentences.

Prefer concrete descriptions: add currents to KCL, construct a four-wire line
model, compare with an independent reference calculation, or check equation
residuals. Reserve terms such as callback, ownership and compilation for sections
that explain the programming interface. Keep public function names unchanged.

State the conditions under which a result holds. Distinguish numerical
convergence, satisfaction of the model equations and validation against physical
data. Give limitations where they affect interpretation, and refer to that
explanation elsewhere rather than repeating a general warning on every page.

## Badge line convention

Each capability page opens with a one-line metadata block so readers can place it
at a glance:

```markdown
> **Kind:** Problem specification · **Maturity:** promotion candidate ·
> **Direction:** inverse · **Temporal:** single-snapshot
```

- **Kind** — Component model / Problem specification / Bespoke algorithm.
- **Maturity** — *proof of concept*, *prototype*, *experimental*, *research prototype*, or
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
- Keep branch plans, PR readiness gates and implementation checklists in
  `dev/notes/`, outside the published documentation. Preserve useful equations,
  implemented scope and validation requirements in the user-facing reference
  guides before relocating a workflow record. Current behavior belongs in the
  reference guide; a development note is not an API contract.

## BMOPFTools compatibility

BMOPFTools is unregistered. The root and docs `Project.toml` source tables pin
commit `df6bc42bfd64fc5475658c438adcdb28a3469777`; package compatibility remains
0.1.0. Update both source pins together, then run the complete tests and
documentation build on the Julia compatibility floor and current stable.

Run `julia --project=. scripts/instantiate_pinned.jl` for local tests and
`julia --project=docs scripts/instantiate_pinned.jl` for documentation. CI uses
these commands too. The script checks agreement between the two declarations
and explicitly installs the commit, including on Julia 1.10 where `[sources]`
is not consumed automatically. A fresh Julia 1.10 documentation environment
is bootstrapped with the local PowerOptLab path before resolving the registered
documentation dependencies.

Documentation CI builds on Julia 1.10 and current stable; only current stable
deploys. To build locally without deployment, use:

```sh
DOCUMENTER_DEPLOY=false julia --project=docs docs/make.jl
```

Both manifest files remain local. A `Pkg.develop` path override means tests
use that checkout's working tree. Re-run the setup script to restore the tested
commit before diagnosing a discrepancy with CI. Record `versioninfo()`,
`Pkg.status()`, the Ipopt version, and explicit solver tolerances with published
results.

Most integrations use BMOPFTools' public staged-model and admittance APIs. HELM
and static operability also need the engine's parsed load decomposition,
which the pinned BMOPFTools revision does not expose publicly. Those imports are isolated in
`src/upstream.jl`; they are the only supported private compatibility adapter.
The upstream API needed to remove it is a public, read-only load-decomposition
function returning each sub-load's terminal pair, complex power, and shunt
admittance. Until BMOPFTools provides that interface, an upstream upgrade must treat
changes to the adapter's imported names as breaking.
