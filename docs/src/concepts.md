# Concepts

PowerOptLab is a **staging ground** for experimental power-system modelling built
on top of the [BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl)
reference optimal-power-flow engine. BMOPFTools ships one small, correct four-wire
rectangular current–voltage OPF and deliberately refuses to become a "model zoo";
PowerOptLab is where the zoo lives. Everything here reuses the engine's device
physics, per-unit handling, and result extraction **through its public extension
seams** — without forking the engine. Anything that matures into accepted practice
can later be folded back into the BMOPF spec.

## The three kinds of contribution

Every capability in PowerOptLab is one of three structurally different things,
distinguished by *which layer of the engine it extends*. The docs and the source
tree (`src/components/`, `src/problems/`, `src/algorithms/`) are organised the
same way.

| Kind | What it contributes | Seam it reuses | Source |
|---|---|---|---|
| **Component model** | a new network *element* (device physics) | `model_hook!` / `solution_hook!` | `src/components/` |
| **Problem specification** | a new *objective + variable/constraint structure* over the whole network | staged `build_opf_model` / `enforce_kcl!` / `generation_cost` / `extract_result` | `src/problems/` |
| **Bespoke algorithm** | a new *solution method* — how you solve, not what | custom solve loops around the staged API | `src/algorithms/` |

- A **component model** answers *"what new thing can sit on the network?"* It is
  stamped into a snapshot as a per-phase current injection added to the engine's
  Kirchhoff-current-law accumulators, optionally with extra variables (e.g. a
  state-of-charge state linking snapshots).
- A **problem specification** answers *"what question do we ask of the same
  physics?"* It builds one or more snapshots into a JuMP model via the staged API,
  swaps in its own objective and constraints, and reads back a bespoke result. It
  changes the *formulation*, not the elements.
- A **bespoke algorithm** answers *"how do we solve it?"* — decomposition,
  sequential linearization, warm-start schemes, or any custom loop around the
  staged API. HELM is the first: it uses the engine's augmented admittance
  matrix in a custom power-series solve.

## Cross-cutting dimensions

The three kinds above are the primary axis, but each capability also carries a few
orthogonal attributes, surfaced as a badge line at the top of its page:

- **Maturity** — *prototype* (experimental, expect churn) vs *promotion
  candidate* (stable enough to fold back into the BMOPF spec). This is the whole
  point of a staging ground: knowing what is ready to graduate.
- **Direction** — *forward* (given parameters, find an operating point:
  dispatch/OPF) vs *inverse* (given measurements, find the state or parameters:
  state and parameter estimation).
- **Temporal structure** — *single-snapshot* vs *inter-temporal* (state-of-charge
  coupling, shared parameters, or per-interval envelopes spanning many snapshots).

- **Uncertainty/security treatment** — deterministic, multi-scenario, robust,
  or explicitly corner-secure. Operating envelopes already exercise this axis.

## Interface conventions

Everything is **SI at the interface** (watts, vars, watt-hours, volts); per-unit
conditioning inside each solve is handled by the engine's `ctx.bases`. Where it
applies, `per_unit=true` and `per_unit=false` give identical results.

Cross-cutting code uses four shared contracts:

- [`AbstractDevice`](@ref) with validation, stamping, temporal linking, and
  extraction methods;
- [`TimeGrid`](@ref), which gives every snapshot an explicit positive duration;
- [`AbstractMeasurement`](@ref) with common value/uncertainty accessors;
- [`build_multi_context`](@ref), which builds ordered snapshots into one JuMP
  model, plus [`solve_status`](@ref) / [`solve_diagnostics`](@ref) for a stable
  result-status view.

JuMP-backed results retain the exact normalized [`SolveStatus`](@ref) produced
at solve time, including whether a primal candidate existed and whether it was
strictly publishable. Compatibility fields such as `termination_status` remain,
but new code should branch on `solve_status(result).publishable`.

See [Contributing](contributing.md) for how to add each kind of contribution and
the path back to the BMOPF spec.
