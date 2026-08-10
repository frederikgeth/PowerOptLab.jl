# Concepts

PowerOptLab studies **four-wire distribution-network decisions under state and
model uncertainty**. It is built on the
[BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl) reference
current–voltage OPF engine and reuses the engine's neutral-explicit device
physics, per-unit handling, and result extraction through public extension
seams, without forking the engine.

The scientific organizing principle is an evidence-to-decision loop:

1. **Model evidence:** compile telemetry, metadata, and exact circuit equations.
2. **Forensics:** retain the states, parameter regions, or discrete models that
   the evidence cannot distinguish.
3. **Experiment design:** seek safe measurements or interventions that separate
   alternatives which matter to a decision.
4. **Decision verification:** test an operating decision across the remaining
   model and forecast alternatives.
5. **Revision:** return new evidence to the inference stage.

Only parts of this loop exist today. The constrained estimator, parameter
estimator, and inverse-Carson solver provide distinct forms of model evidence;
the operating-envelope solver consumes explicitly supplied network scenarios;
the bilevel and HELM work provide local sensitivity and solution diagnostics.
There is not yet a joint model-hypothesis API, an active-probing optimizer, or a
global robust-feasibility certificate. The [Research program](research_program.md)
keeps those proposed capabilities separate from the callable API.

The source tree is organized by implementation layer. That engineering taxonomy
supports the research loop; it is not the project's scientific identity.

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

- **Maturity** — *proof of concept* (demonstrates a narrow idea), *prototype*
  (test-backed but experimental), *research prototype* (intended for declared
  studies rather than operational use), or *promotion candidate* (stable enough
  to propose for the BMOPF spec).
- **Direction** — *forward* (given parameters, find an operating point:
  dispatch/OPF), *inverse* (given measurements, find the state or parameters),
  or *hierarchical* (different decision-makers or response levels).
- **Temporal structure** — *single-snapshot* vs *inter-temporal* (state-of-charge
  coupling, shared parameters, or per-interval envelopes spanning many snapshots).

- **Uncertainty/security treatment** — deterministic, explicit scenarios,
  tested utilization points, or globally certified. These labels describe
  different claim strengths and must not be used interchangeably.

## Interface conventions

Everything is **SI at the interface** (watts, vars, watt-hours, volts); per-unit
conditioning inside each solve is handled through the engine's
`opf_bases(ctx)` accessor. Where it applies, `per_unit=true` and
`per_unit=false` give identical results.

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

See [Contributing](contributing.md) for how to add each kind of contribution,
the research-fit test, and the path back to the BMOPF spec.
