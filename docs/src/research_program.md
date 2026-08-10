# Research program

This page distinguishes a coherent research direction from the functionality
that PowerOptLab exposes today. It is a prioritization guide, not an API promise.

The target workflow is:

```text
telemetry and metadata
        ↓
model forensics → informative, safe interventions
        ↓                    ↓
plausible model ensemble ← new evidence
        ↓
verified operating decisions
```

The defining requirement is to preserve consequential ambiguity. A good result
does not force one network model when several models explain the data. It reports
what is unidentified, whether the alternatives change a downstream decision,
and what additional evidence could separate them.

## Current boundary

The repository contains useful foundations, but not the integrated workflow.

| Area | Implemented now | Important boundary |
|---|---|---|
| Constrained state estimation | Neutral-explicit state, exact equations versus stochastic residuals, tangent-space observability, selected local covariance, constraint multipliers | Network identity and topology are fixed at compile time; there is no topology-hypothesis search |
| Operational parameter estimation | Shared line lengths and transformer taps over multiple snapshots | The line prototype is single-terminal per end; it does not jointly infer full mutually coupled impedances, phase maps, switch states, or grounding models |
| Inverse Carson reconstruction | Enumeration of overhead-construction candidates, compatible-candidate retention, local rank and profile intervals, materialized primitive four-wire models | It uses sequence linecode data, not joint operational telemetry, and does not combine candidates with the state estimator |
| Operating envelopes | Shared allocation across explicit forecast/model scenarios; bound-point, corner, and custom-utilization checks | Scenario models must be supplied by the caller; corner checks are not a global robust-feasibility certificate |
| Differentiated controller experiment | Local DiffOpt sensitivity of smooth lower-level PV/controller responses | Derivatives are conditional on the returned local solution and stable active set; discrete modes and alternative solution branches are not covered |
| HELM power flow | No-load-connected branch with physical residual and series diagnostics | It does not discover every power-flow branch or certify non-existence when its finite series fails |

These boundaries are deliberate. Terms such as *compatible*, *locally
observable*, *corner-secure*, and *converged* must retain the precise meanings
given on the corresponding capability pages.

## Priority 1: joint model forensics

The first major capability should unify the existing inverse tools around an
explicit model-hypothesis representation. The uncertain model may include
switch and neutral continuity, phase labels, transformer terminal maps and taps,
line construction or mutually coupled impedance, source impedance, grounding,
and meter sign, scale, or reference errors.

A useful result would expose:

- ranked discrete hypotheses and compatible continuous parameter regions;
- locally unidentified directions and observational equivalence classes;
- groups that are different physically but equivalent for a declared decision;
- model discrepancies that affect the intended operating decision; and
- candidate evidence that could distinguish material alternatives.

This is substantially more than widening the current parameter-estimation
vector. It needs explicit discrete hypotheses, shared snapshot data, consistent
four-wire circuit compilation, and result semantics that never turn one local
optimum into a claim of unique model identification.

## Priority 2: differentiable experiment design

Once model alternatives are explicit, PowerOptLab should be able to choose safe
additional evidence: meter placement or synchronization, inverter active/reactive
probing, a tap or capacitor action, or a suitable operating interval.

Candidate objectives include increasing the smallest identifiable singular
value, reducing a selected parameter or decision covariance, and separating
discrete hypotheses. Every intervention must also satisfy the same voltage,
current, converter, and disturbance limits used for operation.

The existing residual Jacobians, constrained tangent-space diagnostics, staged
four-wire model, DiffOpt proof of concept, and HELM coefficient recursion are
ingredients. They are not yet an experiment-design solver. In particular,
nonconvex sensitivities need the solution-validity checks described below.

## Priority 3: counterexample-guided decision verification

The operating-envelope implementation already states the central gap correctly:
feasibility at the simultaneous bound, or even at every box corner, does not
prove feasibility throughout a nonconvex AC utilization box.

The next step is an outer allocation and inner counterexample search over both
actual customer utilization and the plausible model ensemble. Each discovered
violation becomes a new scenario for the allocation problem. Results should
distinguish at least:

- a falsification search that found a violating point;
- a search-stable result for which no violation was found under a declared
  multistart and stopping policy; and
- a certificate backed by a valid global relaxation bound.

Only the last case is a robust-feasibility certificate. Naming this hierarchy is
important even before all levels are implemented.

## Foundational track: solution validity and differentiability

The inference and intervention program depends on local nonlinear solves. A
parallel foundational track should therefore add pseudo-arclength continuation,
multiple-solution discovery, active-set transition tracking, KKT regularity and
conditioning diagnostics, and finite-difference or continuation checks of
DiffOpt sensitivities.

The scientific output is not just a derivative. It should state which solution
branch was differentiated, whether the local solution is isolated, where the
derivative becomes ill-conditioned or discontinuous, and whether another local
response changes the decision.

## Later, question-driven extensions

These directions fit after the core workflow has a concrete use for them:

- **Multiconductor relaxations and global bounds:** develop SDP, SOCP, or QC
  relaxations when they can certify a forensics or operating-envelope claim;
  compare them against the neutral-explicit IVR model.
- **Hybrid controller equilibria:** represent discrete taps, switched
  capacitors, deadbands, hysteresis, and interacting inverter controls when the
  question is reachability, multiplicity, or cycling—not merely optimal dispatch.
- **Fundamental-plus-harmonic optimization:** add frequency-indexed primitives
  when harmonic voltage, neutral heating, resonance, or emission uncertainty
  changes a hosting-capacity or operating-envelope decision.

Moving-horizon estimation, decision-focused forecasting, safe learned local
controllers, decomposition, and protection-aware envelopes are also plausible
extensions when they close a demonstrated gap in this loop.

## Deliberate non-priorities

PowerOptLab should not accumulate generic OPF-learning surrogates, generic
scenario wrappers, transmission-market models without a multiconductor research
question, a full EMT simulator, or standalone device models with no inference,
intervention, verification, or solution-validity use case. Data-driven methods
belong when their topology, grounding, uncertainty, and model-mismatch semantics
are explicit.
