# DOE development roadmap

> **Status date:** 5 September 2026 · **Purpose:** implementation boundary and
> next research increments; scientific claims and citations live in the
> [literature evidence register](doe_literature_evidence.md)

PowerOptLab treats DOE quantification as a framework for comparing declared
geometries, network models, uncertainty sets, information structures, and
evidence levels. It does not prescribe one operational policy, and it retains
pointwise perfect recourse for reproducing anticipative formulations.

## Current research capability

| Axis | Implemented | Boundary |
|---|---|---|
| Geometry | One-sided active-power allocation per connection | No paired import/export band or coupled P–Q region |
| Network model | Nonlinear, unbalanced, neutral-explicit AC model | Fundamental-frequency steady state; local nonlinear solve |
| Utilization coverage | Bound point, exact small-N corners, bound-point dropout faces, explicit points, scrambled Halton and margin-directed adaptive search | Finite search; no general continuous nonlinear box certificate. Refinement reach depends on the finite move budget and minimum-step floor; dropout faces are seeded explicitly |
| Exogenous uncertainty | Typed scenarios, Gaussian and box-conditioned draws, residual bootstrap, categorical materialization seam | Model validity and physical transforms remain caller responsibilities |
| Control recourse | Issue, scenario, local-law, and context stages for supported controls | Arbitrary controls need registration; generator P/Q lacks a safe linked handle |
| Validation | Tri-state numerical evidence, policy-compatible context verdicts, provenance-matched fixed-control replay, candidate confirmation, multistart, analytic references | General AC replay shares the formulation/solver family; the analytic containment proof applies only to its idealized fixture |
| Statistical evidence | Held-out coverage, leakage audits, capacity curves, model pairing, shift and probability-calibration diagnostics | Assumptions such as independence, exchangeability, and stationarity are declared, not inferred |
| Reproducibility | Study hashes, manifests, benchmark rows, eight executable tutorials, committed synthetic fixtures | Real datasets and estimator-specific materializers remain study-specific |

The cleanup adds [faithful replay](../tutorials/doe_faithful_replay.md) and an
[analytic reference suite](../tutorials/doe_analytic_reference.md). One physical
IBR may bind only one participant; terminal declarations must match the engine,
and legacy ports are single-phase. Numerical limits remain unresolved rather
than becoming candidate violations. The separate PF adapter rejects droop
profiles whose rating bases the pinned engine would alter.

## Branch-completion gate

The current DOE research branch is ready for review when:

- the focused DOE tests and documentation build pass;
- the committed benchmark writes stable-schema rows from a clean checkout;
- each tutorial states question, model, runtime, data, evidence, and limitations;
- perfect recourse is explicit and remains available for replication;
- finite scenario or utilization evidence is never labelled a certificate; and
- the PR reports known unrelated test failures separately.

## Next research branch 1: published coverage reproduction

Reproduce the Liu–Braslavsky 2.91 kW counterexample on the documented Australian
feeder, subject to data availability and licensing. The contribution is not the
known counterexample; it is a reproducible baseline against which nonlinear
search, recourse policies, and four-wire variants can be compared.

Required outputs:

- exact network/data provenance and transformation;
- reported bound point and the published customer reductions;
- AC replay of voltage magnitude and location;
- solver-start and tolerance sensitivity; and
- a clear comparison between published and reproduced quantities.

## Next research branch 2: four-wire sensitivity filtering

Re-derive or numerically characterize utilization sensitivities with explicit
neutral impedance and grounding. Apply the published nonlinear sensitivity-
filtering procedure, then search omitted vertices and interior points.

The research question is whether a reduced set remains sufficient under
declared network and operating conditions. Do not assume that Kron reduction
simply removes coupling: it can create dense effective impedances while still
omitting neutral voltage and grounding states.

Report missed violations, sensitivity sign changes, neutral/grounding
dependence, filtering threshold sensitivity, and computational scaling. Random
utilization samples are supporting evidence, not a dominance certificate.

## Next research branch 3: certified linear baseline

Implement a fixed-polyhedron robust decoupled-box baseline using support
functions or the Motzkin transposition construction. Its purposes are to:

- provide a genuinely certified result for the declared linear model;
- separate containment conservatism from AC-model fidelity error;
- benchmark corner enumeration and sensitivity screening; and
- quantify the capacity gap between a certificate and nonlinear search.

Any bilinear expansion may affect global maximality, but every returned box
claimed as certified must have independently checked containment constraints.

## Later research branches

### Alternative geometries

Add paired import/export bounds, then coupled P–Q sets and feeder/aggregator
budgets. Compare capacity, communication burden, privacy, disaggregation
responsibility, and ex-post violation evidence.

### Intertemporal envelopes

Start with finite-horizon feasible trajectories for storage energy, EV
departure requirements, ramping, tap wear, and transformer temperature. Use
viability or controlled-invariance language only when the study horizon and
operational question require it.

### Fairness and value

Define the stakeholder, benefit or burden, entitlement, and time horizon before
choosing an objective. Compare offered capacity, realized energy, curtailment,
customer value, and reliability rather than treating one interval's Jain index
as a complete fairness result.

### Model-fidelity ladder

Compare linear, relaxed, and learned candidates against the nonlinear AC oracle
with false-safe and false-conservative rates. Add harmonics, protection,
dynamics, or controller equilibria only when the omitted physics changes the
decision under study.

## Documentation maintenance

Every new research increment should update four places:

1. the [DOE start page](../tutorials/doe_getting_started.md), only if the core
   learning path changes;
2. the [problem specification](operating_envelope.md), for stable semantics;
3. the [literature evidence register](doe_literature_evidence.md), for new
   scientific evidence; and
4. this roadmap, for implementation status and the next acceptance gate.

Detailed historical reasoning remains in the [DOE quantification scientific
audit](doe_quantification_review.md). API additions belong in the [DOE API
reference](../api/operating_envelope.md), not in introductory tutorials.
