# DOE implementation and documentation workflow

> Historical developer workflow record, excluded from the documentation build.
> Extracted from [doe_quantification_review.md](../../docs/src/problems/doe_quantification_review.md).
> Original audit review date: 2026-09-05. Check the published reference guides
> for current capabilities and evidence.

## Control-recourse implementation sequence

1. **Implemented:** add discovery and classification while retaining historical
   context recourse as the compatibility default.
2. **Implemented:** add `:issue` and `:scenario` linking for transformer taps
   and free IBR P/Q handles; prescribed IBR control profiles are classified as
   `:local_law`, not linked setpoints.
3. **Implemented first slice:** a paired STATCOM test shows ideal recourse
   admitting more capacity than one issue-time setting; verification repeats
   each feasible context with optimized free controls fixed.
4. **Partly implemented:** `DOEControlRegistration` lets extensions declare
   stable semantic handles, units, classifications, laws, and provenance.
   Generator P/Q still needs a public upstream power handle; currents are not
   used as a surrogate.
5. After migration warnings and compatibility tests, require an explicit policy
   for multi-context studies with free controls. Keep both operational and
   perfect-recourse presets stable and tested.
6. **Implemented for current finite methods:** the same policy is used in
   corner checks, explicit utilization allocation, custom verification, Halton
   screening, and counterexample-guided allocation. Eventual stochastic/robust
   formulations must retain this path.

## Recommended development program

### Stage 0 — trustworthy experimental record

Deliver this before claiming a new quantification method:

1. **Implemented:** `DOEStudySpec` records interval/scenario network hashes,
   utilization coverage, control information structure, objective,
   normalization, solver/options, random seeds, metadata, and software versions.
   Typed scenario and uncertainty-construction records cover Gaussian,
   box-conditioned Gaussian, and empirical residual-bootstrap inputs; automated
   estimator-specific materializers and external dataset adapters remain.
2. **Implemented first slice:** return per-context status, margins, identifiers,
   snapshots, controls, joint timing, and independent replay timing/evidence.
3. **Partly implemented:** independently evaluate JuMP primal residuals and
   repeat each context with free controls fixed. A second power-flow engine or
   formulation is still needed for stronger independent physics replay.
4. **Implemented first slice:** deterministic registered-variable multistart
   records the start policy, every status, accepted runs, and capacity spread.
5. **Implemented first slice:** `doe_benchmark_rows`,
   `doe_context_benchmark_rows`, and `scripts/run_doe_benchmark.jl` emit stable
   TSV-ready evidence without mixing empirical runs into unit tests.
6. **Partly implemented:** explicit extension controls can register stable
   semantic handles and fail closed. Automatically detecting an undeclared
   custom free control without confusing it with a state variable remains an
   upstream registration problem.

Acceptance gate: a third party can recreate every table row from one committed
configuration and can identify exactly what was and was not tested.

### Stage 1 — search-stable range verification

1. **Adaptive heuristic slice implemented:** deterministic Halton coverage can
   now be followed by coordinate refinement around points with the least
   normalized voltage, current, or unbalance headroom. The accumulated set is
   jointly verified under the selected recourse policy and the search retains
   every score and round. A true continuous inner adversarial optimization and
   a globally valid violation bound remain future work.
2. **Implemented for finite Halton and adaptive screens:** separate outer loops
   solve, screen, append the complete tested set, and repeat. The adaptive loop
   replays the issue/scenario controls actually selected by each allocation and
   retains every intermediate allocation. Replace both screens with a true
   violation-maximizing oracle when available.
3. **Implemented first confirmation slice:** candidate points can be repeated
   from deterministic start perturbations and are labelled repeated,
   not-reproduced, or inconclusive. This improves numerical evidence without
   relabelling repeated local infeasibility as a global certificate.
4. Compare full corners, random/quasi-random points, sensitivity-screened
   points, and adaptive counterexamples on the same feeders.

Acceptance gate: intentionally unsafe bound-point candidates are reliably
falsified, and stopping records are sufficient to reproduce missed/found
counterexamples.

### Stage 2 — uncertainty and time

1. **Typed provenance slice implemented:** `DOEScenario` and `DOEScenarioSet`
   record scenario IDs, roles, optional relative weights, source, construction
   method, seed, timestamp, and metadata; these fields participate in the study
   identity. Typed uncertainty samples, PSD-aware reproducible Gaussian draws,
   and deep-copy materialization with required callback identity now provide a
   generic path from covariance or model candidates to scenarios. A first
   componentwise physical-support policy now draws from a box-conditioned
   Gaussian by rejection with an explicit draw budget, no clipping, and
   recorded acceptance diagnostics. An empirical residual bootstrap now
   supports row resampling and contiguous moving blocks, with full-library
   identity, source/target time separation, and block provenance. Automatic
   selection of DSSE-derived quantities, non-box transforms/support rules, and
   inverse-Carson/topology insertion remains future work.
2. **Held-out evaluation slice implemented:** role-selected coverage reports
   per-scenario outcomes and empirical, weighted, and conservative candidate-
   violation rates. A one-sided Hoeffding bound is emitted only after an
   explicit i.i.d. assertion. A chronological blocked split helper now exposes
   its exclusion gap and unassessed group/site leakage; finite capacity-
   violation curves retain issued controls; descriptive reference/shifted
   comparisons refuse to claim that distribution shift was detected. Group-
   aware split policies and calibration/evaluation audits now expose ID,
   exact-network, timestamp, group, provenance, and effective-sample-size
   diagnostics under declared requirements. Declared numeric metadata can now
   be compared with multivariate energy distance and per-feature effects;
   permutation inference requires an explicit scenario- or group-level
   exchangeability assertion and remains limited to those covariates.
   Separately issued violation-probability forecasts can now be evaluated with
   reliability tables, proper scores, calibration errors, explicit unresolved
   treatment, and independence-conditional concentration intervals. A first
   temporal shift design preserves contiguous serial order through circular
   window shifts, checks ordering/spacing, and requires an explicit invariance
   assertion stronger than stationarity. An uncertainty-model sensitivity
   layer now holds issued capacities and evaluation settings fixed across
   labelled scenario constructions, preserving paired discordance when stable
   identities exist and labelling independent ensembles unpaired. It makes no
   causal model-effect or distribution-shift claim. Add stationary/seasonal
   resampling, block-length selection diagnostics, and estimator-specific
   scenario materializers next.
3. Couple envelopes to storage SOC, EV energy, ramping, tap operations, and
   first-stage/recourse decisions using the package's multi-period machinery.
4. Evaluate offered capacity, realized utilization, curtailed energy, customer
   value, and price of fairness over time.

Acceptance gate: every probabilistic result reports calibration data,
confidence/coverage semantics, held-out performance, and distribution-shift
limitations.

### Stage 3 — alternative envelope geometries and scalable certificates

1. Add paired lower/upper active-power bounds, then coupled P–Q sets.
2. Implement a linear robust inner-envelope baseline and compare its
   conservatism against nonlinear counterexample search.
3. Investigate SOCP, SDP, QC, interval, or monotonicity certificates only where
   their assumptions match the four-wire model and can be checked.
4. Add decomposition and learned/convex surrogates only with false-safe tests
   against the nonlinear oracle.

Acceptance gate: “certified” appears only when a valid mathematical bound covers
the entire declared utilization and uncertainty set.

## Tutorial program

Develop the tutorials in this order:

1. **Implemented synthetic runnable case — What does an envelope guarantee?** A
   balanced three-phase bound point succeeds, while asymmetric partial
   utilization violates a tight negative-sequence limit; the tutorial follows
   the failure through corners, a non-corner point, confirmation, and adaptive
   margin-directed search. A faithful reproduction of the published 2.91 kW
   Australian-feeder counterexample remains a separate licensed-data benchmark.
2. **Implemented formulation workflow — Control recourse changes the DOE.**
   Compare issue-time settings, scenario recourse, prescribed local laws, and
   perfect context-dependent recourse. The control audit and paired STATCOM
   regression quantify capacity attributable to additional information; a
   dedicated plotted case remains useful future presentation work.
3. **Implemented synthetic first slice — From uncertainty to a tested DOE.** A
   runnable tutorial declares calibration/test roles and provenance, allocates
   without test leakage, and reports held-out empirical and optional i.i.d.
   bounds. It now demonstrates reproducible covariance draws and explicit
   box-conditioned physical support, plus paired and unpaired uncertainty-model
   sensitivity semantics. It also demonstrates ordinary and moving-block
   residual bootstrap provenance. Extend it with empirical DSSE residuals from
   a versioned dataset, non-Gaussian/model-uncertainty construction, and
   sensitivity to block length and support rules using real data.
4. **Fairness–efficiency over time.** Compare equal kW, nameplate-normalized,
   request-normalized, proportional, max–min, and historical-curtailment rules;
   report who benefits, the price of fairness, and sensitivity to the chosen
   entitlement.
5. **Reproducible solver validation.** Demonstrate multistart, residual checks,
   independent power-flow replay, solver failures, and why local feasibility is
   not global certification.
6. **Scaling beyond exact corners.** Benchmark enumeration, sampling,
   sensitivity screening, and counterexample generation as ``N`` grows.
7. **Two-sided and multi-period flexibility.** Couple import/export limits to a
   battery or EV and show why independent single-interval solves can be
   physically inconsistent.
8. **Power quality as a fidelity ladder.** Start with fundamental voltage and
   ampacity, add voltage unbalance and neutral current, and finally compare with
   a harmonic-aware external/reference model.

Each tutorial should state its claim, assumptions, expected runtime, dataset
license, random seed, solver/version information, and a “what this does not
prove” box. Prefer runnable scripts that generate figures and tables over copied
REPL fragments.

## Documentation architecture

The documentation now separates a short learning path from durable reference
material:

- [Dynamic operating envelopes in 15
  minutes](../../docs/src/tutorials/doe_getting_started.md) introduces the central range
  pitfall with one progressive synthetic case.
- The shorter scenario-design, uncertainty-model, and statistical-validation
  tutorials each answer one research question; the [advanced uncertainty
  laboratory](../../docs/src/tutorials/doe_uncertainty_coverage.md) remains the complete
  API-dense workflow.
- [Operating envelopes](../../docs/src/problems/operating_envelope.md) defines stable semantics and
  result fields, the [literature evidence register](../../docs/src/problems/doe_literature_evidence.md)
  records scientific support, and the [research program](../../docs/src/research_program.md)
  places future experiments in the wider scientific direction.

The maintenance rules are:

- Keep the implemented glossary for *allocation*, *operating point*,
  *utilization set*, *scenario*, *recourse*, *verification*, *search-stable*,
  and *certificate* synchronized with result diagnostics.
- Document the mathematical formulation, including quantifier order and which
  network controls are shared or context-dependent.
- Maintain the implemented field-by-field result schema with units,
  missing-value semantics, and permitted claims as diagnostics evolve.
- Maintain the dated literature/evidence register, separating primary results,
  review evidence, reanalysis, hypotheses, preprints, and industry reports.
- Keep benchmark case metadata, generated result tables, source data, and
  licensing together; do not commit publisher PDFs as repository assets.
- Teach one distinction at a time, then make the learner change one modelling
  axis and predict the result before solving.
- Keep overview pages conceptual, tutorials runnable, the problem
  specification normative, and the API page exhaustive; avoid copying the same
  code across all four layers.
- Replace absolute phrases such as “ensures network integrity” with claims tied
  to a declared model, uncertainty set, utilization set, solver status, and
  validation procedure.

## Decision

The next scientific contribution should not be another allocation objective.
The highest-value step is **counterexample-guided range verification with
explicit control recourse and a reproducible evidence record**. It directly
addresses the largest semantic gap in the current implementation, makes the
existing nonlinear four-wire model scientifically useful, and provides the
oracle needed to evaluate later robust, stochastic, linear, or learned methods.
