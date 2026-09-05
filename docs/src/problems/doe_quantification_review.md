# DOE quantification: scientific review and roadmap

> **Review date:** 4 September 2026 · **Scope:** the PowerOptLab implementation,
> scientific quantification claims, experiments, tutorials, and documentation.

## Executive assessment

The current implementation is a useful **nonlinear, unbalanced AC allocation
and verification prototype**. It has unusually careful result semantics for an
early research code: it does not publish failed iterates, it distinguishes one
simultaneous bound from box corners, it can share an allocation across a finite
scenario set, it retains prescribed inverter controls, and it exposes several
fairness choices.

It is not yet a robust DOE quantification method in the strong sense that the
term *operating envelope* often implies. In particular, it does not establish
that every realization inside a continuous utilization set is network-safe; it
does not attach a probability or coverage guarantee to its scenario set; it
allows controllable network-device decisions to change independently between
represented contexts; and it relies on one local nonlinear solve. The most
defensible current description is therefore:

> a locally solved active-power allocation that is feasible at a declared
> finite set of network scenarios and participant-utilization points.

This wording is not merely cautious. Research has shown that feasibility at the
upper endpoint need not imply feasibility under partial utilization in an
unbalanced network, because cross-phase sensitivities and controllable network
settings can reverse the intuitive monotonicity argument [Liu and Braslavsky
(2022)](https://doi.org/10.1109/ACCESS.2022.3203062).

## What quantity should be called a DOE?

Let participant utilization be ``u \in \mathcal U`` and let the advertised
capacity vector be ``\bar p``. For the present one-sided export formulation,
realized active power is

```math
p_i(u) = u_i \bar p_i, \qquad 0 \le u_i \le 1.
```

Let ``\xi`` contain the load, generation, source, topology, and parameter
conditions not controlled by the participant. Let ``z`` contain network state
and control variables. A robust box envelope requires a precisely ordered
statement such as

```math
\forall u \in \mathcal U,\; \forall \xi \in \Xi,\;
\exists z \in \mathcal Z(u,\xi):
g(\bar p,u,\xi,z)=0,\quad h(\bar p,u,\xi,z)\le 0.
```

The quantifier on ``z`` is an operational policy, not a mathematical detail.
Allowing a tap, capacitor, or STATCOM setpoint to be chosen after observing
every utilization and uncertainty realization is a full-recourse assumption.
A fixed issuance setpoint instead belongs outside the universal quantifiers.
A local autonomous controller is neither: it must be represented by its own
closed-loop law. The prototype currently gives each AC context independent
controllable-asset variables, while prescribed IBR Q–V laws remain fixed.

The following claim ladder should be used throughout the package:

| Level | Required evidence | Permitted wording |
|---|---|---|
| Bound point | one declared scenario and ``u=\mathbf 1`` | feasible simultaneous bound point |
| Finite tested set | all declared scenario/utilization pairs solve locally | locally feasible at tested points |
| Falsified | an admissible point with a verified violation is found | candidate envelope is unsafe for the declared set |
| Search-stable | no violation is found under a recorded global-search budget | no counterexample found; not certified |
| Probabilistic | calibrated chance/scenario method plus out-of-sample coverage | estimated or bounded violation probability |
| Certified robust | a valid global/convex inner-set certificate covers the full set | robustly feasible for the declared model and set |

Corner feasibility belongs to the second row unless additional structure proves
that all relevant constraints attain their worst cases at corners.

## Implementation audit

### What is already strong

- The network model is nonlinear and neutral-explicit rather than a balanced
  positive-sequence approximation.
- Export and import are separate studies and retain positive SI units at the
  interface.
- A finite ensemble of forecast or model networks can share one allocation.
- `security=:corners` represents all ``2^N`` zero/full participant combinations
  for small ``N`` and states that this is not a global certificate.
- Connection-bound IBRs retain their apparent-power, current, phase-topology,
  DC-coupling, and prescribed reactive-control equations.
- Allocation choices include equal, total-capacity, proportional, alpha,
  max–min, and equal-curtailment policies with explicit normalization.
- Infeasible or non-optimal nonlinear iterates are not published as envelopes.
- Issuance metadata, fallback provenance, fixed-capacity re-solving, and a
  standalone DSSE-to-DOE replay provide a useful beginning for an audit trail.

### Principal weaknesses

| Priority | Weakness in the current prototype | Scientific consequence | Recommended response |
|---|---|---|---|
| P0 | The advertised range is checked only at one point, all corners, or caller-supplied samples. | Interior violations can be missed, so `:corners` is not a box-containment proof. | Add an adversarial utilization search and distinguish falsified, search-stable, and certified results. |
| P0 | Ipopt supplies one local nonlinear result with no multistart, objective bound, or branch analysis. | Capacity and even feasibility can depend on initialization or the returned AC branch. | Record starts and residuals; add multistart and independent fixed-dispatch power-flow checks; later add valid relaxations. |
| P0 | The historical default makes network controls independently adjustable in every scenario/corner context. | An omitted or poorly chosen policy can assume unavailable perfect recourse from taps, STATCOMs, or other flexible assets. | The first typed control-policy slice now links native taps and IBR P/Q at `:issue` or `:scenario`, retains `:local_law`, and audits `:context`; extend registration and replay coverage before treating the gap as closed. |
| P0 | Verification returns one feasibility flag for a monolithic multi-context NLP. | A failure does not identify the offending scenario, utilization, constraint, or violation magnitude. | Return one structured record per context and an explicit worst counterexample. |
| P0 | Scenario lists are unweighted and caller-constructed. | “Three scenarios” has no confidence level, coverage statement, or reproducible provenance. | Add typed uncertainty/scenario metadata, weights, seeds, generation method, and out-of-sample evaluation. |
| P1 | The envelope is a one-sided scalar active-power box ``[0,\bar p_i]``. | It cannot express a simultaneous import/export band, a nonzero baseline, asymmetric reserve, or coupled P–Q flexibility. | Generalize to lower/upper P bounds and later polyhedral/ellipsoidal P–Q sets with device feasibility. |
| P1 | Intervals are physically independent; rolling fairness only carries an allocation history. | Storage SOC, EV energy, ramping, tap wear, thermal memory, and forecast recourse are omitted. | Reuse the multi-period infrastructure and separate first-stage envelopes from device recourse. |
| P1 | Rolling fairness accumulates offered capacity, not realized use, denied request, customer value, or cost. | A high Jain index can coexist with unequal economic outcomes or repeated involuntary curtailment. | Define the stakeholder, benefit/burden variable, entitlement reference, time horizon, and price of fairness before selecting a metric. |
| P1 | The legacy port constrains aggregate P and Q but can choose phase currents within the OPF. | A multi-phase teaching connection may acquire an unrealizable phase allocation. | Restrict it to an explicit phase-sharing rule or mark it single-phase only; use bound IBRs for research cases. |
| P1 | Diagnostic margins cover a useful subset, not every inherited constraint. | A converter, transformer, neutral, phase-to-phase, control, or other constraint can bind without appearing in the summary. | Build diagnostics from named constraint families and independently recompute all claimed margins. |
| P1 | Exact corners create ``S2^N`` full AC contexts in one monolithic model. | Memory and solution time become prohibitive before realistic fleet sizes. | Add sensitivity screening, adaptive constraint generation, decomposition, and benchmarked stopping rules. |
| P2 | Fundamental-frequency RMS physics omits harmonics, protection, dynamics, communications, and compliance. | The electrical DOE may not be safe or implementable for the deployed system. | Add these only as question-driven fidelity layers with matched measurements and validation data. |

The result diagnostics now explicitly report `global_certificate=false`,
`solver_class=:local_nonlinear`, finite-scenario uncertainty semantics, shared
allocation coupling, and the selected typed control policy. Custom utilization
checks are labelled `:explicit_utilization_points` rather than as a bound-point
test. A per-control audit reports native classification, stage, rule source,
automatic law, contexts, equality groups, and link count. These fields make the
current boundary machine-readable, but do not close the methodological gaps.

## Control-recourse contract: implemented first slice

Security coverage and control recourse should remain independent API axes:

- `security` defines the network scenarios and participant-utilization points
  at which feasibility is required;
- `control_policy` defines the information available when each non-state
  decision is made and which decisions must be equal across those contexts.

For example, `security=:corners` with issue-time controls means that every
corner is represented but all manual setpoints are linked. The same security
mode with context recourse is the current idealized formulation. Neither mode
changes the fact that bus voltages, currents, and other physical state variables
must be context-specific.

### Control stages

The typed policy classifies every discovered free control into one of these stages:

| Stage | Linking rule | Intended meaning |
|---|---|---|
| `:issue` | one value across all scenarios and utilization points in an interval | manual or scheduled before uncertainty and customer use are known |
| `:scenario` | one value per scenario, shared across utilization points | scenario is observed before the control action; valid only when this timing is defensible |
| `:local_law` | no equality link; response is fixed by an enumerated controller equation | automatic Volt–VAr, Volt–Watt, droop, or another implemented causal law |
| `:context` | independent free value for every scenario/utilization pair | anticipative/perfect-recourse formulation used by some published methods; also an upper benchmark when that timing is not implementable |

The framework does not impose one information structure on every study.
It provides named presets `PerfectRecourse()`,
`IssueFixedControls()`, and `IssuePlusLocalLaws()`, plus per-device rules for
mixed formulations. `PerfectRecourse()` must remain a first-class option for
replicating published formulations that optimize controls separately at each
represented point. When its timing is not operationally realizable, reports
should label it as `ideal_recourse`; that label describes the claim rather than
invalidating the experiment.

For operational studies, `IssuePlusLocalLaws()` is the recommended conservative
preset. Unclassified free controls cause an error under this preset, but
`PerfectRecourse()` deliberately classifies them as context-adaptive. A
future API should require an explicit policy whenever a multi-context model has
free controls. During migration, an omitted policy retains the historical
perfect-recourse behavior and records `control_policy_source=:legacy_default`;
a deprecation warning or later explicit-policy requirement remains a future API
decision.

The public shape is:

```julia
policy = DOEControlPolicy(
    default_stage=:issue,
    rules=[
        DOEControlRule(component=:transformer, id="tx1",
                       quantity=:tap, stage=:issue),
        DOEControlRule(component=:ibr, id="pv17",
                       quantity=:reactive_power, stage=:local_law),
        DOEControlRule(component=:ibr, id="statcom17",
                       quantity=:reactive_power, stage=:issue),
    ],
    on_unclassified=:error,
)

r = solve_operating_envelope(scenarios, cps;
    security=:corners, control_policy=policy)

# Reproduce a formulation in which every represented OPF point can redispatch
# all otherwise-free network controls independently.
published = solve_operating_envelope(scenarios, cps;
    security=:corners, control_policy=PerfectRecourse())
```

This is a typed asset/quantity registry, not a list of arbitrary JuMP
variables. Transformer taps and IBR active/reactive powers are linked through
BMOPFTools public handles. General
generator current handles are available, but linking currents would be
physically wrong when voltage varies; generator active/reactive setpoint
expressions or variables need a public upstream seam before they are supported
as linked controls.

### Required diagnostics

Every result now returns a first-slice control audit with:

- component, identifier, controlled quantity, stage, and policy source;
- equality groups used for non-anticipativity;
- the named automatic control law (immutable controller-parameter capture is a
  remaining provenance improvement);
- controls detected as free but not safely linkable;
- rules that do not match a discovered control are rejected, and issue-time
  controls absent from one or more contexts are rejected; and
- summary fields including `all_discovered_free_controls_classified`,
  `nonanticipativity_enforced`, and `ideal_recourse_used`.

The operational claim is determined jointly by coverage and recourse. A useful
summary is therefore a pair such as
`(coverage=:all_box_corners, recourse=:issue_plus_local_laws)`, not one overloaded
word such as `robust`.

### Framework and replication requirements

PowerOptLab should treat each control policy as a reproducible formulation, not
as a hidden implementation switch. To support replication and fair comparison:

- policy presets must expand to the same typed per-control rules returned in
  diagnostics;
- custom policies and presets must use the same model-building path;
- every result must record whether a rule came from a preset, explicit override,
  network control law, or backwards-compatible default;
- publications should be reproducible by a committed study manifest that names
  the policy and software version;
- method comparisons should solve the same network/scenario/utilization set
  under multiple policies and report the capacity attributable to additional
  recourse; and
- the framework should permit new information partitions without adding a new
  DOE solver—for example, controls shared across a scenario-tree history rather
  than only `:issue`, `:scenario`, or `:context`.

The most general internal representation is therefore a partition of contexts
for each control: contexts in the same information set share one decision.
`:issue`, `:scenario`, and `:context` are convenient presets of that partition;
`:local_law` replaces a free decision with a declared causal equation. This
keeps the implementation extensible to multi-stage stochastic, robust,
distributed, and literature-replication studies.

### Implementation sequence

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

## Open scientific questions

### 1. What set is the operator actually promising?

The first decision is whether the object is a point, an independent nodal box,
a feeder-level budget, a coupled polytope, or a P–Q capability set. Boxes are
easy to communicate but discard correlations between participants. A coupled
set can expose more flexibility but shifts optimization and compliance burden
to the aggregator. An experiment should compare volume, communication burden,
and ex-post violations for at least a box, feeder budget, and coupled set.

### 2. Which variables are here-and-now, locally controlled, or recourse?

DOE capacity, day-ahead tap positions, real-time autonomous Volt–VAr response,
and emergency network dispatch live on different information timelines. The
prototype needs an explicit information structure. The scientific comparison
should quantify the apparent capacity gain caused by perfect recourse and then
repeat with fixed and controller-realizable settings.

### 3. How should uncertainty be calibrated?

Finite stress scenarios, chance constraints, norm-bounded robustness, and
distributionally robust methods answer different questions. A good study should
report empirical coverage on held-out time periods, violation frequency and
magnitude, envelope conservatism, and sensitivity to distribution shift. The
state estimator and inverse-Carson candidate sets in this repository offer a
natural route from measurement/model uncertainty to a DOE ensemble, but that
connection is not yet automated.

### 4. Can range safety be falsified efficiently?

For a proposed ``\bar p``, formulate an inner problem that maximizes normalized
constraint violation over ``u`` and ``\xi``. Multi-start local search can find
counterexamples without being called a proof. Each discovered point can be
added to the outer allocation problem. The key research outputs are violation
found, search budget exhausted, and relaxation-certified upper bound—not a
single Boolean called “robust.”

### 5. When is a linear or learned surrogate safe enough?

Linear unbalanced OPF, sensitivity models, polynomial chaos, and learned convex
surrogates can scale much better than repeated nonlinear AC contexts. Their
value depends on error control. Compare them on the same cases using an exact
AC oracle, stress conditions outside the training/calibration set, and explicit
false-safe versus false-conservative rates. A surrogate should never silently
replace the model used to define the security claim.

### 6. What does fairness mean over time?

Equal kW, equal fraction of nameplate, equal request satisfaction, equal
curtailment, welfare, and historical redress are not interchangeable. The
decision variable may be access, actual energy, curtailment, bill impact, or
reliability. Studies should publish the efficiency–fairness frontier and the
price of fairness, stratified by electrical location and relevant customer
classes, rather than declaring one scalar objective fair.

### 7. How should flexible devices and markets consume a DOE?

A technically feasible nodal limit can make an aggregator's energy/reserve
commitment infeasible, or a customer may be physically unable to remain inside
it under PV uncertainty. Device-level energy constraints and market products
therefore need to be tested downstream. Network-secure bidding work treats the
customer controller and market recourse as part of the formulation rather than
assuming compliance [Attarha et al.](https://doi.org/10.1109/TSG.2021.3138099).

### 8. Which omitted physical constraints change the answer?

Voltage magnitude and ampacity are not always the limiting quantities.
Negative-sequence limits are already available in the prototype, but neutral
heating, transformer ageing, protection reach, RMS voltage/current under
harmonics, THD, and control interactions require targeted studies. A 2026
harmonic OPF study found that harmonic-aware modelling can reduce export
envelopes and alter import envelopes, through both distortion limits and the
RMS definitions themselves [Antić et al.
(2026)](https://doi.org/10.1109/TSG.2026.3660479).

## State of the art and implications

The literature does not point to one universally superior DOE algorithm. It
shows a progression from scalable nominal allocations toward stronger but more
expensive uncertainty, range, device, and power-quality claims.

| Evidence | Contribution | Implication for PowerOptLab |
|---|---|---|
| [Liu et al. (2022)](https://doi.org/10.1109/TSG.2022.3188927) | Linearized three-phase OPF for meter-level import/export envelopes, controllable OLTCs, large realistic MV–LV studies, and Monte Carlo utilization assessment. | Add realistic-scale benchmarks and make network-control timing explicit; retain nonlinear AC replay as an oracle. |
| [Liu and Braslavsky (2022)](https://doi.org/10.1109/ACCESS.2022.3203062) | Demonstrates sensitivity and robustness problems when actual customer utilization differs from the allocation point in unbalanced networks. | Treat bound-point results as allocations, not range guarantees; prioritize counterexample search. |
| [Liu and Braslavsky (2023)](https://doi.org/10.1109/TPWRS.2023.3308104) | Constructs robust inner operating envelopes for unbalanced networks and analyzes proportional-fair allocation. | Implement a certified linear baseline against which nonlinear search methods can be compared. |
| [Liu, Braslavsky, and Mahdavi (2023)](https://doi.org/10.35833/MPCE.2023.000653) | Robust linear unbalanced OPF with norm-bounded load and parameter uncertainty. | Add typed uncertainty sets and separate parameter from operating uncertainty. |
| [Liu and Braslavsky (2024)](https://doi.org/10.17775/CSEEJPES.2024.03220) | Sensitivity filtering for robust DOE computation with non-convex unbalanced OPF. | Use constraint/scenario screening before building an exponential monolith. |
| [Koirala, Geth, and Van Acker (2024)](https://doi.org/10.1016/j.segan.2024.101528) | Stochastic unbalanced OPF with polynomial chaos for day-ahead DOEs and alpha-fairness comparisons. | Add probabilistic calibration, held-out coverage tests, and a horizon-level experiment. |
| [Moring, Farrell, and Mathieu (2024)](https://doi.org/10.1109/ALLERTON63246.2024.10735316) | Formalizes fairness over time for DER coordination and demonstrates different dynamic fairness choices. | Replace allocation-only history with declared benefit/burden histories and compare intertemporal policies. |
| [Kumarawadu et al. (2025)](https://doi.org/10.1016/j.apenergy.2025.125469) | Uses smart-meter-derived sensitivities for a three-phase four-wire model with mutual and neutral effects. | Benchmark joint model-estimation/DOE workflows and propagate estimation error into security results. |
| [Antić et al. (2026)](https://doi.org/10.1109/TSG.2026.3660479) | Nonlinear three-phase harmonic OPF for export/import DOEs on real and 906-node networks. | Treat harmonic/RMS constraints as a future fidelity track and publish when they materially change allocation. |
| [de Carvalho et al. (2026 preprint)](https://arxiv.org/abs/2607.08578) | Formulates full-range robust DOE problems with nonlinear AC and approximate linear models. | Reproduce the formulation as an emerging benchmark, but label the evidence as preprint until peer reviewed. |
| [Zabihi et al. (2026 preprint)](https://arxiv.org/abs/2606.29351) | Couples P–Q envelopes, unbalanced AC physics, VUF constraints, and lexicographic/proportional fairness. | Use as a benchmark for the P–Q and lexicographic extensions, again with preprint status explicit. |

Australian implementation evidence reinforces that DOE design is
socio-technical: terminology, calculation, communication, compliance,
fallbacks, and customer outcomes form one system. The [ARENA *Navigating
Dynamic Operating Envelopes* report](https://arena.gov.au/knowledge-bank/navigating-dynamic-operating-envelopes/)
provides a current cross-industry terminology framework, while [Project EDGE's
fairness report](https://www.aemo.com.au/-/media/files/initiatives/der/2023/the-fairness-in-dynamic-operating-envelope-objectives-report.pdf)
shows why the fairness population and normalization must be stated rather than
hidden inside an objective.

## Recommended development program

### Stage 0 — trustworthy experimental record

Deliver this before claiming a new quantification method:

1. **Implemented first slice:** `DOEStudySpec` records interval/scenario network
   hashes, utilization coverage, control information structure, objective,
   normalization, solver/options, random seeds, metadata, and software versions.
   Typed issuance and uncertainty-construction records remain future additions.
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

## Minimum benchmark protocol

Every method comparison should include:

- balanced and unbalanced radial cases, a floating-neutral four-wire case, and
  at least one real or utility-derived feeder with documented licensing;
- voltage-, thermal-, unbalance-, converter-, and deliberately infeasible
  regimes, plus a case where partial utilization is worse than the bound point;
- central, stress, and held-out uncertainty cases with fixed random seeds;
- at least one controllable asset under fixed, ideal-recourse, and local-law
  assumptions;
- capacity/volume, curtailed energy, losses, fairness and price-of-fairness;
- violation frequency, maximum violation, false-safe rate, and empirical
  coverage under independent AC replay;
- solve time, memory, model size, iteration count, failure rate, and scaling in
  participants, scenarios, and intervals;
- sensitivity to solver starts, tolerances, per-unit base, curve smoothing, and
  forecast/model perturbations.

Unit tests should verify contracts and small known cases. Research benchmarks
should produce versioned artifacts and statistical summaries; they should not
be weakened merely to make continuous integration fast.

## Tutorial program

Develop the tutorials in this order:

1. **Implemented first runnable case — What does an envelope guarantee?** A
   balanced three-phase bound point succeeds, while asymmetric partial
   utilization violates a tight negative-sequence limit; the tutorial follows
   the failure through corners and adaptive margin-directed search. Extend it
   with a non-corner interior failure as the oracle matures.
2. **Control recourse changes the DOE.** Compare a fixed STATCOM/tap setting,
   perfect context-dependent recourse, and a local controller. Plot both
   capacity and ex-post violations.
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

## Documentation changes

The documentation should be organized around scientific claims rather than
only API entry points:

- Keep [Operating envelopes](operating_envelope.md) as the concise API and
  semantics page; link to this review for limitations and research priorities.
- Split the existing broad tutorial as the tutorial program matures. Its
  current material is a good modelling primer, but it does not yet demonstrate
  a failure mode, probabilistic calibration, scaling, or a reproducible
  fairness frontier.
- Add a glossary for *allocation*, *operating point*, *utilization set*,
  *scenario*, *recourse*, *verification*, *search-stable*, and *certificate*.
- Document the mathematical formulation, including quantifier order and which
  network controls are shared or context-dependent.
- Provide a field-by-field result schema with units, missing-value semantics,
  and permitted claims for every diagnostic status.
- Maintain a dated literature/evidence table, separating peer-reviewed work,
  preprints, and industry reports.
- Add benchmark manifests and generated result tables to the documentation,
  with source data and licensing recorded.
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
