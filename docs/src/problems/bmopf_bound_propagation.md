# BMOPF bound propagation and tightening: design proposal

Status: research and implementation design, 22 September 2026. The first
solver-free slice is now implemented; see the [usage and current coverage](../tutorials/bmopf_bounds.md).
The broader APIs and stages below remain proposals where the usage page does
not document an implementation. The source audit uses the
BMOPFTools revision pinned in `Project.toml`,
`72e6cec22a66cf376c37ec4d64aef350b9f1100d`, rather than the newer sibling checkout.

## Recommendation

Build a **solver-independent, connection-aware bound service**, followed by an
optional optimization-based bound tightening (OBBT) service. Both should update
the same domain store and produce the same provenance and coverage report.
Prototype in PowerOptLab; promote the stable parser adapter, propagation rules,
and serialization contract to BMOPFTools once tested across its device families.

The contract should be: **account for every object and enclose every supported
physical quantity over the declared feasible set**. It cannot promise finite
bounds on every quantity for every valid JSON. Syntactically valid networks can
be infeasible, underconstrained, unsupported by a particular formulation, or
genuinely unbounded. Report these distinctions instead of manufacturing limits.

The initial product is useful without a solver: fill missing implied limits,
expose missing assumptions, detect some contradictions, and supply finite domains
for function envelopes and solver variables. The second stage spends optimization
effort on the remaining influential domains.

```text
BMOPF JSON + explicit model/time/control scope
                     |
       validation and semantic normalization
                     |
      quantities + relations + declared domains
                     |
          propagation to a stopping criterion
                     |
              domain/proof store
               /           \
     report + JSON patch    selected OBBT targets
                                  |
                        valid outer relaxation
                                  |
                       accepted objective bounds
                                  |
                         propagate again
```

## What the state of the art supports

This is a targeted review of methods relevant to this design, not an exhaustive
survey or a claim that one formulation is universally strongest.

| Evidence | Established result or approach | Consequence for this project |
|---|---|---|
| [Pyomo FBBT](https://pyomo.readthedocs.io/en/stable/api/pyomo.contrib.fbbt.fbbt.html) | Forward expression evaluation and backward constraint propagation tighten domains without auxiliary optimization. | Start with specialized electrical contractors and a small expression kernel; use a dependency queue. |
| [Coffrin, Hijazi and Van Hentenryck, 2015/2016](https://arxiv.org/abs/1512.04644) | Bound tightening, convex envelopes and valid cuts work together to strengthen AC power-flow relaxations. | Alternate propagation and relaxation strengthening; neither is a one-time preprocessing step. |
| [Sundar et al., 2018/2019](https://arxiv.org/abs/1809.04565) | Strengthened QC relaxations support OBBT for voltage magnitudes and angle differences; the paper distinguishes feasibility-set tightening from tightening with an objective cutoff. | Make feasible-set scope explicit. Bounds conditional on cost must never become unconditional network limits. |
| [Narimani et al., 2023/2024](https://arxiv.org/abs/2310.14463) | Joint envelopes and coordinate rotations strengthen QC relaxations on transmission benchmarks. | Preserve useful relations, not just scalar boxes. Consider stronger envelopes after the basic service works. |
| [Shao et al., online 2024, volume 2025](https://link.springer.com/article/10.1007/s11081-024-09891-7) | An SOCP-based OBBT method targets a subset of variables to reduce the expense of exhaustive tightening. | Evaluate bound improvement per unit time and use a target scheduler. |
| [Cengil et al., 2025](https://link.springer.com/article/10.1007/s10589-025-00715-7) | Dynamic, learned variable ranking accelerates OBBT on large AC-OPF benchmarks. | Start with deterministic impact scores; later learn which subproblems to solve. Predictions select work, never certify bounds. |
| [Solving Three-phase AC Infeasibility Analysis, 2025 preprint](https://arxiv.org/abs/2508.15937) | Three-phase infeasibility optimization combines McCormick relaxations, spatial branching, filtering and decomposition. | Relevant evidence for rectangular three-phase machinery, but not validation of BMOPF's entire explicit-neutral device catalogue. |
| [Awadalla and Bouffard, 2023](https://arxiv.org/abs/2305.10385) | Objective information affects constraint screening in AC OPF. | Keep cost-aware screening and general-purpose bound completion as separate modes. |

The first architectural inference is to combine cheap feasibility-based bound
tightening (FBBT) with selective OBBT. The second is that the formulation must
match BMOPF: transmission QC results do not establish soundness or performance
for four-wire distribution networks with floating neutrals and detailed controls.

[PowerModels' versioned OBBT documentation](https://lanl-ansi.github.io/PowerModels.jl/v0.15/utilities/)
is a useful software precedent, not an API dependency to assume in current
PowerModels versions. [PowerModelsDistribution](https://arxiv.org/abs/2004.10081)
is a useful comparison framework for multiconductor formulations. Conversion
must demonstrate preservation of the original feasible set before its bounds
can be imported; a numerically successful converted solve is insufficient.

## Existing seams and gaps

The following observations are from the pinned BMOPFTools source, available in
the [revision tree](https://github.com/frederikgeth/BMOPFTools.jl/tree/72e6cec22a66cf376c37ec4d64aef350b9f1100d).

| Existing facility | Reuse and boundary |
|---|---|
| `parse_bmopf`, `schema_check`, `spec_conformance_check` | Parsing performs migration and terminal normalization; parsing alone does not validate. Retain the raw-to-normalized mapping for export. Add integrity and supported-model checks. |
| `src/augmentation/voltage_bounds.jl` | Inserts policy/recipe-based voltage limits into missing fields. This changes the feasible set; it is not implication-based propagation. |
| `ext/BMOPFOpfExt/data_utils.jl`, `_limit_current_box!` | Already derives rectangular current bounds from a magnitude limit in selected builders. Generalize the mathematical rule, including its warning about series versus terminal currents. |
| `OpfModelKey`, public key constructors, `opf_object_keys`, `opf_object` | Stable attachment points for applying domains and selecting OBBT targets. Avoid variable-name parsing and private `ctx` dictionaries. |
| `OpfSemanticBlock`, `opf_diagnostic_schema`, `opf_coordinate_bases` | Existing quantity, component and unit metadata. Useful mapping substrate, but not a complete solver-independent relation graph or interval proof system. |
| `build_opf_model`, `enforce_kcl!`, staged builders | Useful exact-model oracle for differential checks. The reference optimizer is a nonconvex IV model, not an existing convex OBBT backend. |
| PowerOptLab's [bounded function formulations](../formulations/bounds_and_relations.md) | Can consume certified sensed-input domains. Currently does not derive them from arbitrary network constraints. |

The principal missing seam is a public, versioned description of **resolved
device connections and equations without constructing a JuMP model**. Add a
small adapter first, with parity tests against upstream; move shared connection
resolution upstream instead of maintaining a second independent interpretation
of loads, transformers and IBRs. Keep solver-specific lowering separate.

## Validity and scope

Let `F(net, scope)` be the states satisfying the declared physics and active hard
constraints. For each scalar quantity `q`, the output must satisfy

```math
\{q(x):x\in F\}\subseteq [\ell_q,u_q].
```

These are outer bounds, not an operating envelope whose every point is feasible.
The Cartesian product of individually valid domains usually contains many
infeasible states. Propagation reaching a fixed point does not prove feasibility
or tightness of the result.

`scope` must specify schema/semantics version, snapshot or time range, fixed
versus adjustable taps and switches, active control laws, hard-limit policy,
parameter uncertainty, and optional extension devices. Nominal values and starts
are not limits unless the device model explicitly makes them fixed quantities.
For example, `p_nom` of a voltage-dependent load is not its actual constant power.

Default to all feasible states under the supplied hard constraints, with no cost
cutoff and no implicit operational voltage band. An optional policy augmentation
runs before propagation, producing a new input and a visible assumption record.
Do not reinterpret an absent limit, a zero rating, and an inactive device as the
same thing: resolve their meanings for the specific field and schema version.

Cost-limited tightening uses `F ∩ {f(x) ≤ C}` and gets a different scope identity.
It preserves globally optimal solutions when `C` is a valid incumbent upper
bound for the same minimization problem, but can remove otherwise feasible
states. A locally found point can supply that cutoff only after its feasibility
and objective have been accepted under the declared numerical policy.

For alternative snapshots, scenarios, or discrete modes, a bound valid for all
alternatives is the **hull of their domains**, not their intersection. When a
shared decision must satisfy all scenario constraints jointly, compile that
joint problem explicitly. Per-snapshot boxes can omit temporal couplings safely
as an outer relaxation, but cannot claim to exploit storage dynamics.

### Finite bounds are conditional on the available information

Examples the service must explain:

- Bounds on phase-to-neutral voltages alone do not bound phase-to-earth or
  neutral-to-earth voltages in a floating system; common-mode motion may remain.
- Bounding the total injection at a bus does not bound two individual ideal
  sources that can inject equal and opposite unlimited currents.
- Fixed endpoint voltages do not bound circulation in a loop of ideal conductors.
- A constant-power relation cannot yield `I_max = S_max / V_min` if the relevant
  connection voltage can approach zero.
- A loose relaxation with an unbounded subproblem proves no finite bound by that
  relaxation. It does **not** prove the original AC problem unbounded.

Every registered quantity has a domain, including infinite endpoints where
necessary. Separately report whether an endpoint is declared, derived, unresolved,
or supported by a proof of unboundedness. Distinguish unsupported semantics from
an understood quantity for which current rules derive no finite bound.

## Internal representation and public interface

Use a typed internal representation and a JSON-compatible sidecar. Do not put
solver variables, inferred material properties, and operational limits into one
dictionary of guessed `*_min` / `*_max` fields.

| Record | Required content |
|---|---|
| `BoundScope` | Model/input fingerprint, schema and adapter versions, snapshots, discrete modes, controls, assumptions and optional objective cutoff. |
| `QuantityKey` | Object family/subtype/id, port or winding, terminal/pair, quantity, component, snapshot; explicit voltage reference and current orientation. |
| `Domain` | Real interval, complex rectangular enclosure plus magnitude interval, or finite discrete set. Unit and endpoint validity are explicit. |
| `Relation` | Affine row, norm, product, connection map, device law or control graph; ordered input/output quantities, applicability conditions, source pointers. |
| `BoundEvidence` | Separate lower/upper derivations: rule version, parent evidence IDs, input JSON pointers, arithmetic policy, assumptions, optional solver certificate. |
| `BoundReport` | Domains, object coverage, conflicts, stop reason, unsupported features, work counts and exportable patch. |

Store nonconvex information such as an annulus, admissible tap values or an angle
arc as structured domains/relations. A rectangular box alone loses useful
information. Initially support real intervals, complex boxes with magnitude
constraints, and finite tap sets; defer general unions except where needed to
avoid deleting admissible modes. Use a conservative hull when a union is not
represented, and report the loss of precision.

Proposed API (illustrative Julia, deliberately not a runnable example):

```julia
scope = BoundScope(problem=:feasibility, snapshot=1, hard_limits=:declared)
report = propagate_bounds(net; scope, arithmetic=:outward)
bound(report, key)
explain_bound(report, key; endpoint=:upper)
coverage(report)

# Same report contract; an explicitly supported relaxation backend is required.
tighter = tighten_bounds(net, report;
    backend=IVConicRelaxation(), targets=:impact, time_limit=60.0)

# Non-mutating export, only where BMOPF has an equivalent limit field.
filled, patch = materialize_bounds(net, tighter; mode=:fill_missing)
apply_bounds!(ctx, tighter)  # checks scope and coordinate mapping first
```

Default calls do not mutate the input. `:fill_missing` adds equivalent missing
fields; a separate `:intersect` mode may strengthen existing fields. The domain
store always intersects with applicable original limits, regardless of export
mode. Reject stale reports after a change to topology, coefficients, profiles,
controls or scope; semantic equality/reuse optimizations can come later.

Each object gets a coverage entry, even if it contributes only constants,
metadata, or an unsupported extension. Auxiliary linecodes, geometry, wire data
and profiles are data dependencies rather than independent operating states.
An unknown device must retain unknown injections; omitting its contribution
from KCL would impose a false zero injection and invalidate other bounds.

### JSON export

Deliver three artifacts: a valid BMOPF network, a versioned bounds sidecar, and a
patch with original values and reasons. Most rectangular coordinates, internal
winding quantities and conditional bounds belong in the sidecar. Filling a
network limit means adding a mathematically redundant constraint for that scope;
it does not infer an equipment nameplate or allowable rating for other studies.

JSON has no infinity. In the sidecar encode an infinite endpoint as `null` plus
an explicit endpoint status; distinguish it from missing/unsupported information.
Omit such endpoints from BMOPF fields. Do not insert arbitrary new properties into
schema-constrained device records. Conditional scenario/cost bounds are not
eligible for unconditional export.

Round-trip validation must check the exported JSON against its schema and check
the field's meaning, terminal order, units and active-constraint interpretation.
Use normalized JSON as the initial export target, with documented migrations;
source-preserving edits require a verified reverse mapping to original fields.

## Stage 1: propagation

Compile declared constraints into a dependency graph. Seed exact constants,
hard bounds, discrete domains and unconditional facts such as nonnegative
magnitudes. Repeatedly run affected contractors and intersect their candidate
domains with the current ones. A contractor must never remove a feasible value.
Keep KCL, KVL and connection equations explicit so improvements travel through
the network rather than stopping at the device that supplied a limit.

Use deterministic queue order, reverse incidence, per-quantity physical scales,
relative/absolute improvement thresholds and a work/time cap. Record
`fixed_point`, `tolerance_stopped`, `budget_exhausted`, or `conflict`; do not call a
tolerance stop an exact fixed point. Skipping a tiny contraction is safe; rounding
an endpoint inward to force convergence is not. Maintain an immutable derivation
node for each accepted update so cyclic network propagation still produces an
acyclic history of evidence.

### Initial rule catalogue

| Quantity / relation | Sound rule | Essential qualification |
|---|---|---|
| Complex magnitude | `|z| ≤ M ⇒ Re(z), Im(z) ∈ [-M,M]` | `|z| ≥ m` does not imply either component has magnitude at least `m`. |
| Rectangle to magnitude | Upper radius from farthest rectangle corner; lower radius from distance of rectangle to origin | Retain both box and magnitude information. A lower norm bound is nonconvex. |
| Connection voltage | `w = v_a - v_b`; `|w| ≤ U_a + U_b`; `|w| ≥ max(0,L_a-U_b,L_b-U_a)` | Applies to the actual pair, including phase-to-neutral and delta branches. |
| Neutral reference | `|v_phase| ≤ |v_phase-v_n| + |v_n|` | Phase-to-neutral upper bound needs a neutral bound to give a phase-to-earth upper bound. |
| Linear rows | Isolate each variable in `Σ a_j x_j = b` and intersect with interval evaluation | Separate real/imaginary components; do not divide by a coefficient interval containing zero. |
| Power | `S = w conj(i)`; `|S| ≤ U_w U_i`; finite P/Q boxes imply `|S| ≤ hypot(maxabs(P),maxabs(Q))` | The voltage and current must be a conjugate pair for that port. An aggregate rating does not bound each cancelling phase contribution. |
| Current from power | `|i| ≤ U_S/L_w`, if `L_w > 0`; `|i| ≥ L_S/U_w`, if `U_w > 0` | Use connection voltage, not nominal voltage; intersect with the direct current limit. |
| ZIP load | Propagate the actual constant-power, constant-current and impedance submodels, then sum their connection currents | Do not treat all nominal P/Q as constant power or silently regularize zero voltage. |
| Line π model | `Δv = Z i_series`; `i_from = i_series + Y_from v_from`; `i_to = -i_series + Y_to v_to` | Retain full mutual coupling, terminal maps and both shunts. Thermal limits are on total terminal currents. |
| Neutral / ground current | Per-terminal KCL and explicit earth paths | Neutral current is not generally the sum of phase magnitudes; the sum is only a conservative upper bound under the relevant connection equation. |
| Shunt / capacitor | Linear current map `i=Yv`; propagate power using the same port relation | Respect delta/wye incidence and fixed versus switchable state. |
| Ideal switch | Open fixes current to zero; closed equates mapped terminal voltages | Closed does not imply zero current; ideal cycles can leave circulation unresolved. |
| Transformer | Propagate each winding's connection, ratio and leakage equations in both directions | Coil and terminal currents differ; include tap range, phase shift, core shunt and subtype-specific rating semantics. |
| Sequence voltage | Apply the complex linear Fortescue map, with the model's actual phase/reference convention | Phase bounds do not transfer unchanged to positive sequence; VUF division needs a positive lower bound on `|V1|`. |
| Control graph | Bound sensed input, evaluate every reachable PWL segment and breakpoint, propagate output/capability intersections | Match equality versus cap, extrapolation, saturation, deadbands and active-power priority. Smoothing is a different law unless enclosed by an error bound. |
| DC network / IBR DC port | Signed voltage/current rows, connection differences and `P=VI` | Pole-to-earth and pole-to-pole differ; no division through zero or assumption of positive voltage without evidence. |

For lines, affine KVL row propagation works even if `Z` is singular. A verified
linear solve can additionally bound `i` from `Δv` when appropriate. Do not invert
a singular matrix or assume a floating-point inverse is exact. Fixed transformer
ratios give affine maps; free ratios introduce interval products and, later,
McCormick relaxations. Discrete tap propagation retains every admissible value.

The coverage roadmap includes all transformer subtypes (including multi-winding,
autotransformers and regulators), voltage sources, load configurations,
generators, IBR AC/DC ports and controls, grounding, and time-series references.
Until a subtype has been implemented, return explicit partial coverage. A report
must not say “complete” because it visited every top-level dictionary.

### Small example and a useful counterexample

Suppose a single load connection has `210 ≤ |v_p-v_n| ≤ 250 V`,
`|v_n| ≤ 10 V`, and fixed consumption `S = 3000 + j4000 VA`.
Then `|v_p| ≤ 260 V`, `|v_p| ≥ 200 V`, and
`20 ≤ |i_load| ≤ 5000/210 ≈ 23.80952381 A`. Its real and imaginary current
components are each enclosed by `[-5000/210,5000/210] A`. The lower current
magnitude does not exclude either component being zero. Production endpoints
must be rounded outwards; the displayed decimal is illustrative.

If a line terminal is limited to 100 A and its shunt current is bounded by 5 A,
then its series current is bounded by 105 A. A bound of 100 A on the series
variable is invalid: a 105 A series current and opposing 5 A shunt current meet
the terminal limit. This is already recognized in BMOPFTools' current-box helper.

### Numerical evidence

Use outward interval arithmetic for certified propagation, including input
conversion, matrix coefficients, square roots, trigonometric constants and
serialization. [IntervalArithmetic.jl](https://juliaintervals.github.io/IntervalArithmetic.jl/dev/manual/api/)
is a candidate implementation; pin and test the operations and configuration
used. Decide explicitly whether the mathematical input is parsed binary numbers
or exact decimal JSON values; the latter requires enclosing decimal conversion.

Padding a final answer with one `nextfloat` cannot certify an arbitrary sequence
of ordinary floating-point operations. Numerical solver results likewise need a
separate evidence class. Distinguish `validated_interval`,
`solver_bound_with_tolerances`, and `heuristic_candidate`; the last must never
enter the accepted domain store. Evidence inherits the weakest premise, so an
outward contractor applied to a solver-tolerance bound is not newly certified.

An empty interval from validated rules gives an infeasibility proof for the
declared scope. Near-conflicts involving solver-tolerance evidence require
revalidation and cannot silently become certified infeasibility. Never repair
an empty intersection by widening an original constraint. Return the conflicting
derivation chain; a minimal conflict set can be a later feature.

## Stage 2: optimization-based tightening

For the current domain box `B`, construct a relaxation `R(B)` with a documented
embedding of every original feasible state. For each selected scalar quantity:

```math
L_q = \inf_{x\in R(B)} q(x),\qquad
U_q = \sup_{x\in R(B)} q(x).
```

In practice solve **two minimizations**, `min q` and `min -q`. If accepted global
lower objective bounds are `d_plus` and `d_minus`, update

```math
\ell_q\leftarrow\max(\ell_q,d_{plus}),\qquad
u_q\leftarrow\min(u_q,-d_{minus}).
```

Apply any justified outward error correction before intersection. Use objective
**bounds**, not the returned feasible objectives: a primal point of `min q`
usually supplies an upper bound on the minimum, which is the wrong direction.
Ipopt extrema, even from many starts, do not establish global ranges. They are
useful witnesses, warm starts and candidate targets.

The adapter should check termination, dual/global-bound availability and solver
semantics. A time-limited run can still contribute an accepted global bound;
absence of one means no update. `LOCALLY_SOLVED`, numerical failure and
`INFEASIBLE_OR_UNBOUNDED` cannot certify a bound or original-model infeasibility.
An infeasible outer relaxation supports original infeasibility only with the
appropriate certificate/numerical evidence classification.

For strict numerical certification, verify dual feasibility or a conservative
Lagrangian bound, accounting for residuals over known variable domains. Some
residual corrections need finite domains; if those are absent, refuse the
certificate rather than assuming a tolerance is harmless. A generic safety
margin on `objective_bound` is an engineering policy, not a mathematical proof.

### First relaxation backend

I recommend a native rectangular IV LP/SOCP backend for the first supported
device slice. This is an engineering recommendation based on BMOPF semantics,
not a claim that it will dominate transmission QC or SDP relaxations.

1. Preserve affine KCL, full matrix KVL, shunts, fixed-ratio winding connections
   and explicit neutral/ground equations.
2. Preserve convex upper magnitude balls as cones, including terminal-current
   limits and compatible P/Q capability disks.
3. Lift products in `P=w_r i_r+w_i i_i` and
   `Q=w_i i_r-w_r i_i`, using all four McCormick inequalities over finite input
   domains. Share repeated lifts.
4. For `y=x²`, use the convex lower epigraph and the upper secant on finite
   `[l,u]`: `x² ≤ y ≤ (l+u)x-lu`. Squared voltage lower limits can then constrain
   sums of lifted squares without incorrectly replacing an annulus by a ball.
5. Relax free tap products and supported PWL control graphs with valid outer
   envelopes; keep admissible discrete alternatives or a justified convex hull.
6. Omit unsupported nonlinear relations only through an explicit safe-relaxation
   rule retaining their free quantities/injections. Never replace them with an
   approximate equality. Report omissions and their effect on coverage.

Finite product domains are a prerequisite for finite McCormick envelopes. If
stage 1 cannot establish them, use a weaker relaxation without those envelopes,
solve a linear/conic preliminary range problem if useful, or return unresolved.
Artificial large boxes are not a valid escape from this dependency.

Do not directly maximize a norm in a convex solver and call it convex OBBT.
Use a valid lifted squared-magnitude representation with a linear target, or
bound coordinates and combine them conservatively. Nonlinear query expressions
need a lifting/enclosure contract with the same original-state embedding.

Clarabel is already a project dependency and is a candidate LP/SOCP solver;
backend acceptance depends on the required bound evidence. Stronger lifted
multiconductor SOCP/SDP constraints, QC variants and verified global solvers can
be added behind the same interface. Benchmark them by cost and domain improvement
on matching physical models. A DC or LinDistFlow approximation is not generally
an outer relaxation of AC physics and cannot supply unconditional valid bounds.

### Scheduling and reuse

Start with deterministic scores based on finite-domain gaps, product-envelope
widths, the number of dependent constraints, controller breakpoints and recent
tightening benefit. For `z=xy`, `(u_x-l_x)(u_y-l_y)` is a useful envelope-width
proxy. Prioritize unbounded prerequisites, neutral and connection voltages,
current coordinates, free taps, and controller sensing domains. Derive cheap
reporting quantities after these improve instead of optimizing every alias.

Reuse models and warm starts, but refresh bound-dependent envelopes after each
round. Propagate accepted updates before selecting more targets. Worker-local
solver models can process a common domain snapshot in parallel; their accepted
bounds can be intersected safely if all came from the same scope. Avoid shared
mutable solver state. Stop by time, solve count, or marginal improvement.

Cache by normalized model, scope, domain version, relaxation, arithmetic and
solver settings. Bounds for a changed load/control/topology are not reusable
without a proof of feasible-set inclusion. Local subnetworks are usable only
with outward boundary-injection domains; fixing boundary values to a previous
solution can cut off valid network states.

## Delivery sequence and acceptance criteria

| Increment | Deliverable | Gate before claiming support |
|---|---|---|
| 1. Contracts and inventory | Quantity catalogue, scope, sidecar schema, declared-bound ingestion, normalization trace, complete object accounting | Every input object is accounted for; no nominal/start value becomes a constraint accidentally; export round-trips. |
| 2. Useful solver-free slice | Bus/connection domains, line π models, KCL, fixed sources, basic loads/generators, shunts and fixed switches | Analytic and adversarial cases; independent evaluation of each rule; input immutability and deterministic stopping. |
| 3. Full device expansion | Transformer subtype adapters, taps, sequence quantities, IBR controls, DC and profile resolution | Per-subtype parity with pinned BMOPFTools equations and explicit status for all remaining gaps. |
| 4. OBBT vertical slice | LP/SOCP IV relaxation for the supported fixed-topology subset, objective-bound adapter, iterative propagation | Prove lifting/containment for each relation; validate bound direction and solver-status handling; no local optimum accepted as a global bound. |
| 5. Performance and promotion | Scheduling, batching, stronger relaxations, stable BMOPFTools integration | Reproducible end-to-end benefit on representative cases; stable semantic/versioning contract. |

Suggested implementation layout: `src/bounds/{types,scope,normalize,registry,
arithmetic,propagate,export}.jl`, `src/bounds/rules/` by relation/device family,
and `src/bounds/obbt/` for relaxations, certificate adapters and scheduling.
Keep the stage-1 kernel independent of JuMP even while prototyped in this package.
Do not start by building a general-purpose symbolic optimizer or a monolithic
AC-OPF global solver.

Tests should establish soundness, not merely shrinking intervals:

- Analytic examples for power/current conversion, affine rows, magnitude/box
  conversion, grounding, taps, PWL breakpoints and both power-flow branches.
- Adversarial floating neutral, zero-voltage constant-power load, phase-to-phase
  single-phase load, delta circulation, parallel ideal sources, cancelling phase
  powers, line-shunt cancellation, singular impedance and contradictory limits.
- Nonstandard terminal labels, terminal permutations, phase permutations with
  consistent sequence convention, SI/per-unit equivalence, and profile/time scope.
- Zero current limits, absent limits, unchanged known limits, unsupported objects,
  arbitrary JSON object order, no input mutation and reproducible provenance.
- For OBBT, deliberate timeouts, unavailable bounds, local-only solves,
  infeasible-or-unbounded statuses, incorrect bound direction and stale reports.
- Cost-cutoff bounds cannot be exported as general bounds; alternative-scenario
  aggregation uses hulls; smoothing error and unsupported relation omissions are
  visible in evidence.

Use exact/rational or high-precision independent references on tiny cases,
interval proofs for contractors, and global solutions where available. Feasible
sample containment is a valuable regression check, but never a proof that all
feasible states are retained. Differential checks against the exact BMOPFTools
model should cover active physical constraints and extension scope, not only
matching objective values.

Benchmark small analytic networks, the checked-in LV3/MV21 snapshots, and ENWL
cases when available, stratified by neutral treatment, transformer/control mix,
and missing-bound pattern. Do not silently add ratings to make a test pass.
Measure finite-endpoint coverage by quantity family, normalized width reduction,
unresolved reasons, propagation work, auxiliary solves, build time, memory,
downstream solve behavior and relaxation gap where a valid objective comparison
exists. Report ratios for initially finite widths separately from transitions
from infinite to finite domains.

The first implementation milestone should be increments 1–2: a sound,
explainable fill-in pass with truthful partial coverage. That establishes the
semantic and numerical foundation on which every later tightening method depends.
