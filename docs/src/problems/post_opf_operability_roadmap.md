# Post-OPF voltage operability roadmap

> **Status:** implementation underway · **Implementation:** first static
> checker slice and opt-in HELM cross-check landed · **Scope:** static,
> fundamental-frequency, unbalanced four-wire networks
>
> **Owner:** PowerOptLab maintainers · **Last reviewed:** 2026-08-21

This page preserves the scientific and implementation plan for a processor that
checks the voltage operability of an OPF-determined equilibrium. It is not an API
promise and it does not claim that a static post-processor proves dynamic voltage
stability.

## Implementation checkpoint

The current branch contains the first bounded implementation slice in
`src/problems/operability.jl`:

- independent static-load endpoint residuals with explicit scope exclusions;
- audited-policy-scaled Jacobian evidence, singular values, and critical-mode
  participation;
- connection-level and sequence-level observables plus uniform load-scale and
  named P/Q sensitivities, with requested-versus-realized power records and
  path-qualified `dP/dV`; opt-in independent finite-difference re-solves now
  validate the implicit uniform-load and named P/Q state derivatives, and the
  report now classifies the path slope as a high-side/near-nose/low-side
  indicator under its declared sign convention, with complementary positive-/
  negative-sequence magnitude and VUF derivatives;
- an opt-in HELM endpoint cross-check that records the energized no-load
  homotopy, supported-physics preflight, HELM status, and endpoint mismatch; and
- a natural-parameter load-scale continuation trace with damped Newton
  correctors, residual/singular-value history, and explicit near-singular or
  corrector-failure events; and
- a first static pseudo-arclength predictor/corrector slice that records
  target-parameter crossings, fixed-λ target refinements, tangent-based fold
  candidates with critical left/right mode participation, voltage-normalized
  arclength provenance, and the same endpoint evidence. An analytic two-bus
  high-/low-voltage branch regression now checks
  that a locally feasible low branch is distinguished from the no-load-connected
  high branch; and
- an explicit `stop_at_target=false` stress mode now continues past the audited
  λ=1 crossing and exercises bordered fold localization on the analytic nose;
- an opt-in bordered-equation fold localizer solving ``F=0``, ``Jv=0``, and
  ``‖v‖₂=1`` from a declared approximate fold state; pseudo-arclength fold
  candidates now retain its localized status and residual evidence when invoked.

These additions are evidence-producing diagnostics, not a claim of globally
complete branch discovery or a production-grade fold certificate. The current
pseudo-arclength slice is limited to the native static ybus load scope; target
refinement improves endpoint matching but does not replace fold refinement or
global branch discovery. The bordered localizer is a local candidate solver and
does not certify uniqueness or global reachability.

## Research question and intended claims

Given a network model and a candidate OPF solution, determine separately:

1. whether the published point satisfies the declared network and device
   equations when evaluated independently;
2. whether its device-terminal, conductor, neutral, and sequence voltages meet
   declared operating criteria;
3. whether the equilibrium is locally regular and how it responds to declared
   active- and reactive-power perturbations;
4. whether it lies on the real solution sheet connected to a specified energized
   base state without crossing a fold or an incompatible control event; and
5. how far that sheet can be followed along declared stress directions before a
   saddle-node or an operational limit is encountered.

The result must retain a claim hierarchy. Each check reports `pass`, `fail`,
`inconclusive`, or `not_applicable`, with its model, path, scaling, tolerance,
and numerical evidence. There will be no aggregate `operable = true` that hides
which of these different claims was established.

The four per-check outcomes are initially **operability-local vocabulary**, not
a repository-wide replacement for result-specific solver statuses or
[`SolveStatus`](@ref). Milestone 0 will define one local type and aggregation
rules before result structs ship. Promotion to a shared repository vocabulary
requires a separate review against the existing result contracts.

## Terminology

Use **energized-base-connected branch** or **no-load-connected branch**, rather
than assuming that the desired solution is componentwise the "highest-voltage"
solution. Under unbalance, reverse power flow, distributed controls, and
voltage-dependent devices, componentwise voltage ordering need not define a
unique branch.

In the initial implementation, *reachable* will mean reachable through a
specific, documented homotopy. It will not mean dynamically reachable under
arbitrary disturbances or controls. Reachability is path-dependent.

The processor studies a static equilibrium. Induction-motor dynamics, load
restoration, protection, tap timing, controller dynamics, and short- or
long-term dynamic voltage stability are outside its initial claims.

## Equilibrium model

Write the post-OPF equilibrium in real coordinates as

```math
F(z, \lambda; m,c)=0,
```

where ``z`` contains every free rectangular voltage component and any algebraic
controller states retained by the chosen closure, ``\lambda`` is a declared
loading/stress parameter, ``m`` identifies load laws, and ``c`` identifies
connection and control modes.

Two closures answer different questions and must not be mixed:

- **frozen dispatch:** optimized setpoints, taps, and other independent controls
  are fixed at the OPF values; and
- **operational equilibrium:** physical algebraic controls such as droop or PV
  regulation remain in ``F``, with their active limit mode recorded.

The OPF KKT matrix is not the power-flow Jacobian and does not answer which
power-flow branch contains the primal voltage point.

### Scaling provenance is part of the audited point

The checker must inherit the effective scaling policy used to build the audited
OPF, rather than choose a new per-unit convention. BMOPFTools exposes
`AbstractOpfScalingPolicy`, `SIUnitsScaling`, `ClassicPerUnitScaling`,
`ConsistentPerUnitScaling`, `ZonePerUnitScaling`, and
`opf_scaling_policy(ctx)`; `opf_research_provenance(ctx)` retains a serialisable
description of that effective policy.

The preferred input is therefore an OPF context or explicit research
provenance alongside the SI-valued solution. If only a detached solution
dictionary is available and its scaling provenance was not retained, the
caller must supply the policy explicitly. The checker must not infer one from
nominal voltage or silently fall back to `ClassicPerUnitScaling`.

Endpoint mismatch is evaluated in physical SI coordinates and also normalized
with the audited policy's coordinate bases for tolerances and attribution. The
Jacobian factorization, singular values, and rank thresholds use that same
policy. A mismatched or incomplete policy is an input/provenance failure, not a
new operability conclusion.

At a regular endpoint,

```math
J = F_z,
\qquad
\frac{dz}{d\lambda}=-J^{-1}F_\lambda.
```

For a derived observable ``g(z,\lambda)``, report the directional sensitivity

```math
\frac{dg}{d\lambda}=g_z\frac{dz}{d\lambda}+g_\lambda.
```

Near a generic fold, ``J`` loses rank and natural-parameter sensitivities
become unbounded. Pseudo-arclength continuation augments the equations so that
the solution curve can still be followed and the fold can be bracketed.

## Connection-aware device derivatives

All voltage checks and derivatives start in the actual device-terminal basis.
For a real terminal-incidence matrix ``C``, let

```math
u=Cv,
\qquad
i_{\mathrm{node}}=C^{\mathsf T}i_{\mathrm{branch}}(u).
```

Rows of ``C`` form phase-to-neutral voltages for wye devices and line-to-line
voltages for delta devices. With

```math
C_R=\begin{bmatrix}C&0\\0&C\end{bmatrix},
```

the rectangular contribution is

```math
J_{\mathrm{device}}=C_R^{\mathsf T}J_{\mathrm{branch}}C_R.
```

This representation is required: it captures cross-phase delta coupling,
floating-neutral displacement, and the distinction among phase-to-earth,
phase-to-neutral, and phase-to-phase voltages without special sensitivity
formulas for each topology.

Load laws must contribute their exact incremental characteristic. For example,
a ZIP active-power branch may use

```math
P(|u|)=P_0\left(a_P+a_I\frac{|u|}{u_0}
                    +a_Z\frac{|u|^2}{u_0^2}\right),
```

with a corresponding reactive-power law. Constant-power, constant-current,
constant-impedance, ZIP, and supported exponential laws therefore produce
different equilibrium Jacobians even if they consume equal power at the OPF
point.

The implementation should obtain this physics through a public BMOPFTools
equilibrium/residual seam if one exists or is added. The current upstream audit
is already recorded in `src/upstream.jl`: the pinned BMOPFTools release lacks a
public read-only seam for the load decomposition, so HELM isolates
`_load_subloads`, `_subload_S`, `_subload_yz`, `_stamp_pair!`,
`_neutral_terminal`, `_neutral_labels`, and `_DEFAULT_CONFIG` in that single
compatibility adapter. Milestone 0 must begin from that finding, extend the same
adapter only if unavoidable, and propose the smallest upstream public seam. It
must not create a second private-import boundary.

The processor should not create a second, silently divergent copy of the OPF
device equations. An independent evaluation path is still required for
validation; independence should come from the residual assembly and numerical
method, not from changing the physics.

## Interpreting ``dP/dV``

``dP/dV`` is retained as a path-specific reported diagnostic, not the primary
classifier.

For the scalar constant-power two-bus relation

```math
P(V)=\frac{V(E-V)}{R},
\qquad
\frac{dP}{dV}=\frac{E-2V}{R},
```

the consumption-sign convention gives a negative slope on the high-voltage
side, zero at the nose, and a positive slope on the low-voltage side. That sign
does not generalize without declaring:

- which requested power or stress direction is being varied;
- which device-terminal or derived voltage is observed;
- whether ``P`` means requested power, realized voltage-dependent load power,
  or network injection; and
- which controls are fixed or responding.

For voltage-dependent loads, the derivative of realized power contains the
load law's own voltage dependence. The tip of a plotted realized-P versus V
curve need not identify the equilibrium fold. The implementation will therefore
store requested/stress power and realized device power as different quantities.

For a one-dimensional continued curve with arclength ``s``, a displayed
``dP/d|u|`` may be formed as

```math
\frac{dP/ds}{d|u|/ds}
```

where the denominator is nonzero. The continuation tangent, fold history, and
full Jacobian evidence remain the branch evidence.

## Phase, terminal, and sequence observables

The reporting hierarchy is:

1. device-terminal quantities;
2. conductor and per-phase quantities; and
3. sequence aggregates at complete three-phase buses.

The endpoint report should include, where defined:

- phase-to-neutral magnitude and angle for wye devices;
- line-to-line magnitude and angle for delta devices;
- phase-to-earth voltage when grounding or insulation makes it relevant;
- neutral displacement and neutral/reference integrity;
- phase-angle separation;
- zero-, positive-, and negative-sequence phasors;
- voltage-unbalance factor ``|V_2|/|V_1|``; and
- all corresponding directional sensitivities requested by the study.

For a sequence transform ``T``, both state and tangent are transformed:

```math
V_{012}=T V_{abc},
\qquad
\frac{dV_{012}}{d\lambda}=T\frac{dV_{abc}}{d\lambda}.
```

For example,

```math
\frac{d|V_1|}{d\lambda}
=\Re\left\{\frac{V_1^*}{|V_1|}\frac{dV_1}{d\lambda}\right\}.
```

Positive-sequence quantities are useful system-level summaries, but they can
hide a weak phase, negative/zero sequence, and neutral displacement. They are
also not intrinsically defined at one- or two-phase buses. They will never
replace terminal-level checks.

## Local regularity and sensitivity evidence

Use the full rectangular current-balance Jacobian, including whichever
controller equations belong to the selected closure. Transform state and
residual coordinates using the audited OPF's effective scaling policy before
factorization; do not introduce a processor-local per-unit system.

Report at least:

- smallest scaled singular value ``\sigma_{\min}(J)``;
- a reciprocal condition estimate;
- left and right critical singular vectors;
- bus, terminal, device, phase, and sequence participation in the critical
  mode; and
- declared directional ``d|V|/d\lambda``, ``d\theta/d\lambda``,
  ``d|V_1|/d\lambda``, and selected ``d|V|/dP_k`` / ``d|V|/dQ_k`` values.

The determinant is not a headline metric because it is badly scaled and its
sign depends on coordinate conventions. A reduced Q-V Jacobian is not the
default because high-R/X and unbalanced distribution systems do not justify the
usual decoupling assumptions.

Local nonsingularity is not a branch classifier: high- and low-voltage sheets
are both locally regular away from a fold.

PowerOptLab already exports `finite_difference_jacobian`,
`fixed_point_oracle`, and `screen_fixed_point_gain`, together with network
sensitivity result types, from `src/problems/closed_loop_evidence.jl`. Milestone
2 should reuse their input validation, scale-aware step conventions, and result
semantics where applicable. The current `finite_difference_jacobian` requires a
square real map (`length(y) == length(x)`), so it directly checks ``F_z`` for a
square equilibrium residual but does not cover ``F_\lambda`` or rectangular
directional maps without a deliberate generalization. That limitation must not
be hidden by adding a near-duplicate helper.

## Branch reachability and margin

### HELM cross-check

PowerOptLab's existing [four-wire HELM implementation](../algorithms/helm.md)
starts from the energized no-load germ and analytically continues its selected
branch. When the HELM physics covers the case and its endpoint converges to the
OPF voltage vector, that agreement is strong independent evidence for
no-load-connected reachability. The [HELM versus nonlinear power-flow
tutorial](../tutorials/helm_vs_nonlinear_power_flow.md) states the existing
cross-method evidence contract.

The current HELM scope is constant-power and constant-impedance load parts; it
does not yet reproduce arbitrary current/exponential loads, generators, IBRs,
or OPF controls. Because `helm_series` is the proposed programmatic cross-check,
its relevant status symbols are `:series_diverged` and `:max_order_reached`;
they remain numerical outcomes, not non-existence certificates.

Before calling `helm_series`, the operability processor must run an explicit
supported-physics preflight. Unsupported physics returns `not_applicable` with
the named devices and reasons. The expected `ArgumentError` currently raised by
HELM for unsupported loads must not escape as an uncaught processor failure.
Unexpected exceptions remain errors rather than being relabelled as scientific
inconclusiveness.

### General continuation

Harden and extend the first predictor-corrector pseudo-arclength slice on the
same connection-aware equilibrium model:

1. construct and verify the energized base equilibrium;
2. orient a homotopy from that base toward the OPF endpoint;
3. adapt the arclength step using corrector convergence and curvature;
4. monitor the scaled Jacobian, tangent, and device/control events;
5. detect and refine folds;
6. detect every crossing of the target parameter ``\lambda=1``; and
7. compare the first pre-fold target crossing with the OPF voltage and retained
   controller states.

The current implementation covers the static native-ybus slice and records
items 1, 2, 4, and 6 provisionally, with fixed-λ refinement when a target
crossing is bracketed. Fold refinement, deflation, curvature-based step control,
and controller/generator closures remain Milestone 3 work.

The initial default homotopy proposed for review is:

- keep passive network elements, physical shunts, grounding, and the source
  boundary energized;
- scale all demand nameplate terms, including P/I/Z fractions, from zero to the
  endpoint values while retaining their fractions and connection;
- scale fixed non-slack generator/IBR injections from zero to their OPF values;
  and
- keep optimized independent settings fixed, while offering a separate
  operational-control continuation when their physical laws are available.

A successful endpoint match establishes reachability only along that recorded
homotopy. A first fold below ``\lambda=1`` rejects reachability on the pre-fold
sheet along that path. A corrector failure without a refined singular point is
inconclusive.

Continue beyond the endpoint, when requested, to find the first of:

- a saddle-node/fold;
- a voltage or current limit;
- reactive-power or inverter-capability saturation;
- a tap, deadband, or other discrete-control event; or
- a model-domain failure such as a device terminal approaching zero voltage.

Margins are path-specific. Initial stress directions should include uniform
fixed-power-factor demand, feeder/area demand, weakest-phase demand, DER trip or
export increase, and reactive-power loss. A later closest-boundary calculation
in a multidimensional injection space is a different and harder research task.

## Relationship to operating envelopes

[`solve_operating_envelope`](@ref) allocates feasible import/export capacities
over declared forecast/model scenarios, while [`verify_operating_envelope`](@ref)
checks fixed allocations at selected utilization points. The operability
processor answers an adjacent question about one returned equilibrium: endpoint
validity, local regularity, branch reachability, and path-specific margin.

The initial APIs remain deliberately separate. An operating-envelope solve does
not automatically acquire a branch or voltage-stability certificate, and the
operability processor does not allocate capacity. Later,
`verify_operating_envelope` or a counterexample-guided envelope study may invoke
the processor for selected snapshots and retain its evidence, but it must report
the finite scenario/utilization/path scope rather than upgrade it to a global
robust-feasibility claim.

## Sufficient certificates

Where its assumptions match the compiled model, add the multiphase fixed-point
conditions of Bernstein et al. as a sufficient existence/uniqueness and
Jacobian-nonsingularity certificate. A passed condition is valuable evidence;
a failed condition is `inconclusive`, not evidence of multiplicity or collapse.

HELM convergence, fixed-point certificates, local Jacobian regularity, and CPF
each establish different claims. The report will not silently substitute one
for another.

## Proposed result structure

The eventual structured result should separate:

| Section | Evidence |
|---|---|
| provenance | network/solution identity, units, grounding, load laws, connections, controls, tolerances |
| endpoint | independent current/power mismatch and model-domain checks |
| voltage quality | terminal, phase, neutral, angle, and sequence observables versus caller-supplied limits |
| local regularity | scaled singular values, conditioning, rank tolerance, critical modes |
| sensitivities | named perturbation directions and phase/terminal/sequence responses |
| reachability | base state, homotopy, continuation history, folds/events, endpoint distance |
| margin | first limiting mechanism and path-specific parameter/arclength margin |
| certificates | sufficient conditions, assumptions, pass or inconclusive outcome |
| cross-checks | HELM/nonlinear/OpenDSS or other independent-oracle comparisons |

Every numerical array must carry an ordering map back to buses, terminals,
devices, and controller variables. Raw singular vectors without that provenance
are not publishable research output.

## Implementation roadmap

### Milestone 0 — contracts and analytic fixtures

- Finalize terminology, sign conventions, tolerances, and operability-local
  result-status semantics.
- Require the effective upstream `AbstractOpfScalingPolicy` through a live
  context, retained research provenance, or an explicit caller argument; define
  how missing or mismatched provenance fails before analysis.
- Specify the frozen-dispatch and operational-equilibrium closures.
- Start from the existing `src/upstream.jl` compatibility audit; determine which
  additional connection-aware residual/device derivative data are missing,
  extend only that adapter if unavoidable, and propose the smallest upstream
  public seam.
- Add closed-form two-node fixtures for constant P, I, Z, ZIP/exponential cases,
  including both sides of a known saddle-node where applicable.
- Add balanced wye/delta equivalence, unbalanced delta, floating-neutral, and
  deliberately perturbed endpoint fixtures.

**Exit criterion:** the scientific contract and expected fixture results are
reviewed before a public API is fixed.

### Milestone 1 — endpoint and voltage audit

- Compile the free-voltage ordering and device-terminal incidence maps.
- Reconstruct the candidate OPF phasor state without using solver start values.
- Evaluate independent KCL/device mismatch and model-domain checks.
- Use the audited OPF scaling policy for normalized residual tolerances and
  record its serialised provenance in the report.
- Publish terminal/phase/neutral/sequence observables with configurable limits.
- Record unsupported devices explicitly; never omit them silently.

**Exit criterion:** a valid endpoint passes, a deliberately perturbed endpoint
fails, and all reported quantities trace to their terminal definitions.

### Milestone 2 — local Jacobian and directional sensitivities

- Assemble the full rectangular equilibrium Jacobian.
- Reuse or deliberately generalize the existing `finite_difference_jacobian`
  for square ``F_z`` checks; validate ``F_\lambda`` and rectangular directional
  maps with separately named tests rather than a duplicate square-map helper.
- Validate analytic/automatic derivatives against scale-aware central
  differences and directional residual tests.
- Add audited-policy-scaled SVD/conditioning and critical-mode attribution.
- Solve named P/Q stress directions and transform their responses into device,
  phase, and sequence coordinates.
- Report path-qualified ``dP/dV`` only after requested and realized powers are
  separated.

**Exit criterion:** derivative checks pass across P/I/Z, wye/delta, unbalanced,
and floating-neutral fixtures, and conclusions are invariant to equivalent unit
scaling.

### Milestone 3 — continuation and branch matching

- Implement adaptive pseudo-arclength predictor/corrector continuation.
- Refine folds and record left/right critical modes.
- Match the first pre-fold ``\lambda=1`` crossing to the OPF endpoint.
- Preflight HELM support, map unsupported physics to a reasoned
  `not_applicable`, and cross-check compatible fixtures against `helm_series`.
- Expose failed correctors as inconclusive unless a limiting mechanism is
  independently established.

**Exit criterion:** analytic saddle-node locations and high/low endpoint labels
are recovered, and the checker distinguishes endpoint feasibility from
base-connected reachability.

### Milestone 4 — controls and limit-induced events

- Retain smooth droop/PV/controller equations in the operational closure.
- Add PV-PQ, inverter capability, and other continuous limit transitions.
- Add explicit piecewise-smooth event handling for taps, deadbands, and switched
  devices where the network model supports it.
- Distinguish a saddle-node from a limit-induced loss of the current operating
  mode.

**Exit criterion:** control events are reproducible and carry pre/post model
provenance; no derivative is reported through an unacknowledged discontinuity.

### Milestone 5 — certificates, robustness, and studies

- Add the sufficient multiphase fixed-point certificate on its supported scope.
- Run multiple named stress directions and model ensembles for ZIP/exponential
  uncertainty.
- Add selected contingency workflows without treating a finite list as a global
  guarantee.
- Provide table-ready rows and a documented research-study tutorial.

**Exit criterion:** every positive claim is tied to a theorem, continuation
trace, or declared finite study; conservative certificate failures remain
inconclusive.

## Proposed repository layout

The processor and the continuation algorithm have different roles. A proposed
layout, subject to review, is:

```text
src/problems/operability.jl          # public types and orchestration
src/problems/operability/equilibrium.jl # residual, ordering, incidence
src/problems/operability/metrics.jl  # terminal/phase/sequence observables
src/problems/operability/sensitivity.jl # Jacobian and directional solves
src/algorithms/continuation.jl       # reusable pseudo-arclength engine
test/operability_tests.jl
test/continuation_tests.jl
docs/src/problems/operability.md
docs/src/tutorials/post_opf_operability.md
```

The public entry point might ultimately resemble
`check_opf_operability(network, solution; specification=...)`, but naming and
types should be chosen only after Milestone 0 establishes the upstream seams and
result contract.

## Verification matrix

The maintained test matrix should include:

| Feature | Required oracle |
|---|---|
| two-node constant-P fold | closed-form voltage, slope, and collapse loading |
| constant-I/Z and ZIP/exponential | closed-form or independently evaluated device law and finite differences |
| balanced wye/delta | symmetry and analytically equivalent terminal loading |
| unbalanced delta | cross-phase residual/Jacobian blocks and OpenDSS comparison where available |
| four-wire wye | neutral rise and phase-to-neutral versus phase-to-earth distinction |
| high and low solutions | continuation history from the energized base, not solver initialization |
| scaling-policy covariance | physically identical SI/classic/consistent/zone-policy solves, each audited under its own retained upstream policy, give invariant physical and dimensionless conclusions; missing/mismatched provenance is rejected |
| control saturation | known first event and active-mode transition |
| HELM-compatible cases | HELM/CPF/nonlinear endpoint agreement on common physics |
| unsupported physics | explicit `not_applicable` or `inconclusive`, never silent omission |

ForwardDiff, analytic derivatives, and central finite differences should be
compared on small fixtures, reusing `finite_difference_jacobian` for its square
map scope. Sparse factorizations may replace dense reference calculations only
after numerical parity is established.

## Review decisions before implementation

1. Confirm the proposed default homotopy, especially whether fixed non-slack
   injections scale from zero and whether every load ZIP fraction scales by the
   same nameplate multiplier.
2. Confirm whether Milestone 1 initially covers only native BMOPFTools static
   devices or must also compile PowerOptLab advanced IBR equations.
3. Decide which optimized controls are frozen by default and which physical
   controls remain endogenous.
4. Confirm that voltage/unbalance thresholds are caller-supplied study data;
   the library should not embed one jurisdiction's limits as universal facts.
5. Decide whether the public result should retain complete continuation states
   or a compact trace plus an opt-in full trajectory.
6. Agree on the smallest public BMOPFTools seam needed to prevent device-physics
   duplication, beginning from and preserving the single adapter in
   `src/upstream.jl`.
7. Confirm that the audited OPF's effective `AbstractOpfScalingPolicy` is a
   required provenance input and that detached solution dictionaries without
   that provenance must supply it explicitly.

## Scientific references

- V. Ajjarapu and C. Christy, [“The continuation power flow: a tool for steady
  state voltage stability analysis”](https://doi.org/10.1109/59.141737), IEEE
  Transactions on Power Systems, 1992.
- I. Dobson and H.-D. Chiang, [“Towards a theory of voltage collapse in electric
  power systems”](https://doi.org/10.1016/0167-6911(89)90072-8), Systems &
  Control Letters, 1989.
- T. J. Overbye, [“Effects of load modelling on analysis of power system voltage
  stability”](https://doi.org/10.1016/0142-0615(94)90037-X), International
  Journal of Electrical Power & Energy Systems, 1994.
- V. Ajjarapu and S. Battula, [“Effect of load modeling on steady state voltage
  stability”](https://doi.org/10.1080/07313569508955639), Electric Machines &
  Power Systems, 1995.
- A. Bernstein et al., [“Load-Flow in Multiphase Distribution Networks:
  Existence, Uniqueness, Non-Singularity and Linear
  Models”](https://doi.org/10.1109/TPWRS.2018.2823277), IEEE Transactions on
  Power Systems, 2018; [open manuscript](https://arxiv.org/abs/1702.03310).
- A. Trias, [“Fundamentals of the Holomorphic Embedding Load-Flow
  Method”](https://arxiv.org/abs/1509.02421), 2015.
- P. S. Nirbhavane et al., [“TPCPF: Three-Phase Continuation Power Flow Tool for
  Voltage Stability Assessment of Distribution Networks With Distributed
  Energy Resources”](https://doi.org/10.1109/TIA.2021.3088384), IEEE
  Transactions on Industry Applications, 2021.
- P. Juanuwattanakul and M. A. S. Masoum, [“Analysis and Comparison of Bus
  Ranking Indices for Balanced and Unbalanced Three-Phase Distribution
  Networks”](https://www.researchgate.net/publication/254012429_Analysis_and_comparison_of_bus_ranking_indices_for_balanced_and_unbalanced_three-phase_distribution_networks),
  AUPEC, 2011.
- M. Jereminov et al., [“Impact of Load Models on Power Flow
  Optimization”](https://arxiv.org/abs/1902.04154), IEEE PES General Meeting,
  2019.
