# Advanced IBR review: resolution plan

## Green integration and coordinated Australian capability reference

At `2ead3b4`, the complete offline suite passes 20,431/20,431 assertions
(7m42.5s). The LCL current-scale regression is fixed and the independent-write
serialization fixture has deterministic provenance. All 77 extended Australian
cases also pass again at the matching physical-model source hashes; their
reported powers, currents and voltages are unchanged from the preceding source.

The coordinated Australian driver is validated separately with 163 assertions
covering balanced/free sequence allocation, optional customer voltage bounds,
three numerical bases, input preservation and rejection of a raw-success point
with failed native SI equations. Its 24 predeclared candidate cases pass physical,
policy, native SI, voltage and fixed-PQ replay gates. Saved-point reconstruction
passes 1,148 checks and recomputes 5,472 physical inequalities across 288 devices.
Eight Julia record tests and four Python CSV tests reject modified currents,
missing records, stale source hashes, changed metrics and false success flags.

The benchmark releases autonomous VV+VW Q laws, and distinguishes balanced
converter currents from additional free I₂ actuation on unchanged three-leg
hardware. A 253 V upper phase-neutral bound is an illustrative extra planning
criterion, not regional conformity. Source/table/figure publication follows
the committed-driver rerun; local feasible solutions are not global upper bounds.
The full goal remains active; four-leg control and representative EMT follow
completion and commit of this benchmark.

The full suite at `2b6d4de` clears the LCL regression (20,422 passing assertions)
but exposed an unrelated clock-sensitive Kron serialization fixture. Two writes
crossed a UTC second and differed only in automatically generated `meta.created`.
The fixture now supplies a fixed timestamp and explicitly checks its preservation;
the complete electrical/provenance dictionary comparison remains unchanged.

## LCL numerical-scale regression

The first full run after `0874954` passed 20,418 assertions with one failure:
the watt-priority, converter-target LCL test returned `OTHER_ERROR`. The new
converter-current characteristic had used `i_grid_max` when `i_max` was omitted.
Those arm currents differ across the LCL shunt. The characteristic now uses a
declared converter rating or physical voltage/power data, independently of the
grid-arm rating. The grid current constraint itself remains unchanged.

All six priority/target reproductions now publish and pass physical checks,
along with the engineering and complete PCC power-boundary regressions. A
metamorphic regression checks that changing only a grid-arm rating cannot
change the converter-current characteristic. The final integration and refreshed
study provenance are recorded above. Sizing and
physical acceptance tolerances are unchanged.

## Optional-current scaling, loss error budget and operating-policy acceptance

The engineering audit reproduced raw `LOCALLY_SOLVED` DC demand understated by
about 0.30 W for an unrated three-leg plant with `a_loss=100 W/A`. The prior
physical gate correctly rejected it. The shifted-norm error is now capped at
0.0001 W in total across semiconductor legs, while retaining any smaller
rating-relative width. Current characteristics derive from declared currents or
physical voltage starts/ratings, rather than treating a unity SI coordinate base
as a physical voltage. LCL/three-leg/balanced-arm current rows retain engineering
scaling without introducing an omitted current bound. Australian sizing is fixed.

Constrained-device extraction now independently reconstructs operating P/Q from
terminal phasors and checks current balance, smooth VV/VW policy, fixed Q/PF and
power boxes. Standalone fixed-dispatch requests are audited too. A deliberately
coarse raw-success solve at 60 W / 0 var for a 2,500 W / −600 var request passes
all hardware checks but fails the new Q/dispatch checks; its raw flags remain
visible and publication is rejected. The unchanged physical floors apply.

`dc_source_line_frequency_diagnostics` reports ω/2ω source currents with explicit
RMS/peak conventions, resistance dissipation and a normalized cancellation
screen. Separate capacitor-node KCL and period integration validate the source
response. The mean availability/bounds apply at the regulated DC-link input,
including converter/capacitor ESR loss; upstream source/DC–DC mean-path and
harmonic dissipation are outside this boundary. The diagnostic adds no upstream
availability bound, thermal rating, regulator or dynamic-stability model.

All 304 new engineering assertions pass, including SI/10 kVA/1 MVA/100 MVA
paired constrained solves, independent policy counterexamples and actual raw
solver-success regressions for violated Q/dispatch, PF and volt-watt relations.
The full frozen-source suite and refreshed Australian campaigns remain validation
steps at this source checkpoint. Coordinated OPF,
four-leg control and representative EMT remain part of the active goal.

## Validated coupling and expanded Australian campaigns

At core checkpoint `11a9723`, the complete offline suite passes
20,119/20,119 assertions (7m46.0s), including all 102 new fleet-coupling assertions.
The study setup at `e1f7327` passes 186 separate source-sequence/profile/analytic
energy-accounting assertions. All 41 predeclared snapshots and all 36 synthetic
chronology intervals are accepted without per-case retries. Each environment
records the same 55 executed source/setup hashes. The saved-record audit checks
22,356 physical inequalities, terminal phasors, case completeness and independently
recomputes duration-weighted per-device energy. Four record-regression cases
cover incomplete plans, nonfinite failed candidates and false physical success.

The new `voltage_aware_snapshots` and `voltage_aware_chronology` campaigns are
distinct from the earlier fixed-structure scaling controller. They preserve the
same hardware and terminal reactive-priority boundary settings. Source U₂
magnitude/angle and load phase allocation vary independently. The largest
fixed-PQ cross-formulation phasor difference is 65.1 nV; this replay still uses
the same network engine and is not an independent external circuit oracle.
Standalone PNG/SVG figures and generated tables explain numerical iteration,
unbalance-angle dependence and operating-mode energy costs.

The synthetic day uses essentially all 405 kWh of DC availability within the
physical tolerance. Its POC energy shortfall is modeled converter/filter loss,
not material DC curtailment. Energized night operation imports 6.18–6.93 kWh
across the fixed 12-inverter fleet while retaining reactive absorption. Opening
the grid-side contactors removes AC services/losses; disconnected DC standby
remains undefined and excluded. These are hypothetical quasi-static profiles,
not measured yield or dynamic trajectories. Mean and sequence sensing nearly
coincide here; sequence sensing does not itself request I₂ compensation.

The goal remains active. Outstanding work includes the optional-arm SI and
upstream/bridge DC-source boundary audit, coordinated Australian OPF, deliberate
four-leg control and representative EMT validation. The retained earlier
1 MVA/300-iteration fixed-structure mean-controller failure remains diagnostic;
the separate common 1,000-iteration comparison is documented below.

## Voltage-aware feeder coupling and contactor modes

`solve_voltage_aware_inverter_fleet` now couples the local Q→P→I₂ allocator to
native feeder power flows. Grid-arm currents and the POC shunt enter KCL exactly
once. Numerical feeder probes omit inverter hardware constraints so a saturated
rail does not prevent evaluating a perturbed voltage map. Returned physical
states are reconstructed and independently audited at the actual feeder voltage;
acceptance additionally requires local allocation success, current/power/policy
closure, voltage-map convergence and native SI equality residual checks. The
native residual audit evaluates registered equations, not an independent circuit
implementation or all native operational inequalities.

Explicit energized and grid-side-disconnected modes support complete per-device
mode maps. Energized zero-PV operation retains Q support and modeled AC loss
imports. Disconnection removes the complete converter/filter/POC-shunt AC model;
DC standby consumption is outside that model. Unselected native devices remain
native. Failures retain explicit reasons and cannot publish; a failed local NLP
does not establish physical or policy infeasibility. Sizing is unchanged.

The 102 new assertions pass: three numerical bases, independent diagonal-feeder
Ohm/KCL and loss identities, native-device retention, warm starts/damping, mixed
operating modes, LCL passive I₀/neutral sensing, controller-domain failures,
allocation failures and explicit iteration exhaustion. Full-suite validation and
expanded campaigns are recorded at the validated checkpoints above. This is
a quasi-static local fixed point, not a dynamic stability or uniqueness proof.

The broader goal remains active: controlled-unbalance/chronological studies,
coordinated Australian OPF, four-leg control, representative EMT and the remaining
optional-arm SI/DC-source boundary audit are not yet complete.

## Validated source checkpoint e19cbae

The frozen source passes the full offline package suite: 20,017/20,017 assertions
(7m53.2s), including the 56 replay and 117 service-allocation regressions. All 152
canonical Australian main/sensitivity/terminal cases pass again, as does the
380-check independent replay audit. The 14 campaign environments match the same
52 source-file fingerprints at source commit `e19cbae`.

The retained 1 MVA diagnostic still reaches its 300-iteration limit for mean
VV+VW. An additional predeclared common 1,000-iteration profile accepts all three
laws at that base. Against the 100 kVA reference, maximum per-device differences
are 4.29 µW P, 23.31 µvar Q and 6.95 nA maximum phase-current magnitude; maximum
phase-neutral voltage difference is 0.275 µV (neutral-to-ground differs by up to
5.34 µV). Both the failed 300-iteration and accepted 1,000-iteration profiles are
preserved. This demonstrates a stopping-budget failure for that original case,
not a physical infeasibility or hardware-sizing requirement. It does not claim
robustness at every base/operating point or change the default 100 kVA campaign.

The goal remains active: local allocation is implemented, but its network
coupling and the unbalance/disconnection/chronological/OPF/four-leg/EMT program
are not complete.

## Local voltage-aware Q→P→I₂ allocation

`allocate_inverter_services` now supplies staged physical allocation at a fixed
measured voltage. Q tracking runs first, then achieved Q is fixed while P is
maximized, then both powers are fixed while requested I₂ ray support is maximized.
Every stage independently audits fixed voltage, currents/sequence requests,
preserved powers, source and volt-watt conditions, thermal limits and continuous
rails. No finite priority weight or per-case retry is used. Ipopt uses one common
adaptive/no-scaling profile with explicit caller overrides preserved.

Within-request mode maximizes a bounded Q-service coordinate and fixes that
coordinate for P allocation. The initial equality-plus-Q-interval representation
produced a false P-stage infeasibility at the Australian worst-phase request;
removing that redundancy recovers it without changing hardware or the service law.
Protective Q override is an explicit different mode; the strict failed request
never silently acquires it. A failed stage issues no command. An independent
unrestricted plant projection can provide a feasible witness, but a failed local
policy NLP still establishes no global policy infeasibility.

`evaluate_voltage_aware` supplies exact mean/worst/positive-sequence VV+VW requests
from separately supplied sensed and physical POC phasors. It supports reactive
priority and voltage-admittance negative-sequence requests, and explicitly rejects
current-dependent ripple blending and balance-first ordering. The existing
fixed-structure stamp remains unchanged. This is a solver-based local reference;
feeder feedback, response times and dynamic stability remain unverified.

Validation: 117 new assertions cover analytic Q/P capability at three bases,
saturated Q, the Australian active-Q boundary regression, rail protection,
explicit override versus unresolved strict request, PV ceilings/night loss import,
preserved P/Q under I₂ service, separate sensing boundaries, LCL passive I₀ and
admittance support. The seven-case reference records six accepted commands and
one explicit unresolved strict policy with an accepted physical diagnostic.
Australian hardware is imported unchanged from the study builder.

Remaining goal work includes network coupling of this allocator, the optional-arm
SI/DC-source audit, full controlled-unbalance/disconnection/chronological studies,
coordinated Australian OPF, four-leg control and representative EMT. The active
goal is not complete; inverter sizing is unchanged.

## Fixed-current replay acceptance

The replay helpers underpinning the fixed-point and network-sensitivity oracles
still used raw solver acceptance. A direct reproduction returned publishable
`LOCALLY_SOLVED` at zero current with a full-cycle rail margin near -35.38 V.
Both single and fleet paths now refine rails, normalize command equalities,
use a base-invariant watt objective/default zero Ipopt bound relaxation, and gate
publication on each independent physical/current report. Raw MOI fields and solver
messages remain visible, and unused selected native IBR outputs are removed.

Fixed-point iterations use virtual imposed-current capability probes. Applying
PV availability to a deliberately perturbed starting voltage prematurely prevents
evaluation of otherwise convergent maps, so the returned point independently
audits request availability instead. Its Pdc upper enclosure uses the reported
value plus the independent auditor's absolute reported-Pdc residual. An overdrawn
returned point cannot publish or report convergence. Individual internal replays
also accept explicit signed DC source bounds. Sensitivity probes remain capability
experiments and do not imply available PV generation. These are quasi-static
calculations, not transient response simulations.

Validation: all 115 existing closed-loop assertions and 56 new acceptance/source/
cross-base assertions pass. The new cases cover single/fleet false rail successes,
one rejected device in a mixed fleet, source-limited and non-charging replays,
rejected sensitivity bases, and overdrawn returned points after a virtual iterate.
The full goal remains active; the priority-aware allocation layer and the
Australian scenario/OPF/four-leg/EMT work remain pending.

## Australian source freshness check, 31 October 2026

The [Energex installer standards index](https://www.energex.com.au/our-services/connections/residential-and-commercial-connections/solar-connections-and-other-technologies/connection-information-for-installers)
still lists STNW1170 dated 21 February 2025. The
[release-10 PDF](https://www.energex.com.au/__data/assets/pdf_file/0010/1072549/STNW1170-Standard-for-Small-IES-Connections.pdf)
retains the VV/VW settings in Tables 8–9. Section 4.3.3 requires balanced AC
output for multiphase IES, and the three-phase-inverter definition addresses
simultaneous connection/disconnection. Preserve the distinction between sourced
curve parameters and research sequence/neutral-control capability. Future
disconnection and four-leg cases must state their connection assumptions;
the present physical-model checks do not establish connection conformity.

## Offline voltage-aware allocation reference

`project_inverter_current` now provides the reusable fixed-POC-voltage reference.
It minimizes squared converter-arm or grid-arm current error over the declared
plant with free sequence allocation, non-charging PV bounds by default and
adaptive full-cycle rail refinement. A separate SI reconstruction checks the
requested point through the full primitive/LCL filter. Split-link balancing
charge is explicitly fixed in both calculations. That conditional scope is
part of the API, not an implicit assumption about an optimized balancing actuator.

Request status is feasible, violated or unresolved. Plant status is feasible
only with a direct requested-point witness or a physically accepted allocation;
otherwise it is unresolved. A negative rail lower enclosure alone cannot prove
a violation: a sampled violating witness is required. Singular filter inversion
also remains unresolved. An unsuccessful local solve never becomes a physical
infeasibility assertion. The result preserves allocation solver/physical diagnostics.

The 91-assertion regression suite includes an analytic radial VA projection at
three numerical bases, an omitted-current-rating case at 100 MVA, voltage-rail
conflict and reference-angle invariance, primitive mutual/neutral LCL transfer in
both arm directions, fixed split-link charge, zero-PV loss import, solver failures
with/without a separate feasible witness, singular transfer and input contracts.
The allocated point independently verifies fixed POC voltage and midpoint charge;
a physically feasible point at a shifted boundary cannot publish as a projection.
`scripts/formulations/inverter_current_projection.jl` reproduces three accepted
operating points on one fixed synthetic plant. Its committed JSON snapshot
contains complete hardware, source fingerprints and individual physical checks.
The final frozen source passes the complete offline package suite:
19,844/19,844 assertions (7m28.4s). An earlier mixed-source test run was discarded;
this complete run includes the final boundary acceptance helper and its tests.
The companion Matplotlib figure shows the balanced-current capability slice
and the requested/allocated instantaneous line-voltage envelopes. It identifies
the synthetic fixture separately from the unchanged Australian sizing.

An additional optional-current-rating probe with an unbalanced LCL plant accepts
SI, 1 MVA and 100 MVA runs, but returns `ALMOST_LOCALLY_SOLVED` at 10 kVA under
the same tolerance profile. The relaxed result remains unpublishable; this probe
does not resolve the outstanding optional-arm conditioning audit.

This completes the offline reference step, not the full voltage-aware control
requirement. Remaining work includes online reactive-priority voltage allocation,
controller-domain handling, optional-arm SI conditioning/DC-source boundary audit,
the controlled-unbalance/operating-mode/chronological Australian studies,
coordinated Australian dispatch, four-leg control and representative EMT.
The goal remains active and Australian inverter sizing is unchanged.

## Refreshed physical-acceptance checkpoint

Source commit `c0226a3` passes the complete offline package suite: 19,753/19,753
assertions (7m30.8s). Fresh main (54), sensitivity (80) and matched terminal (18)
snapshots all pass. All 14 campaign environments, including the separate diagnostic,
match the same 50 source hashes. Per-attempt physical reports/checks are retained.
The independently replayed modified-network audit passes 380/380 checks.

`experimental_terminal_scaling_current` separately records the current 1 MVA,
300-iteration profile: mean VV+VW reaches `ITERATION_LIMIT`, while worst-phase
and sequence VV+VW pass. No selected-subset scientific comparison is published.
Historical diagnostics retain their original source metadata; an older sequence
failure must not be described as the current failure. Neither implies physical
infeasibility.

A preliminary fixed-terminal offline allocation experiment also found a feasible
plant behind an infeasible current request. At balanced 253 V, a synthetic 620 V
DC-link plant cannot deliver the requested 4.8 kW unity-PF point (rail margin
-26.907 V), but a local physical projection delivers about 1.886 kW / -5.956 kvar
within the same hardware limits. This is a feasible witness, not a globally nearest
projection or a production allocator. The Australian 750 V hardware is unchanged.
The next implementation step is a reusable offline voltage-aware allocation
reference with explicit request-versus-plant diagnostics.

## PWM publication gate follow-up

The fundamental/continuous-rail checkpoint is commit `6eb7ee9`. Its complete
package suite passes 19,723 tests. Its matched terminal data and 380-check replay
audit retain their exact source fingerprints.

The next audit reproduced a distinct carrier-closure bug: with a one-iteration
budget, a 20 kVA three-leg plant declared a 0.5 A capacitor rating, reserved
0 A carrier current and predicted 16.763 A switching RMS, yet still published
its successful fundamental NLP. The outer loop returned its last point without
certifying the carrier model.

PWM-enabled publication now requires a closed reserve, finite predictions and
modulation headroom, and the actual predicted capacitor-thermal and AC total-RMS
bounds. Capacitor checks replace the allocated carrier channel with the predicted
one while preserving fixed allowances, source-dependent 2omega current and
split-bank fundamental shares. Nonfinite harmonic networks and iteration-budget
exhaustion retain the candidate/MOI diagnostics and reject publication. If the
numerical closure tolerance is met but a physical rating is still violated,
reserve updates continue; the tolerance cannot relax an engineering bound.

Regressions cover the one-iteration counterexample, converged dispatch with the
same hardware, a lossless harmonic antiresonance, and a tight thermal rating
below the numerical closure floor. Existing four-leg, split-bank, finite-source
and AC-current reserve tests now assert physical publication. The affected
plant/certificate/PWM suite passes 741 assertions.

Still pending: optional-arm SI conditioning and the DC-source average-power
boundary; voltage-aware command allocation and its offline reference; refreshed
complete campaign artifacts; controlled unbalance/operating-mode/chronological
studies; coordinated dispatch; four-leg control; and representative EMT. The
active goal remains open and sizing parameters remain unchanged.


## Active goal checkpoint: physical acceptance and continuous rails

Work remains on `codex/ibr-pcc-apparent-power`; inverter sizing values are
unchanged. The active objective still includes the full engineering audit,
voltage-aware allocation, stronger case studies, coordinated dispatch, four-leg
control and representative EMT validation. This checkpoint does not complete
that objective.

Implemented in this batch:

- Independent SI reports recompute current/sequence limits, conductor-domain
  filter equations, selected-boundary power, internal/GFM voltage, DC ripple,
  split-bank midpoint/thermal budgets and reported-power consistency. Wrappers
  add DC availability, POC P/Q bounds and actual command-current residuals.
  Solver `feasible`/`optimal` evidence remains separate from the publication
  gate. Failed candidates can retain diagnostic reports.
- Capability and linking rows default to engineering normalization. Single-phase
  ripple/modulation rows are included. Unsupported single-phase DC-bus/source
  and neutral bounds now reject input rather than being silently ignored.
  Ipopt wrappers default to zero bound relaxation, with explicit caller options
  preserved and other optimizers' option vocabularies respected.
- Full-cycle rail headroom is enclosed by adaptive Taylor bounds for each signed
  degree-two trigonometric inequality. Standalone solves add separating angles
  and repeat, for at most 12 refinements. The physical gate rejects an unresolved
  rail violation. This is not yet a voltage-aware command allocator.
- Filter starts satisfy linear KVL/KCL and the three-leg zero-current sum.
  Constant-service and zero-availability PCC graphs initialize auxiliary
  definitions and quadratic roots in dependency order. The balanced test fixture
  now uses exact 120-degree angles; a separate common-mode case tests the strict
  LCL balanced-grid incompatibility rather than accepting it through tolerance.
- At zero DC availability, total-PCC loss accounting implies mandatory loss
  import and identically zero active PV service. The smooth law now uses that
  exact identity and leaves reactive/unbalance support active. The positive-
  sequence legacy offset can contain signed I2 power and is not simplified by
  the PCC identity.
- New study runs save per-attempt `physical_reports.csv` and `physical_checks.csv`,
  including rejected attempts and raw solver status. Source fingerprints include
  new/untracked Julia files. Generated CSVs are cleared before a campaign so an
  early rejection cannot leave stale accepted metrics.

Validation at this checkpoint: 189 new certificate/enclosure/initialization
assertions, the affected inverter suites, and 380 Australian replay-audit checks.
The matched 18-point terminal sweep accepts all points with the fixed 100 kVA
base/profile and no retries: 216 device reports and 6,264 individual checks pass.
Publication artifacts will retain the exact source fingerprints used at runtime.
Intermediate rejected runs remain diagnostics; none establishes physical
infeasibility on its own.

Remaining work, in dependency order:

1. Complete the engineering ledger beyond this fundamental-model certificate:
   test nonconverged/nonfinite PWM reserve closure publication, actual carrier
   RMS versus declared bounds, optional-arm SI conditioning, and the DC-source
   impedance/average-power boundary. Preserve the existing source-model oracles.
2. Build the offline closest-feasible converter allocation reference, then add
   voltage/modulation-aware allocation and explicit plant-versus-policy domain
   diagnostics. Adaptive rail cuts only enforce the plant; they do not redesign
   the local command when that command is infeasible.
3. Refresh the complete main/sensitivity artifacts from a committed source
   checkpoint. Extend matched studies with VUF angle/magnitude/common mode and
   phase permutations, disconnected versus energized operation, and chronological
   profiles before energy claims. Keep hardware fixed.
4. Add coordinated dispatch using the same plant and engineering checks. A local
   nonconvex OPF solution is a feasible benchmark, not a certified global upper
   bound; any claimed upper bound needs separate relaxation/global evidence.
5. Develop four-leg per-phase control and split-bank protection, then validate
   representative binding cases with parameter-identical EMT. Retain the
   distinctions between the fundamental certificate, carrier replay and EMT.


## 31 October 2026: engineering limits and next priorities

The validated thermal/PCC study checkpoint is commit `80f6d8d` on
`codex/ibr-pcc-apparent-power`. Subsequent work deliberately leaves inverter
sizing values unchanged.

The next audit found a real hard-bound counterexample: a declared 0.4 A
negative-sequence rating returned about 2.19 A with `LOCALLY_SOLVED`. Raw
per-unit squared-current rows and their auxiliary links were not sufficient
physical acceptance evidence. Sequence auxiliaries and linking rows now use
their engineering ratings as scales. Independent fixed-current tests reject
contradictory I₀/I₁/I₂ ratings and accept feasible counterparts at 10 kVA,
1 MVA and 100 MVA numerical bases.

Both controller backends now allocate converter-arm I₁/I₂ limits, including
loss import and LCL shunt offsets. Common scaling retains its ray; legacy
sequence priority clips the independent sequence services before coupled
protection. PCC allocation retains its service order and avoids an artificial
I₂ intersection on analytically balanced P/Q directions after balance-first
saturation. Equal-valued ratings are tested by sequence identity, not value.
Fleet rows expose the declared ratings and remaining margins in amperes.

The new teaching figures distinguish exact local-law illustrations from
accepted feeder equilibria, separate losses from unused PV power, show physical
headroom utilization, and present the three matched reactive-priority terminal
controls together. Hardware sizing and benchmark inputs remain unchanged.

Validation: 225 new engineering-limit assertions; 2,352 assertions across 106
affected model/controller testsets in total. The Australian replay audit passes
128 checks. Earlier rejected solver/base probes remain diagnostic evidence.

Recommended order for further work:

1. **Complete the bound/acceptance ledger.** Audit every engineering quantity,
   its terminal and RMS/peak convention, stamping scale, allocator coverage and
   independent post-solve residual. Extend the sequence-bound scaling lesson to
   all optional current, voltage, capacitor-bank and DC-source constraints.
   Treat a solver status as necessary evidence rather than a physical certificate.
2. **Close voltage-protection and controller-domain gaps.** Certify continuous
   rail feasibility with an analytical peak/separation oracle or adaptive cuts.
   Build an offline closest-feasible allocation reference before extending
   voltage/modulation protection: reducing Q is not always protective because
   reactive absorption can reduce internal voltage. Distinguish a physically
   infeasible plant from an infeasible strict service or mandatory-offset anchor.
3. **Strengthen numerical and feedback screening.** Guard auxiliary-scale
   domains, test independent starts and unit bases, and retain failed probes.
   Extend the reduced loop screen to include U₀ voltage response: the three-wire
   actuator has four real current coordinates, but phase-magnitude sensors see
   all six voltage coordinates. This remains an equilibrium/iteration screen,
   not proof of dynamic converter stability.
4. **Make the study more causal.** Keep hardware fixed while crossing source
   VUF magnitude/angle, common-mode voltage and load-phase permutations. Add
   explicit disconnected versus energized zero-PV operating modes; availability
   alone does not encode inverter status. Then add chronological quasi-static
   profiles before making energy claims. Add
   a matched single-phase baseline with consistent losses, filter and terminal
   semantics before using mixed-fleet adoption cohorts as a realistic comparison.
   Then test additional real LV feeders and loaded MV+LV coupling.
5. **Add useful capability comparators.** Use the same physical plant to compute
   a coordinated OPF upper bound, with a clearly stated objective and service
   constraints. Add four-leg per-phase control for independent I₀ actuation,
   followed by split-link midpoint/bank protection. Derive LCL/asymmetric passive
   sequence admittances from the full circuit rather than assuming a scalar Z.
6. **Validate higher-fidelity mechanisms.** Use parameter-identical EMT cases
   for representative binding limits. Two-stage DC/DC source coupling, DC input
   current/voltage envelopes and explicit PV/no-charging semantics are useful
   model extensions after the preceding contracts are clear. Parameter fitting
   and hardware sizing optimization are outside this batch.

## 31 October 2026: terminal-profile fleet conditioning

Work continues on `codex/ibr-pcc-apparent-power`. The three stressed positive
voltage-control cases now pass, together with an 18-case sweep covering all
three sensing policies, both source levels and zero/half/full PV. Every case
uses the same 100 kVA base, adaptive/no-automatic-scaling Ipopt profile and
300-iteration cap, without retries.

The scalar capability allocator now keeps squared norms as expressions and
normalizes the root of larger lifted curve/thermal graphs by service magnitude.
It rescales the residual and smooth-maximum width together, preserving the
quadratic, requested smoothing and service order. The maximum is evaluated
from a strictly positive anchor to prevent cancellation on trial points.
Constant-service graphs retain their established root coordinates. No curves,
ratings, protection bounds, P/Q priority or acceptance tolerances were relaxed.

Validation: 2,127 assertions across 102 model/controller testsets, including the
18 operating points and stressed selection-objective invariance. The network
audit passes 128 checks, including independent fixed-PQ replay of all three
terminal equilibria (maximum voltage difference below 0.2 µV).

Earlier rejected positive-profile runs are retained under
`results/experimental_terminal_var_priority_before_fix`. Conditioning remains
sensitive to numerical configuration: the 1 MVA adaptive-profile sequence case
still fails and is retained in `results/experimental_terminal_numerical_scaling`.
The full negative-sequence-controller PCC fleet has not been revalidated by
this positive-voltage-control campaign. Neither historical nor current failed
probes establish physical infeasibility.

## 30 October 2026 follow-up: terminal semantics and thermal allocation

Historical status before the 31 October conditioning fix:

Mean and worst-phase equivalence fixtures now share reactive priority,
converter-current actuation and PCC apparent-power limiting. This is a named
Australian voltage-control comparison. The established proportional-PQ/grid-
target/converter-VA full-controller baseline is retained: attempted promotion
of PCC semantics to negative-sequence fleet cases exposed unresolved numerical
failures, preserved in `results/experimental_pcc_fleet`. They are not evidence
of physical infeasibility. The stressed terminal-profile sensitivity also fails
for mean VV+VW; sequence VV+VW converges with the default solver but fails the
selected adaptive profile. Both profiles are retained as numerical diagnostics,
and no partial matched comparison is published. Under PCC reactive priority, common and power-priority
constructors share Q→P→I₂ ordering; names alone are not distinct policies.

Three-leg controller backoff now includes the plant's capacitor thermal budget:
source-dependent 2ω sharing, frequency weights and fixed/current-dependent
carrier allowances. An implied looser thermal bound needs no additional
controller intersection, while the plant continues enforcing both bounds.
PCC curve outputs and service directions are lifted for the larger nonlinear
graphs to prevent repeated expression-tree expansion in fleet solves.

The campaign adds independent capacitor ratings, balanced-source/load controls,
half-width smoothing and a matched-gain/angle passive-admittance comparison.
It records runtime provenance and checks dense switching, VA and thermal
margins. Australia A curve values now cite public Energex tables; the capacitor
paper's author list and public two-stage-model reference are corrected.

Still open: manufacturer-calibrated losses/filter/DC-source data, full
gain/angle/source-VUF sweeps, phase permutations, mixed fleets, loaded MV+LV
integration, coordinated OPF bounds and parameter-identical EMT validation.
These updates do not establish certification, dynamic GFM performance or
representative Australian daily energy outcomes.

Validation: 2,013 targeted model/controller assertions across 101 testsets;
80 network/provenance audit checks. The main 54-case campaign and 77 cases in
ten completed sensitivity campaigns pass their physical gates. The separate
three-case terminal-profile diagnostic accepts one case and rejects two.

Branch: `ibr-review-fixes` (from `main` @ `2c49436`). This note tracks every
finding of the October 2026 critical review of the advanced IBR model, the
phase-aware controllers and the two network studies. Each item stays here until
it has code, a regression test and documentation, or an explicit decision to
defer it. Update the status column in the same commit as the change.

Status: `todo`, `doing`, `done`, `deferred` (with reason).

## A. Correctness

| ID | Finding | Plan | Status |
|---|---|---|---|
| A1 | DC power availability is never enforced. `P_dc` is unbounded, so a controlled PV supplies its own losses from nothing: +57 W at 5 kW available with typical loss coefficients, and 20 W at zero PV with POC P = 0. The docs list DC availability as a protection priority. | The controller limits the positive-sequence active-power command to `p_available - L`, where `L = P_dc - 3Re(U1 conj(I1_target))` is the solved plant offset (all losses plus negative-sequence active power). The exact evaluator receives the same offset from the solved plant, as the existing capability backoff does. Report `dc_availability_margin_W = p_available - p_dc` in the device rows. At zero PV the inverter then imports its losses from the grid. Regression: lossy plant at full and zero PV. Done: `dc_power_offset`, `p_dc_W` and `dc_availability_margin_W` rows; the oracles refresh the offset from each plant re-solve. Follow-up: the study gate caught a grid-forming night case overdrawing 4–8 W because the backoff scaled the loss import. The import is now an unscaled offset through every limiter (1e-3-of-rating selector, one-sided bias), and the `dv2_max` backoff uses the line-frequency source factor. Regression added. | done |
| A2 | The common-scale limiter scales down I1 when I2 saturates a phase, so negative-sequence support curtails active power. This confounds the sequence-droop comparison. | Add `SequencePriorityLimiter(priority=:power Done: `current_scale` reports the I1 scale; rows add `negative_sequence_scale`. \| :balance)`. The reserved sequence is limited alone; the other sequence gets the largest exact per-phase scale `t` with `\|a + t b\| <= I_max`, using the rationalised quadratic root. The smooth model overestimates the root, so it stays conservative. Plant-aware backoff remains the final common protection layer. Exact and smooth tests, plus a network test showing power priority keeps P where common scaling curtails it. Follow-up: the corrected study showed the binding 2ω ripple limit in the plant-aware backoff still common-scaled both sequences (power priority kept only 25.5 of 58.5 kW). The backoff now honours the priority in the stamp and the exact oracle; regression added. | done |
| A3 | The split-DC model is internally inconsistent: neutral current divides by `C_u/C_Σ` (stiff source at ω) while the 2ω and carrier models assume an open source. Only unequal banks are affected. | Use one DC-source model at every frequency. Derivation: `i_s = ½J(C_l-C_u)/(C_Σ + jωC_uC_l Z_s)`, `i_cu = i_s - ½J`, `i_cl = i_s + ½J`. This also gives a fundamental total-bus ripple that enters the half-rail hull. Equal banks are unchanged for every `Z_s`. Decide on parameter naming against the existing `pwm_dc_source_*` fields during implementation. Done: new `dc_source_r`/`dc_source_l` (line-frequency branch, default open) govern both the 2ω ripple and the split fundamental; the carrier branch stays `pwm_dc_source_*`. Added `dv_bus_fundamental` and the 2ω midpoint term to the hull and all audits. Behaviour change: unequal banks now default to equal neutral shares (open source). The old stiff-ω/open-2ω pairing is removed as inconsistent; `dc_source_r=0` gives the fully stiff limit. A nodal-circuit oracle covers four source models. | done |
| A4 | Sensing reference is tied to wire count: `neutral=nothing` forces phase-to-ground sensing, although many three-wire products sense phase-to-neutral. | Add `sensing_neutral` to `SequenceController`: `:device` (current behaviour), `:ground`, or a bus terminal name. Route the smooth stamp, the exact evaluation, the fixed-point oracles and the sensitivity helpers through one helper. Test the common-mode metamorphic contract with PN sensing on a three-wire plant. Done as `sensing_reference` (`:device`, `:ground` or a terminal name) plus the exported `sensing_terminal`. Also documented: the mean-phase law responds to U0 only at second order. | done |
| A5 | The `:dominant` conflict rule changes the Volt-var curve without a conflict: near the deadband edge q ∝ s², up to about 0.15·ε. | Scale the transition width by a conflict measure `ρ = min(l,h)/(min(l,h)+ε)`. The firmware law then reproduces the standard curve exactly when only one branch is active, while keeping continuity and zero at the tie. The smooth model uses the zero-preserving `nonnegative_min`. | done |
| A6 | Docs say the positive-sequence policy is invariant to U0. Only its Volt-var is; the default worst-phase Volt-watt guard is not. | Correct the wording in the control-law page, the `SequenceController` docstring and the Australian README. | done |
| A7 | `screen_fixed_point_gain` reports only Picard criteria (spectral radius, ∞-norm). Those are stricter than needed for lagged/incremental droop. | Also report `max_real_eigenvalue`. The criterion for a small-step relaxed iteration is `Re λ < 1`; document both. Review finding withdrawn: `maximum_real_eigenvalue`, `continuous_time_margin` and `alpha_max` already existed. Re-reading found `alpha_max = 2/\|1-λ\|`, which is wrong for complex eigenvalues; fixed to `2Re(1-λ)/\|1-λ\|²` with a rotation regression. | done |
| A8 | Found while verifying C2: `grid_forming=true` on `:THREE_LEG` imposed a ground-referred balanced EMF. With I0 = 0 that forces the POC zero-sequence voltage to zero, so any network with zero-sequence unbalance was infeasible. | For a three-leg bridge, constrain only the internal negative-sequence EMF to zero and bound the positive-sequence magnitude; U0 follows the network. Regression: plant I2 = -U2/Z_f on a source with V0 ≠ 0. Docstring and foundations §10 updated. | done |
| A9 | Found while preparing B5: network voltage sensitivity used 6×6 phase coordinates for a three-leg grid-current target. Independent phase perturbations violate ΣI = 0, so the perturbed plant solves cannot succeed. | Use sequence coordinates for every three-leg target. Regression: Z2 self-impedance of the resistive fixture at both targets. | done |
| A10 | Upstream (BMOPFTools `solve_pf`): `_strip_operational_limits!` deletes IBR `s_max`, which zeroes S_MAX/VAR_MAX droop bases. Native Volt-var/Volt-watt then gives zero output in a power flow, contrary to the docstring's claim that profile-governed IBRs remain determined. | Not fixed here (other repository). Report upstream. The status-quo study case (C3) uses the OPF builder with bus voltage bounds removed and a maximum-export objective. Filed as [BMOPFTools.jl#447](https://github.com/frederikgeth/BMOPFTools.jl/issues/447), with a minimal reproduction and the related Volt-watt determinacy gap. | deferred (upstream) |

## B. Parameters and studies

| ID | Finding | Plan | Status |
|---|---|---|---|
| B1 | The Australian Volt-watt curve (250→1, 260→0 of available) is stricter than AS/NZS 4777.2:2020 Australia A (253→1, 260→0.2 of rated). It drives the headline curtailment result. | Default to the Australia A curve on a rated basis; keep the pedagogical curve as a named sensitivity variant. Re-run and update the README tables. Australia A (BMOPFTools `Aus_A` preset) on a rated basis is the default; the pedagogical curve is the `sensitivity_pedagogical_volt_watt` run. | done |
| B2 | 2 mF for 6 kVA with no `dv2_max` means DC-link stress never binds. | Use a representative film link (about 50–80 µF/kVA, stated as an assumption) and declare `dv2_max` and `i_cap_max`. Report 2ω ripple and capacitor current. 300 µF (50 µF/kVA) with `dv2_max` = 7.5 V; 150 µF and 2 mF sensitivity runs. | done |
| B3 | The filter is fixed at 0.05 + j0.15 Ω (about 0.6 % per unit for 6 kVA, roughly 10× below typical) and is not scaled with rating in ENWL. | Specify filters in per unit of each device's rating (stated assumption, e.g. X = 5 %, X/R = 10) in both scripts. 5% series reactance at X/R = 10 in per unit of each rating (both scripts); per-unit loss curve added. | done |
| B4 | The negative-sequence gain κ is in A/V and not rating-scaled. ENWL's κ is a no-op (g ≈ 0.02). | Both scripts declare κ in per unit of `I_rated/V_nom` and convert per device. Calibrate the ENWL gain to the Australian per-unit value. Australian κ is derived for 20% aligned fleet attenuation and recorded in per unit. The first design divided the fleet-aggregate impedance by N twice, making κ 12× too small; this was caught in the results and corrected to κ = 1.63 A/V (43 pu). ENWL uses the same per-unit value, marked as an assumption. | done |
| B5 | `impedance_angle = atan(0.3)` is not linked to the feeder. | Measure the local negative-sequence Thévenin angle with `inverter_control_network_voltage_sensitivity` (or a direct Z2 estimate). Record it per customer and choose `phi2` from it. Keep the fixed angle only as a sensitivity. Measured with the fleet sequence sensitivity. Aggregate Z2 angles are 10–21° (median 19.4°): the feeder is resistive, so the earlier 16.7° angle was reasonable. The review was wrong to expect a transformer-dominated angle; what it lacked was the measurement. | done |
| B6 | ENWL variants use different smoothing widths, which confounds a matched comparison. | Use one smoothing configuration across all variants. ENWL variants now use default smoothing throughout. | done |

## C. Comparison set

| ID | Missing comparator | Plan | Status |
|---|---|---|---|
| C1 | Power-priority sequence allocation | Covered by A2; add as an Australian study mode. Mode `sequence_droop_power`. | done |
| C2 | Grid-forming passive negative-sequence admittance (balanced EMF behind `Z_v` gives `I2 = -U2/Z_v`) | Constructor `VirtualImpedanceNegativeSequence(r, x)` for a flat admittance at `arg Z_v`, documented as the steady-state negative-sequence equivalent of a balanced-EMF voltage source. Add a study mode. Done as `negative_sequence_virtual_impedance(r, x)`; verified against a grid-forming plant solve. | done |
| C3 | Status quo: single-phase PV with per-phase Volt-var | Australian study variant keeping native single-phase BMOPF PV (PN_PER_PHASE) at the original customer phases. Mode `single_phase_native`: native BMOPFTools single-phase PV on the original phases with the same loads, solved around upstream bug A10. | done |
| C4 | Four-leg per-phase Volt-var | Needs a FOUR_LEG controller with I0 enabled (new policy and stamp). Larger API change; planned as a follow-up PR unless time allows. Deferred to a follow-up PR: needs a FOUR_LEG controller with I0 enabled. | deferred (follow-up PR) |

## D. Literature

| ID | Gap | Plan | Status |
|---|---|---|---|
| D1 | No BPSC/PNSC/AARC/IARC taxonomy | Cite Wang, Duarte & Hendrix (TPEL 2011) and map every implemented and candidate law onto it. | done |
| D2 | Peak-current-constrained voltage support in R–X grids | Cite Camacho et al. (TIE 2013; TPEL 2018) and Afshari et al. (TIE 2017); position A2 against them. New references are numbered 54–59 in a separate section so existing item numbers stay stable. | done |
| D3 | LV three-phase inverter unbalance predecessor | Cite Meersman et al. (EPSR 2011) as the predecessor of the conductance (`phi2 = 0`) law. | done |
| D4 | Grid-code negative-sequence injection | Add k-factor context (VDE-AR-N 4110/4120, NC RfG), marked unverified like the other standards. | done |
| D5 | Grid-forming natural unbalance response | Add a paragraph in the foundations and control-law pages connecting `grid_forming=true` and C2. | done |

All new bibliography entries carry the existing † "not checked against live
sources" marker until verified.

## E. Documentation: parameter justification and code sync

No page justifies the hardware or control parameters. The studies label values
"illustrative" or "pedagogical" without sources, ranges or sensitivities, and
the docs-to-code check is manual.

| ID | Gap | Plan | Status |
|---|---|---|---|
| E1 | No parameter-selection guidance | New page `docs/src/ibr/parameter_selection.md`. Cover: filter impedance in per unit (L and LCL ranges, X/R); DC-link voltage by topology and the modulation limit; DC-link capacitance (film vs electrolytic, µF/kVA, ripple-current rating); current headroom over rated current; loss-curve coefficients from a datasheet efficiency curve; AS/NZS 4777.2 Volt-var/Volt-watt settings by region and basis; κ from the target `g = κ\|Z2\|` and rating; `phi2` from a measured Z2 angle; smoothing widths from stated error budgets. Each value gets a source (marked † until verified) or is explicitly a stated assumption. Done as `docs/src/ibr/parameter_selection.md`; the Australia A values are cross-checked against the BMOPFTools `Aus_A` preset. | done |
| E2 | Study parameters are not justified per study | Both study READMEs get an assumptions table: value, source or rationale, and the sensitivity run that covers it (links B1–B6). Assumptions table in the Australian README; per-unit assumptions listed in the ENWL README. | done |
| E3 | Docs made stale by this branch | Update as each item lands: A1 changes the `:rated` basis text (`inverter_controls.md` ≈141), the study boundary (≈596), the protection-priority text in `phase_aware_control_laws.md` (≈192, ≈492) and the Australian README "drawing active power for losses" (now true and stated precisely). Done alongside each item; stale statements updated in the control, foundations, verification and study pages. | done |
| E4 | No systematic docs/code sync check | Audit every IBR page against the code: API names, defaults, keyword lists and "implemented/future" status tables. Hard-coded numeric result tables (study READMEs) must be regenerated from the committed result CSVs, not hand-copied. Record the audit result here. Done: a script checked every backticked identifier on the IBR pages against package names, fields and keywords (no drift found); verification rows were updated for the new oracles. Five table cells whose inline pipes broke rendering were fixed and confirmed in the docs build. Study README tables are generated by `summarize.jl`. | done |

## Order of work

1. A1, A2 (largest effect on study conclusions), then A5–A7.
2. A4, A3.
3. C2, then D1–D5.
4. E1 (parameter guide), then B1–B6 plus C1–C3 in the study scripts with E2;
   re-run the Australian study, and the
   ENWL dry run (its data checkout is external).
5. E3 continuously with each item; E4 audit before opening the PR.
6. C4 as a follow-up unless scope allows.


## Outcome

Study results after all fixes (Australian LV3, full PV, 1.03 pu source, 1%
source VUF) are generated in `studies/australian_pv/README.md`. All 54 main
and 27 sensitivity cases are accepted and the replay audit passes 80/80. The
review's main concerns are borne out:

- the headline curtailment came from the pedagogical Volt-watt curve;
- unbalance attenuation is limited by the DC link (capacitance sensitivity);
- common scaling gives up most of the active power for that attenuation,
  while power priority (extended to the plant backoff) keeps it;
- the single-phase status quo has the worst unbalance and neutral displacement.

One review claim was withdrawn: the measured negative-sequence impedance angle
is resistive (about 19°), so the original fixed angle was reasonable.

## PR #61 merge-review follow-up

- **Stiff DC source:** skip the voltage-ripple capability backoff when the
  source factor is zero. An ideal stiff source then has zero bus ripple without
  introducing an infinite limit into either the smooth or exact controller.
- **Physical loss import:** use the unregularized positive-sequence voltage
  denominator for the mandatory import, keeping the declared voltage floor on
  the service command. The actual imported power now covers the DC shortfall
  even with a 20 V service floor. A nonzero required import at zero
  positive-sequence voltage is rejected by the numeric evaluators.
  Controller stamps use four sequence-current equalities when the plant already
  enforces zero current sum, avoiding two dependent phase-current equations.
  LCL grid targets retain all phase equations because their shunt can carry I0.
  Controlled solves with four independent commands omit unused lossless
  current-magnitude epigraphs; standalone and LCL grid-target solves retain
  their legacy treatment. This removes free barrier directions at zero dispatch.
- **Reduced loop gain:** reconstruct voltage perturbations with the measured
  zero-sequence offset held fixed, for both single-device and fleet screens.
  The reduced screen still omits changes in that channel; this limitation is
  stated in the documentation.
- **Julia 1.10 docs:** include the registered documentation dependencies in the
  bootstrap resolution before adding the local package. This allows Documenter
  to select compatible transitive dependencies instead of inheriting an
  incompatible DocStringExtensions version from an unrelated seed environment.

`test/inverter_control_merge_review_tests.jl` covers the three controller
reproductions, both current targets, common and power-priority limiting, and
independent finite differences at nonzero sensed zero-sequence voltage.
It also checks the physical DC budget through the import transition under
every P/Q priority and binding current limits.
