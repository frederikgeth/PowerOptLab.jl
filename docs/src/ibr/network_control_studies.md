# Network-scale control studies

`solve_controlled_inverter_fleet` is the single-snapshot bridge between the
phase-aware controller/advanced-inverter model and prepared BMOPFTools network
datasets. It is intended for comparative studies of local control laws,
per-leg current requirements, and DC-link stress across many independent LV
network snapshots.

## Scientific scope

The solver answers a deliberately narrow question: given one network state,
one physical inverter design per selected IBR, one locally measured
voltage-phasor control law, and one available-power request, what simultaneous
network and inverter operating point satisfies all algebraic equations and
physical limits?

It is a controlled power flow, not model-predictive control and not a local
optimization. Every controller is a fixed smooth computation graph. Ipopt is
used to solve the coupled network/device equations at scale; it does not choose
the control policy. Independent time, penetration, and uncertainty snapshots
belong in an outer study loop so failures and limit activation remain visible.

The initial fleet implementation supports selected dataset IBR records with
`THREE_LEG` or `FOUR_LEG` AC port descriptions, each exposing three ordered
phase conductors. Their PowerOptLab replacements are three-leg
`AdvancedInverter`s. Native DC-network-coupled IBRs are rejected because
replacing their AC formulation would otherwise remove native DC balance
physics.

Selecting a native `FOUR_LEG` record therefore means deliberately studying a
three-leg replacement at the same four-wire connection; it is not a claim that
the hardware topologies are equivalent. A native `THREE_LEG` record is accepted
only when the replacement has `neutral=nothing` and has neither an LCL midpoint
shunt nor a POC shunt. These conditions preserve the native line-to-line port:
the replacement has zero summed phase current and introduces no phase-to-ground
or neutral current path.

## Replacement, not parallel stamping

The dataset continues to identify an inverter by `network["ibr"][id]`. A
`ControlledInverterFleetSpec` maps that same identifier to a
`ControlledDevice` and an `InverterControlRequest`:

```julia
spec = ControlledInverterFleetSpec(
    Dict("pv_17" => ControlledDevice(advanced_inverter, controller)),
    Dict("pv_17" => InverterControlRequest(
        p_available=6.4e3, p_rated=8e3, q_scale=8e3)))

result = solve_controlled_inverter_fleet(network, spec)
```

PowerOptLab assigns each selected `(:ibr, id)` to a per-component
`OpfDeviceBuilder`. BMOPFTools therefore omits the native IBR equations and KCL
injection for that identifier while retaining native physics for every
unselected IBR. The advanced plant adds the only live terminal injection.

BMOPFTools declares its semantic native IBR current variables before ownership
is resolved. PowerOptLab fixes those unused placeholders to zero. Published
results then follow one unambiguous ownership rule:

- `result.devices[id]` contains each replaced, controlled advanced inverter;
- `result.network["ibr"]` contains only unselected native IBRs; and
- `result.build_manifest.component_owners[(:ibr, id)]` records `:PowerOptLab`
  or `:BMOPFTools`.

This contract is stronger than inferring replacement from a coincident bus or
from a near-zero duplicate injection. A selected native record with tight
native power bounds cannot constrain the replacement plant. The number of
placeholders fixed is derived from BMOPFTools' resolved dataset neutral labels,
not assumed to equal the replacement's three phases. This also closes every
unused variable when a non-conventional neutral spelling makes BMOPFTools
declare an additional semantic IBR-current pair.

One upstream qualification remains: if BMOPFTools does not recognise that
neutral spelling, its native `:network_cost` construction also treats the
terminal as a fourth phase. A costed record must consequently provide four cost
entries or BMOPFTools rejects it while indexing the objective. This affects
native and replaced IBRs alike; PowerOptLab closes the extra unused current
pair, but does not reinterpret native cost data. Prefer BMOPFTools-recognised
neutral metadata/spelling for production datasets.

## Input contract

For each selected identifier, construction validates before model mutation:

1. the device and request dictionaries have identical nonempty key sets;
2. the dictionary key equals the wrapped advanced-inverter identifier;
3. the identifier exists in the network's native `ibr` collection;
4. native and replacement bus identifiers agree;
5. phase order agrees exactly, and a `FOUR_LEG` record places the matching
   neutral after the three phases; and
6. a `THREE_LEG` record is paired only with a neutral-free, shunt-free
   replacement; and
7. the native IBR has no DC-network coupling.

Phase order is a physical convention, not presentation metadata: it determines
the symmetrical-component transform and the reconstructed phase currents.
Silently sorting phase labels would change the control law.

The controller formulation is currently supported only with `per_unit=true`.
Configuration and extracted results nevertheless remain in SI units. The
advanced inverter must use `pwm_strategy=:NONE`; calibrated PWM reserve terms
may be supplied explicitly, but the separate PWM outer iteration is not part of
the fleet solve.

## Objective semantics

Local controller equalities determine their current commands. The objective
only selects freedom left elsewhere in the algebraic model:

- `selection_objective=:loss` (default) minimizes selected inverter
  semiconductor and capacitor losses;
- `:zero` is an objective-invariance diagnostic; and
- `:network_cost` applies BMOPFTools' native generation-cost expression to
  dispatchable unselected network devices.

There is intentionally no silent weighted sum of economic cost and electrical
loss: the quantities have no universal commensurate scale. For comparative
control-law studies, pin exogenous loads, generation, and unselected IBR
setpoints in each prepared snapshot, then use `:loss`. If dispatch of native
devices is part of the experiment, state and use the economic objective
explicitly. Results from different objective semantics are not interchangeable.

## Result tables and aggregation

`ControlledInverterFleetResult` retains both the device-local results and the
ordinary BMOPFTools network result, plus a defensive copy of the study
specification used for provenance. `solve_status(result)` is the authoritative
publication gate. Do not drop infeasible or non-publishable snapshots; retain
their scenario identifiers and classify the failure mechanism.

The bus snapshot is canonicalized once per fleet result and shared by
`result.network["bus"]` and every `result.devices[id].bus`. Treat extracted
results as immutable: mutating that shared dictionary changes every view of the
same snapshot.

Two dependency-free extraction helpers produce table-ready named tuples:

```julia
device_rows = controlled_inverter_rows(result)
phase_rows = controlled_inverter_phase_rows(result)
```

`device_rows` has one row per inverter and scalar SI columns for POC and
converter power, losses, voltage extrema and sequences, VUF, requested power,
limiter scales, maximum converter/grid current, sequence current, 2ω ripple,
DC-link ripple voltage, capacitor current, and the exact-versus-smooth current
residual. `phase_rows` has one row per inverter and phase, retaining rectangular
and magnitude voltage, converter-current, and grid-current values. Both helpers
sort identifiers lexically, independent of dictionary insertion order.

The conventional `voltage_unbalance_factor` is ``|U_2|/|U_1|``. It is distinct
from `regularized_voltage_unbalance`, the controller's finite-denominator
quantity. Keeping both prevents the controller's numerical regularization from
being mistaken for a reported power-quality metric. It is also distinct from the
NEMA/IEEE maximum-deviation approximations that share the name "voltage
unbalance"; state which definition is reported (see Pillay and Manyage, and
Girigoudar and Roald on unbalance metrics inside optimization models, in
[IBR references](@ref ibr-references)), and note that neither this quantity nor
its regularized companion is measured over a standardised power-quality window.

`solve_diagnostics(result)` supplies lightweight fleet aggregates, including
total controlled POC P/Q, total selected-inverter losses, the largest
converter/grid phase current, and the largest exact-versus-smooth current
residual. Study scripts should add scenario keys—time, feeder, penetration,
random seed, and hardware/control configuration—before writing rows.

## Verification obligations

The unit tests exercise discriminative rather than merely feasible cases:

- selected native records have near-zero native ratings while their replacement
  plants deliver kilowatts, detecting accidental parallel native stamping;
- selected and unselected IBRs coexist, with ownership and result channels
  checked explicitly;
- average-voltage and worst-phase Volt-watt laws share one unbalanced network
  state and produce materially different requests;
- output ordering and aggregate identities are deterministic; and
- each solution satisfies three-leg zero-sum current, per-leg current limits,
  converter-terminal complex-power reconstruction, filter active-power balance,
  and DC-link power balance.

These are software and model-consistency checks, not validation of a hardware
controller. The study methodology still requires a stratified subset of
balanced, unbalanced, current-limited, and ripple-limited cases against an
independent averaged model and a switched EMT model. OpenDSS is suitable as a
network-voltage oracle for cases whose inverter injection can be frozen, but it
cannot independently validate the advanced converter, DC capacitor, or the
simultaneous smooth controller equations. Use it to isolate network stamping;
use averaged/EMT models for converter and control-law validation.

## Study workflow

For a publishable comparison, freeze a factorial design before inspecting
results: feeder/scenario, PV penetration and phase allocation, irradiance,
source-voltage unbalance, load realization, hardware rating, current target,
positive-sequence policy, negative-sequence policy, and limiter settings. Use
the same snapshot and hardware across controller variants.

Report at least:

- voltage compliance and phase extrema;
- conventional VUF and negative-sequence current;
- energy curtailed relative to the declared available-power trajectory;
- empirical quantiles and worst cases of converter- and grid-side per-leg
  current;
- frequency and identity of binding physical limits;
- 2ω ripple power, DC voltage ripple, capacitor RMS/thermal current, and
  capacitor stored energy ``C_{dc}V_{dc}^2/2``; and
- solver status, exact-versus-smooth residual, and sensitivity to smoothing
  widths on a representative subset.

Current and capacitor sizing conclusions require outer parameter sweeps. A
feasible point at one installed rating does not identify the minimum rating,
and a smooth-NLP failure is not by itself a physical infeasibility certificate.

## Matched batch experiments

`InverterControlStudyCase` and `run_inverter_control_study` implement the outer
single-snapshot experiment loop without coupling unrelated operating points into
one NLP. A case combines:

- `scenario_id`: the feeder, placement, load/PV, weather, and uncertainty
  realization that must remain common across variants;
- `variant_id`: the controller or hardware configuration being compared;
- one prepared BMOPF network and `ControlledInverterFleetSpec`;
- a positive interval duration in hours and positive scenario weight; and
- copied dataset metadata such as feeder, penetration, seed, and timestamp.

```julia
cases = [
    InverterControlStudyCase(
        scenario_id="feeder_01/time_0042/seed_7",
        variant_id="average_voltage",
        network=network,
        fleet=average_fleet,
        duration_h=0.5,
        weight=1.0,
        metadata=Dict("penetration" => 0.75, "feeder" => "feeder_01")),
    InverterControlStudyCase(
        scenario_id="feeder_01/time_0042/seed_7",
        variant_id="worst_phase",
        network=network,
        fleet=worst_phase_fleet,
        duration_h=0.5,
        weight=1.0,
        metadata=Dict("penetration" => 0.75, "feeder" => "feeder_01")),
]

study = run_inverter_control_study(
    cases; solver_options=("max_iter" => 500, "tol" => 1e-8))
```

The returned `study.settings` records the common per-unit base, optimizer type,
objective semantics, verbosity, and solver attributes. Retain this with exported
tables; it is part of the computational provenance, not an incidental runtime
detail.

Cases are solved serially in deterministic `(scenario_id, variant_id)` order.
Serial execution is the reference reproducibility path. Large campaigns should
partition the flat case list across worker processes and combine table rows;
there is no mutable cross-snapshot state in this layer.

By default (`continue_on_error=true`, equivalently `:validation`), a case-level
`ArgumentError` from dataset/model validation is retained with
`error_class=:validation`, its exception type/message, and the study continues.
Set `continue_on_error=false` to rethrow every exception. Programming,
configuration, and unexpected solver exceptions are otherwise rethrown: they
must not become apparent physical-failure statistics. For an unattended
campaign, the explicit `continue_on_error=:all` mode retains nonfatal unexpected
exceptions as `error_class=:unexpected`; these remain software/data failures
and are counted separately from validation failures and physical solve status.
Fatal process-level exceptions are always rethrown. A solve that returns a
non-publishable status is retained as a fleet result with masked numerical
data, distinct from a thrown error.

Extraction is split by purpose:

- `inverter_control_study_case_rows` contains every case, including errors and
  non-publishable solves;
- `inverter_control_study_device_rows` and
  `inverter_control_study_phase_rows` prefix fleet rows with the matched case
  keys, duration, weight, and metadata;
- `inverter_control_study_summary_rows` always groups by control/hardware
  variant, optionally adds metadata cohorts, reports failure fractions over all
  cases, and computes energy/tail metrics only from publishable cases; and
- `inverter_control_paired_rows` reports device-level `variant - baseline`
  differences for every scenario. Missing, errored, non-publishable, or
  device-mismatched pairs remain as rows with `publishable_pair=false` and NaN
  differences. Duration, weight, metadata, device identifiers, physical
  hardware, requests, and the shared network must agree across a pair;
  `matched_case_definition` exposes that check and prevents unmatched exposure
  or confounding from entering a paired estimate. A separately materialized
  network may be accepted only after external equivalence validation. The
  `baseline_*` and `variant_*` columns retain each arm's own achieved value
  whenever that arm published, even when the pair is refused, with
  `baseline_published`/`variant_published` naming which arms those are, so a
  dropped pair can be characterized rather than only counted. Only `delta_*`
  is gated on `publishable_pair`; do not difference the two arms of an
  unmatched row by hand, because the refusal is the finding; and
- `inverter_control_paired_summary_rows` reports candidate, definition-matched,
  publishable, and dropped pair counts plus exposure-weighted mean and
  p05/p50/p95 deltas. Negative, positive, and numerically indistinguishable
  fractions use a declared magnitude-relative tolerance with a one-unit floor
  in each metric's reported unit. Per-metric finite-pair counts expose
  unavailable quantities. `dropped_pair_known_baseline_count` states how many
  dropped pairs have an inspectable baseline; comparing those baselines with
  the retained ones is what turns an attrition count into a bound on selection
  bias. Unmatched rows retain the baseline metadata, so they are attributed to
  the baseline cohort when grouped.

Summary energy is the weighted sum of power times `duration_h` and is reported
in kWh. Current, VUF, capacitor-current, and ripple-voltage p50/p95/p99 values
are duration-and-weight inverse empirical CDFs over finite, publishable
inverter-device points, not unweighted quantiles or quantiles of scenario
maxima. Every metric has its own `*_finite_points` count; optional ratings can
therefore never produce a confident-looking utilization tail from an
undisclosed subset. The summary reports both ordinary and duration-times-weight
exposure-adjusted publication fractions. Do not
interpret metrics from different publication fractions as a fair control-law
comparison; examine retained failure rows and paired results first.

A publication fraction on its own does not say why cases failed, and the two
reasons are scientifically opposite. `termination_status_counts` gives the
per-cohort histogram of case termination statuses (with `"ERROR"` for thrown
cases), and `iteration_limit_case_count`,
`locally_infeasible_case_count`, and `numerical_error_case_count` are named
conveniences drawn from it. A `LOCALLY_INFEASIBLE` case is a candidate physical
statement about the controller and hardware under test. An `ITERATION_LIMIT` or
`NUMERICAL_ERROR` case is a statement about the smooth NLP's conditioning, and
it is not portable: the same source, Julia version, and pinned dependencies can
converge on one platform's linear-solver build and not on another's. Report the
histogram beside any failure fraction, and never attribute unconverged cases to
a control law without first re-solving them from a different start.

Fleet tails include both absolute current and rating-normalized utilization.
The total converter/grid current companions combine declared carrier reserves
with fundamental current in RMS quadrature; capacitor utilization uses the
thermally weighted equivalent current. This prevents a heterogeneous fleet's
largest nameplate from being mistaken for its most highly stressed device.

Case construction defensively copies metadata once. Extracted case, device, and
phase rows then share that copied dictionary to avoid one metadata allocation
per inverter or phase. Treat row metadata as immutable during analysis.

The library accepts already-prepared cases. Dataset ingestion, scenario
generation, random customer placement, persistence formats, plots, and cluster
scheduling remain in `scripts/studies/inverter_controls/` so the core API does
not encode one dataset's schema.

## Outer hardware grids and sizing diagnostics

Hardware is varied outside the nonlinear network problem. First call
`validate_inverter_control_campaign(cases)` to enforce a complete controller
matrix and identical snapshot, hardware, request, duration, weight, and
metadata across variants. By default it requires variants to share the same
network object; this is a cheap drift guard, while a dataset adapter may opt out
only after independently proving network equivalence.

This keeps every
solve interpretable as one realizable nameplate design and avoids treating a
current or capacitance rating as a freely optimized operating variable.
`InverterHardwareSweepPoint` scales installed ratings relative to each
inverter's base data, so heterogeneous devices retain their relative sizes:

```julia
hardware = [
    InverterHardwareSweepPoint(
        id="reference"),
    InverterHardwareSweepPoint(
        id="more_silicon_and_capacitance",
        converter_current_scale=1.20,
        grid_current_scale=1.20,
        dc_capacitance_scale=1.50,
        capacitor_current_scale=1.20),
]

sweep_cases = expand_inverter_hardware_cases(cases, hardware)
sweep = run_inverter_control_study(sweep_cases)

requirements = inverter_control_hardware_requirement_rows(
    sweep; allowed_dc_ripple_fraction=0.02, current_margin=1.10)
summaries = inverter_control_study_summary_rows(
    sweep; group_by=["penetration", "hardware_point_id"])
```

The expansion appends the hardware identifier to `scenario_id` but leaves
`variant_id` as the control-law identifier. Paired rows therefore compare laws
only at the same feeder snapshot and hardware point. The expansion records the
base scenario, hardware identifier, and every effective scale in metadata
(`1.0` for an unchanged axis) and rejects reserved-key collisions. The
sweep-point struct still uses `nothing` to mean "do not require or alter this
optional rating." The network is still shared by reference; each fleet and
advanced-inverter specification is new.

Converter-leg current, grid-side current, DC capacitance, and capacitor thermal
current are distinct axes. A rating is changed only when its scale is explicitly
supplied, and the all-`nothing` default is a true identity point. Asking to scale
an absent `i_max`, `i_grid_max`, `c_dc`, or `i_cap_max`
throws instead of silently introducing a physical constraint. Apparent-power
rating and DC voltage remain fixed in this sweep; study them as separately
declared hardware factors if required.

`inverter_control_hardware_requirement_rows` reports achieved converter/grid
RMS current, capacitor thermally equivalent current, and the closed-form
monolithic-link requirement

```math
I_{2\omega,rms} = \frac{|\widetilde S|}{\sqrt{2}V_{dc}},\qquad
C_{2\omega,req} =
\frac{|\widetilde S|}{2\omega V_{dc}\Delta V_{allow}},\qquad
E_{2\omega,req}=\frac12 C_{2\omega,req}V_{dc}^2,\qquad
\Delta V_{allow}=\epsilon_{dc}V_{dc}.
```

The explicit 2ω current separates the low-frequency ripple mechanism from the
total thermally equivalent capacitor current, which may also contain carrier
and neutral-current contributions. Required and installed stored energy make
capacitance alternatives comparable at their declared DC voltage. The rows
also contain absolute ratings, explicitly named `*_requirement_utilization`
ratios, and nullable compliance flags. These ratios include the declared
current margin; the similarly named fleet operating utilizations do not. A
current margin of at least one may be applied explicitly.
Failed/non-publishable solves retain the known installed hardware but make no
performance claim: requirements and utilizations are NaN, and compliance is
`missing`.

Every scaled absolute rating is checked after multiplication. A scale that
overflows to infinity or underflows to zero is rejected before model
construction rather than becoming an invalid nameplate silently.

The allowed ripple is a zero-to-peak sinusoidal amplitude: a fraction of 0.02
means ±2% and therefore 4% peak-to-peak. Rows expose both the achieved `dv2`
and the model's enforced `dv2_max`, plus binding flags for converter current,
grid current, capacitor thermal current, DC ripple, and converter apparent
power. Controller-command curtailment and power/current allocation activity are
reported separately. A `missing` binding flag means that no corresponding
installed limit was declared or the solve was not publishable. A `false` flag
only rules out that named limiter: it is not evidence that the operating point
is unconstrained. In particular, the table does not exhaustively classify
sampled modulation rails, internal-voltage bounds, sequence/neutral limits, or
every optional plant constraint.

The independent `dc_capacitance_scale` changes ripple headroom only. It does not
change `esr_dc`, capacitor-current rating, cost, tolerance, or lifetime, so it
must not be interpreted as a realizable parallel-bank or capacitor-loss axis.
A physical bank study must construct coupled plant alternatives—for an ideal
parallel multiplication by `n`, typically `C → nC`, current rating `→ nI`, and
`ESR → ESR/n`—and state its component assumptions explicitly.

There are two scientifically different uses:

1. Run every law at one common, demonstrably non-binding oversized hardware
   point. The achieved-current and closed-form capacitance columns are useful
   first-pass requirements.
2. When a rating changes dispatch, limiter activation, switching feasibility,
   or network voltage, run the explicit hardware grid and judge each point
   using publication status plus the declared voltage, unbalance, curtailment,
   loss, and ripple service criteria.

Do not infer a minimum rating from the current achieved at a binding rating;
that is circular because the controller may curtail to satisfy the installed
limit. Do not apply bisection unless the chosen service predicate has been
shown monotone for that controller and scenario. Retain the full grid when
active regimes change.

Hardware points are counterfactual alternatives, not additional time samples.
Study summaries automatically add `"hardware_point_id"` to their grouping when
all cases carry that metadata; summing energy across hardware points would
multiply the represented exposure. The 2ω formula assumes the
modeled monolithic three-leg link supplies the ripple. It does not size split
half-banks, hold-up energy, transients, tolerance/ageing, thermal lifetime, or a
frequency-dependent DC source, and it does not replace averaged-model and EMT
validation.

## API

```@docs
ControlledInverterFleetSpec
ControlledInverterFleetResult
solve_controlled_inverter_fleet
controlled_inverter_rows
controlled_inverter_phase_rows
InverterControlStudyCase
InverterControlStudyCaseResult
InverterControlStudyResult
run_inverter_control_study
inverter_control_study_case_rows
inverter_control_study_device_rows
inverter_control_study_phase_rows
inverter_control_study_summary_rows
inverter_control_paired_rows
inverter_control_paired_summary_rows
InverterHardwareSweepPoint
resize_controlled_inverter_fleet
validate_inverter_control_campaign
expand_inverter_hardware_cases
inverter_control_hardware_requirement_rows
```
