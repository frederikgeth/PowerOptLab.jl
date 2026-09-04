# DOE uncertainty provenance and held-out coverage

> **Claim:** an envelope calibrated on declared scenarios can be evaluated on a
> disjoint typed test set without presenting in-sample feasibility as
> probabilistic coverage.
>
> **Model:** synthetic nonlinear AC feeder; finite scenario ensemble
>
> **Expected runtime:** under one minute after Julia compilation
>
> **Data license:** no external data; the complete case builder is below
>
> **Seeds:** recorded per scenario; the fixture itself is deterministic

This tutorial separates scenario construction, allocation, and held-out
evaluation. The example deliberately calibrates under higher demand, which
absorbs local export, then tests the resulting envelope under a more restrictive
low-demand realization.

## Build reproducible network realizations

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf
using Dates

function coverage_feeder(p1, p2)
    parse_bmopf("""
    {"bus":{
        "src":{"terminal_names":["1","n"],
               "perfectly_grounded_terminals":["n"]},
        "bus1":{"terminal_names":["1","n"],
                "perfectly_grounded_terminals":["n"],
                "v_min":[216.0],"v_max":[245.0]},
        "bus2":{"terminal_names":["1","n"],
                "perfectly_grounded_terminals":["n"],
                "v_min":[216.0],"v_max":[245.0]}},
     "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],
         "v_magnitude":[230.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.4}},
     "line":{
        "l1":{"bus_from":"src","bus_to":"bus1",
              "terminal_map_from":["1"],"terminal_map_to":["1"],
              "linecode":"lc","length":1.0},
        "l2":{"bus_from":"bus1","bus_to":"bus2",
              "terminal_map_from":["1"],"terminal_map_to":["1"],
              "linecode":"lc","length":1.0}},
     "load":{
        "d1":{"bus":"bus1","terminal_map":["1","n"],
              "configuration":"SINGLE_PHASE","p_nom":[$p1],"q_nom":[0.0]},
        "d2":{"bus":"bus2","terminal_map":["1","n"],
              "configuration":"SINGLE_PHASE","p_nom":[$p2],"q_nom":[0.0]}}}
    """; from_string=true)
end

cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
       ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
```

## Declare roles and provenance

```julia
calibration = DOEScenario(
    id="calibration-high-load",
    network=coverage_feeder(5000.0, 5000.0),
    role=:calibration,
    weight=1.0,
    source="synthetic tutorial fixture",
    generation_method=:deterministic_fixture,
    seed=11,
    timestamp=DateTime(2026, 1, 1, 12),
)

test_low = DOEScenario(
    id="test-low-load",
    network=coverage_feeder(200.0, 200.0),
    role=:stress,
    weight=0.7,
    source="synthetic tutorial fixture",
    generation_method=:held_out_fixture,
    seed=12,
    timestamp=DateTime(2026, 1, 15, 12),
)

test_high = DOEScenario(
    id="test-high-load",
    network=coverage_feeder(5000.0, 5000.0),
    role=:test,
    weight=0.3,
    source="synthetic tutorial fixture",
    generation_method=:held_out_fixture,
    seed=13,
    timestamp=DateTime(2026, 1, 15, 12),
)

scenarios = DOEScenarioSet(
    [calibration, test_low, test_high];
    dataset_id="doe-coverage-tutorial-v1",
    metadata=Dict("license" => "synthetic", "split_method" => "declared"),
)

calibration_set = select_doe_scenarios(scenarios; roles=:calibration)
test_set = select_doe_scenarios(scenarios; roles=(:test, :stress))
```

IDs must be unique within an interval. Weights are optional, but if one scenario
in an interval is weighted, all must be weighted. They are relative weights;
the coverage calculation normalizes them within each interval.

For historical samples, prefer a chronological split over manual role labels:

```julia
time_split = split_doe_scenarios_by_time(
    scenarios;
    calibration_end=DateTime(2026, 1, 2),
    test_start=DateTime(2026, 1, 15),
    split_name="blocked-holdout")

time_split.calibration
time_split.test
time_split.excluded_scenario_ids
time_split.diagnostics
```

Samples before `calibration_end` become calibration cases, samples on or after
`test_start` become test cases, and samples in between are explicitly excluded.
The helper retains each original role in scenario metadata. It guarantees no
timestamp overlap, but it does not assess leakage through feeders, sites,
customers, weather events, or scenario-generator fitting. The current helper
accepts one forecast interval at a time so it cannot silently confuse forecast
intervals with statistical sample time.

## Allocate without looking at the test scenarios

```julia
allocation = solve_operating_envelope(
    calibration_set, cps;
    security=:bound_point,
    control_policy=PerfectRecourse())

spec = DOEStudySpec(
    calibration_set, cps;
    security=:bound_point,
    control_policy=PerfectRecourse(),
    metadata=Dict("experiment" => "held-out tutorial"))

manifest = doe_study_manifest(spec)
manifest["scenario_provenance"]
manifest["scenario_set_metadata"]
```

The study identity covers scenario roles, weights, sources, construction
methods, seeds, timestamps, metadata, and network hashes. Changing a split or
provenance field therefore changes the study ID even when the network numbers
are identical.

## Evaluate held-out outcomes

```julia
coverage = evaluate_operating_envelope_coverage(
    scenarios, cps, allocation;
    roles=(:test, :stress),
    utilizations=:bound_point,
    control_policy=PerfectRecourse(),
    iid_assumption=false)

coverage.outcome
coverage.scenario_rows
coverage.metrics["empirical_candidate_scenario_frequency"]
coverage.metrics["weighted_candidate_scenario_frequency"]
coverage.metrics["one_sided_hoeffding_upper_bound"]  # missing
```

The low-demand stress scenario is more restrictive and produces a candidate
violation; the high-demand test realization passes. The result reports context
and scenario frequencies separately. A scenario is a candidate violation if
any tested utilization point fails. Unresolved cases are excluded from the
ordinary empirical rate and counted as violations in the conservative rate.

Issue-time controls recorded in `allocation` are fixed during held-out replay.
Scenario-stage controls are allowed to adapt to each newly observed scenario,
consistent with their declared information stage. Passing a capacity dictionary
instead loses the issuance record and is labelled `:capacity_values_only`.

## Trace capacity against held-out violations

A single operating point hides how quickly empirical performance deteriorates
as offered capacity rises. Evaluate a declared scale grid:

```julia
curve = evaluate_operating_envelope_coverage_curve(
    scenarios, cps, allocation;
    scales=(0.25, 0.5, 0.75, 1.0),
    roles=(:test, :stress),
    utilizations=:bound_point,
    control_policy=PerfectRecourse(),
    iid_assumption=false)

curve.rows
curve.diagnostics["first_candidate_scale"]
curve.diagnostics["candidate_count_reversals"]
```

Scales are sorted and deduplicated. Issue-time controls from `allocation` are
retained rather than re-optimized against the held-out set. The first candidate
scale is only a threshold on the tested grid. It is not a continuous safe
capacity limit. A decrease in candidate count at a larger scale is preserved:
it can arise from genuine physical non-monotonicity, discrete or control
effects, or irregular local nonlinear solves, and should be investigated rather
than silently smoothed.

## Compare reference and shifted conditions

Here the high-demand test case acts as a reference and the deliberately
restrictive low-demand stress case as the shifted evaluation set:

```julia
reference_coverage = evaluate_operating_envelope_coverage(
    scenarios, cps, allocation;
    roles=:test,
    control_policy=PerfectRecourse())

stress_coverage = evaluate_operating_envelope_coverage(
    scenarios, cps, allocation;
    roles=:stress,
    control_policy=PerfectRecourse())

shift_comparison = compare_doe_coverage_shift(
    reference_coverage, stress_coverage;
    reference_label="held-out reference",
    shifted_label="low-demand stress")

shift_comparison.outcome
shift_comparison.metric_deltas
shift_comparison.diagnostics["distribution_shift_detected"]  # false
```

This comparison reports changed held-out performance. It intentionally does not
claim that a distribution shift has been detected: that requires a separately
chosen two-sample, covariate, residual, or concept-shift test with assumptions
matched to the data-generating process. A named stress set is also not evidence
of prevalence in the operating population.

## Add a statistical bound only when justified

If the selected held-out scenarios are genuinely independent draws from the
target distribution, make that assumption explicit:

```julia
iid_coverage = evaluate_operating_envelope_coverage(
    scenarios, cps, allocation;
    roles=(:test, :stress),
    control_policy=PerfectRecourse(),
    iid_assumption=true,
    confidence=0.95)

iid_coverage.metrics["one_sided_hoeffding_upper_bound"]
```

The one-sided Hoeffding bound uses the unweighted scenario count and treats
unresolved scenarios as candidate violations. Declared scenario weights affect
the weighted empirical estimate, not this concentration bound. With only two
held-out scenarios the bound is intentionally weak.

Do not set `iid_assumption=true` for hand-picked stress cases, overlapping time
windows, scenarios reused during method development, or a split affected by
temporal leakage. In those cases, publish the empirical and weighted results
without a confidence claim.

## Connecting a DSSE workflow

For measurement-derived studies, construct one `DOEScenario` per materialized
network realization and record at least:

- the DSSE dataset/version and time index in `source` and `timestamp`;
- covariance/profile sampling method and code revision in
  `generation_method` and `metadata`;
- the exact random seed;
- a time-block or site-level split that prevents leakage; and
- whether impedance candidates, topology alternatives, and forecast errors are
  sampled jointly or independently.

PowerOptLab does not yet turn a DSSE covariance into network scenarios
automatically. The typed layer preserves the experiment once those networks are
constructed; it does not validate the scenario generator.

## What this does not prove

- Empirical test performance is not a robust box certificate.
- An i.i.d. concentration bound does not protect against distribution shift.
- Scenario weights are not automatically calibrated probabilities.
- Local nonlinear infeasibility remains candidate evidence rather than a global
  infeasibility proof.
- Testing only `:bound_point` says nothing about untested partial utilization;
  combine held-out scenarios with corners or adaptive utilization search when
  the operational promise is a range.
