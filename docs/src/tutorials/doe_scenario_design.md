# DOE scenario design and held-out evaluation

> **Question:** how can scenarios support a reproducible empirical statement
> without being mistaken for a probability model or robust set?
>
> **Prerequisite:** complete [Dynamic operating envelopes in 15
> minutes](doe_getting_started.md)
>
> **Model:** synthetic nonlinear AC feeder · **Runtime:** under one minute after
> compilation · **Data:** CC0 synthetic · **Seeds:** 11–13 in the fixture

An uncertainty study has three distinct objects:

1. an uncertainty description, such as residual data or an impedance candidate
   set;
2. materialized network scenarios; and
3. an issued envelope evaluated against selected scenarios.

Keeping these objects separate prevents test cases from leaking into allocation
and prevents scenario counts from being presented as confidence levels.

## 1. Load one declared scenario set

```julia
using PowerOptLab

include("scripts/cases/doe_uncertainty_tutorial.jl")
case = doe_uncertainty_tutorial_case()

scenarios = case.scenarios
cps = case.connection_points
```

The case contains one calibration scenario and two held-out scenarios. The
low-demand stress scenario is more restrictive for export because less local
demand absorbs the participant injections.

```julia
[(scenario.id, scenario.role, scenario.weight)
 for scenario in scenarios.intervals[1]]
```

The role answers *how the scenario may be used*. The weight records a declared
relative importance. Neither field proves that a scenario is an independent
draw or that its weight is a calibrated probability.

## 2. Allocate without seeing the test cases

```julia
calibration = select_doe_scenarios(scenarios; roles=:calibration)
held_out = select_doe_scenarios(scenarios; roles=(:test, :stress))

allocation = solve_operating_envelope(
    calibration, cps;
    security=:bound_point,
    control_policy=PerfectRecourse())
```

This deliberately simple allocation uses only the high-demand calibration
case. Its capacity is not yet evidence about the held-out conditions.

For time-series or multi-site data, use [`split_doe_scenarios_by_time`](@ref)
and [`audit_doe_scenario_calibration`](@ref) to inspect timestamp, exact-network,
site/group, and provenance overlap. A role label alone cannot prevent leakage.

## 3. Evaluate the fixed issued envelope

```julia
coverage = evaluate_operating_envelope_coverage(
    held_out, cps, allocation;
    roles=(:test, :stress),
    utilizations=:bound_point,
    control_policy=PerfectRecourse())

coverage.outcome
coverage.scenario_results
coverage.metrics
```

The capacity and any issued controls are held fixed. The evaluation reports
which finite scenarios pass, fail, or remain unresolved, together with empirical
and weighted summaries.

The permitted statement is:

> The issued allocation was locally evaluated on the two declared held-out
> network scenarios, with the reported finite outcomes.

It is not “70% reliable,” even if the failing scenario has weight 0.7. A
probability statement needs a sampling model, independence unit, calibration
argument, and appropriate statistical estimator.

## 4. Combine exogenous scenarios with participant utilization

Network uncertainty and customer utilization are different axes:

```julia
corner_coverage = evaluate_operating_envelope_coverage(
    held_out, cps, allocation;
    roles=(:test, :stress),
    utilizations=:corners,
    control_policy=PerfectRecourse())
```

With two scenarios and two participants, this represents eight AC contexts.
Passing both bound points would not imply passing all corners; passing all
corners would still not certify every nonlinear interior point.

## 5. Preserve provenance

```julia
spec = DOEStudySpec(
    calibration, cps;
    security=:bound_point,
    control_policy=PerfectRecourse(),
    seeds=Dict("scenario_fixture" => 20260905),
    metadata=case.metadata)

manifest = doe_study_manifest(spec)
manifest["study_id"]
manifest["scenario_provenance"]
```

A useful scenario record identifies its physical quantity and unit, data source
or model, construction method, time, random seed, role, materializer revision,
and relationship to other scenarios.

## Research exercise

Add a second calibration scenario with 200 W demand at each bus. Predict how
the allocated export changes, then explain why using that same scenario as a
test case would make the held-out result uninformative.

## What this does not prove

- The three scenarios are not a calibrated probability distribution.
- Scenario weights are not automatically event probabilities.
- Local nonlinear infeasibility is candidate evidence, not global
  infeasibility.
- Bound-point evaluation does not test independent participant use.
- Passing a held-out dataset does not protect against later distribution shift.

Continue with [Constructing and comparing DOE uncertainty
models](doe_uncertainty_models.md), then [Statistical validation of DOE
claims](doe_statistical_validation.md). The [advanced uncertainty
laboratory](doe_uncertainty_coverage.md) contains the complete API examples.
