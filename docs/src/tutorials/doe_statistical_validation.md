# Statistical validation of DOE claims

> **Question:** when may finite out-of-sample evaluations support a probability
> statement, and what must remain descriptive?
>
> **Prerequisite:** [DOE scenario design and held-out
> evaluation](doe_scenario_design.md)
>
> **Scope:** empirical coverage, concentration bounds, distribution-shift
> diagnostics, and probability calibration

Statistical validation does not turn a nonlinear local solve into a robust
certificate. It answers a different question: how often and how severely does
a fixed issued envelope fail under a declared data-generating and evaluation
procedure?

## 1. Name the estimand and observation unit

Before computing a rate, state:

- the population or operating period to which it refers;
- whether one observation is a timestamp, household, feeder, day, or event;
- whether observations are independent, grouped, or serially dependent;
- what counts as a violation, including tolerance and unresolved solves;
- whether controls remain fixed at their issued values; and
- whether participant utilization is observed, sampled, or adversarially
  selected.

“Two of 100 scenarios failed” is an empirical sample result. It is not a 2%
physical failure probability until the scenario sampling and observation unit
justify that interpretation.

## 2. Evaluate a fixed envelope first

```julia
using PowerOptLab

include("scripts/cases/doe_uncertainty_tutorial.jl")
case = doe_uncertainty_tutorial_case()

calibration = select_doe_scenarios(case.scenarios; roles=:calibration)
held_out = select_doe_scenarios(case.scenarios; roles=(:test, :stress))

issued = solve_operating_envelope(
    calibration, case.connection_points;
    security=:bound_point,
    control_policy=PerfectRecourse())

coverage = evaluate_operating_envelope_coverage(
    held_out, case.connection_points, issued;
    roles=(:test, :stress),
    utilizations=:corners,
    control_policy=PerfectRecourse())
```

Report the scenario count, context count, candidate violations, unresolved
cases, maximum violation magnitude, and both weighted and unweighted summaries.
Treat unresolved cases conservatively in a safety summary, while retaining them
as a separate numerical category for diagnosis.

## 3. Add a confidence bound only when justified

```julia
iid_coverage = evaluate_operating_envelope_coverage(
    held_out, case.connection_points, issued;
    roles=:test,
    utilizations=:corners,
    control_policy=PerfectRecourse(),
    iid_assumption=true,
    confidence=0.95)

iid_coverage.metrics["one_sided_hoeffding_upper_bound"]
```

The flag is a scientific assertion by the caller, not an independence test. Do
not use it for hand-picked stress cases, overlapping windows, repeated
measurements from one weather event, or scenarios consulted during method
development. Scenario weights do not increase the number of independent
observations.

For serial data, preserve ordering and use an evaluation design appropriate to
the dependence structure. A moving-block bootstrap can represent a declared
short-range dependence model, but choosing a block length after seeing the DOE
outcomes invalidates a simple confirmatory interpretation.

## 4. Diagnose distribution shift without overclaiming

[`test_doe_covariate_shift`](@ref) compares explicitly named numeric covariates.
It does not test every aspect of two unknown distributions. Its permutation
inference requires an exchangeability assertion at the scenario or declared
group level.

[`test_doe_time_series_covariate_shift`](@ref) uses circular shifts of ordered,
regularly spaced windows. Its inference requires circular-shift invariance,
which is stronger than ordinary stationarity and can fail under seasonality,
trends, interventions, or artificial wraparound.

Use shift diagnostics to qualify transportability:

> Held-out performance was measured on this period and these recorded
> covariates; material differences were or were not observed under the declared
> test.

Do not conclude that “no distribution shift exists” from a nonsignificant
finite diagnostic.

## 5. Evaluate probability forecasts as forecasts

If a separate method issues a predicted violation probability for each case,
evaluate those predictions with [`evaluate_doe_probability_calibration`](@ref):

- reliability tables;
- Brier and logarithmic scores;
- calibration error;
- resolution and uncertainty components where applicable; and
- explicit handling of unresolved network evaluations.

A calibrated probability model and a network-secure envelope are different
artifacts. Good average calibration can coexist with unacceptable rare-event
severity or poor subgroup performance. Always report violation magnitude and
stratified results alongside probability scores.

## 6. Use an evidence table in publications

| Claim | Minimum supporting record |
|---|---|
| Feasible on test set | fixed envelope and controls, scenario identities, outcomes, and solver evidence |
| Estimated violation frequency | sampling frame, observation unit, empirical numerator and denominator |
| Confidence bound | explicit independence/group assumption and stated estimator |
| Transported performance | reference/target periods, shift diagnostics, and limitations |
| Calibrated probability | independently issued predictions, proper scores, reliability, and unresolved handling |
| Certified robustness | a complete-set mathematical certificate; none of the rows above substitutes for it |

## Research exercises

1. Duplicate every held-out row ten times and recompute the naïve i.i.d. bound.
   Explain why the apparent sample size changes while the evidence does not.
2. Evaluate the same envelope by timestamp and by day. Compare the estimands and
   uncertainty intervals.
3. Construct a perfectly calibrated predictor with rare, very large voltage
   violations. Explain why calibration alone is inadequate for operation.

## What this does not prove

- A finite sample cannot establish complete nonlinear set containment.
- A caller assertion does not verify independence or exchangeability.
- Failure to detect covariate shift is not evidence of identical distributions.
- Calibration does not limit the physical magnitude of a violation.
- Reusing allocation data for evaluation produces optimistic evidence unless
  the dependence is explicitly accounted for.

The [advanced uncertainty laboratory](doe_uncertainty_coverage.md) contains the
complete coverage-curve, pairing, shift, temporal, calibration, covariance, and
bootstrap APIs in one reproducible reference workflow.
