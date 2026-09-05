# Replay an issued DOE faithfully

<!-- doe-executable -->

All Julia blocks run in order from the repository root and are executed by CI.
The experiment separates scenario information, participant utilization and
numerical uncertainty. It uses a small synthetic feeder; the control policy is
a research assumption, not an inferred property of a deployed controller.

## Issue scenario-dependent controls

One STATCOM setting may depend on which network scenario occurs, but must be
shared by every represented utilization point in that scenario. The participant
nameplates are deliberately small so identity replay is tested away from a
binding capacity optimum.

```julia
using PowerOptLab
using BMOPFTools: add_statcom!
include("scripts/cases/doe_uncertainty_tutorial.jl")

net = doe_uncertainty_feeder(500, 600)
add_statcom!(net, "bus2"; s_max=2000.0)
stat = net["ibr"]["statcom_bus2"]
stat["p_min"] = [0.0]; stat["p_max"] = [0.0]
stat["q_min"] = [-2000.0]; stat["q_max"] = [2000.0]
other = deepcopy(net)
other["load"]["d1"]["p_nom"] = [800.0]
cps = [ConnectionPoint(id="d1", bus="bus1", export_max=2000.0),
       ConnectionPoint(id="d2", bus="bus2", export_max=2000.0)]
a = DOEScenario(id="forecast-a", network=net, role=:test)
b = DOEScenario(id="forecast-b", network=other, role=:test)
scenarios = DOEScenarioSet([a, b]; dataset_id="two-forecasts")
policy = PerfectRecourse(rules=[DOEControlRule(component=:ibr,
    id="statcom_bus2", quantity=:reactive_power, stage=:scenario)])
issued = solve_operating_envelope(scenarios, cps;
    security=:corners, control_policy=policy)
@assert length(issued.diagnostics[1]["issued_control_values"]) == 2
```

## Reorder or select recorded information states

The result supplies direction, policy, per-unit base and smoothing defaults to
verification. Scenario settings follow IDs **and** network hashes, not vector
indices. Connection bindings must match the original physical devices.

```julia
reordered = DOEScenarioSet([b, a]; dataset_id="reordered")
replay = verify_operating_envelope(reordered, cps, issued; utilizations=:corners)
subset = verify_operating_envelope(
    DOEScenarioSet([b]; dataset_id="subset"), cps, issued; utilizations=:corners)
@assert all(replay.feasible) && all(subset.feasible)
@assert replay.diagnostics[1]["control_policy_source"] == :issued_result
@assert replay.diagnostics[1]["issued_control_replay_count"] == 2
original_q = issued.diagnostics[1]["issued_control_values"][2]["value"]
replayed_q = subset.diagnostics[1]["issued_control_values"][1]["value"]
@assert isapprox(original_q, replayed_q; atol=1e-3)
```

For a new held-out information state, the recorded scenario setting is not
available. Held-out coverage deliberately retains issue-time settings and lets
scenario settings adapt. That is a different declared information experiment.

```julia
new_state = DOEScenarioSet([
    DOEScenario(id="new-state", network=other, role=:test)]; dataset_id="new")
@assert try
    verify_operating_envelope(new_state, cps, issued)
    false
catch error_
    error_ isa ArgumentError
end
held_out = evaluate_operating_envelope_coverage(new_state, cps, issued;
    utilizations=:corners)
@assert held_out.verification.diagnostics[1]["replay_control_stages"] == [:issue]
@assert ismissing(held_out.metrics["one_sided_hoeffding_upper_bound"])
```

Using a capacity dictionary discards issuance metadata. Explicitly changing the
policy or disabling issued-control replay permits control re-optimization and
is labelled accordingly. These are useful comparisons, but they must not be
described as replaying the original policy.

## Keep numerical non-results visible

An exhausted iteration budget cannot identify a physically violating scenario.
Its conservative coverage treatment remains adverse, while its classification
remains unresolved.

```julia
limited = evaluate_operating_envelope_coverage(new_state, cps, issued;
    solver_options=(max_iter=0,))
@assert limited.outcome == :inconclusive
@assert limited.metrics["candidate_scenario_count"] == 0
@assert limited.metrics["unresolved_scenario_count"] == 1
@assert limited.metrics["conservative_scenario_frequency"] == 1
```

An individual successful solve also cannot establish a common shared-control
policy. Inspect `verification_verdict` and interval evidence alongside the raw
`context.feasible`. A failed or incomplete requested fixed-control replay is
unresolved. For physical violation magnitudes, use a valid power-flow witness
or an [analytic reference](doe_analytic_reference.md), not an unpublished iterate.

## Connect DSSE through an explicit adapter

The committed DSSE example materializes estimated injections into the
operational network, then compares the DOE point with a separate fixed-dispatch
power flow. The synthetic measurement state has zero DER injection, so each
estimated net injection identifies one load. This assumption would need an
explicit disaggregation model on a real feeder.

```julia
include("scripts/validate_doe_from_dsse.jl")
include("scripts/cases/doe_dsse_validation_demo.jl")
validation = run_doe_validation(doe_validation_case())
@assert validation.max_dsse_voltage_error_V < 1e-3
@assert validation.statcom.max_doe_pf_voltage_difference_V < 1e-3
```

The separate PF adapter freezes supported free IBR P/Q, respects export/import
sign, and retains constant power-factor laws. The pinned PF engine strips
ratings used to normalize droop, so the adapter rejects Volt-VAr/Volt-Watt
profiles; the main DOE verifier preserves those laws. It also rejects unsupported
free controls, custom context hooks and custom smoothing.
This is a separate build using the same network engine, not independent software
validation or a proof of continuous envelope containment.
