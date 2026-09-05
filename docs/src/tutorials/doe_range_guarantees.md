# What does a dynamic operating envelope guarantee?

<!-- doe-executable -->

All Julia blocks on this page run in order from the repository root. CI executes
the exact blocks with `julia --project=. scripts/run_doe_tutorials.jl`.

> **Claim:** a feasible simultaneous upper point need not imply that partial
> utilization of the advertised box is feasible.
>
> **Model:** synthetic, unbalanced nonlinear AC, fundamental-frequency RMS
>
> **Expected runtime:** under one minute after Julia compilation
>
> **Data license:** no external data; the complete synthetic case is below
>
> **Random seed:** none; the allocation, corners, and refinement schedule are deterministic
>
> **Solver evidence:** local nonlinear Ipopt solutions and same-formulation AC replay

This example isolates a range-safety failure. Three equal single-phase export
limits are connected to a balanced three-phase feeder with a tight
negative-sequence-voltage limit. Simultaneous equal export remains balanced,
but unequal customer utilization creates voltage unbalance.

The motivating field result is the published counterexample of [Liu and
Braslavsky (2022)](https://doi.org/10.1109/ACCESS.2022.3203062), in which a
customer moving strictly inside an equal 2.91 kW envelope produces a larger
cross-phase voltage violation. This tutorial demonstrates the same logical
pitfall on a self-contained case. It is **not** a reproduction of the CSIRO
feeder, its parameters, or its numerical result.

## Build the complete case

```julia
using PowerOptLab
include("scripts/cases/doe_range_benchmark.jl")
include("scripts/cases/doe_analytic_reference.jl")
using .DOEAnalyticReference
case = doe_benchmark_case()
net, cps = case.nets, case.connection_points
```

The case is intentionally small enough to inspect. It is a methodological
fixture, not a calibrated distribution feeder.

## Allocate only at the simultaneous bound

```julia
bound = solve_operating_envelope(
    net, cps;
    security=:bound_point,
    control_policy=PerfectRecourse())

bound.total_capacity
bound.diagnostics[1]["security_scope"]

bound_check = verify_operating_envelope(
    net, cps, bound;
    utilizations=:bound_point,
    control_policy=PerfectRecourse())

@assert bound_check.feasible == [true]
```

At ``u=(1,1,1)``, equal phase exports preserve symmetry. The returned result
supports only the claim “locally feasible at the simultaneous bound point.”
It does not establish box containment.

## Ask the range question

For three participants, exact corner enumeration is still inexpensive:

```julia
corner_check = verify_operating_envelope(
    net, cps, bound;
    utilizations=:corners,
    control_policy=PerfectRecourse())

corner_check.feasible
corner_check.context_results[1]
```

Several asymmetric corners violate the negative-sequence constraint.

## Trace the violation against the number of participants who back off

Corner enumeration answers “is some corner unsafe?” but not “how much
asymmetry does it take?”. The operationally meaningful question is what
happens when participants simply do not use what they were issued — a customer
whose PV is clouded, whose battery is charging, or whose inverter is offline.
[`doe_dropout_utilizations`](@ref) generates exactly that family: every way of
taking `1..depth` participants to zero while the remainder stay at the issued
limit.

```julia
dropouts = doe_dropout_utilizations(length(cps), length(cps))

dropout_check = verify_operating_envelope(
    net, cps, bound;
    utilizations=dropouts,
    control_policy=PerfectRecourse())

for (point, context) in zip(dropouts, dropout_check.context_results[1])
    dropped = count(iszero, point)
    println(dropped, " dropped → feasible=", context.feasible,
            "  worst normalized margin=",
            get(context.diagnostics, "minimum_normalized_margins", missing))
    p = point .* [bound.envelope[cp.id][1] for cp in cps]
    println("  analytic negative sequence [V]=", negative_sequence(phase_voltages(p)))
end
```

The number of participants that back off is the natural x-axis for this
network: a single dropout already breaks the phase symmetry that made the
bound point safe, while dropping every participant returns to the trivially
feasible zero point. A failed NLP has no published physical margin. The independent analytic
reference supplies the violation magnitude: 10 V negative sequence against a
1 V limit for the nontrivial corners. See the [analytic reference
tutorial](doe_analytic_reference.md) for the assumptions and derivation.

This is `length(cps)` points at depth 1 and
`sum(binomial(n, k) for k in 1:depth)` overall. For fixed small depth this is
polynomial in `n`; taking `depth=n`, as in this three-participant illustration,
is exponential.

## Search without enumerating

The same failure can be found by enumerating only the single-dropout seeds.
If those pass, the algorithm can refine around the most stressed point:

```julia
search = search_operating_envelope_adversarial(
    net, cps, bound;
    seed_samples=0,
    refinement_rounds=1,
    restarts=2,
    initial_step=0.5,
    control_policy=PerfectRecourse())

@assert search.outcome == :candidate_counterexample

search.diagnostics["candidate_locations"]
search.worst_interval
search.worst_context.utilization
```

The seeded dropout faces ``(0,1,1)``, ``(1,0,1)`` and ``(1,1,0)`` are evaluated
in the first round. The search stops immediately when it finds a candidate;
no interior refinement occurs in this run. A separate run with
`dropout_depth=0` isolates interior refinement. The search retains every verification round and every
normalized-headroom score.

!!! note "Read refinement reach against the finite budget"
    Starting from the bound, an optimistic minimum reachable coordinate after
    `r` decreasing moves is

    ```
    max(0, 1 - sum(max(initial_step * step_decay^k, minimum_step) for k in 0:r-1))
    ```

    `refinement_reachable_depth_from_bound` reports this budget bound, not an
    assurance about which paths the search takes. A positive `minimum_step`
    invalidates an infinite geometric-series argument. Increasing the step,
    decay factor, or number of rounds can extend reach; seeded dropout faces
    test complete backing-off immediately.

    The Halton seeds have a second budget limit: for a dimension whose prime
    base exceeds `seed_samples`, the radical inverse is just `index / base`, so
    those coordinates collapse into a narrow band near zero and correlate with
    each other. Beyond a handful of participants, pass
    `halton_scramble_seed` and check
    `diagnostics["halton_seed_stratified"]`.

Repeat the candidate from several deterministic OPF starts before presenting it
as numerical evidence:

```julia
confirmation = confirm_operating_envelope_counterexample(
    net, cps, bound, search.worst_context.utilization;
    start_scales=(1.0, 0.95, 1.05),
    control_policy=PerfectRecourse())

confirmation.outcome            # :repeated_candidate in this fixture
confirmation.diagnostics["candidate_locations_by_run"]
```

Finally, add the discovered utilization points to the allocation and repeat the
search:

```julia
adaptive = solve_adversarial_search_stable_operating_envelope(
    net, cps;
    max_rounds=2,
    control_policy=PerfectRecourse(),
    search_keywords=(
        seed_samples=0,
        refinement_rounds=1,
        restarts=2,
        initial_step=0.5,
    ))

adaptive.outcome
[allocation.total_capacity for allocation in adaptive.allocations]
adaptive.diagnostics["final_allocation_screened"]
```

The second allocation is lower because it represents the asymmetric points
that falsified the original bound-point allocation. `:search_stable` still
means only that no further counterexample was found within this recorded finite
budget.

## Run the committed benchmark

The same synthetic network is packaged as a versioned case with explicit
method labels, seed provenance, claim text, and CC0 licensing metadata:

```sh
julia --project=. scripts/run_doe_benchmark.jl \
    scripts/cases/doe_range_benchmark.jl \
    results/doe-range.tsv
```

The output contains one row for a simultaneous bound-point calculation with
explicit ideal recourse and one for finite corner security with issue-time
controls plus local laws. Each row has its own study identity and records
`global_certificate=false`. The smoke test checks the row schema and these
claim fields without fixing solver-dependent timings or capacities.

## Interpret the result

This experiment distinguishes four statements:

1. The simultaneous upper point is locally feasible.
2. The advertised independent box is falsified by a partial-utilization point.
3. Exact corner verification also finds unsafe points in this particular case.
4. Neither the local infeasibility status nor the heuristic search is a global
   mathematical certificate for a general nonlinear AC model.

For a publication-quality study, rerun the candidate from multiple starts,
independently recompute the binding network quantity, preserve the complete
study manifest, and repeat with operationally defensible control recourse.

## What this does not prove

- It does not show that all worst cases occur at corners.
- It does not assign a violation probability or uncertainty coverage level.
- It does not certify global AC infeasibility at the candidate point.
- It does not validate the synthetic impedance or limit values against field
  data.
- `PerfectRecourse()` is used because this case has no free network controller;
  studies with taps, STATCOMs, or other controls must state their information
  structure explicitly.
