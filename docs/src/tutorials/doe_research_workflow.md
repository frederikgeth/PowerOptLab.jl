# Reproducible DOE recourse, verification, and search

> **Audience:** researchers comparing dynamic-operating-envelope formulations
>
> **Claim level:** finite tested-point evidence from local nonlinear AC models
>
> **Model:** caller-supplied nonlinear unbalanced AC cases; examples use the
> synthetic cases introduced by the modelling and range-guarantee tutorials
>
> **Expected runtime:** case dependent; the small examples run in under one
> minute after Julia compilation
>
> **Data and seeds:** record both in `DOEStudySpec`; no external dataset is
> implied by the schematic fragments below

This tutorial turns the DOE framework into a reproducible experiment. It keeps
three questions separate:

1. Which participant-utilization and uncertainty points are represented?
2. Which controls know what information when they act?
3. What numerical and replay evidence supports the reported allocation?

None of the workflows below produces a global robust-feasibility certificate.
They are intentionally formulation-neutral: pointwise recourse remains a
first-class replication mode, while issue-time and local-law policies represent
different information structures rather than a replacement definition of DOE.

## 1. Name the information structure

Use pointwise recourse deliberately when reproducing a paper that re-optimizes
network controls at every represented operating point:

```julia
ideal = solve_operating_envelope(scenarios, cps;
    security=:corners,
    control_policy=PerfectRecourse())
```

For an operational comparison, issue manual/free settings before scenario and
participant utilization are known while retaining implemented automatic laws:

```julia
issued = solve_operating_envelope(scenarios, cps;
    security=:corners,
    control_policy=IssuePlusLocalLaws())

recourse_gain = ideal.total_capacity .- issued.total_capacity
```

This difference is the capacity attributable to additional information and
control flexibility. It is not automatically deployable capacity.

A mixed literature model can override one device family:

```julia
policy = PerfectRecourse(rules=[DOEControlRule(
    component=:ibr,
    id="statcom17",
    quantity=:reactive_power,
    stage=:scenario,
)])
```

Here the STATCOM obtains one setting after the scenario is known, but that
setting is shared by all utilization points in the scenario.

## 2. Inspect the control audit

Never infer the information structure from the function call alone. Preserve
the resolved audit returned with the result:

```julia
audit = issued.diagnostics[1]["control_audit"]

for control in audit
    println((control["component"], control["id"], control["quantity"]),
            " => ", control["stage"],
            ", links = ", control["link_constraints"])
end
```

Each record includes the native classification, canonical unit, automatic law,
policy source, contexts in which the control exists, equality groups, and
provenance metadata. `all_discovered_free_controls_classified=true` refers only
to the declared discovery scope; it is not evidence that an arbitrary extension
contains no hidden control variables.

## 3. Register a research extension

Extensions should register stable semantic handles rather than relying on JuMP
variable names. The model hook creates and registers the variable; the policy
registration describes its physical meaning:

```julia
using BMOPFTools: OpfModelKey, opf_model, register_opf_object!
using JuMP

support_key = OpfModelKey(:variable, :research_support, "setting")

function support_hook!(ctx)
    setting = @variable(opf_model(ctx), lower_bound=-1.0, upper_bound=1.0)
    register_opf_object!(ctx, support_key, setting)
end

support = DOEControlRegistration(
    component=:research_device,
    id="support1",
    quantity=:setting,
    handle=support_key,
    canonical_unit=:per_unit_setting,
    metadata=Dict("implementation" => "src/my_support_model.jl"),
)

policy = IssuePlusLocalLaws(registrations=[support])
r = solve_operating_envelope(net, cps;
    security=:corners,
    control_policy=policy,
    context_hook! = support_hook!)
```

Free generator P/Q remains deliberately unsupported for linked stages until the
upstream model exposes voltage-independent P/Q handles. Linking generator
currents would impose a different physical policy when voltage changes.

## 4. Verify every represented context

Verification now returns one `OperatingEnvelopeContextResult` per scenario and
utilization point:

```julia
check = verify_operating_envelope(scenarios, cps, issued;
    utilizations=:corners,
    control_policy=IssuePlusLocalLaws())

for context in check.context_results[1]
    println("scenario=", context.scenario,
            " utilization=", context.utilization,
            " feasible=", context.feasible)
end
```

When the joint policy model is feasible, each context is independently solved
again with its optimized free control values fixed. Inspect:

```julia
replay = check.context_results[1][1].diagnostics["independent_replay"]
replay["feasible"]
replay["control_replay_complete"]
replay["maximum_voltage_difference_V"]
```

If the joint model fails, the verifier solves every context separately. If all
individual contexts pass, the result is reported as
`shared_control_incompatibility_or_joint_nlp_failure`: the evidence does not
silently blame one utilization point when the conflict may arise from a common
setpoint or local solver behavior.

When `issued` is an `OperatingEnvelopeResult`, verification also reuses its
recorded issue/scenario control values by default:

```julia
check.diagnostics[1]["issued_control_replay_source"]
check.diagnostics[1]["issued_control_replay_count"]
```

This prevents the validation set from receiving a newly optimized manual
setting. A bare capacity dictionary has no issuance record and is explicitly
reported as `:capacity_values_only`.

## 5. Screen the interior of the advertised box

Corners do not prove that a non-convex AC feasible set contains the box
interior. A deterministic Halton screen supplies reproducible continuous
utilization points:

```julia
screen = search_operating_envelope_utilizations(
    scenarios, cps, issued;
    samples=64,
    sequence_offset=0,
    control_policy=IssuePlusLocalLaws())

screen.outcome                 # :search_stable, :candidate_counterexample, or :inconclusive
screen.utilization_points
screen.candidate_contexts
```

`search_stable` means local AC feasibility at the generated points only.
`candidate_counterexample` requires confirmation; `inconclusive` commonly
means a shared-control conflict could not be assigned to one point. Every
outcome records `global_certificate=false`.

Use margin-directed refinement when uniform coverage is less useful than
concentrating the same finite budget near stressed regions:

```julia
adversarial = search_operating_envelope_adversarial(
    scenarios, cps, issued;
    seed_samples=16,
    refinement_rounds=4,
    restarts=3,
    initial_step=0.25,
    control_policy=IssuePlusLocalLaws())

adversarial.outcome
adversarial.point_scores
adversarial.worst_context
adversarial.diagnostics["worst_score"]
```

The score is the negative minimum normalized headroom across declared voltage,
ampacity, and negative-sequence limits. Larger scores are more stressed; zero
denotes a binding limit. Each refinement round adds coordinate neighbours of
the highest-scoring points and jointly re-verifies the accumulated set, so
issue-time controls remain shared. This is a deterministic black-box
falsification heuristic—not a proof that the returned point globally maximizes
violation. Record all `verifications`, point scores, step settings, and the
candidate context.

Confirm a candidate across deterministic starts, then place the adaptive search
inside the allocation loop:

```julia
candidate = adversarial.worst_context
confirmation = confirm_operating_envelope_counterexample(
    scenarios, cps, issued, candidate.utilization;
    start_scales=(1.0, 0.9, 1.1),
    control_policy=IssuePlusLocalLaws())

adaptive_margin = solve_adversarial_search_stable_operating_envelope(
    scenarios, cps;
    max_rounds=4,
    control_policy=IssuePlusLocalLaws(),
    solve_keywords=(fairness=FairnessPolicy(
        kind=:max_min, normalization=:capacity),),
    search_keywords=(seed_samples=16, refinement_rounds=3),
)
```

Candidate confirmation records repeated local failure, successful feasibility
reproduction, or an inconclusive mix; none is global proof. The outer loop
replays the control values selected by each allocation while searching, then
re-optimizes them only when the enlarged utilization set is allocated again.

For a counterexample-guided allocation loop:

```julia
adaptive = solve_search_stable_operating_envelope(
    scenarios, cps;
    samples_per_round=32,
    max_rounds=4,
    control_policy=IssuePlusLocalLaws(),
    solve_keywords=(fairness=FairnessPolicy(
        kind=:max_min, normalization=:capacity),),
)
```

The complete screened set is added to the next allocation when a candidate or
policy conflict is found. `budget_exhausted` is not renamed as robust.

## 6. Measure local-solver branch sensitivity

Use deterministic multistart before treating a small objective difference as a
methodological result:

```julia
multi = solve_operating_envelope_multistart(
    scenarios, cps;
    start_scales=(1.0, 0.9, 1.1),
    security=:corners,
    control_policy=IssuePlusLocalLaws())

multi.selected
multi.diagnostics["capacity_spread_W"]
multi.diagnostics["maximum_primal_constraint_violation"]
```

Start scaling perturbs registered native variables and clamps starts to their
bounds. It is a diagnostic, not a global method. Publish all starts, statuses,
accepted runs, and capacity spread rather than only the most favorable run.

## 7. Preserve the study identity

Build the manifest from the exact networks and declarations used in the study:

```julia
spec = DOEStudySpec(scenarios, cps;
    security=:corners,
    control_policy=IssuePlusLocalLaws(),
    fairness=FairnessPolicy(kind=:max_min, normalization=:capacity),
    solver_options=(tol=1e-8,),
    seeds=Dict("scenario_generation" => 42),
    metadata=Dict(
        "dataset" => "feeder-export-2026-08",
        "case_builder_commit" => "<git commit>",
    ))

manifest = doe_study_manifest(spec)
manifest["study_id"]
```

The study identity covers SHA-256 network hashes, connection points, coverage,
control/fairness policies, solver options, seeds, versions, and supplied
provenance. Custom hook functions cannot be serialized, so record their file and
commit explicitly in metadata.

## Minimum evidence bundle

Archive these objects for each reported table row:

- `doe_study_manifest(spec)`;
- the complete `OperatingEnvelopeResult`, including control and primal-residual
  diagnostics;
- `OperatingEnvelopeVerification.context_results`;
- utilization-search points, normalized scores, refinement rounds, outcomes,
  and candidates;
- every multistart run and capacity spread; and
- an explicit statement that the result is locally solved, finitely tested, and
  either operationally non-anticipative or ideal-recourse by construction.

For repeatable batch comparisons, define a case file containing
`doe_benchmark_case()` with `nets`, `connection_points`, and optional
`label => NamedTuple` methods, then run:

```sh
julia --project=. scripts/run_doe_benchmark.jl case.jl results.tsv
```

The runner builds a separate study identity for each method and writes stable
interval rows. Use `doe_context_benchmark_rows` when a second table should retain
one row per verification context.

The committed `scripts/cases/doe_range_benchmark.jl` fixture is a complete
example of this interface. It records its synthetic-data licence, method labels,
seed provenance, scientific claim, and claim limitation in the study metadata.

See the [DOE development roadmap](../problems/doe_development_roadmap.md) for the
remaining certification, uncertainty-calibration, temporal, fairness, and
scaling agenda, and the [detailed scientific
audit](../problems/doe_quantification_review.md) for its rationale.

For typed calibration/test splits and empirical coverage reporting, continue to
[DOE scenario design and held-out evaluation](doe_scenario_design.md). The
[advanced uncertainty laboratory](doe_uncertainty_coverage.md) collects every
implemented uncertainty and statistical API in one longer workflow.

## What this does not prove

- A control appearing in the audit does not establish that its timing is
  implementable; the selected policy must be justified operationally.
- Perfect recourse is a legitimate formulation and upper benchmark, not
  automatically deployable capacity.
- Fixed-control replay with the same formulation and solver is not independent
  model validation.
- Multistart and adaptive search provide stronger numerical evidence but no
  global optimality, branch-completeness, or box-containment certificate.
