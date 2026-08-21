# Single-snapshot post-OPF operability study

> **Audience:** power-system researchers · **Scope:** one solved native-static
> equilibrium, with optional finite stress and load-model comparisons.

This tutorial shows how to turn one SI-valued OPF/PF solution into a compact,
auditable operability record. The workflow is deliberately snapshot-first:
contingencies and operating envelopes are not implied. If several cases are
available, evaluate each case independently and attach the case label outside
the processor.

The repository includes a runnable companion study at
`scripts/run_post_opf_operability_study.jl`. It builds a small balanced
four-wire feeder, evaluates constant-power and ZIP load models, and prints the
snapshot rows plus finite stress summaries. It also includes an unbalanced
DELTA case and a floating-neutral WYE case so the connection-level voltage,
sequence, and neutral-displacement evidence are exercised. A final
single-phase transformer case checks a separate low-voltage zone and its
local voltage bounds:

```sh
julia --project=. scripts/run_post_opf_operability_study.jl
```

The script asserts that all declared snapshots and finite rows pass. Its final
`PASS` line means only that this finite, declared study completed; it is not a
global solvability or operating-envelope claim.

## 1. Define the claim before running the checker

`check_opf_operability` separates several claims that are often conflated:

| Claim | Evidence | Typical status |
|---|---|---|
| The reported point satisfies the audited static equations | endpoint residual and domain checks | `:pass` or `:fail` |
| Declared voltage or VUF limits hold | terminal and sequence checks | `:pass` or `:fail` |
| The local equilibrium is regular | scaled Jacobian and critical mode | `:pass` or `:inconclusive` |
| The point is connected to the energized no-load branch | HELM or continuation evidence | `:pass` or `:inconclusive` |
| A sufficient local uniqueness region contains the point | fixed-point certificate | `:pass` or `:inconclusive` |

The aggregate result is not a hidden stability verdict. A certificate or
reachability failure is conservative evidence, not proof of non-existence or
dynamic voltage instability.

## 2. Audit one solved point

For a workflow that may receive heterogeneous network dictionaries, run the
cheap scope audit before solving. It distinguishes native frozen-dispatch
readiness from an unsupported control seam and catches missing or dangling
voltage-source references without treating them as infeasibility:

```julia
scope = operability_scope_audit(net)
scope["status"]
scope["topology"]
scope["unsupported_reasons"]
```

The checker consumes the SI-valued solution returned by BMOPFTools. Pass the
audited scaling policy explicitly when the solution is detached from its OPF
context. For a non-SI policy, provide voltage and current bases for every bus;
the checker rejects incomplete base maps rather than inventing defaults.

```julia
using BMOPFTools
using PowerOptLab

pf = solve_pf(net; per_unit=false)

spec = OperabilitySpec(
    scaling_policy = SIUnitsScaling(),
    voltage_min = 0.95 * vbase,
    voltage_max = 1.05 * vbase,
    vuf_max = 0.02,
    compute_helm = true,
    compute_fixed_point_certificate = true,
)

report = check_opf_operability(net, pf; spec)
report.status
```

The detailed result keeps the ordering and physical interpretation needed for
research review:

```julia
report.checks["endpoint"]
report.checks["terminal_voltage_bounds"]
report.checks["sequence_unbalance"]
report.branch_evidence["dP_dV"]
report.branch_evidence["critical_mode"]
report.branch_evidence["reachability"]
report.branch_evidence["fixed_point_certificate"]
report.provenance["operability"]["model_inventory"]
```

The `dP_dV` records distinguish requested from realized power and retain the
path-qualified sign convention. Positive- and negative-sequence derivatives
are complementary summaries; they must not replace phase- or connection-level
evidence in an unbalanced system.

## 3. Produce a study-table row

Use `operability_snapshot_row` when the study needs one deterministic row per
case. It includes endpoint and conditioning evidence, voltage/VUF extrema,
branch-indicator counts, HELM/certificate statuses, and unsupported-scope
counts. It is intentionally not a contingency or envelope row.

```julia
row = operability_snapshot_row(report; snapshot_id="base")

row.status
row.minimum_terminal_voltage
row.maximum_vuf
row.high_side_indicator_count
row.near_nose_indicator_count
row.fixed_point_certificate_status
row.helm_reachability_status
row.scope                    # "single_snapshot_static_ybus"
```

For a collection of already solved cases, map the stable checker over the
cases and retain the case identifier in `snapshot_id`:

```julia
rows = [
    operability_snapshot_row(
        check_opf_operability(case.network, case.solution; spec),
        snapshot_id = case.id,
    )
    for case in solved_cases
]
```

This is the appropriate place to add a contingency label later. The processor
itself still evaluates one network/solution pair at a time.

## 4. Add finite, declared stress directions

Stress campaigns are useful for sensitivity studies, but their evidence is
finite and path-specific. Supply the solve callback explicitly and keep the
same checker specification for each point:

```julia
directions = [
    OperabilityStressDirection(:uniform),
    OperabilityStressDirection(:reactive_removed;
        p_scale=1.0, q_scale=0.0),
]

stress_rows = operability_stress_rows(net, pf;
    spec = OperabilitySpec(
        scaling_policy = SIUnitsScaling(),
        compute_fixed_point_certificate = true,
    ),
    directions = directions,
    lambdas = [0.0, 0.5, 1.0],
    solve = network -> solve_pf(network; per_unit=false),
)

stress_summary = operability_stress_summary(stress_rows)
```

Each row retains the direction, loading parameter, endpoint status, certificate
status, condition margin, and any solve/checking error. A missing non-pass
boundary is reported as `:not_observed`; it is not an unlimited margin.

For ZIP or exponential uncertainty, label each finite campaign and combine the
rows with `operability_stress_ensemble_rows`:

```julia
ensemble = operability_stress_ensemble_rows(Dict(
    "nominal" => stress_rows,
    "alternative_load_model" => alternative_rows,
))
```

The resulting table supports comparisons across declared models. It does not
promote a finite ensemble into a global operating envelope.

## 5. Interpret unsupported scope explicitly

The current residual scope is the native static `ybus_linearized` model:
constant-P/I/Z, ZIP, and finite exponential loads with supported WYE,
SINGLE_PHASE, and DELTA connections. Generators, IBR controls, unsupported
load laws, and unsupported connections are returned as `:not_applicable` with
reasons in three places:

```julia
report.status
report.unsupported
report.checks["scope"]
report.provenance["operability"]
```

Do not reinterpret `:not_applicable` as a pass or a failed power-flow solve.
It means that the requested claim lies outside this processor's declared
equation seam.

## 6. Minimum publication record

For each snapshot, retain:

1. the network/model identifier and the SI-valued solution provenance;
2. the scaling policy and complete coordinate bases when non-SI coordinates
   were used;
3. residual tolerances, voltage/VUF limits, and closure (`:frozen_dispatch` in
   the current native slice);
4. the snapshot row plus detailed `dP_dV`, sequence, Jacobian, and certificate
   evidence when those checks were requested; and
5. any finite stress direction, lambda grid, solver callback version, and
   non-pass boundary observed.

This record makes the positive claim auditable without implying contingency
coverage, global branch discovery, or dynamic voltage stability.

See [Post-OPF operability](../problems/operability.md) for the full result
contract and [the roadmap](../problems/post_opf_operability_roadmap.md) for the
remaining certificate and study extensions.
