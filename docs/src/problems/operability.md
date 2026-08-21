# Post-OPF voltage operability

`check_opf_operability` evaluates an already extracted, SI-valued power-flow
solution as a diagnostic point. The first implementation is intentionally
scoped to the native `BMOPFTools.ybus_linearized` static network residual:
constant-P, constant-I, constant-Z, ZIP, and exponential loads, with wye,
single-phase, and delta connection records.

The checker reports:

* the independent non-source current-balance residual;
* a coordinate-scaled voltage Jacobian, singular values, condition number, and
  a rank/regularity check with critical left/right singular-mode participation;
* connection-level terminal voltages, requested-versus-realized powers, and
  derivatives under uniform load scaling plus named P/Q perturbations;
* optional independent finite-difference re-solves that validate the implicit
  uniform-load and named P/Q sensitivities; and
* path-qualified `dP/dV` branch indicators, with explicit sign convention and
  high-side/near-nose/low-side classifications; and
* positive-/negative-sequence magnitude and VUF path derivatives, retained as
  complementary evidence rather than a replacement for per-phase checks; and
* positive-, negative-, and zero-sequence voltages and VUF on complete
  three-phase buses; and
* an opt-in HELM cross-check of the no-load-connected branch for supported
  constant-power/constant-impedance cases; and
* an explicit natural-parameter load-scale continuation trace with corrector,
  residual, singular-value, and event history; and
* an opt-in pseudo-arclength predictor/corrector trace for the same static
  load-scale path, including voltage-normalized arclength provenance, target
  crossings, fixed-λ target refinement, and fold-candidate events with critical
  left/right mode participation and optional bordered fold localization; and
* a path-specific continuation-margin summary that identifies the first
  declared voltage-limit or localized fold boundary relative to λ=1; and
* an opt-in Bernstein-style Z-bus fixed-point certificate on the native
  constant-P/I/Z, ZIP, and finite exponential-law scope, with connection-aware
  wye/delta incidence, conservative voltage-domain Lipschitz bounds, an
  invariant contraction region around the energized no-load solution, and an
  independent fixed-point oracle trace; and
* deterministic one-row-per-point continuation records for table-ready
  residual, conditioning, curvature, and event summaries; and
* bounded named stress-direction snapshots with explicit P/Q and
  connection-weight transformations, caller-supplied solve callbacks, and
  table-ready residual/certificate rows; and
* deterministic direction-level and model-ensemble summaries with first
  observed non-pass boundaries and finite path-specific margins; and
* `operability_snapshot_row`, a compact single-snapshot projection for study
  tables that does not imply contingency or envelope coverage; and
* explicit `:not_applicable` scope evidence when generator or IBR equations are
  outside the residual seam.

Scaling is part of the evidence contract. Supply the staged OPF `context` so
the audited policy, AC coordinate bases, and research provenance are inherited,
or pass an explicit `OperabilitySpec(scaling_policy=...)`. Non-SI policies also
require `scaling_bases` keyed by every bus; the checker rejects an incomplete
base map and never invents a per-unit base. Equivalent SI and per-unit audits
should preserve the physical conclusion (endpoint claim, voltage/VUF limits,
branch-indicator classification, and certificate status), while normalized
residuals and singular values remain coordinate-dependent evidence.
The current equilibrium closure is explicitly `:frozen_dispatch` and is
recorded in provenance; other closures are rejected until their equations are
available through a public upstream seam.

The aggregate result is `:pass` only when a requested primary operational claim
(terminal-voltage bounds, sequence-unbalance bounds, HELM reachability, or the
fixed-point certificate) has passed. If no such claim was requested, the result
remains `:not_applicable` even when endpoint and Jacobian evidence are good.

If the network contains a generator, IBR, unsupported load law, or unsupported
connection, the checker returns `:not_applicable` before evaluating a partial
residual. The reasons are retained in `result.unsupported`, the `scope` check,
and `result.provenance["operability"]`; this is explicit scope evidence, not a
pass and not a claim that the omitted equations are feasible.

The fixed-point certificate is a sufficient, local-in-model-region result. It
uses the source-eliminated implicit Z-bus map and conservative connection-aware
componentwise complex-voltage polydisc bounds around the no-load solution,
retaining a uniform-radius fallback. A `:pass` means the
candidate lies inside a contraction region that contains a unique equilibrium
and has a nonsingular equilibrium Jacobian on that region. A failed sufficient
condition, unsupported load model, or candidate outside the region is
`:inconclusive`; it is not evidence of non-existence, multiplicity, or a
low-voltage branch.

```julia
using BMOPFTools
using PowerOptLab

pf = solve_pf(net; per_unit=false)
report = check_opf_operability(net, pf;
    spec = OperabilitySpec(
        scaling_policy = SIUnitsScaling(),
        voltage_min = 0.95 * vbase,
        voltage_max = 1.05 * vbase,
        vuf_max = 0.02,
        compute_helm = true,
        compute_fixed_point_certificate = true,
    ))

report.status                    # :pass, :fail, :inconclusive, or :not_applicable
report.checks["jacobian_regular"]
report.sensitivities["load_scale"]
report.sensitivities["directions"]["P"]
report.branch_evidence["critical_mode"]
report.branch_evidence["dP_dV"]
report.branch_evidence["sequence_sensitivity"]
report.branch_evidence["reachability"]
report.branch_evidence["fixed_point_certificate"]
# Includes `connection_law_terms`, the certified radius, lower connection
# voltage, contraction factor, finite-difference `law_bound_validation`, and
# the independent Z-bus oracle summary.
report.provenance["operability"]["model_inventory"]

trace = continue_opf_operability(net, pf;
    spec = OperabilitySpec(scaling_policy = SIUnitsScaling()),
    continuation = OperabilityContinuationSpec(initial_step = 0.1))
trace.status
trace.events

pseudo = continue_opf_operability_pseudo_arclength(net, pf;
    spec = OperabilitySpec(scaling_policy = SIUnitsScaling()),
    continuation = OperabilityPseudoArclengthSpec(initial_step = 0.05))
pseudo.status
pseudo.provenance["continuation"]["pseudo_arclength"]
pseudo.provenance["continuation"]["arclength_state_scale"]
pseudo.provenance["continuation"]["curvature_history"]

# To continue the declared stress path beyond the audited endpoint:
stress = continue_opf_operability_pseudo_arclength(net, pf;
    spec = OperabilitySpec(scaling_policy = SIUnitsScaling()),
    continuation = OperabilityPseudoArclengthSpec(max_steps = 200),
    stop_at_target = false)
stress.events
stress.provenance["continuation"]["lambda_min"]

# Optionally terminate when a declared terminal-voltage limit is first crossed:
limited = continue_opf_operability_pseudo_arclength(net, pf;
    spec = OperabilitySpec(scaling_policy = SIUnitsScaling(), voltage_min = 0.95vbase),
    continuation = OperabilityPseudoArclengthSpec(),
    stop_on_voltage_limit = true)
limited.status
limited.events
limited.provenance["continuation"]["margin"]
rows = operability_continuation_rows(pseudo)
rows[1].lambda
rows[end].event_kinds

directions = [
    OperabilityStressDirection(:uniform),
    OperabilityStressDirection(:reactive_removed; p_scale=1.0, q_scale=0.0),
]
stress_rows = operability_stress_rows(net, pf;
    spec = OperabilitySpec(scaling_policy = SIUnitsScaling(),
                           compute_fixed_point_certificate = true),
    directions = directions,
    lambdas = [0.0, 0.5, 1.0],
    solve = network -> solve_pf(network; per_unit=false))
stress_rows
operability_stress_summary(stress_rows)
operability_stress_ensemble_rows(Dict("base" => stress_rows))
operability_snapshot_row(report; snapshot_id = "base")

validated = check_opf_operability(net, pf;
    spec = OperabilitySpec(
        scaling_policy = SIUnitsScaling(),
        compute_sensitivity_validation = true))
validated.checks["load_scale_sensitivity_validation"]
validated.checks["directional_sensitivity_validation"]
validated.sensitivities["validation"]["load_scale"]
validated.sensitivities["validation"]["directions"]

# Starting from a voltage/state near a suspected nose, solve the bordered
# equations F=0, J*v=0, ||v||₂=1 at a declared load scale.
fold = locate_opf_operability_fold(net, approximate_fold_solution;
    lambda = 1.04,
    spec = OperabilitySpec(scaling_policy = SIUnitsScaling()))
fold.status
fold.lambda
fold.critical_mode
```

This is a local post-solve screen, not a proof of global solvability or dynamic
voltage stability. HELM agreement is path-qualified evidence for its energized
no-load homotopy; HELM non-convergence is inconclusive, not a non-existence
certificate. The current pseudo-arclength implementation is a first static
slice: it does not yet provide refined fold points, deflation/global branch
discovery, controller closures, or fuller generator/IBR equilibrium seams. A
low-voltage equilibrium can still pass the local endpoint residual check while
failing the no-load-connected target comparison; that distinction is deliberate
branch evidence, not a claim that the low branch is infeasible.
`locate_opf_operability_fold` is likewise local: convergence establishes a
bordered-equation fold candidate for the declared static model and initial
guess, not uniqueness or global reachability.

```@docs
OperabilitySpec
OperabilityCheck
OperabilityResult
check_opf_operability
OperabilityContinuationSpec
OperabilityContinuationResult
continue_opf_operability
OperabilityPseudoArclengthSpec
continue_opf_operability_pseudo_arclength
OperabilityFoldResult
locate_opf_operability_fold
operability_continuation_margin
operability_continuation_rows
OperabilityStressDirection
operability_stress_network
operability_stress_rows
operability_stress_summary
operability_stress_ensemble_rows
operability_snapshot_row
```
