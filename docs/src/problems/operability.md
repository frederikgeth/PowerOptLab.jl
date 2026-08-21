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
* positive-, negative-, and zero-sequence voltages and VUF on complete
  three-phase buses; and
* an opt-in HELM cross-check of the no-load-connected branch for supported
  constant-power/constant-impedance cases; and
* an explicit natural-parameter load-scale continuation trace with corrector,
  residual, singular-value, and event history; and
* an opt-in pseudo-arclength predictor/corrector trace for the same static
  load-scale path, including voltage-normalized arclength provenance, target
  crossings, fixed-λ target refinement, and fold-candidate events with critical
  left/right mode participation; and
* explicit `:not_applicable` scope evidence when generator or IBR equations are
  outside the residual seam.

Scaling is part of the evidence contract. Supply the staged OPF `context` so
the audited policy, AC coordinate bases, and research provenance are inherited,
or pass an explicit `OperabilitySpec(scaling_policy=...)`. Non-SI policies also
require `scaling_bases` keyed by bus; the checker never invents a per-unit base.

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
    ))

report.status                    # :pass, :fail, :inconclusive, or :not_applicable
report.checks["jacobian_regular"]
report.sensitivities["load_scale"]
report.sensitivities["directions"]["P"]
report.branch_evidence["critical_mode"]
report.branch_evidence["reachability"]

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
```
