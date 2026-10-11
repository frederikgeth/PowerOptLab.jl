# Power flow and voltage operability

The methods in this chapter calculate power-flow solutions, examine static
voltage behaviour and reduce the network model under specified grounding
assumptions.

| Study | Guide | Function reference |
|:--|:--|:--|
| Compare an iterative power flow with the HELM solution connected to the no-load state | [HELM comparison](../tutorials/helm_vs_nonlinear_power_flow.md) | [HELM API](../api/helm.md) |
| Check a solved equilibrium and trace voltage response as loading changes | [Post-OPF operability](../problems/operability.md) and [worked study](../tutorials/post_opf_operability_study.md) | [Operability API](../api/operability.md) |
| Compare an explicit neutral model with an ideally grounded neutral reduction | [Neutral reduction](../api/network_reduction.md) | [Reduction API](../api/network_reduction.md) |

The HELM and operability guides list the load models, connections and controls
they support. Check the equation residuals and operating limits at the calculated point.
Continuation gives information along a specified loading path; sufficient
fixed-point conditions establish existence or uniqueness only in the region
covered by those conditions.

## Available methods

- [**HELM power flow**](helm.md) — the Holomorphic Embedding Load-flow Method: a
  non-iterative 4-wire power flow that expands each voltage as a power series in
  a load-scaling parameter and evaluates it by Padé analytic continuation, with
  equation residuals, agreement between Padé approximants, coefficient decay
  and a heuristic singularity estimate for assessing numerical convergence.
- **Static operability continuation** — `continue_opf_operability` provides the
  natural-parameter trace, while `continue_opf_operability_pseudo_arclength`
  provides a pseudo-arclength predictor–corrector method for the supported
  static load models, including refinement when the trace crosses a target
  loading level; `stop_at_target=false` continues a declared stress path beyond
  λ=1. `locate_opf_operability_fold` separately refines a declared fold
  candidate with bordered equations, and `stop_on_voltage_limit=true` can
  terminate at the first declared terminal-voltage violation. The
  `operability_continuation_margin` helper summarizes the first declared
  voltage-limit or localized fold boundary relative to λ=1; pseudo-arclength
  traces also retain tangent-turning curvature and accepted-step histories.
  `operability_continuation_rows` returns the trace as records for study tables.

## Further research

The [research programme](../research_program.md) identifies several extensions
needed for studies with uncertain network models and operating conditions:

- **Continuation and branch discovery** — the implemented path and margin
  methods are described in [Post-OPF operability](../problems/operability.md);
  the current pseudo-arclength slice still needs refined singular-point
  handling, deflation, and systematic multistart to identify which nonlinear
  solution a sensitivity or decision uses.
- **Sensitivity validation** — KKT regularity, active-set transition, and
  finite-difference/continuation comparisons for DiffOpt and, where useful,
  direct sensitivities through the HELM coefficient recursion.
- **Counterexample generation** — adversarial utilization/model searches for
  falsifying operating-envelope candidates, with explicit local-search stopping
  semantics.
- **Relaxations and valid bounds** — multiconductor convex relaxations when a
  valid bound can strengthen network identification or robust operating decisions.
- **Decomposition** — spatial, temporal, or hypothesis decomposition only when
  the integrated inference or verification studies require it for scale.

## Adding a solution method

A solution method changes the numerical procedure used to solve the network
problem. Examples include continuation, decomposition and successive
linearisation. These methods can reuse the BMOPFTools network equations or
admittance matrices. [Contributing](../contributing.md) describes the programming
interfaces and how to document a new method.
