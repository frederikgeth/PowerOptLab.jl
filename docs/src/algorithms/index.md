# Power flow and operability

Choose a method from the question about the operating point, its branch or the
network model. A nonlinear solve, a HELM calculation, a static operability check
and a neutral reduction provide different evidence.

| Your task | Guide | API |
|:--|:--|:--|
| Compare nonlinear power flow with the no-load-connected HELM branch | [HELM versus nonlinear power flow](../tutorials/helm_vs_nonlinear_power_flow.md) | [HELM API](../api/helm.md) |
| Audit a solved native static equilibrium and follow a declared stress path | [Post-OPF operability](../problems/operability.md) and [study tutorial](../tutorials/post_opf_operability_study.md) | [Operability API](../api/operability.md) |
| Study an ideal-ground neutral projection and its retained/lost physics | [Neutral reduction reference](../api/network_reduction.md) | [Reduction API](../api/network_reduction.md) |

HELM and operability have explicit supported-physics boundaries. Follow the
linked guides before applying their results to a network containing generators,
IBRs or controls. A fixed-point certificate, a path-specific continuation margin
and an endpoint residual are distinct claims.

## Available

- [**HELM power flow**](helm.md) — the Holomorphic Embedding Load-flow Method: a
  non-iterative 4-wire power flow that expands each voltage as a power series in
  a load-scaling parameter and evaluates it by Padé analytic continuation, with
  physical-residual, Padé-spread, coefficient-tail, and heuristic singularity
  diagnostics exposed for independent validation.
- **Static operability continuation** — `continue_opf_operability` provides the
  natural-parameter trace, while `continue_opf_operability_pseudo_arclength`
  provides an opt-in first fold-capable predictor/corrector slice for the same
  native static-load scope, including fixed-λ refinement when a target crossing
  is bracketed; `stop_at_target=false` continues a declared stress path beyond
  λ=1. `locate_opf_operability_fold` separately refines a declared fold
  candidate with bordered equations, and `stop_on_voltage_limit=true` can
  terminate at the first declared terminal-voltage violation. The
  `operability_continuation_margin` helper summarizes the first declared
  voltage-limit or localized fold boundary relative to λ=1; pseudo-arclength
  traces also retain tangent-turning curvature and accepted-step histories.
  `operability_continuation_rows` provides deterministic table-ready records.

## Question-driven candidates not yet built

The next algorithms should support the [research program](../research_program.md),
rather than form an unconnected catalogue:

- **Continuation and branch discovery** — the implemented path and margin
  contracts are documented in [Post-OPF operability](../problems/operability.md);
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
  valid bound can strengthen a forensics or robust-decision claim.
- **Decomposition** — spatial, temporal, or hypothesis decomposition only when
  the integrated inference or verification studies require it for scale.

## Implementation layers

A **bespoke algorithm** is a new *solution method*: it changes *how* a problem is
solved rather than *what* is solved. Where a [problem specification](../problems/multiperiod.md)
hands one monolithic model to Ipopt, an algorithm wraps the staged API (or the
engine's admittance-matrix primitives) in a custom loop.

See [Contributing](../contributing.md) for extension interfaces and source-tree
conventions. Register a new guide with the study topic it supports.
