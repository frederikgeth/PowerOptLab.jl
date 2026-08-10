# Bespoke algorithms

A **bespoke algorithm** is a new *solution method*: it changes *how* a problem is
solved rather than *what* is solved. Where a [problem specification](../problems/multiperiod.md)
hands one monolithic model to Ipopt, an algorithm wraps the staged API (or the
engine's admittance-matrix primitives) in a custom loop.

## Available

- [**HELM power flow**](helm.md) — the Holomorphic Embedding Load-flow Method: a
  non-iterative 4-wire power flow that expands each voltage as a power series in
  a load-scaling parameter and evaluates it by Padé analytic continuation, with
  physical-residual, Padé-spread, coefficient-tail, and heuristic singularity
  diagnostics exposed for independent validation.

## Question-driven candidates not yet built

The next algorithms should support the [research program](../research_program.md),
rather than form an unconnected catalogue:

- **Continuation and branch discovery** — pseudo-arclength continuation,
  singular-point handling, deflation, and systematic multistart to identify which
  nonlinear solution a sensitivity or decision uses.
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

When one lands, add it under `src/algorithms/`, give it a page here, and document
the validity of its claims. See [Contributing](../contributing.md) for the pattern.
