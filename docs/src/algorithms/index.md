# Bespoke algorithms

A **bespoke algorithm** is a new *solution method*: it changes *how* a problem is
solved rather than *what* is solved. Where a [problem specification](../problems/multiperiod.md)
hands one monolithic model to Ipopt, an algorithm wraps the staged API (or the
engine's admittance-matrix primitives) in a custom loop.

## Available

- [**HELM power flow**](helm.md) — the Holomorphic Embedding Load-flow Method: a
  non-iterative 4-wire power flow that expands each voltage as a power series in
  a load-scaling parameter and evaluates it by Padé analytic continuation, so
  voltage collapse is a *certified* outcome (Stahl's theorem) and the loading
  margin falls out of every solve.

## Candidates not yet built

- **Decomposition** — spatial (per-feeder) or temporal (per-snapshot) splitting
  with a coordinating master problem.
- **Sequential linearization / SLP-SQP** — solve a sequence of tractable
  approximations to convergence.
- **ADMM / operator splitting** — distributed consensus across sub-networks.
- **Warm-start schemes** — reuse one solve's state to accelerate the next
  (e.g. across the snapshots of a multi-period run).

When the next one lands, add it under `src/algorithms/` and give it a page here.
See [Contributing](../contributing.md) for the pattern.
