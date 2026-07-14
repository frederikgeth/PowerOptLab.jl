# Bespoke algorithms

*This section is a reserved slot — there are no bespoke algorithms yet.*

A **bespoke algorithm** is a new *solution method*: it changes *how* a problem is
solved rather than *what* is solved. Where a [problem specification](../problems/multiperiod.md)
hands one monolithic model to Ipopt, an algorithm wraps the staged API in a custom
loop. Candidates include:

- **Decomposition** — spatial (per-feeder) or temporal (per-snapshot) splitting
  with a coordinating master problem.
- **Sequential linearization / SLP-SQP** — solve a sequence of tractable
  approximations to convergence.
- **ADMM / operator splitting** — distributed consensus across sub-networks.
- **Warm-start schemes** — reuse one solve's state to accelerate the next
  (e.g. across the snapshots of a multi-period run).

When the first one lands, add it under `src/algorithms/` and give it a page here.
See [Contributing](../contributing.md) for the pattern.
