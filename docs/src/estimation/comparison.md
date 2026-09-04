# Choosing a formulation

Both estimators fit the same physics to the same telemetry. The choice between
them is a modelling choice, not a performance one.

## The one-sentence version

Use the **[legacy WLS](../problems/state_estimation.md)** estimator when every
piece of information you have is uncertain and the network is
grounded-neutral. Use the **[constrained NLLS](../problems/constrained_state_estimation.md)**
estimator when something in your model is *exactly* true, when neutral
displacement is part of the answer, or when you need branch telemetry,
covariance, or a time series.

## Side by side

| | Legacy WLS | Constrained NLLS |
|---|---|---|
| **Problem** | one weighted sum of squares | ``\min \tfrac12\lVert r\rVert^2`` subject to ``c(x)=0`` |
| **Solver** | Ipopt (general interior-point NLP) | Byrd–Omojokun composite step (dense) or Hachtel augmented system with sparse QR |
| **State** | rectangular ``V`` at non-source nodes, plus a free injection current per measured terminal | rectangular ``V`` of free conductors only |
| **Neutrals** | per-phase contract; a modelled neutral is not a first-class state | explicit conductor; an ungrounded neutral stays in the state |
| **Node measurements** | `:vmag`, `:pinj`, `:qinj` | `:vr`, `:vi`, `:vmag`, `:pinj`, `:qinj` |
| **Branch measurements** | none | `:ire`, `:iim`, `:imag`, `:pflow`, `:qflow` on a named line |
| **Exact device laws** | none | `ConstantPowerDevice`, `ConstantCurrentDevice`, `ZIPDevice` |
| **Zero injection** | hard (the engine's KCL), declared per phase | hard `c(x)=0`, declared per conductor |
| **Input contract** | strict — rejects loads, limits, and uncovered buses | none; rank diagnostics instead |
| **Observability** | rank of the measurement + KCL Jacobian at the solution | rank of ``HZ`` on the feasible tangent space, where ``CZ=0`` |
| **Covariance** | not provided | `selected_state_covariance`, `derived_covariance` |
| **Time series / priors** | no | `solve_time_series_state_estimator`, `StatePrior` |
| **Warm start** | Ipopt-internal, seeded from the readings | explicit `x0`, reused across snapshots |
| **Error model** | Gaussian today; the JuMP objective could express any smooth ``-\log f``, and WLAV/Huber exactly by reformulation | Gaussian **structurally** — Gauss–Newton and the Hachtel step require a sum of squares |
| **Bad data** | neither. `standardized` is a raw σ-normalised residual; multipliers are leads | |

## The differences that actually decide it

### Exact information has nowhere to go in WLS

The WLS statement has exactly one channel for information: a residual with a
σ. If a quantity is known exactly, the only way to express it is a very small
σ — which is not the same statement. It makes the problem ill-conditioned
rather than constrained, and it lets a conflicting measurement trade against a
physical law that should not be negotiable.

The constrained formulation keeps three categories apart, and the separation
is the point:

| Category | Where it lives | Examples |
|---|---|---|
| Network identity | compiled into `SEStructure` | Ybus, conductor incidence, closed-switch aliases |
| Exact equation | ``c(x)=0`` | true zero injection, an exact device law |
| Stochastic information | whitened residual ``r(x)`` | meters, forecasts, nominal loads, priors |

A small variance does not make a reading exact. Exact constraints shrink the
feasible tangent space and *can make a model inconsistent* — which is
informative, and much better than a soft compromise that hides the conflict.

### Four-wire networks

The WLS estimator states zero injection per phase, against a return terminal.
The constrained estimator states KCL per conductor against earth, so an
ungrounded neutral is a state and neutral displacement is estimated rather
than assumed away. If your question involves neutral-to-earth voltage,
unbalance, or a broken neutral, only the second formulation can express it.

Note the matching difference in the `zero_injection` keyword: a bare bus id
expands over *phase* terminals for WLS and over *all* conductors, neutral
included, for the compiled estimator. Both are correct for their own
formulation.

### Error models run the other way

Almost everything on this page favours the constrained estimator. Error models
do not. Minimising ``\tfrac12\lVert r\rVert^2`` *is* the Gaussian
log-likelihood, and the Gauss–Newton step, the Hachtel system and the
``(H^\top H)^{-1}`` covariance all depend on that structure — so Laplace,
Huber, Beta and mixtures are not expressible there at all. An NLP formulation
has no such restriction and gets WLAV exactly, for almost no cost, via an
epigraph reformulation. See [estimators, not curve fits](maximum_likelihood.md).

### The input contract

`solve_state_estimation` refuses a network that carries loads, generators,
IBRs, or operational limits, and refuses a bus that is neither measured nor
declared zero-injection. That contract catches a real class of silent error:
a retained load adds its injection *on top of* the estimated one, and a
retained voltage limit pins the estimate to the limit — both while still
reporting `LOCALLY_SOLVED`.

`compile_state_estimator` has no equivalent contract. It is the more
expressive interface and the sharper instrument; it will happily compile an
under-determined problem and tell you so afterwards through
`observability_diagnostics`. Check the diagnostics; nothing else will — and see
[observability](observability.md) for how to read them and what to do when a
set turns out to be under-observed.

## Reading the status

Neither estimator's success status is a uniqueness claim.

* WLS: trust `primal_status == "FEASIBLE_POINT"`, then read
  `observability.observable`. Voltages are `NaN` when the solve did not
  converge — an unconverged iterate is not published as an estimate.
* Constrained: trust `:converged_unique` and `:converged_underobserved`; the
  latter is feasible but not identified. `:converged_unique` reports that the
  reduced Jacobian has full rank **at the returned point**. Two different
  points can both earn it — see
  [current-magnitude measurements](current_magnitude.md), where a mirror pair
  of states both fit the data exactly and both come back `:converged_unique`.

## Agreement between them

Where both can express the same problem, they agree. `test/` checks the
compiled estimator against a BMOPFTools `solve_pf` state recovered from exact
telemetry, and checks the dense and sparse solvers against each other on that
problem; the WLS estimator is checked against a power flow of the same feeder.
The two formulations are validated against the same physics, not against each
other's output.

## Known limitations of both

* Diagonal covariance only; no correlated whitening.
* No χ² test, no largest-normalised-residual bad-data identification, no
  topology-error hypothesis search.
* No multistart or global search on a nonconvex problem.
* Both observability diagnostics use dense SVDs and are sized for small and
  medium cases. The sparse solver's *step* is sparse; its rank diagnostic runs
  only on exit paths, and is dense.
* Jacobian assembly in the compiled estimator is still dense. See
  [the roadmap](state_of_the_art.md#roadmap).
