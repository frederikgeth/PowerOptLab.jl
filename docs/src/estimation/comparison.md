# Choosing a formulation

The native WLS and compiled estimators use the passive network physics and a Gaussian measurement
objective. Their APIs differ in injection coverage, supported exact device
models, diagnostics, and computational cost.

## The one-sentence version

For new work, use the compiled API with dense, sparse, or classic Ipopt solution
over one shared evaluator. See [measurement and electrical compatibility](compatibility.md).

Use the **[legacy WLS](../problems/state_estimation.md)** estimator for a
JuMP/Ipopt model with complete paired P/Q injection coverage and explicit zero
injections. Use the **[constrained NLLS](../problems/constrained_state_estimation.md)**
estimator for exact nonlinear device laws, branch telemetry, covariance, or a
time series. Both retain ungrounded neutral voltages and impose exact network
equations; uncertain readings remain weighted residuals in both.

## Side by side

| | Legacy WLS | Constrained NLLS |
|---|---|---|
| **Problem** | one weighted sum of squares | ``\min \tfrac12\lVert r\rVert^2`` subject to ``c(x)=0`` |
| **Solver** | Ipopt (general interior-point NLP) | Byrd–Omojokun composite step (dense) or Hachtel augmented system with sparse QR, or shared Ipopt |
| **State** | rectangular ``V`` at non-source nodes, plus a free injection current per measured terminal | rectangular ``V`` of free conductors only |
| **Neutrals** | explicit BMOPFTools voltage variables with terminal-to-reference injection currents | explicit conductor voltages and compiled return-current constraints |
| **Node measurements** | `:vr`, `:vi`, `:vmag`, `:pinj`, `:qinj` | `:vr`, `:vi`, `:vmag`, `:pinj`, `:qinj` |
| **Branch measurements** | none | `:ire`, `:iim`, `:imag`, `:pflow`, `:qflow` on finite lines, transformer ports and shunts; net-injection and custom scalar functions also available |
| **Exact device laws** | none | `ConstantPowerDevice`, `ConstantCurrentDevice`, `ZIPDevice` |
| **Zero injection** | hard (the engine's KCL), declared per phase | hard `c(x)=0`, declared per conductor |
| **Input contract** | rejects loads, limits, and uncovered phase injections | electrical preflight rejects unsupported models and reports omissions; no complete injection-coverage requirement |
| **Observability** | measurement rank on the exact-KCL tangent space at the solution | rank of ``HZ`` on the feasible tangent space, where ``CZ=0`` |
| **Covariance** | not provided | `selected_state_covariance`, `derived_covariance` |
| **Time series / priors** | no | `solve_time_series_state_estimator`, `StatePrior` |
| **Warm start** | phase-aware voltage and injection-current starts, seeded from the readings | explicit `x0` or load-free initialiser, reused across snapshots |
| **Error model** | Gaussian today; the JuMP model can express solver-compatible ``-\log f`` objectives and exact reformulations of WLAV/Huber | Gaussian today; robust/general likelihoods require changes to step weights, globalisation, and uncertainty calculations |
| **Bad data** | neither. `standardized` is a raw σ-normalised residual; multipliers are leads | |

## The differences that actually decide it

### Exact information and the available APIs

Weighted least squares specifies an objective, not whether constraints are
allowed. Equality-constrained WLS is classical, and both implementations impose
hard zero injections. The API distinction here is support for explicit exact
nonlinear device equations in the compiled estimator. Treating an uncertain
reading as a very small-variance residual is different from imposing an equality.

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

Both estimators preserve explicitly modelled, ungrounded neutral voltages.
The legacy estimator attaches free current to a measured terminal and its
return with opposite signs. The compiled estimator eliminates those currents
and imposes a combined current-sum equation for each connected injection return
component. Several measured phases sharing one neutral therefore contribute one
combined balance, independently of the uncertain P/Q readings.
Use `injection_connections` in the compiled API to preserve those paths when
power telemetry drops out; the default infers paths from the available P/Q rows.
Distinct original customer terminals retain separate powers at merged buses,
whose powers and independent variances add. Repeated readings of the same
original customer instead share its power and combine by precision.

The compiled formulation preserves floating neutrals. Native WLS implicitly
grounds source neutrals, so its wrapper requires an explicit source earth bond. Supply the physical grounding and return paths, then inspect the
rank diagnostics. The compiled interface also supports explicit nonlinear
device laws and selective covariance for neutral-related quantities.

Note the matching difference in the `zero_injection` keyword: a bare bus id
expands over *phase* terminals for WLS and over *all* conductors, neutral
included, for the compiled estimator. Use explicit conductor tuples when moving between APIs.

### Error models run the other way

Almost everything on this page favours the constrained estimator. Error models
do not. Minimising ``\tfrac12\lVert r\rVert^2`` is the Gaussian
log-likelihood for standardised measurement errors, and the implemented
Gauss--Newton step, Hachtel system and ``(H^\top H)^{-1}`` covariance all use
that structure. Robust losses are possible in specialised least-squares
algorithms, but would require new weights, convergence logic and uncertainty
semantics here. A general NLP is the easier experimental platform and can
represent WLAV exactly through an epigraph, at the cost of additional variables
and inequalities. See [likelihood, loss, and priors](maximum_likelihood.md).

### The input contract

`solve_state_estimation` refuses a network that carries loads, generators,
IBRs, or operational limits, and refuses a bus that is neither measured nor
declared zero-injection. That contract catches a real class of silent error:
a retained load adds its injection *on top of* the estimated one, and a
retained voltage limit pins the estimate to the limit — both while still
reporting `LOCALLY_SOLVED`.

`compile_state_estimator` runs an electrical support preflight and records
network injections and operating limits that it omits. It does not require
complete injection telemetry or reject every underdetermined problem: rank
diagnostics describe the information available at the evaluation point.
Inspect both `structure.preflight` and `observability_diagnostics`, and see
[observability](observability.md) for how to read them and what to do when a
set turns out to be under-observed.

## Reading the status

Neither estimator's success status is a uniqueness claim.

* WLS: check `solve_status(result).publishable` (termination and primal feasibility), then read
  `observability.observable`. Voltages are `NaN` when the solve did not
  converge — an unconverged iterate is not published as an estimate.
* Constrained: read `:converged_unique` and `:converged_underobserved` as first-order termination; the
  latter is feasible but not identified. `:converged_unique` reports that the
  reduced Jacobian has full rank **at the returned point**. Two different
  points can both earn it — see
  [current-magnitude measurements](current_magnitude.md), where a mirror pair
  of states both fit the data exactly and both come back `:converged_unique`.

## Agreement between them

Tests establish agreement on specific compatible problems. `test/` checks the
compiled estimator against a BMOPFTools `solve_pf` state recovered from exact
telemetry, and checks the dense and sparse solvers against each other on that
problem; the WLS estimator is checked against a power flow of the same feeder.
The regression suite also compares both compiled solvers directly with the
legacy JuMP/Ipopt fit on noisy phase-to-neutral measurements while checking
exact return-current balance. Merged-customer tests check summed power,
unequal uncertainties, repeated readings, missing components and row permutations
against the JuMP formulation, including profiled voltage information and
covariance against aggregate telemetry. These small cases establish agreement within
their tested models and tolerances, rather than equivalence on every network.

## Known limitations of both

* Diagonal covariance only; no correlated whitening.
* No χ² test, no largest-normalised-residual bad-data identification, no
  topology-error hypothesis search.
* No multistart or global search on a nonconvex problem.
* Both observability diagnostics use dense SVDs and are sized for small and
  medium cases. The sparse solver's *step* is sparse; its rank diagnostic runs
  only on exit paths, and is dense.
* Jacobian assembly is sparse, but rank/covariance decompositions remain dense. See
  [the open questions](state_of_the_art.md#Open-questions).
