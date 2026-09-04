# State estimation

> **Kind:** Problem specification · **Direction:** inverse · **Maturity:** prototype

State estimation is the inverse of power flow. Power flow asks *given the
injections, what is the state?* Estimation asks *given noisy, incomplete,
partly-wrong telemetry, what state best explains it — and how much should I
trust the answer?*

PowerOptLab carries **two independent estimators** over the same network
physics. They are not two solvers for one problem; they are two different
problem statements, and the difference matters more than the numerics:

| | [Legacy WLS](../problems/state_estimation.md) | [Constrained NLLS](../problems/constrained_state_estimation.md) |
|---|---|---|
| Entry point | `solve_state_estimation` | `compile_state_estimator` + a solver |
| Statement | ``\min \sum_i (z_i - h_i(x))^2/\sigma_i^2`` | ``\min \tfrac12\lVert r(x)\rVert^2`` s.t. ``c(x)=0`` |
| Engine | JuMP + Ipopt | own Gauss–Newton / Hachtel solvers |

[Choosing a formulation](comparison.md) is the head-to-head, with the
measured behaviour of both on the same problem.

## What both of them assume

Read this before either page. These assumptions are shared, they are not all
enforceable, and every one of them is a way an estimate can be confidently
wrong:

* **Errors are independent, zero-mean, and Gaussian**, with a diagonal
  covariance. Both estimators whiten by ``1/\sigma``; neither supports a
  correlated covariance. A meter whose error is a calibration offset rather
  than noise violates this, and no amount of redundancy repairs it. This is a
  *likelihood* assumption, not a weighting convention — see
  [likelihood, loss, and priors](maximum_likelihood.md).
* **Topology and line parameters are exact.** Neither estimator hypothesises a
  different switch status or a different impedance. A topology error surfaces
  as a large residual that looks exactly like bad data. See
  [parameter estimation](../problems/parameter_estimation.md) and
  [inverse Carson](../problems/inverse_carson.md) for the complementary
  problems.
* **The estimate is a local stationary point.** Both problems are nonconvex.
  Convergence is not uniqueness, and neither estimator does multistart. The
  [current-magnitude tutorial](current_magnitude.md) shows two different
  states that both fit the same data exactly, each reported as converged.
* **Redundancy is not observability, and observability is not a binary.**
  Placement decides what is identified, and the smallest singular value
  decides whether the identification is worth anything. See
  [observability](observability.md).
* **A σ is a modelling claim, not a formality.** It sets the weight, and it is
  the only channel through which you say how much a reading should move the
  answer. A pseudo-measurement given a meter's σ *is* a meter.

## Where to go next

* [Choosing a formulation](comparison.md) — the comparison, and which to reach for.
* [Likelihood, loss, and priors](maximum_likelihood.md) — when these estimators
  coincide with Gaussian MLE, when a prior makes the result MAP, and which error
  models each engine can represent.
* [Legacy WLS state estimation](../problems/state_estimation.md) — reference.
* [Constrained NLLS state estimation](../problems/constrained_state_estimation.md) — reference.
* [Modelling tutorial](../tutorials/constrained_nlls_state_estimation.md) — exact
  vs uncertain information, grounding, and diagnostics.
* [Observability and under-observed estimation](observability.md) — which
  parts of an estimate are supported by measurements and which are decoration.
* [Current-magnitude measurements](current_magnitude.md) — why ampere readings
  are genuinely hard, demonstrated rather than asserted.
* [Background and roadmap](state_of_the_art.md) — selected literature, and what
  these prototypes do not yet do.
