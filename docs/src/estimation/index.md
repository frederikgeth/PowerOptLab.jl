# State estimation

> **Kind:** Problem specification · **Direction:** inverse · **Maturity:** prototype

For calibration or line-construction questions, begin with the
[estimation overview](studies.md). This page compares the two
state estimators; their functions and arguments are listed in the
[state estimation API](../api/state_estimation.md).

State estimation determines bus voltages from measurements and an assumed
network model. PowerOptLab provides a JuMP/Ipopt WLS estimator and a constrained
nonlinear least-squares estimator. Their treatment of exact device equations,
measurements and numerical solution differs.

| | [Legacy WLS](../problems/state_estimation.md) | [Constrained NLLS](../problems/constrained_state_estimation.md) |
|---|---|---|
| Entry point | `solve_state_estimation` | `compile_state_estimator` + a solver |
| Statement | ``\min \sum_i (z_i - h_i(x))^2/\sigma_i^2`` | ``\min \tfrac12\lVert r(x)\rVert^2`` s.t. ``c(x)=0`` |
| Engine | JuMP + Ipopt | own Gauss–Newton / Hachtel solvers |

[Choosing a formulation](comparison.md) compares the two estimators on the
same network and measurements.

## What both of them assume

Both methods rely on the following assumptions. Errors in those assumptions
can bias the estimate even when the numerical solver converges:

* **Measurement weights are diagonal.** Both estimators scale residuals by
  ``1/\sigma`` and do not support correlated measurement errors. Independent,
  zero-mean Gaussian errors give the weighted least-squares objective its
  maximum-likelihood interpretation. Least-squares fitting can also be used
  under other error distributions, but the likelihood, covariance and residual
  tests then need separate justification. A systematic meter offset requires
  analysis of bias as well as random noise. See
  [likelihood, loss, and priors](maximum_likelihood.md).
* **Topology and line parameters are exact.** Neither estimator hypothesises a
  different switch status or a different impedance. A topology error can
  produce large residuals or be partly absorbed into the estimated voltages. See
  [parameter estimation](../problems/parameter_estimation.md) and
  [inverse Carson](../problems/inverse_carson.md) for the complementary
  problems.
* **The estimate is a local stationary point.** Both problems are nonconvex.
  Convergence is not uniqueness, and neither estimator does multistart. The
  [current-magnitude tutorial](current_magnitude.md) shows two different
  states that both fit the same data exactly, each reported as converged.
* **Observability requires rank; precision also depends on conditioning.**
  Measurement placement determines which state directions are identified. Small
  singular values indicate directions that are weakly supported by the data. See
  [observability](observability.md).
* **Standard deviations determine the relative measurement weights.** Assign
  uncertainties according to the evidence for each reading or forecast. Giving
  a load forecast the precision of a calibrated meter can overstate the
  information it contributes.

## Where to go next

* [Worked state-estimation study](../tutorials/constrained_nlls_state_estimation.md) —
  estimate a displaced neutral, inspect residuals and remove measurements to test
  observability. The complete code and assertions run in CI.
* [Choosing a formulation](comparison.md) — the comparison, and which to reach for.
* [Likelihood, loss, and priors](maximum_likelihood.md) — when these estimators
  coincide with Gaussian MLE, when a prior makes the result MAP, and which error
  models each engine can represent.
* [Legacy WLS state estimation](../problems/state_estimation.md) — reference.
* [Constrained NLLS state estimation](../problems/constrained_state_estimation.md) — reference.
* [Observability and under-observed estimation](observability.md) — which
  state directions are identifiable and how to report partial observability.
* [Current-magnitude measurements](current_magnitude.md) — how distinct states
  can fit the same current magnitudes.
* [Background and open questions](state_of_the_art.md) — selected literature, and what
  these prototypes do not yet do.
