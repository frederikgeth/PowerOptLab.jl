# Likelihood, loss, and priors

With measurement rows only, both estimators in this package implement the
**constrained maximum-likelihood estimate for independent Gaussian additive
errors**. They do not use squares merely because squares are convenient. This
statement is conditional on the probability model, however; adding a genuine
state prior changes the statistical interpretation to maximum a posteriori
(MAP), and using the same objective under a misspecified likelihood is better
described as weighted nonlinear least squares or quasi-likelihood.

The distinction is not pedantic. It determines what ``\sigma`` means, what the
objective is allowed to be, what the residuals can be tested against, and which
of the two engines can represent a given error model at all.

## Weighted least squares is a corollary, not a definition

Given measurements ``z_j`` with independent errors distributed as ``f_j``, the
most likely state maximises the likelihood, or equivalently minimises the
negative log-likelihood:

```math
\hat x \;=\; \arg\min_x \; \sum_j -\log f_j\bigl(z_j - h_j(x)\bigr)
\qquad \text{subject to}\qquad c(x) = 0 .
```

Substitute a zero-mean Gaussian, ``f_j = \mathcal N(0, \sigma_j^2)``:

```math
-\log f_j(r_j) = \frac{r_j^2}{2\sigma_j^2} + \log\bigl(\sigma_j\sqrt{2\pi}\bigr)
```

The additive term does not depend on ``x``, so it drops, and what remains is

```math
\min_x \sum_j \frac{\bigl(z_j - h_j(x)\bigr)^2}{\sigma_j^2}
```

— weighted least squares. Thus WLS is the exact MLE under the independent
Gaussian additive-error model above. That does not make WLS invalid whenever
errors are non-Gaussian: under weaker moment assumptions it can still be used as
a generalized/quasi-least-squares estimator, but its finite-sample likelihood,
covariance, and residual-test interpretations no longer follow automatically.
When a Gaussian model *is* asserted, the weight ``1/\sigma_j^2`` is an inverse
variance rather than a preference score.

This is why the [shared assumptions](index.md#What-both-of-them-assume) put the
error model first. Choosing `sigma` *is* choosing a likelihood.

## What each distribution turns into

| Error distribution | ``-\log f`` (up to constants) | Objective form | Classical name |
|---|---|---|---|
| Gaussian ``\mathcal N(0,\sigma^2)`` | ``r^2/2\sigma^2`` | sum of squares | WLS |
| Laplace ``(0,b)`` | ``\lvert r\rvert / b`` | sum of absolute values | WLAV / LAV |
| Huber density ``f\propto e^{-\rho_H}`` | Huber loss ``\rho_H(r)`` | quadratic core, linear tails | Huber M-estimator |
| Student ``t_\nu(0,s)`` | ``\tfrac{\nu+1}{2}\log(1 + r^2/(\nu s^2))`` | smooth, generally nonconvex | heavy-tailed robust MLE |
| Beta, Gamma, Weibull, LogNormal | general ``-\log f(x)`` on bounded or half-line support | smooth on the interior of support | non-Gaussian model for a pseudo-measured quantity |
| Gaussian mixture | ``-\log \sum_i \omega_i \mathcal N_i`` | smooth and generally nonconvex; may be uni- or multimodal | GMM pseudo-measurement |

The first two rows are additive residual models. Positive- or bounded-support
families are usually models for the pseudo-measured quantity itself, not for a
zero-centred additive residual; writing every row as ``f(z-h(x))`` would be the
wrong model for them. Reading the first two rows the other way round: an
``\ell_2`` likelihood fit asserts Gaussian residuals, while an ``\ell_1``
likelihood fit asserts Laplace residuals. Neither is a neutral choice of norm.

## What the contract on distributions actually is

It is tempting to say the only requirement is "a maximum-likelihood estimate
exists". That is close to right for a general nonlinear-programming
formulation, and **not** right for this package's second estimator. The two
engines have genuinely different contracts.

### For an NLP formulation (the JuMP/Ipopt path)

The objective is whatever expression you write, so the requirements are those
of the solver and of the estimator being well posed:

1. **Support.** ``f_j(r_j) > 0`` wherever the solver may step, or ``-\log f_j``
   is ``+\infty`` there. Bounded-support densities (Beta especially) are a real
   hazard: a reading that falls outside the assumed support makes the problem
   infeasible rather than merely unlikely.
2. **Smoothness.** Interior-point methods want ``-\log f`` twice continuously
   differentiable in the residual. Non-smooth likelihoods are admissible only
   through an exact reformulation — see the next section.
3. **Attainment.** The minimum must exist. Likelihoods that are unbounded above
   (mixtures with a free variance collapsing onto a point) have no MLE to find.
4. **Log-concavity is optional but valuable for optimisation.** It makes
   ``-\log f`` convex *in the modelled random variable* and would give a convex
   objective if the measurement map and feasible set were convex. It does not
   imply robustness: the Gaussian is log-concave, yet its linear score gives
   least squares an unbounded influence function.

So: smooth on the represented domain, with a finite attained optimum.
Log-concavity improves objective geometry; bounded influence and resistance to
contamination are separate statistical properties.

### For the Gauss–Newton / Hachtel engine

Much stricter, and this is structural rather than a missing feature.
[`solve_compiled_state_estimator`](@ref) and
[`solve_sparse_state_estimator`](@ref) minimise ``\tfrac12\lVert r(x)\rVert_2^2``.
The Gauss–Newton step, the null-space composite step, the Hachtel augmented
system and the ``(H^\top H)^{-1}`` covariance **all** exploit the sum-of-squares
structure — that is what makes ``J^\top J`` a Hessian approximation at all.

A sum of squared *standardised measurement residuals* is the Gaussian negative
log-likelihood. The implemented engine and its covariance are therefore
Gaussian-only. That is an implementation boundary, not a theorem that
Gauss--Newton ideas can never support robust losses: Huber or Student objectives
can use iteratively reweighted/generalised Gauss--Newton steps, but the weights,
globalisation, convergence tests, and covariance interpretation would all need
to change. Laplace requires a nonsmooth or reformulated step; general density
models require still broader objective support.

That is the sharpest architectural difference between the two estimators, and
it cuts the opposite way to most of the [comparison](comparison.md): the
compiled estimator wins on exact constraints, four-wire modelling, covariance
and time series, while the JuMP path is presently the easier place to experiment
with general error models.

## WLAV exactly, by reformulation

The absolute-value objective is not differentiable at zero, which is usually
where the discussion stops. It should not: ``\ell_1`` has an **exact** smooth
reformulation, so an NLP-based estimator can deliver the Laplace MLE without
approximating anything.

Introduce one epigraph variable per residual:

```math
\min_{x,\,t}\; \sum_j \frac{t_j}{b_j}
\qquad\text{s.t.}\qquad
-t_j \le r_j(x) \le t_j,\quad c(x)=0 .
```

At the optimum ``t_j = \lvert r_j\rvert``, because nothing else pushes ``t_j``
down. The equivalent split-variable form, ``r_j = u_j - v_j`` with
``u_j, v_j \ge 0`` minimising ``\sum (u_j+v_j)/b_j``, is the same trick.

The reformulation adds linear constraints and linear objective terms, while the
sole nonlinear expression remains ``r_j(x)``. It is exact, but it is not
computationally free: it adds one variable and two inequalities per residual,
enlarges the interior-point system, and can change convergence behaviour. That
overhead should be benchmarked rather than inferred from algebra alone. It is
nevertheless a natural next experiment for the JuMP-based estimator.

Huber also has exact auxiliary-variable formulations. More general convex losses
need a solver-compatible representation of their epigraph; piecewise smoothness
alone does not guarantee a finite smooth NLP formulation.

## Priors change the estimand

A likelihood describes data conditional on a state, ``p(z\mid x)``. A state
prior describes ``p(x)``. Combining them gives

```math
\hat x_{\mathrm{MAP}} = \arg\min_x
  \{-\log p(z\mid x)-\log p(x)\}\quad\text{s.t.}\quad c(x)=0,
```

not an MLE. A `StatePrior` row can instead be interpreted as an additional
pseudo-measurement of the state; then it belongs to an augmented likelihood,
but its independence from the meter rows must be defensible.

The previous-snapshot option is deliberately modest. It inserts the previous
point estimate with a user-chosen movement sigma; it does not propagate the
previous covariance, specify a process model, or account for correlation caused
by reusing earlier measurements. It is useful temporal regularisation, not a
Kalman filter or a complete sequential Bayesian estimator. Covariance reported
after adding such rows is conditional on that modelling choice.

## What this package does today

Be clear about the gap between architecture and implementation:

| | error model implemented today | architecture admits |
|---|---|---|
| [`solve_state_estimation`](@ref) (JuMP/Ipopt) | Gaussian only — the hook writes a hard-coded quadratic objective | general solver-compatible ``-\log f`` models; WLAV/Huber through reformulation |
| compiled constrained NLLS | Gaussian only | robust/general likelihoods require changes to the step, globalisation and uncertainty calculation |

Neither estimator currently offers a non-Gaussian option. For the JuMP path
that is a missing feature with a clear route; for the compiled path it is a
property of the method. Both are recorded on the
[roadmap](state_of_the_art.md#Roadmap).

## Related work: Vanin et al. (2023)

A directly relevant implementation in the same JuMP ecosystem is Vanin, Van
Acker, D'hulst and Van Hertem,
*Exact Modeling of Non-Gaussian Measurement Uncertainty in Distribution System
State Estimation*, IEEE Trans. Instrumentation and Measurement, 72:9002911,
2023, implemented in the open-source
[PowerModelsDistributionStateEstimation.jl](https://github.com/Electa-Git/PowerModelsDistributionStateEstimation.jl).

It takes exactly the position argued above — SE is MLE, WLS and WLAV are its
Gaussian and Laplacian special cases — and builds an "optimization-first"
unbalanced DSSE in JuMP that admits any smooth pdf, reporting that the tool
supports Laplacian, Beta, Gamma, Weibull, Gaussian, LogNormal, GMM and
polynomial log-pdf fits. It also corroborates two choices made here
independently: zero injections enter as equality constraints rather than
high-weight pseudo-measurements, and the estimator is posed over an extended
OPF-like variable space rather than as a reduced normal-equation system.

It is the right paper to read next. A few things are worth weighing before
adopting its conclusions wholesale:

* **The headline accuracy gain is narrower than the abstract suggests.** For
  the Beta pseudo-measurements the authors describe as realistic for Belgian
  residential demand (Cases I and II), the paper says the differences between
  the exact model and its Gaussian approximation are "minimal", with the
  approximation showing mainly a few more outliers. The large gain appears in
  Case III, whose polynomial pdf was chosen for being strongly non-Gaussian. So
  the demonstrated result is closer to "approximating a badly non-Gaussian
  density by a Gaussian hurts" than to "exact modelling helps in practice".
* **It assumes the pdf is known, and does not test what happens when it is
  not.** The entire benefit comes from using the correct density, and the
  authors concede that both estimators "work better when the provided pdfs are
  correct and known in advance". In deployment the pseudo-measurement density
  is itself estimated from limited, non-stationary data. The missing experiment
  is sensitivity to pdf misspecification — plausibly the dominant error term,
  and possibly larger than the Gaussian-approximation error being removed.
* **Bounded supports need explicit domain handling.** A Beta density is zero
  outside ``[x^{\min}, x^{\max}]``. The corresponding state or pseudo-measurement
  variable must remain inside that support; inconsistent physics or other data
  can therefore make the model infeasible. A Gaussian approximation has no hard
  support boundary, although that is not by itself evidence that it is the
  better statistical model.
* **The PV correlation constraint is very strong.** Correlated irradiance is
  handled by ``x_i = x_j`` — a *hard equality* asserting identical irradiance
  across the selected users. The numerical study separately assumes equal PV
  capacities, while noting that unequal capacities can scale the common
  irradiance. The defensible criticism is therefore narrower: strong spatial
  correlation is not exact equality, and a joint stochastic model would expose
  rather than suppress imperfect correlation.
* **Nonconvexity compounds.** The formulation keeps exact nonconvex ACR power
  flow and can add objectives that are themselves nonconvex; a GMM can also be
  multimodal. Ipopt returns a local stationary point. The paper reports no
  multistart or uniqueness check, and its own remark that log-concave pdfs
  would permit convex methods *if the power flow were convex* concedes the
  point.
* **It is explicitly not robust estimation.** The authors are careful about
  this: the aim is accuracy under known statistics, not resistance to gross
  errors. Bad data therefore remains unhandled, exactly as it does here.

Their acknowledgement that the fixed reference-angle assignment "might
introduce estimation errors in unbalanced networks", with pointers to virtual
and improved angular references, connects directly to the gauge-freedom
discussion in [observability](observability.md#Four-wire-observability-and-gauge-freedoms).

## References

* Vanin, M., Van Acker, T., D'hulst, R. and Van Hertem, D., "Exact modeling of
  non-Gaussian measurement uncertainty in distribution system state
  estimation", *IEEE Trans. Instrumentation and Measurement*, 72:9002911, 2023.
  DOI: 10.1109/TIM.2023.3287253.
* Vanin, M., Van Acker, T., D'hulst, R. and Van Hertem, D., "A framework for
  constrained static state estimation in unbalanced distribution networks",
  *IEEE Trans. Power Systems*, 37(3):2075–2085, May 2022.
  DOI: 10.1109/TPWRS.2021.3116291.
* Mínguez, R., Conejo, A. J. and Hadi, A. S., "Non Gaussian State Estimation in
  Power Systems", in B. C. Arnold, N. Balakrishnan, J. M. Sarabia and R.
  Mínguez (eds.), *Advances in Mathematical and Statistical Modeling*,
  Birkhäuser Boston, pp. 141–156, 2008. DOI: 10.1007/978-0-8176-4626-4_10.
* Merrill, H. M. and Schweppe, F. C., "Bad data suppression in power system
  static state estimation", *IEEE Trans. Power Apparatus and Systems*,
  PAS-90(6):2718–2725, 1971. DOI: 10.1109/TPAS.1971.292925.
* Zhao, J., Mili, L. and Pires, R. C., "Statistical and numerical robust state
  estimator for heavily loaded power systems", *IEEE Trans. Power Systems*,
  33(6):6904–6914, 2018. DOI: 10.1109/TPWRS.2018.2849325.
* Singh, R., Pal, B. C. and Jabr, R. A., "Distribution system state estimation
  through Gaussian mixture model of the load as pseudo-measurement", *IET
  Generation, Transmission & Distribution*, 4(1):50–59, 2010.
  DOI: 10.1049/iet-gtd.2009.0167.
