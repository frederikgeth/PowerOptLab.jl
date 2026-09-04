# Estimators, not curve fits

Both estimators in this package are **maximum-likelihood estimators**. Neither
is an ``\ell_p`` fitting engine that happens to use squares because squares are
convenient.

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

— weighted least squares. **WLS is the maximum-likelihood estimator for
Gaussian errors, and only for Gaussian errors.** The weight ``1/\sigma_j^2`` is
not a tuning knob expressing how much you like a reading; it is the inverse
variance of an error distribution you are asserting.

This is why the [shared assumptions](index.md#What-both-of-them-assume) put the
error model first. Choosing `sigma` *is* choosing a likelihood.

## What each distribution turns into

| Error distribution | ``-\log f`` (up to constants) | Objective form | Classical name |
|---|---|---|---|
| Gaussian ``\mathcal N(0,\sigma^2)`` | ``r^2/2\sigma^2`` | sum of squares | WLS |
| Laplace ``(0,b)`` | ``\lvert r\rvert / b`` | sum of absolute values | WLAV / LAV |
| Huber | quadratic core, linear tails | piecewise smooth | Huber M-estimator |
| Student ``t_\nu`` | ``\tfrac{\nu+1}{2}\log(1 + r^2/\nu)`` | smooth, nonconvex | heavy-tailed robust MLE |
| Beta, Gamma, Weibull, LogNormal | general ``-\log f`` on a bounded or half-line support | smooth on support | non-Gaussian pseudo-measurement |
| Gaussian mixture | ``-\log \sum_i \omega_i \mathcal N_i`` | smooth, nonconvex, multimodal | GMM pseudo-measurement |

Reading the table the other way round: an ``\ell_2`` fit *assumes* Gaussian, an
``\ell_1`` fit *assumes* Laplace. Neither is a neutral choice of norm.

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
4. **Log-concavity is optional but valuable.** It makes ``-\log f`` convex *in
   the residual*, which bounds the damage a single reading can do and would
   give convexity outright if the power-flow model were linear. It is not,
   so the composite problem stays nonconvex either way.

So: smooth, positive on the search region, attained. Log-concavity buys
well-behavedness, not admissibility.

### For the Gauss–Newton / Hachtel engine

Much stricter, and this is structural rather than a missing feature.
[`solve_compiled_state_estimator`](@ref) and
[`solve_sparse_state_estimator`](@ref) minimise ``\tfrac12\lVert r(x)\rVert_2^2``.
The Gauss–Newton step, the null-space composite step, the Hachtel augmented
system and the ``(H^\top H)^{-1}`` covariance **all** exploit the sum-of-squares
structure — that is what makes ``J^\top J`` a Hessian approximation at all.

A sum of squares is the Gaussian negative log-likelihood. So this engine is
Gaussian-only by construction. Laplace, Huber, Beta and mixtures are not
"unimplemented" here; they are not expressible without changing the algorithm
(to IRLS, to a general NLP, or to a smooth surrogate that is no longer exact).

That is the sharpest architectural difference between the two estimators, and
it cuts the opposite way to most of the [comparison](comparison.md): the
compiled estimator wins on exact constraints, four-wire modelling, covariance
and time series, and loses decisively on error models.

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

Both add linear constraints and linear objective terms only. The sole
nonlinearity is ``r_j(x)``, which was already there. So **WLAV costs
essentially nothing over WLS in an NLP formulation** — a genuine strength of
the modelling-language approach, and the reason robust estimation is a natural
next step for the JuMP-based estimator rather than for the compiled one.

Huber admits the same treatment with one extra split, and any piecewise-smooth
convex ``-\log f`` can be handled by its epigraph.

## What this package does today

Be clear about the gap between architecture and implementation:

| | error model implemented today | architecture admits |
|---|---|---|
| [`solve_state_estimation`](@ref) (JuMP/Ipopt) | Gaussian only — the hook writes a hard-coded quadratic objective | any smooth ``-\log f``, plus WLAV/Huber by exact reformulation |
| compiled constrained NLLS | Gaussian only | Gaussian only, structurally |

Neither estimator currently offers a non-Gaussian option. For the JuMP path
that is a missing feature with a clear route; for the compiled path it is a
property of the method. Both are recorded on the
[roadmap](state_of_the_art.md#roadmap).

## Related work: Vanin et al. (2023)

The closest published treatment is Vanin, Van Acker, D'hulst and Van Hertem,
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
* **Bounded supports are operationally risky.** A Beta density is zero outside
  ``[x^{\min}, x^{\max}]``, so a reading outside that window does not merely
  become unlikely: the log-likelihood is undefined and the problem infeasible.
  A Gaussian approximation degrades gracefully where the exact model does not.
* **The PV correlation constraint is very strong.** Correlated irradiance is
  handled by ``x_i = x_j`` — a *hard equality* asserting identical irradiance
  and identical installed capacity across unmetered users. That is the same
  trap the [comparison page](comparison.md) warns about from the other
  direction: strong correlation is not exact equality, and a correlated
  covariance is the more defensible instrument.
* **Nonconvexity compounds.** The formulation keeps exact nonconvex ACR power
  flow and adds objectives that are themselves nonconvex and, for a GMM,
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
* Vanin, M., Van Acker, T., D'hulst, R. and Van Hertem, D., "A framework for
  constrained static state estimation in unbalanced distribution networks",
  *IEEE Trans. Power Systems*, 37(3):2075–2085, May 2022.
* Mínguez, R., Conejo, A. J. and Hadi, A. S., *Non-Gaussian State Estimation in
  Power Systems*, Birkhäuser, 2008.
* Merrill, H. M. and Schweppe, F. C., "Bad data suppression in power system
  static state estimation", *IEEE Trans. Power Apparatus and Systems*,
  PAS-90(6):2718–2725, 1971.
* Zhao, J., Mili, L. and Pires, R. C., "Statistical and numerical robust state
  estimator for heavily loaded power systems", *IEEE Trans. Power Systems*,
  33(6):6904–6914, 2018.
* Singh, R., Pal, B. C. and Jabr, R. A., "Distribution system state estimation
  through Gaussian mixture model of the load as pseudo-measurement", *IET
  Generation, Transmission & Distribution*, 4(1):50–59, 2010.
