# State of the art and roadmap

What the literature settles, what these prototypes implement, and what they do
not. This page is deliberately blunt about the gaps.

## What the field settled long ago

**The formulation.** Schweppe and Wildes (1970) posed static state estimation
as weighted least squares over a nonlinear measurement model, and that is
still the frame. Both PowerOptLab estimators are instances of it.

**Zero injections are constraints, not measurements.** Modelling a known zero
injection as a pseudo-measurement with a tiny σ is the classic source of
ill-conditioning: it inflates the gain matrix's condition number by the square
of the weight ratio. Treating them as equality constraints and solving the
resulting KKT system is the standard cure. The
[constrained NLLS estimator](../problems/constrained_state_estimation.md) does
this natively; the [WLS estimator](../problems/state_estimation.md) gets the
same effect from the engine's KCL rather than from a weighted row.

**The augmented (Hachtel) system beats the normal equations.** Gjelsvik, Aam
and Holten (1985) showed that solving the augmented system instead of forming
``H^\top W H`` greatly improves numerical stability on ill-conditioned
networks, at some extra storage. `solve_sparse_state_estimator` uses exactly
this: it factors the blocked ``[I\;H\;0;\;H^\top\;\gamma I\;C^\top;\;0\;C\;0]``
system with sparse QR and never forms ``H^\top H``.

**Bad data needs the residual covariance, not the raw residual.** Monticelli
and Garcia (1983) established largest-normalised-residual processing:
``r^N_i = r_i/\sqrt{S_{ii}}`` with ``S`` the residual sensitivity matrix.
Dividing by σ alone is not this, and does not have the χ² distribution the
test needs. **Neither estimator implements it** — see the gaps below.

**Ampere measurements are a modelling hazard, not a conditioning one.** Abur
and Expósito (1997) showed that current magnitudes admit multiple solutions.
The [current-magnitude tutorial](current_magnitude.md) reproduces this and
also shows the conditioning story is a myth.

**Distribution estimation is measurement-poor.** Modern DSSE surveys
(Dehghanpour et al., 2019; Primadianto and Lu, 2017) agree that the binding
constraint is telemetry, not algorithms: real feeders are observable only via
pseudo-measurements from billing or load allocation, whose errors are large,
non-Gaussian, and correlated. That last word is the one that matters here —
both estimators assume a diagonal covariance.

## Where these prototypes stand

| Capability | WLS | Constrained NLLS | Notes |
|---|:--:|:--:|---|
| Nonlinear WLS fit | ✅ | ✅ | |
| Exact zero injection | ✅ | ✅ | constraints, not weighted rows |
| Exact device laws | ❌ | ✅ | constant P, constant I, ZIP |
| Four-wire / floating neutral | ❌ | ✅ | |
| Branch telemetry | ❌ | ✅ | lines only; no transformers |
| Local observability | ✅ | ✅ | dense SVD, both |
| State covariance | ❌ | ✅ | refuses rank-deficient rather than inventing |
| Time series / priors | ❌ | ✅ | |
| Hachtel sparse step | ❌ | ✅ | |
| χ² detection | ❌ | ❌ | |
| Largest-``r^N`` identification | ❌ | ❌ | |
| Robust estimator (LAV / Huber) | ❌ | ❌ | |
| Correlated covariance | ❌ | ❌ | |
| Topology-error hypotheses | ❌ | ❌ | |
| Multistart / global search | ❌ | ❌ | |
| PMU / synchrophasor model | ❌ | partial | `:vr`/`:vi` exist; no reference-angle or PMU error model |

## Roadmap

Ordered by how much each would change what a user can defensibly claim.

### 1. Bad-data processing (largest gap)

Neither estimator can tell a bad reading from a bad model. This needs the
residual sensitivity ``S = I - H(H^\top W H)^{-1}H^\top W``, the χ² test on
the objective, and normalised-residual identification. The constrained
estimator already has the pieces — `_derived_covariance` builds
``(H_{red}^\top H_{red})^{-1}`` on the tangent space — so ``S`` is reachable
without new linear algebra. Constraint multipliers are exposed today as
*leads* only, and the documentation is careful to say so.

### 2. Genuinely sparse Jacobian assembly

`residual_jacobian` and `constraint_jacobian` return **dense** matrices, and
build them with `Vector(E[i, :])` per row and `M[i, :]` slices of a CSC matrix
(which scan every column). At 1600 states the residual Jacobian is already a
10 MB dense array. The sparse *step* is sparse; the assembly around it is not,
so `solve_sparse_state_estimator` does not yet scale to a real feeder. This is
a contained, well-understood refactor: assemble triplets directly.

The dense SVD rank diagnostic used to run every iteration and dominated the
solve; it now runs only on exit paths, but it is still ``O(n^3)`` and should
become a sparse rank-revealing factorisation before large cases are attempted.

### 3. Robust estimation

Least squares has an unbounded influence function: one bad reading moves the
estimate without limit. LAV (Abur and Çelik) or a Huber M-estimator would give
the prototypes a defensible answer under the non-Gaussian pseudo-measurement
errors that DSSE actually faces. This changes the objective, not the
structure, so the compiled evaluator can be reused.

### 4. Correlated covariance

Both estimators whiten with a diagonal ``1/\sigma``. Pseudo-measurements
derived from a common load allocation are strongly correlated, and treating
them as independent overstates the information content — the estimate looks
more certain than it is. A Cholesky factor of a block covariance would slot
into the existing whitening step.

### 5. Multistart

Both problems are nonconvex, and the current-magnitude tutorial exhibits two
exactly-fitting solutions. A cheap multistart with a spread comparison would
turn "converged" into something closer to a claim.

### 6. Observable-island decomposition

The unobservable directions are computed and can be inspected
([observability](observability.md)); turning them into maximal observable
islands, and into a minimal set of pseudo-measurements that merges them, is a
combinatorial step this codebase does not take.

### 7. Transformer and source telemetry

`BranchMeasurement` covers lines via `line_yprim`. Transformers, switches,
neutral currents and source phasors are not modelled. Note also that
`:pflow`/`:qflow` are referenced to **earth** at the measured conductor, with
no per-measurement `reference`, unlike node `Measurement`s — a real meter
reports phase-to-neutral power, so on a four-wire line with a displaced
neutral these are not the same quantity.

## References

* Schweppe, F. C. and Wildes, J., "Power system static-state estimation, Part
  I: Exact model", *IEEE Trans. Power Apparatus and Systems*, PAS-89(1),
  1970.
* Gjelsvik, A., Aam, S. and Holten, L., "Hachtel's augmented matrix method — a
  rapid method improving numerical stability in power system static state
  estimation", *IEEE Trans. Power Apparatus and Systems*, PAS-104(11):2987–2993,
  November 1985.
* Monticelli, A. and Garcia, A., "Reliable bad data processing for real-time
  state estimation", *IEEE Trans. Power Apparatus and Systems*, PAS-102(5),
  1983.
* Abur, A. and Expósito, A. G., "Detecting multiple solutions in state
  estimation in the presence of current magnitude measurements", *IEEE Trans.
  Power Systems*, 12(1):370–375, February 1997.
* Abur, A. and Expósito, A. G., *Power System State Estimation: Theory and
  Implementation*, Marcel Dekker, 2004.
* Monticelli, A., "Electric power system state estimation", *Proceedings of
  the IEEE*, 88(2):262–282, 2000.
* Baran, M. E. and Kelley, A. W., "State estimation for real-time monitoring
  of distribution systems", *IEEE Trans. Power Systems*, 9(3), 1994.
* Primadianto, A. and Lu, C.-N., "A review on distribution system state
  estimation", *IEEE Trans. Power Systems*, 32(5):3875–3883, 2017.
* Dehghanpour, K., Wang, Z., Wang, J., Yuan, Y. and Bu, F., "A survey on state
  estimation techniques and challenges in smart distribution systems", *IEEE
  Trans. Smart Grid*, 10(2):2312–2322, 2019.
