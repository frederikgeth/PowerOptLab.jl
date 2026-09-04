# Observability and under-observed estimation

> **Audience:** anyone who has a state estimate and needs to know which parts
> of it are supported by measurements and which parts are decoration.

An estimator always returns something. Observability analysis is what tells
you whether that something means anything. This page covers what PowerOptLab
computes, how to read it, what it deliberately refuses to do, and which parts
of the problem are still open research.

## What "observable" means here

Both estimators report **local numerical observability**: the rank of the
measurement Jacobian at the returned point.

For the [constrained estimator](../problems/constrained_state_estimation.md)
the test is on the *feasible tangent space*. With ``C x = 0`` linearising the
exact equations and ``Z`` a basis for ``\ker C``, the relevant matrix is
``HZ``, not ``H``:

```math
\text{observable} \iff \operatorname{rank}(HZ) = \dim Z .
```

This distinction is not cosmetic. Exact constraints *remove directions the
measurements no longer have to explain*, so a set that cannot identify the
free state can still identify the constrained one. [Below](#Constraints-change-the-question)
shows a two-row set going from 2-of-4 to 2-of-2 by adding one zero-injection
constraint.

Three properties of this test are worth internalising before trusting it:

* **It is local.** It is a statement about one point, computed from a
  linearisation there. It is not a global uniqueness certificate, and no local
  method can supply one.
* **It is numerical.** It asks whether a matrix has full rank in floating
  point, not whether the network graph admits a spanning structure. See
  [numerical versus topological](#Numerical-versus-topological-observability).
* **It is not binary in practice.** Rank is a yes/no, but the smallest
  singular value is a continuum, and it is the number that actually predicts
  whether your estimate is useful.

## Reading the diagnostic

```julia
d = observability_diagnostics(structure, parameters, x)
```

| field | meaning |
|---|---|
| `constraint_rank` | independent exact equations at this point |
| `tangent_dimension` | free directions left after the constraints, ``\dim Z`` |
| `observable_dimension` | ``\operatorname{rank}(HZ)`` — directions the measurements pin down |
| `unobservable_dimension` | the difference; zero means locally identified |
| `min_singular` | ``\sigma_{\min}(HZ)``; the inverse of the worst-case standard deviation |
| `condition_number` | ``\sigma_{\max}/\sigma_{\min}`` |

The WLS estimator's `result.observability` reports the analogous
`(observable, n_states, rank, redundancy, min_singular, cond,
tangent_dimension)`. Exact KCL equations are ranked in their own units, then the
**whitened** measurement Jacobian is ranked on their tangent space. `rank` is
the sum of those two independent ranks; `redundancy` counts measurement rows
beyond the tangent directions they identify. `min_singular` and `cond` come
from the reduced measurement Jacobian, so they have the same interpretation as
the constrained estimator's.

!!! note "Row scaling and what it does and does not affect"
    Exact rank is invariant to any nonzero row scaling; *numerical* rank is not,
    because a tolerance is compared with singular values. An earlier diagnostic
    stacked exact KCL rows in amperes with measurement rows whitened from volts,
    watts and vars. Its singular values were meaningless, and sufficiently
    extreme relative scaling could also change its numerical rank. The current
    implementation ranks the exact block first and the measurement block on its
    tangent space, each in its own scale.

## Redundancy is not observability

The single most common misreading. More measurements than states does not
mean the states are determined — *placement* decides, not count. On a two-bus
feeder (4 states: ``V_{re}``, ``V_{im}`` at buses `b` and `c`):

| measurement set | rows | tangent | observable | unobservable | ``\sigma_{\min}`` |
|---|---:|---:|---:|---:|---:|
| three ``\|V_b\|`` readings | 3 | 4 | **1** | 3 | 0 |
| ``\|V_b\|``, ``\|V_c\|`` | 2 | 4 | 2 | 2 | 0 |
| ``V_{re}``, ``V_{im}`` at both | 4 | 4 | **4** | 0 | 10.0 |

Three readings for four states, and rank one: stacking copies of the same
measurement adds redundancy against *noise* without adding information about
*direction*. Rows two and three have the same count and opposite verdicts.

## What you cannot see

Rank tells you *how much* is missing. [`unobservable_directions`](@ref) tells
you **what**, and this is the part worth using:

```julia
s = compile_state_estimator(net, [vmag_b, vmag_c])
p = SEParameters(s, [vmag_b, vmag_c])
U = unobservable_directions(s, p, x)     # 4x2, columns span the null space
```

With the state ordered ``[V_{re,b}, V_{re,c}, V_{im,b}, V_{im,c}]`` and
``V_b = 211.58 + 7.79j``:

```
dir 1: -0.0368  +0.0000  +0.9993  +0.0000
dir 2: +0.0000  -0.0539  +0.0000  +0.9985
```

Each direction touches one bus, and direction 1 is exactly

```math
\frac{1}{|V_b|}\bigl(-\operatorname{Im} V_b,\; \operatorname{Re} V_b\bigr)
= (-0.0368,\; 0.9993),
```

the tangent to the circle ``|V_b| = \text{const}`` — a **pure phase
perturbation**. The code has rediscovered, numerically and without being told,
the textbook fact that a magnitude-only measurement set leaves the phase
angles free. ``\|Hu\| \approx 3\times10^{-17}``.

This is the diagnostic to reach for when a result looks wrong. "Rank 6 of 8"
is not actionable; "the unobservable direction is the neutral common mode at
these three buses" is.

No basis is stored on the structure — directions are generated on request, so
asking the question costs nothing until you ask it.

## Constraints change the question

Two measurements, four states, and one exact zero injection:

| set | rows | tangent | observable | unobservable |
|---|---:|---:|---:|---:|
| ``V_{re}, V_{im}`` at `b` | 2 | 4 | 2 | 2 |
| the same + `zero_injection=["c"]` | 2 | **2** | 2 | **0** |

The measurements did not change. The constraint removed two directions from
the space they had to explain, and the set became sufficient. Testing
``\operatorname{rank}(H)`` alone would have called both cases unobservable and
been wrong about the second.

This is why the classical literature treats zero injections separately, and
why this implementation generalises the idea: any exact device law — constant
power, constant current, ZIP, a delta connection — enters ``c(x)=0`` and
contributes to ``Z`` in the same way.

## Observability is a continuum

The binary is the least useful thing the diagnostic reports. Here is a set
approaching criticality for a concrete physical reason: a ``|V_b|`` and a
``|I_\ell|`` measurement on a resistive two-bus feeder, as the current angle
relative to the source falls toward zero and the two rows rotate onto each
other ([why](current_magnitude.md#Real-problem-1:-sensitivity-rows-can-collapse-onto-each-other)).

| source-referenced ``\angle I`` (rad) | ``\sigma_{\min}(HZ)`` | ``\mathrm{cond}(HZ)`` | sd of estimated ``\angle V_b`` |
|---:|---:|---:|---:|
| −0.6 | 5.61e+01 | 3.86 | 7.2e−05 rad |
| −0.3 | 2.91e+01 | 7.63 | 1.6e−04 rad |
| −0.1 | 9.78e+00 | 22.8 | 4.8e−04 rad |
| −0.03 | 2.94e+00 | 76.1 | 1.6e−03 rad |
| −0.01 | 9.80e−01 | 228 | 4.9e−03 rad |
| −0.003 | 2.94e−01 | 761 | 1.6e−02 rad |
| 0 | 0 | ``\infty`` | undefined |

Every row except the last reports `observable`, and the current is 40 A
throughout. The rank test passes each of them without comment while the
uncertainty in the estimated angle grows by more than two orders of magnitude.

``\sigma_{\min}`` is the quantity to watch. Under this page's standing
assumptions — independent Gaussian errors, a diagonal covariance, and the model
linearised at the returned point — the estimate covariance on the tangent space
is ``(H_{red}^\top H_{red})^{-1}``, so ``1/\sigma_{\min}`` is the standard
deviation in the worst-determined direction. That is the Cramér–Rao bound *of
the linearised Gaussian problem*, not of the nonlinear one: it is a local
first-order quantity and it inherits every assumption in
[the shared list](index.md#What-both-of-them-assume). Report it next to the
estimate, with those caveats.

## Critical measurements: the bad-data blind spot

A **critical** measurement is one whose removal destroys observability. It has
a property that should worry anyone relying on residuals: its residual is
identically zero, whatever it says. In an exactly determined set, *every*
measurement is critical.

Four rectangular readings on four states, then one of them corrupted by 50σ:

| | status | ``\|r\|`` | estimated ``V_b`` |
|---|---|---:|---|
| clean | `:converged_unique` | 8.9e−15 | 211.5788 + 7.7884j |
| one reading +50σ | `:converged_unique` | **8.9e−15** | **216.5788** + 7.7884j |

The estimate moved by exactly the injected error, 5.0 V, and the residual did
not move at all. Nothing in the output distinguishes these two runs. Add a
single redundant row — ``|V_b|``, one extra measurement — and the same
corruption is unmistakable:

```
per-row whitened residuals:  -24.97,  -0.90,  0.00,  0.00,  +24.98
```

This is the concrete argument for local redundancy, and it is also why
[bad-data processing is the top roadmap item](state_of_the_art.md#Roadmap).
Turning per-row residuals into a test needs the residual **sensitivity** matrix
``S = I - H(H^\top W H)^{-1}H^\top W`` and the residual **covariance**
``\Omega = S R`` (with ``R = W^{-1}`` the measurement covariance); the
normalised residual is ``r^N_i = r_i/\sqrt{\Omega_{ii}}``. A critical
measurement has ``\Omega_{ii} = 0``, which is the formal statement of the row
of zeros above. Neither estimator computes either matrix today.

!!! warning "Redundancy must be local"
    Global redundancy does not protect a locally critical measurement. A
    feeder can have twice as many readings as states and still have a pocket
    where one meter is unchecked.

## Restoring observability, and being honest about it

The standard distribution practice is to fill gaps with pseudo-measurements —
load allocation from billing data, forecasts, nominal profiles. In this
codebase that is a [`StatePrior`](@ref): an ordinary whitened residual row, not
a constraint.

It restores rank. The question is what the answer then means, and the
covariance answers it:

| prior σ | observable | sd(``V_{re,b}``) | sd(``V_{re,c}``) |
|---:|---:|---:|---:|
| — (no prior) | 2 / 4 | *covariance refused* | *refused* |
| 1.0 V | 4 / 4 | 0.1000 V | **1.0000 V** |
| 10.0 V | 4 / 4 | 0.1000 V | **10.0000 V** |
| 100.0 V | 4 / 4 | 0.1000 V | **100.0000 V** |

Bus `b` has a meter, and its standard deviation is the meter's σ (0.1 V)
regardless of the prior. Bus `c` has no meter, and its standard deviation is
*exactly the prior σ*, every time. The covariance is reporting, with no
ambiguity, that the estimate at `c` is 100% assumption and 0% data.

That is the honest way to run an under-observed estimator: restore rank if you
must, then publish the covariance so the reader can see which buses were
measured and which were asserted.

### The estimator refuses to invent a covariance

```julia
selected_state_covariance(s, p, x, [1])
# ArgumentError: requested covariance is not finite:
#   2 tangent directions are unobservable
```

An unobservable direction has infinite variance. Returning a finite matrix
anyway — as a pseudo-inverse silently would — reports confidence that does not
exist. Both [`selected_state_covariance`](@ref) and
[`derived_covariance`](@ref) throw instead.

`derived_covariance` also answers the question users actually ask, which is
rarely about rectangular components. Given ``\partial|V|/\partial x``:

```
|V_b| = 211.7234 V   sd = 0.1016 V
|V_c| = 205.1897 V   sd = 0.1033 V
```

## Four-wire observability and gauge freedoms

Distribution estimation has a reference problem that transmission does not.
Transmission fixes one angle. A four-wire feeder additionally needs a
*potential* reference for the neutral system, and if nothing bonds it to
earth, the neutral voltages are determined only up to a common shift. That is
a **gauge freedom**: a genuine unobservable direction, not a measurement gap,
and no amount of extra telemetry on phase-to-neutral quantities removes it.

Because the constrained estimator keeps ungrounded neutrals as explicit
states, the freedom shows up in the diagnostic instead of being assumed away
by grounding the neutral for convenience:

| network | ``n_{\text{states}}`` | rank | tangent | redundancy | ``\sigma_{\min}`` |
|---|---:|---:|---:|---:|---:|
| bonded source neutral | 8 | **8** | 4 | +2 | 1.35e+00 |
| floating neutral system | 10 | **9** | 4 | +3 | 2.06e−15 |

The second network has *more* redundancy and *less* observability. The missing
direction is the gauge — and `redundancy` counting measurement rows beyond the
directions they identify is exactly why it can rise while observability falls.

!!! warning "Evaluate at a meaningful operating point"
    A local diagnostic depends on where it is evaluated. The same well-posed
    set on the bonded network reports rank 8 at a power-flow point and
    **rank 4 at a flat start**, because a flat profile is itself degenerate —
    every line carries zero current, so every injection row loses its
    sensitivity. Evaluate at a solved state, not at the initial guess, or you
    will diagnose the start rather than the measurement set.

## Numerical versus topological observability

Two traditions answer this question, and they answer slightly different ones.

**Topological** methods (Krumpholz, Clements and Davis, 1980; Monticelli and
Wu, 1985) work on the network graph, looking for a spanning tree that the
measurements can cover. They are combinatorial, immune to round-off, and give
observable islands directly — but they assume *generic* parameter values, so
they can call a network observable when a particular numerical coincidence
makes it not.

**Numerical** methods factor the Jacobian or gain matrix and count rank. They
see the actual numbers, handle mixed measurement types and exact constraints
without special cases, and give ``\sigma_{\min}`` — the continuum above. Their
weakness is that rank in floating point needs a threshold, and with many
injection measurements the disparity of magnitudes in ``H`` makes that
threshold delicate (Expósito and Abur, 1998).

PowerOptLab is entirely numerical: dense SVD with a relative tolerance
(`rank_rtol`, default ``\sqrt{\varepsilon}``), on ``HZ`` for the constrained
estimator and on the measurement-plus-KCL Jacobian for WLS. That choice buys
the constraint handling and the singular values; it costs ``O(n^3)`` and a
threshold. Treat a ``\sigma_{\min}`` near the tolerance as "ask again", not as
an answer.

## Academic questions suggested by the prototype

These are research leads, not novelty claims. A systematic review and sharper
problem statement are required before presenting any of them as open:

1. **Observability of four-wire networks with explicit neutrals.** The
   literature already includes complete four-wire LV state estimators and
   neutral-to-earth-voltage estimation (Bandara et al., 2021; Melo, Teixeira and
   Mingorança, 2023). The narrower question exposed here is how bonding,
   multi-grounding, broken neutrals, and the chosen measurement references alter
   the gauge/null-space structure, and how that relates to existing neutral
   reductions and observability results.
2. **Observability under general exact device equations.** Classical treatments
   special-case zero injections. The ``HZ`` formulation admits any exact law —
   ZIP, delta connections, an inverter control characteristic. A useful study
   would classify when these laws contribute usable rank and compare the result
   with established equality-constrained and augmented-state formulations. A
   constant-power device is a nonlinear constraint whose contribution to ``Z``
   depends on the operating point, so local observability can change with
   loading.
3. **Quantitative observability for meter placement.** With
   `derived_covariance` the design objective is available directly: maximise
   ``\sigma_{\min}`` (E-optimality) or minimise the variance of a quantity the
   operator actually cares about, such as the worst phase-to-neutral voltage.
   The continuum table above is the argument that binary placement criteria
   are the wrong objective.
4. **Critical measurements and bad data in the constrained setting.** Criticality
   interacts with exact constraints: a measurement can be critical on the
   tangent space while looking redundant in the raw count. Characterising
   critical ``k``-tuples of ``HZ`` — and the corresponding "critical
   constraints", whose multipliers the sparse solver already exposes — is a
   clean extension of Monticelli and Garcia's framework.
5. **Local versus global identifiability with magnitude measurements.** The
   [mirror solutions](current_magnitude.md) are locally unique and globally
   not. The research question is which existing tools from polynomial-system
   solving, interval analysis, relaxations, or certified continuation can give
   useful uniqueness certificates for bounded distribution-network models.
6. **Observable-island decomposition from the null space.** The unobservable
   directions are computed; turning them into maximal observable islands and a
   minimal set of pseudo-measurements that merges them is a combinatorial step
   this codebase does not yet take.

## What to do with an under-observed estimate

1. Read `unobservable_dimension` and `min_singular` **before** the state.
2. Call [`unobservable_directions`](@ref) and interpret them physically. They
   are usually a gauge freedom, an unmetered pocket, or a magnitude-only
   phase ambiguity, and the three have different fixes.
3. Prefer a constraint to a pseudo-measurement when the physics is genuinely
   exact — it costs no rank and adds no fictitious information.
4. If you add a prior, publish the covariance alongside the estimate.
5. Never report a state whose covariance the estimator refused to compute
   without saying which directions were unsupported.
6. Evaluate the diagnostic at a solved operating point.

## References

* Krumpholz, G. R., Clements, K. A. and Davis, P. W., "Power system
  observability: a practical algorithm using network topology", *IEEE Trans.
  Power Apparatus and Systems*, PAS-99(4), 1980.
* Monticelli, A. and Wu, F. F., "Network observability: theory", *IEEE Trans.
  Power Apparatus and Systems*, PAS-104(5), 1985.
* Clements, K. A., Krumpholz, G. R. and Davis, P. W., "Power system state
  estimation with measurement deficiency: an observability/measurement
  placement algorithm", *IEEE Trans. Power Apparatus and Systems*, PAS-102(7),
  1983.
* Expósito, A. G. and Abur, A., "Generalized observability analysis and
  measurement classification", *IEEE Trans. Power Systems*, 13(3), 1998.
* Monticelli, A. and Garcia, A., "Reliable bad data processing for real-time
  state estimation", *IEEE Trans. Power Apparatus and Systems*, PAS-102(5),
  1983.
* Abur, A. and Expósito, A. G., *Power System State Estimation: Theory and
  Implementation*, Marcel Dekker, 2004, chapters on observability and bad data.
* Bandara, W. G. C., Almeida, D., Godaliyadda, R. I., Ekanayake, M. P. and
  Ekanayake, J., "A complete state estimation algorithm for a three-phase
  four-wire low voltage distribution system with high penetration of solar PV",
  *International Journal of Electrical Power & Energy Systems*, 124:106332,
  2021. DOI: 10.1016/j.ijepes.2020.106332.
* Melo, I. D., Teixeira, M. O. N. and Mingorança, J. S., "Neutral-to-Earth
  Voltage (NEV) and state estimation for unbalanced multiphase distribution
  systems based on an optimization model", *Electric Power Systems Research*,
  217:109123, 2023. DOI: 10.1016/j.epsr.2023.109123.
