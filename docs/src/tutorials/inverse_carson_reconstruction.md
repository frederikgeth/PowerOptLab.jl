# Inverse Carson reconstruction: a modelling tutorial

> **Audience:** power-system researchers · **Scope:** reconstructing plausible
> overhead-line constructions from diagonal sequence impedance data.

Inverse Carson is an inverse *construction* problem, not a magic conversion
from sequence impedance to a unique geometry. Diagonal zero- and
positive-sequence data discard coupling information; several conductor catalogs,
spacings, and neutral arrangements can fit the same reported values. This
tutorial shows how to frame the problem honestly, fit a bounded catalog, assess
ambiguity, and materialize models without erasing their assumptions.

## 1. Start from the observation model

[`SequenceLineObservation`](@ref) represents the reported diagonal sequence
data. Series observations are `(R0, X0, R1, X1)`; optional shunt observations
are `(B0, B1)`. State the units, frequency, soil resistivity, earth-return
approximation, transposition convention, neutral treatment, and sequence
convention before fitting anything.

```julia
using PowerOptLab

obs = SequenceLineObservation(
    z0 = 0.5952 + 1.5873im,  # Ω/km
    z1 = 0.4472 + 0.3692im,  # Ω/km
    frequency = 50.0,
    sigma = [0.006, 0.016, 0.0045, 0.0037],
)
```

`sigma` is in the declared input units and defines both residual weighting and
candidate compatibility. If measurement errors are correlated, provide a
positive-definite full `covariance` in the ordering `(R0, X0, R1, X1[, B0, B1])`.
The objective then uses a Mahalanobis residual; compatibility remains reported
as readable marginal standardized residuals.

### Pitfall: interpreting published rounding as metrology

When a paper reports only rounded sequence values, its last digit is not a
calibrated measurement uncertainty. A nominal one-percent `sigma` is useful for
descriptive screening, but it does not give likelihood-ratio intervals a
frequentist interpretation. Include catalog tolerance, earth-model mismatch,
and rounding/source uncertainty in the uncertainty model—or explicitly label
the result as a tolerance study.

## 2. Represent a construction catalog, not an unconstrained line

The method enumerates discrete conductor/geometry families outside the NLP and
optimizes only bounded continuous geometry and temperature variables inside each
candidate. For example:

```julia
candidate = OverheadCarsonCandidate(
    id = "mars-horizontal-3",
    geometry = :horizontal_3,
    r_ac_ref = fill(4.5e-4, 3),  # Ω/m at 20 °C
    gmr = fill(4.0e-3, 3),
    radius = fill(5.0e-3, 3),
    lower = [0.30, 6.0, 20.0],   # half-span [m], height [m], temperature [°C]
    upper = [0.80, 12.0, 90.0],
)
```

The supported canonical families are `:horizontal_3`, `:triangle_3`,
`:horizontal_4`, and `:neutral_under_4`. Physical bounds matter: they keep
conductor separation positive, prevent a nonlinear solver from exploring
meaningless geometry, and encode the actual construction domain under study.

### Pitfall: treating bounds as harmless optimizer settings

Candidate bounds are prior physical knowledge. A fit at a span, height, or
temperature bound says the data would prefer to leave the candidate domain; it
does not validate that construction. Report bound-active fits and test whether
a different catalog family explains the observation more credibly.

## 3. Fit every credible candidate and retain ambiguity

```julia
result = solve_inverse_carson(obs, [candidate_a, candidate_b];
                              starts=16, acceptance_sigma=3.0)

result.compatible_candidates
result.fits[1].parameters
result.fits[1].standardized_residual
```

The solver uses deterministic multistart and bounded smooth weighted least
squares. `acceptance_sigma` is a componentwise compatibility rule: every
standardized residual must lie within the threshold. It is intentionally not a
selection rule that collapses all plausible candidates to a single winner.

### Pitfall: reporting only the minimum objective candidate

The numerically best candidate can be practically indistinguishable from other
compatible constructions. Choosing it alone turns catalog ambiguity into false
precision. Report all compatible ids, their residual vectors, parameter bounds,
and local rank diagnostics. If a downstream study needs one model, describe the
external engineering criterion used to select it.

## 4. Know what sequence data discard

The forward model constructs a primitive phase/neutral impedance and capacitance
matrix, applies the stated neutral-elimination convention for comparison, then
forms the observed diagonal sequence quantities. The returned fit includes a
full reconstructed `Z_sequence` (and `B_sequence` when applicable), including
off-diagonal sequence coupling absent from the observation.

Those off-diagonal entries are *model-derived completions*, not new
measurements. Likewise, for four-wire candidates a neutral may be Kron-reduced
only for the sequence comparison; the materialized operational model keeps the
neutral explicit.

### Pitfall: confusing a data transformation with an operational neutral model

Imposing zero neutral voltage to form a Schur complement is an observation/data
convention. It does not say that a downstream LV feeder has a perfectly grounded
neutral. Preserve the primitive four-wire representation for unbalanced power
flow and state estimation, then use the appropriate neutral termination in
those operational models.

## 5. Assess local identifiability before quoting error bars

Each [`InverseCarsonFit`](@ref) returns singular values and rank of the
ForwardDiff prediction Jacobian at the fitted point. A rank-deficient Jacobian
means one or more local parameter directions are not identified. In that case
the local covariance and normal confidence intervals are deliberately absent:

```julia
fit = result.fits[1]
fit.jacobian_rank
fit.jacobian_singular_values
fit.local_parameter_covariance
```

For series impedance alone, vertically translating all conductors does not
change the series impedance, so absolute height is expected to be weak or
unidentified. Shunt data can help, but height sensitivity to shunt uncertainty
can make the practical interval broad.

### Pitfall: using a pseudoinverse covariance as evidence of certainty

A pseudoinverse can assign a convenient-looking number to an unidentified
direction. This implementation instead refuses local covariance when rank is
deficient. Treat that as an inference result: additional data, a narrower
catalog, or a different question is required.

## 6. Profile nonlinear and bound-limited directions

Local covariance is only a linear approximation around one fit. Use profile
likelihood to fix one parameter, reoptimize the others, and trace the connected
region below the one-degree-of-freedom threshold:

```julia
profile = profile_inverse_carson(fit, candidate, obs;
                                 parameter=:height,
                                 confidence_level=0.95)

profile.lower_status
profile.upper_status
```

Profile endpoints distinguish a threshold crossing, a candidate bound, and a
solver-limited endpoint. Continuation warm-starts neighboring profile points and
retries failures from deterministic fitted/candidate starts.

### Pitfall: calling a profile interval global

The profile follows the connected region reached from the selected fit. It can
miss disconnected feasible regions or another candidate family. Use
deterministic multistart, broader catalog enumeration, or a global method before
claiming global confidence regions.

## 7. Keep earth-return and forward-model mismatch visible

The inverse fit uses power-frequency modified Carson. OpenDSS commonly uses the
Deri complex-depth approximation; they are related but not interchangeable.
The included OpenDSS benchmark is intentionally rejected under tight uncertainty
when the earth-return mismatch is visible, rather than inventing a geometry to
absorb it.

### Pitfall: calibrating geometry against the wrong forward model

If an observation was produced with a different earth-return model, frequency,
soil resistivity, transposition, conductor resistance convention, or neutral
elimination than the inverse model, a low residual can be physically misleading
and a high residual can be the correct result. Treat forward-model provenance as
part of the data, not a hidden software detail.

## 8. Materialize a candidate without mutating the study network

After selecting or retaining a fit, create BMOPF-ready construction blocks:

```julia
construction = materialize_inverse_carson(fit, candidate)
construction["wire_data"]
construction["line_geometry"]
```

Materialization does not mutate a network. It permits a transparent next step:
build one operational network per compatible construction and propagate the
structural ambiguity into power flow, DSSE, or DOE scenarios.

### Pitfall: converting ambiguity into one deterministic linecode too early

Operational studies often care about voltage/current sensitivity that was not
observed in diagonal sequence data. Use all compatible constructions as model
scenarios where the decision is sensitive. A single selected linecode is a
decision, not an inverse-problem discovery.

## 9. Minimum research-report checklist

For every study, report:

1. input values, units, frequency, covariance/sigma source, and rounding model;
2. earth-return formulation, soil resistivity, transposition, and sequence
   convention;
3. neutral-elimination and shunt observation conventions;
4. conductor catalog, geometry family, parameter bounds, temperature/resistance
   convention, multistart count, and compatibility threshold;
5. all compatible candidates and their standardized residuals, not only the
   nominal winner;
6. Jacobian rank/singular values, local covariance availability, and profile
   endpoint statuses;
7. treatment of forward-model mismatch and every downstream use of the
   materialized primitive model.

With this record, inverse Carson reconstruction is a reproducible constrained
inference exercise rather than an overconfident reverse-engineering claim. See
[Inverse Carson reconstruction](../problems/inverse_carson.md) for the full API
and benchmark conventions.
