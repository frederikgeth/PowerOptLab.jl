# Inverse Carson reconstruction

> **Kind:** Problem specification · **Maturity:** prototype · **Direction:** inverse · **Scope:** overhead lines

[`solve_inverse_carson`](@ref) reconstructs plausible overhead-line construction
data from diagonal zero- and positive-sequence impedance, with optional shunt
susceptance. It follows Tam, Geth & Mithulananthan, *The Inverse Carson's
Equations Problem: Definition, Implementation and Numerical Experiments*
([arXiv:2404.08210](https://arxiv.org/abs/2404.08210)), but deliberately reports
all compatible candidates instead of forcing a unique construction.

## What is—and is not—observed

Diagonal sequence data discard the off-diagonal sequence coupling of an
untransposed line. Many physical constructions can therefore produce the same
four series values `(R0, X0, R1, X1)`. The result is a physically plausible
completion of missing construction data, not a uniquely observed primitive
matrix.

The prototype uses a catalog-first workflow:

1. Supply discrete conductor and geometry candidates.
2. For each candidate, optimize bounded geometry and temperature variables.
3. Rank by standardized measurement residual.
4. Retain every candidate within the stated uncertainty.
5. Report the singular values and rank of the local prediction Jacobian.

The forward primitive `Z` and `C` matrices come from
`BMOPFTools.overhead_line_constants`. A neutral is Kron-reduced only for the
sequence-data comparison; the returned and materialized model keeps that neutral
explicit for downstream power flow and state estimation.

## Example

```julia
using PowerOptLab

mars_horizontal = OverheadCarsonCandidate(
    id = "mars-horizontal-3",
    geometry = :horizontal_3,
    r_ac_ref = fill(4.5e-4, 3),       # Ω/m at 20 °C
    gmr = fill(4.0e-3, 3),            # m
    radius = fill(5.0e-3, 3),         # m
    lower = [0.30, 6.0, 20.0],        # half-span [m], height [m], temp [°C]
    upper = [0.80, 12.0, 90.0],
)

obs = SequenceLineObservation(
    z0 = 0.5952 + 1.5873im,           # Ω/km
    z1 = 0.4472 + 0.3692im,           # Ω/km
    frequency = 50.0,
    sigma = [0.006, 0.016, 0.0045, 0.0037],
)

result = solve_inverse_carson(obs, [mars_horizontal]; starts=16)
result.compatible_candidates
result.fits[1].parameters
result.fits[1].jacobian_singular_values

# Produce BMOPF-ready wire_data/line_geometry blocks without mutating a network.
construction = materialize_inverse_carson(result.fits[1], mars_horizontal)
```

## Interpreting identifiability

With modified Carson series data alone, translating every overhead conductor
vertically does not change the series impedance. Absolute height therefore has
zero sensitivity and should be expected to be rank-deficient. Reliable shunt
data can constrain height, but the paper shows it is highly sensitive to small
susceptance errors.

Similarly, three-wire horizontal and triangular families can match the same
diagonal sequence series components by changing their span. In that case both
candidates appear in `compatible_candidates`; this is the intended result.

Each fit also returns the full reconstructed `Z_sequence`, including the
off-diagonal sequence coupling absent from the observation, plus `B_sequence`
when shunt data were supplied. These entries are model-derived completions and
must not be interpreted as separately measured quantities.

`sigma` controls both weighting and compatibility. When it is omitted, the
constructor supplies a descriptive 1% tolerance with an absolute floor in the
declared input units. Scores should only be treated as statistical quantities
when `sigma` genuinely represents measurement and source-model uncertainty.

A full `covariance` matrix may be supplied instead of `sigma`, in the ordering
`(R0, X0, R1, X1[, B0, B1])` and in the declared input units. The objective then
uses the Mahalanobis residual obtained from a Cholesky factor. Candidate
compatibility still uses marginal standardized residuals, so it remains readable
and does not depend on the arbitrary ordering used by a triangular whitening
factor. Covariance must be symmetric positive definite; singular covariance is
not silently regularized.

## Confidence and profile intervals

For a locally full-rank fit, `local_parameter_covariance` is the inverse Fisher
approximation `(J'J)^-1`, mapped back from normalized variables to physical
units. It is evaluated from the Jacobian SVD rather than by forming and inverting
the more ill-conditioned normal matrix. `local_confidence_intervals` contains
its two-sided normal intervals at `confidence_level`. Both fields are `nothing`
when the Jacobian is rank
deficient: using a pseudoinverse there would misleadingly assign zero variance
to an unidentified direction.

[`profile_inverse_carson`](@ref) provides the more defensible nonlinear check.
It fixes one parameter at a time, reoptimizes all others, and traces the connected
set below the one-degree-of-freedom likelihood-ratio threshold. Endpoint status
distinguishes a threshold crossing from a candidate bound or solver failure.
Profiles can reveal asymmetric, bound-limited, and flat directions that the
local covariance cannot. They remain local connected profiles; deterministic
multistart or a global method is needed to rule out disconnected regions.

These intervals only have coverage meaning when the supplied covariance is a
credible total error model. Paper rounding, uncertainty in conductor catalogs,
and a mismatch between modified Carson and Deri are source-model errors, not
meter noise, and should not be disguised as very precise measurement sigma.

## Numerical formulation

- Discrete candidates are enumerated outside the NLP.
- Continuous variables are affinely scaled to `[0,1]`.
- The default objective is smooth weighted least squares.
- Ipopt uses limited-memory Hessians, exact (unrelaxed) variable bounds, and
  deterministic multistart. This prevents trial points from escaping the
  candidate domain validated before optimization.
- Alternative optimizers do not receive Ipopt-specific raw attributes and must
  honor the JuMP variable bounds during nonlinear evaluations.
- Distances remain strictly positive through physical candidate bounds.
- Candidate compatibility requires every standardized residual to lie within
  `acceptance_sigma` (three by default).
- A ForwardDiff Jacobian is evaluated at the best local solution for rank and
  conditioning diagnostics.

The prototype supports modified Carson and four canonical overhead families:
`:horizontal_3`, `:triangle_3`, `:horizontal_4`, and `:neutral_under_4`.
Concentric-neutral and tape-shield cables, multi-frequency observations, global
confidence regions, and joint operational-data estimation remain future work.

## Validation datasets and conventions

The test suite contains two machine-readable benchmark families under
`test/data/inverse_carson`:

- `paper_table_iv.toml` reconstructs all five overhead cases in Table IV of Tam,
  Geth & Mithulananthan from Tables I, II and V. The published sequence values
  are reproduced within two units of their last printed decimal. The paper does
  not explicitly state frequency or earth resistivity; 50 Hz and 100 Ω·m are
  labelled as *inferred* from the constants in Table I, not as reported facts.
- `opendss_mars_triangle.toml` and its `.dss` input retain the primitive matrices
  from DSS C-API 0.14.3/OpenDSS SVN 3723. It explicitly records Deri earth return,
  no transposition, no neutral, units, soil resistivity, and the choice
  `Rdc=Rac` that prevents an unrelated frequency interpolation.

The independent OpenDSS case is intentionally not forced to agree with modified
Carson. With tight uncertainty it is rejected as model mismatch, primarily in
zero-sequence reactance. This is a validation success: changing earth-return
formulation must not be absorbed into a fictitious conductor geometry.

Every benchmark records frequency, soil resistivity, earth model, transposition,
neutral elimination, symmetrical-component convention, units, and the status of
measurement uncertainty. Published cases distinguish unavailable measurement
uncertainty from a rounding-only surrogate; deterministic OpenDSS output records
reproducibility tolerance rather than pretending it is measurement noise.

The experiments additionally cover correlated noise, rounded inputs, a wrong
Libra-vs-Mars conductor catalog, incorrect soil resistivity, and loose shunt
susceptance. The latter visibly widens the height confidence interval, agreeing
with the paper's warning that height is highly sensitive to shunt error.

## Literature and modelling context

The electromagnetic starting point is Carson's homogeneous half-space ground
return formulation [Carson (1926)](https://doi.org/10.1002/j.1538-7305.1926.tb00122.x).
The implementation uses the familiar power-frequency modified-Carson
approximation for inverse fitting. Deri et al.'s complex-depth image plane
provides a different, efficient approximation
([Deri et al. (1981)](https://doi.org/10.1109/TPAS.1981.317011)); the OpenDSS
benchmark makes this model choice observable rather than treating all
"Carson-like" calculations as interchangeable. Keshtkar, Solanki & Solanki
compare these approximation families and their frequency/soil sensitivity
([IEEE TPWRD 2014](https://doi.org/10.1109/TPWRD.2013.2276061)).

For a four-wire construction, the paper comparison uses the impedance Schur
complement obtained by imposing zero neutral voltage. That is a data convention,
not a claim that the downstream operational network has a perfectly grounded
neutral. The returned primitive remains four-wire. The modelling consequences
of impedance transformations under sparse neutral grounding are discussed by
[Geth, Heidari & Koirala (2022)](https://doi.org/10.1145/3538637.3538844), while
[Low (2024)](https://arxiv.org/abs/2403.17391) treats reverse Kron reduction as a
separate network-identification problem. Shunt comparison uses the phase
principal block corresponding to a grounded eliminated neutral; another neutral
termination requires another observation model.

Finally, Jacobian rank addresses local structural sensitivity, whereas finite
noise creates *practical* identifiability limits. The distinction and the use of
profile likelihood are standard in inverse modelling; see
[Raue et al. (2009)](https://doi.org/10.1093/bioinformatics/btp358). Operational
time-series line estimation is complementary rather than equivalent: for
example, [Vanin et al.](https://arxiv.org/abs/2209.10938) estimate impedance
matrices jointly with unbalanced network state, whereas inverse Carson maps
linecode data to a constrained construction catalog.

## Relationship to parameter estimation

[Parameter estimation](parameter_estimation.md) fits line lengths and transformer
taps from multiple operational snapshots. Inverse Carson instead starts from a
per-length sequence linecode and infers plausible construction data. Keeping the
two formulations separate makes their distinct data assumptions and
identifiability limits explicit. A future joint estimator can use compatible
inverse-Carson candidates as structural priors.
