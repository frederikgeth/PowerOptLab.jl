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
constructor supplies a descriptive 1% tolerance. Scores should only be treated
as statistical quantities when `sigma` genuinely represents measurement and
source-model uncertainty.

## Numerical formulation

- Discrete candidates are enumerated outside the NLP.
- Continuous variables are affinely scaled to `[0,1]`.
- The default objective is smooth weighted least squares.
- Ipopt uses limited-memory Hessians and deterministic multistart.
- Distances remain strictly positive through physical candidate bounds.
- Candidate compatibility requires every standardized residual to lie within
  `acceptance_sigma` (three by default).
- A ForwardDiff Jacobian is evaluated at the best local solution for rank and
  conditioning diagnostics.

The prototype supports modified Carson and four canonical overhead families:
`:horizontal_3`, `:triangle_3`, `:horizontal_4`, and `:neutral_under_4`.
Concentric-neutral and tape-shield cables, profile intervals, correlated
measurement covariance, and joint operational-data estimation remain future
work.

## Relationship to parameter estimation

[Parameter estimation](parameter_estimation.md) fits line lengths and transformer
taps from multiple operational snapshots. Inverse Carson instead starts from a
per-length sequence linecode and infers plausible construction data. Keeping the
two formulations separate makes their distinct data assumptions and
identifiability limits explicit. A future joint estimator can use compatible
inverse-Carson candidates as structural priors.
