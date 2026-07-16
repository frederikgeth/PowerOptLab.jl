# Constrained nonlinear least-squares state estimation

> **Kind:** compiled problem specification · **Maturity:** prototype · **Direction:** inverse · **Temporal:** single snapshot or sequential time series

This is the new four-wire estimator.  It solves the equality-constrained
nonlinear least-squares problem

```math
\min_x \frac12\lVert r(x)\rVert_2^2 \qquad \text{subject to}\qquad c(x)=0.
```

It is separate from the [legacy WLS/Ipopt formulation](state_estimation.md).
Use this API when the distinction between **uncertain information** and a
**genuinely exact electrical equation** matters, when neutral displacement must
remain in the state, or when branch telemetry/exact nonlinear devices are
required.

!!! warning "Prototype scope"
    The dense and sparse solvers are reference/prototype implementations.  The
    sparse Hachtel step is available, but large-feeder benchmarks, bad-data
    processing, topology hypotheses, transformer branch telemetry, and sparse
    rank/covariance backends remain future work.

## Modelling decisions

### Three kinds of information

The formulation deliberately keeps three categories separate.

| Category | Representation | Examples |
|---|---|---|
| Network identity | compiled into `SEStructure` | passive Ybus, conductor incidence, closed-switch aliases |
| Exact equation | `c(x)=0` | true zero injection, exact device law, ideal source relationship |
| Stochastic information | whitened residual `r(x)` | meters, forecasts, nominal loads, state priors |

A small meter variance does **not** make a reading exact.  Exact constraints
reduce the feasible tangent space and can make a model inconsistent; use them
only for physics known without uncertainty.

### State, grounding, and references

The state is rectangular conductor-to-ground voltage,
`[real(V_free); imag(V_free)]`.  A source terminal with an imposed phasor is
eliminated into the fixed-voltage parameter vector.  Perfectly grounded BMOPF
terminals are earth (`V=0`); an ungrounded neutral remains an explicit state.
Thus a global phasor reference is not silently confused with a neutral ground.
Floating/common-mode freedoms appear in the tangent-space observability result.

`compile_state_estimator` imports BMOPFTools' passive `I = YV` in SI units.
Closed ideal switches use BMOPFTools' node aliases.  The compiled structure owns
ordering and sparsity; [`SEParameters`](@ref) owns readings, standard deviations,
source phasors, device values, and priors.  Update parameters between snapshots
instead of recompiling whenever topology is unchanged.

### Sign convention

Passive and branch currents are positive **from a bus into an element**.  An
exact device branch current is positive from `TerminalConnection.positive` to
`.negative`; positive constant power therefore denotes consumption.  Use signed
negative power for generation.  Branch `:pflow`/`:qflow` have the same
into-the-line convention at the requested `side`.

## Build and solve

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

# `net` is a BMOPF network.  It may retain ungrounded neutral conductors.
measurements = [
    Measurement(kind=:vr, bus="b1", terminal="1", reference=nothing,
                value=230.0, sigma=1.0),
    Measurement(kind=:vi, bus="b1", terminal="1", reference=nothing,
                value=0.0, sigma=1.0),
    Measurement(kind=:vmag, bus="b1", terminal="2", value=230.0, sigma=1.0),
]

structure = compile_state_estimator(net, measurements;
                                    neutral="n",
                                    zero_injection=[("b2", "1")])
parameters = SEParameters(structure, measurements)
x0 = zeros(2length(structure.free_state_map))

# Dense: transparent small-system reference solver.
dense = solve_compiled_state_estimator(structure, parameters, x0)

# Sparse: augmented Hachtel/SuiteSparse-QR step solver.
sparse = solve_sparse_state_estimator(structure, parameters, x0)
```

Trust an estimate only for `:converged_unique` or
`:converged_underobserved`.  The latter is feasible but non-unique.  Other
statuses distinguish, among other conditions, constraint restoration failure,
invalid device domain, trust-region stall, and numerical failure.

## Measurements

All values are SI and every scalar has an independent standard deviation
`sigma`; residual rows are `(prediction - value) / sigma`.

### Terminal measurements

[`Measurement`](@ref) supports:

- `:vr`, `:vi`, `:vmag` — voltage component/magnitude across terminal to its
  `reference` (default: `neutral`; `nothing`: earth);
- `:pinj`, `:qinj` — terminal injection power into the passive network.

Magnitude derivatives are undefined at zero.  Set
`SEParameters(...; magnitude_epsilon=...)` only when smoothing below the
instrument's meaningful resolution is appropriate.

### Line telemetry

[`BranchMeasurement`](@ref) attaches to a named BMOPF line and `side=:from` or
`:to`:

```julia
BranchMeasurement(kind=:imag,  line="l1", side=:from, terminal="1",
                  value=12.3, sigma=0.2)      # amperes
BranchMeasurement(kind=:pflow, line="l1", side=:to, terminal="1",
                  value=-2_000.0, sigma=50.0) # watts into the line
```

Supported kinds are `:ire`, `:iim`, `:imag`, `:pflow`, and `:qflow`.  They use
BMOPFTools' public `line_yprim` primitive, so linecode truncation and shunt
stamping match the passive network model exactly.  This requires the
BMOPFTools release containing `line_yprim` (introduced by PR #348).

## Exact devices and continuation

Use [`ExactDeviceEquation`](@ref) only when a device law is truly exact.
`TerminalConnection` is a general oriented branch: phase-neutral connections
model wye/single-phase devices; phase-phase connections model delta devices.

```julia
load = ExactDeviceEquation(ConstantPowerDevice(
    [TerminalConnection(("b1", "1"), ("b1", "n"))],
    ComplexF64[4_000 + 1_000im],
))

structure = compile_state_estimator(net, measurements; exact_devices=[load])
parameters = SEParameters(structure, measurements; exact_devices=[load],
                           voltage_min_model=1.0)
result = solve_with_continuation(structure, parameters, x0;
                                 alphas=[0.0, 0.5, 1.0])
```

Available models are [`ConstantPowerDevice`](@ref),
[`ConstantCurrentDevice`](@ref), and [`ZIPDevice`](@ref).  Constant-power
evaluation rejects trial states below `voltage_min_model`.  Continuation starts
with a regularised internal law and **must finish at `α=1`** before accepting a
physical result.

## Diagnostics, multipliers, and uncertainty

[`observability_diagnostics`](@ref) evaluates rank on the tangent space
`H*Z`, where `C*Z=0`; rank of the raw measurement Jacobian alone is not the
relevant test.  [`unobservable_directions`](@ref) generates local directions
only on request.

Use [`selected_state_covariance`](@ref) or [`derived_covariance`](@ref) for
requested covariance blocks/derived quantities.  They throw when the tangent
space is rank deficient rather than returning a fictitious finite covariance.

The sparse result exposes `constraint_multipliers` in the same order as
`evaluation.constraints`.  Large, physically scaled multipliers are a useful
lead for a bad zero-injection label or an incorrect exact-device specification;
they are diagnostics, not automatic proof of a bad constraint.

## Time series

For unchanged topology, provide one parameter object per snapshot:

```julia
p1 = SEParameters(structure, measurements_t1; exact_devices=[load])
p2 = SEParameters(structure, measurements_t2; exact_devices=[load])
series = solve_time_series_state_estimator(structure, [p1, p2], x0;
                                            previous_state_sigma=10.0,
                                            solver=:sparse)
```

The prior from snapshot `t-1` is a whitened residual, not a hard voltage
constraint.  The driver warm-starts each snapshot, stops on the first failed
one, and returns `:time_series_stalled` without publishing later estimates.

## Current limitations

- Diagonal measurement covariance only; no correlated whitening yet.
- Line telemetry only; transformer, source-phasor, neutral-current, and angle
  measurement models are pending.
- No robust bad-data test, automated bad-constraint ranking, topology-error
  hypothesis search, or full filter globalisation.
- Sparse rank/uncertainty methods and cached symbolic factorisations are not yet
  production-scale implementations.

See [`SEStructure`](@ref), [`SEParameters`](@ref), and the solver/result types
in the API reference for the full callable interface.
