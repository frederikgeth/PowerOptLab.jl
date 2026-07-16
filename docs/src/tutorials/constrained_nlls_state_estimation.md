# Constrained nonlinear least-squares state estimation: a modelling tutorial

> **Audience:** power-system researchers · **Scope:** four-wire distribution
> state estimation with exact network/device equations and uncertain telemetry.

Convergence is not automatically a defensible state estimate. A plausible LV
voltage profile can arise from an incorrectly grounded neutral, an uncertain
forecast treated as an equality, or an unobservable measurement set. This
tutorial explains how PowerOptLab's constrained nonlinear least-squares (NLLS)
estimator makes those choices explicit.

It solves

```math
\min_x \frac12\lVert r(x) \rVert_2^2 \qquad \text{subject to}\qquad c(x)=0,
```

where `r(x)` contains whitened *uncertain* information and `c(x)` only physics
regarded as exact. It is distinct from the legacy [`solve_state_estimation`](@ref)
WLS/Ipopt formulation.

## 1. Compile topology; update numerical information

The state is rectangular conductor-to-ground voltage. An ungrounded neutral is
therefore a state; fixed source phasors are parameters, and perfectly grounded
terminals are earth. A source reference must not be confused with a grounded
downstream neutral.

Compile the invariant electrical structure once and create a parameter object
for each snapshot's readings and uncertainty:

```julia
using PowerOptLab

structure = compile_state_estimator(net, measurements;
    neutral="n", zero_injection=[("b2", "1")])
parameters = SEParameters(structure, measurements)
x0 = zeros(2length(structure.free_state_map))
result = solve_sparse_state_estimator(structure, parameters, x0)
```

`SEStructure` owns terminal ordering, passive Ybus, source/ground treatment,
measurement incidence, and sparsity. `SEParameters` owns readings, standard
deviations, source phasors, exact-device values, and priors. Reuse the compiled
structure while topology is unchanged.

### Pitfall: silently collapsing the neutral

Grounding a neutral merely for convenience removes neutral displacement from the
state and can distort phase-to-neutral voltage and current estimates. Ground
only physically grounded terminals. A floating/common-mode freedom should show
up in observability diagnostics, not be hidden by an arbitrary reference.

## 2. Put each datum in the correct category

| Category | Representation | Examples |
|---|---|---|
| Network identity | compiled in `SEStructure` | topology, Ybus, conductor incidence, closed-switch aliases |
| Exact equality | `c(x)=0` | genuine zero injection, truly exact device law |
| Stochastic information | whitened `r(x)` | meters, forecasts, nominal P/Q, source uncertainty, state prior |

Every scalar measurement uses an independent `sigma`; its residual is
`(prediction - value) / sigma`.

```julia
measurements = [
    Measurement(kind=:vmag, bus="b1", terminal="1", value=229.4, sigma=0.8),
    Measurement(kind=:pinj, bus="b1", terminal="1", value=-4_200.0, sigma=200.0),
    Measurement(kind=:qinj, bus="b1", terminal="1", value=-1_100.0, sigma=120.0),
]
```

Terminal injection is positive into the passive network, so a load is normally
negative and generation positive.

### Pitfall: zero injection for an uncertain small load

`zero_injection` adds an exact KCL equation. It is appropriate only when the
injection is actually known to be zero. A small, missing, or forecast load is a
P/Q residual with an honest uncertainty. Replacing uncertainty by exactness can
make the model inconsistent and create a falsely sharp state estimate.

## 3. Model what the instrument actually measures

Terminal measurements support `:vr`, `:vi`, `:vmag`, `:pinj`, and `:qinj`.
Specify the physical meter return: by default it is the compilation `neutral`;
`reference=nothing` is terminal-to-earth.

```julia
measurements = [
    Measurement(kind=:vr, bus="b1", terminal="1", reference="n",
                value=228.0, sigma=0.5),
    Measurement(kind=:vi, bus="b1", terminal="1", reference="n",
                value=-2.0, sigma=0.5),
    Measurement(kind=:vmag, bus="b2", terminal="2", reference="n",
                value=231.0, sigma=1.0),
]
```

Line telemetry uses the same primitive admittance as the passive network:

```julia
telemetry = [
    BranchMeasurement(kind=:imag, line="l1", side=:from, terminal="1",
                      value=18.2, sigma=0.3),
    BranchMeasurement(kind=:pflow, line="l1", side=:to, terminal="1",
                      value=-3_800.0, sigma=80.0),
]
structure = compile_state_estimator(net, vcat(measurements, telemetry))
```

Branch power/current is positive from the requested bus *into* the line. Do not
assume a feeder-downstream sign convention without checking `side`.

### Pitfall: promoting ``|V|`` into a phasor

A magnitude meter has no angle information. Replacing it with invented real and
imaginary readings adds information and can make an unobservable system look
observable. Use `:vmag` for a magnitude meter; use `:vr`/`:vi` only for actual
rectangular phasor information. `magnitude_epsilon` is a near-zero numerical
smoothing option, not extra measurement information.

## 4. Use exact nonlinear devices sparingly

`ExactDeviceEquation` adds a law to `c(x)`. Use it for an exact benchmark or a
device law known without uncertainty:

```julia
connection = TerminalConnection(("b1", "1"), ("b1", "n"))
device = ExactDeviceEquation(ConstantPowerDevice(
    [connection], ComplexF64[4_000 + 1_000im]))

structure = compile_state_estimator(net, measurements; exact_devices=[device])
parameters = SEParameters(structure, measurements;
                           exact_devices=[device], voltage_min_model=1.0)
```

`ConstantPowerDevice`, `ConstantCurrentDevice`, and `ZIPDevice` support
phase-neutral and phase-phase `TerminalConnection`s. Positive device power
means consumption. Constant-power trials below `voltage_min_model` are rejected
rather than allowing an undefined `S/V` evaluation.

### Pitfall: making a customer forecast an exact constant-power device

A precise forecast is still uncertain. An exact forecast forces its error into
other measurements and voltages. Use a P/Q residual for forecasts; reserve
exact-device equations for a physical law that is genuinely enforced.

## 5. Diagnose feasibility and identifiability separately

The dense composite-step solver is a transparent small-system reference; the
sparse Hachtel/SuiteSparse-QR solver is the sparse pathway:

```julia
dense = solve_compiled_state_estimator(structure, parameters, x0)
sparse = solve_sparse_state_estimator(structure, parameters, x0)
```

Trust only `:converged_unique` or `:converged_underobserved`. The latter is
feasible but non-unique; restoration failure, invalid device domain,
trust-region stall, and numerical failure are not published estimates.

Observability belongs on the feasible tangent space:

```julia
obs = observability_diagnostics(structure, parameters, sparse.state)
obs.unobservable_dimension
obs.condition_number
```

If `C*Z = 0`, the relevant rank is that of `H*Z`, not the raw measurement
Jacobian. Request local ambiguous directions only for diagnosis:

```julia
directions = unobservable_directions(structure, parameters, sparse.state; count=1)
```

### Pitfall: declaring observability from solver convergence

An optimizer can converge to one point on an underdetermined manifold. An
incorrect exact constraint can also remove a degree of freedom and improve rank
artificially. Report tangent-space observable/unobservable dimensions,
conditioning, and the exact-constraint set—not only the exit status.

## 6. Quantify uncertainty only when it is finite

The residual Jacobian is whitened, so local Gauss--Newton covariance is
available for selected state blocks or derived quantities:

```julia
indices = [1, length(structure.free_state_map) + 1]
Σ_voltage = selected_state_covariance(structure, parameters, sparse.state, indices)
Σ_g = derived_covariance(structure, parameters, sparse.state, J_g)
```

The routines throw if tangent directions are unobservable, rather than returning
a fictitious finite covariance from a pseudoinverse.

### Pitfall: treating covariance as complete physical uncertainty

The available covariance is local Gauss--Newton uncertainty under *diagonal*
measurement covariance. It excludes correlated meter error, topology and
parameter uncertainty, gross bad data, and model discrepancy. State the
covariance model and use explicit scenarios for uncertainty outside it.

## 7. Use continuation and temporal priors as stated assumptions

Constant-power devices can be hard to solve from a poor initial state. Continue
an internal regularisation from `α=0` to the physical model at `α=1`:

```julia
continued = solve_with_continuation(structure, parameters, x0;
                                    alphas=[0.0, 0.25, 0.5, 0.75, 1.0])
```

Accept the result only if the sequence reaches `α=1`. For unchanged topology,
reuse the compiled structure across time:

```julia
p1 = SEParameters(structure, measurements_t1)
p2 = SEParameters(structure, measurements_t2)
series = solve_time_series_state_estimator(structure, [p1, p2], x0;
    previous_state_sigma=10.0, solver=:sparse)
```

The earlier state is a whitened prior residual, not a hard equality. The driver
warm-starts each snapshot and stops at the first failure.

### Pitfall: hard-constraining successive states to be equal

Load, PV, switching, and control actions move distribution states. A hard
`x_t=x_{t-1}` equality can conceal a real change or force conflict into meter
residuals. Select and report a prior sigma representing plausible movement.

## 8. Investigate inconsistency before excluding a datum

Residuals show tension with uncertain measurements. Sparse constraint
multipliers identify exact equations that are expensive to enforce. Large,
physically scaled multipliers can point to a bad zero-injection label or exact
device specification, but are not an automated bad-data test. The prototype has
no robust M-estimator, leverage-adjusted residual test, topology-hypothesis
search, or automatic bad-constraint ranking.

Use a documented sequence: check terminal/reference and signs; verify timestamp
and sigma; inspect residuals/multipliers; test credible topology or parameter
alternatives; only then exclude or downweight data with a recorded rationale.

## 9. Research-report checklist

Record: (1) four-wire topology, grounding, sources, and parameters; (2) meter
type/reference/sign/timestamp/sigma; (3) every exact constraint and its
rationale; (4) device model, voltage guard, solver, initialization, and
continuation; (5) residual/constraint/multiplier diagnostics; (6) tangent-space
observability and conditioning; (7) covariance scope and excluded uncertainty;
and (8) the temporal-prior model. This turns a numerical voltage solution into
a reproducible inference claim.

For the full API, see [Constrained nonlinear least-squares state estimation](../problems/constrained_state_estimation.md).
