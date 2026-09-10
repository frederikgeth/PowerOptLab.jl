# Constrained nonlinear least-squares state estimation

> **Kind:** Problem specification (compiled) · **Maturity:** prototype · **Direction:** inverse · **Temporal:** single snapshot or sequential time series

This is the new four-wire estimator.  It solves the equality-constrained
nonlinear least-squares problem

```math
\min_x \frac12\lVert r(x)\rVert_2^2 \qquad \text{subject to}\qquad c(x)=0.
```

It is separate from the [legacy WLS/Ipopt formulation](state_estimation.md); see
[Choosing a formulation](../estimation/comparison.md) for the head-to-head.
Use this API when the distinction between **uncertain information** and a
**genuinely exact electrical equation** matters, when neutral displacement must
remain in the state, or when branch telemetry/exact nonlinear devices are
required.

!!! warning "Prototype scope"
    The dense and sparse solvers are reference/prototype implementations.  The
    sparse Hachtel step is available, but large-feeder benchmarks, bad-data
    processing, topology hypotheses, transformer branch telemetry, and sparse
    rank/covariance backends remain future work.

## Electrical support preflight and initialization

Call `state_estimator_preflight(net, measurements; zero_injection, exact_devices)`
before compiling to inspect the snapshot interpretation. The same report is
retained as `structure.preflight`; compilation throws `SEUnsupportedNetwork`
with its report if any error is present. `supported` means representable by this
voltage-only model. It does **not** certify observability, a correct network
calibration, or convergence.

```julia
report = state_estimator_preflight(net, measurements)
report.findings                 # severity, code, element, explanation
report.elements                 # passive elements, aliases, fixed taps
report.omitted                  # injections, controls, operating bounds
report.sources                  # prescribed phasors
report.grounded_terminals       # declared perfect grounds
report.measurements             # locations and voltage references
report.exact_equations          # requested zero injections and devices

s = compile_state_estimator(net, measurements)
p = SEParameters(s, measurements)
x0 = initial_state_estimator(s, p)
result = solve_sparse_state_estimator(s, p, x0)
```

The initializer solves the free-conductor, load-free current equations with
sparse QR using the current source phasors. Transformer ratios and phase shifts
therefore enter the initial guess. This is a numerical starting point, not a
load-flow operating point, a prior, or an extra measurement. A source-free or
reference-deficient network may require a physical warm start; floating modes
are not repaired by initialization. Trust-region voltage scaling uses local
load-free bus levels (with a 1 V floor), so a low-voltage secondary does not
inherit the primary's voltage scale. Neutral components inherit their bus's
phase scale without being grounded.

| Electrical feature | Contract |
|---|---|
| Finite lines, shunts and capacitors | Valid finite primitives and declared terminal maps required |
| Closed ideal switches / negligible-impedance line aliases | Imported as node aliases; open switches remain open |
| Single-phase, center-tap, wye–delta, delta–wye transformers | Finite leakage, valid winding maps and positive fixed ratios |
| Single-phase and open-delta regulators | Finite leakage and known fixed taps; no automatic control or tap estimation |
| General n-winding transformers | Finite primitive, positive winding voltages, WYE/DELTA connections, equal nonzero phase counts and a/b/c/n names required; not covered by the component oracle matrix |
| Ideal or degenerate transformers | Rejected; auxiliary current states and exact winding relations require a future augmented formulation |
| Loads, generators and inverters in network data | Not automatically imported as injections; listed in `omitted` |
| Operating limits and controls | Inventoried but not imposed on the estimated state |

Do not insert an arbitrary small transformer impedance to bypass rejection.
Use physical leakage data or an appropriate augmented formulation. The guard
matches the pinned BMOPFTools exporter: ordinary transformer leakage components
all at or below `1e-6` ohm are classified as ideal, and series-line norm at or
below `1e-4` ohm follows its alias rule. These are exporter thresholds, not
physical accuracy guarantees. Primitive warnings become preflight errors,
rather than allowing a skipped or shunt-only substitute to pass silently.

For n-winding data, fixed ratios are encoded in the winding `v_nom` values;
separate tap fields are rejected because the exporter does not apply them.

Recompile after changing taps, topology, impedances or terminal maps.
`SEParameters` updates measurements, source phasors and supported exact-device
parameters; it cannot change compiled transformer ratios. A forecast remains an
uncertain measurement with an appropriate covariance, rather than an exact
zero-injection constraint. A successful preflight does not turn nominal load
data into pseudomeasurements.

## Transformer and regulator validation

`test/state_estimation_network_tests.jl` exercises both compiled solvers on
single-phase (including a 1.03 primary tap), wye–delta, delta–wye and center-tap transformers, single-phase and
open-delta fixed regulators, and **SE-IEEE13-finite**. The component fixtures and
OpenDSS decks are versioned under `test/`; their provenance and modifications
are recorded alongside the decks.

SE-IEEE13-finite is a deliberately modified IEEE 13-node network: finite winding
resistance of 0.5% per winding, regulator XHL of 1%, transformer XHL of 2%, fixed
unequal phase taps, constant-power loads, omitted line shunts, and a source at
bus 650. It includes mixed phases, delta loads, mutual line coupling and a
4.16 kV / 480 V voltage transition. It is **not** a reproduction of the original
IEEE operating case. The open-delta component includes physical 1 kΩ
phase-to-earth shunts to establish its common-mode reference in both models.

OpenDSS supplies conductor voltages. Injection readings are constructed from
specified load currents and those voltages, independently of the estimator's
Ybus evaluator. Actual OpenDSS source-terminal voltages are used as the fixed
boundary, avoiding an artificial disagreement from source impedance. Tests use
redundant rectangular-voltage and active/reactive injection readings, exact
zero-injection equations, and reproducible Gaussian perturbations. Both solvers
start from the load-free initializer. They must converge uniquely, satisfy the
exact equations, agree with each other, and yield finite selected covariance.

Voltage-component errors are divided by the corresponding oracle phasor
magnitude (100 V floor), rather than by an imaginary component near zero. The
regression bounds are 0.2% for clean data and 0.6% for the noisy realization.
Voltage standard deviations are 0.2% of that phasor scale; power standard
deviations are 1% of apparent injection (10 W/var floor). Noise is one standard
deviation with seed 349. The estimator stationarity tolerance is `1e-6` in its
residual-equivalent norm; exact-current constraints retain the default solver
tolerance. JuMP objective agreement is checked at `atol=rtol=1e-5`, and fitted
state components agree within 0.01 V. These tolerances are regression acceptance
criteria, not confidence bounds.

A separate JuMP/Ipopt WLS formulation differentiates scalar measurement
equations independently and checks the fitted state and objective. This is an
independent optimization check, **not** an independent admittance model: it
shares the compiled electrical matrices. OpenDSS provides the separate
physical-model comparison. Noisy fits need not interpolate every reading.

This is a regression matrix, not a Monte Carlo uncertainty-calibration study
or a proof of coverage for general distribution networks. It uses generous
redundant instrumentation; sparse field telemetry, correlated load forecasts,
bad data, unknown taps, switching hypotheses, automatic regulator controls,
ideal transformers and large-feeder runtime/memory scaling need separate
validation. Transformer branch telemetry is still outside this API.

Run the complete reproducible matrix with `julia --project -e 'using Pkg;
Pkg.test()'`; OpenDSSDirect is a declared test dependency.

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

!!! note "`zero_injection` covers the neutral here"
    Four-wire KCL is stated per conductor against earth, so a bare bus id in
    `zero_injection` expands over **all** of the bus's non-grounded terminals,
    its neutral included: nothing attached to a bus means nothing injected into
    any of its conductors.  The [WLS estimator](state_estimation.md) states zero
    injection per phase against a return terminal and deliberately expands over
    phase terminals only.  Passing `(bus, terminal)` tuples is explicit either
    way.

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

# Minimal runnable feeder, with all quantities in SI units.
net = parse_bmopf("""
{"bus":{"src":{"terminal_names":["1"]},
        "b1":{"terminal_names":["1"]},"b2":{"terminal_names":["1"]}},
 "voltage_source":{"s":{"bus":"src","terminal_map":["1"],
     "v_magnitude":[230.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.1}},
 "line":{"l1":{"bus_from":"src","bus_to":"b1","terminal_map_from":["1"],
     "terminal_map_to":["1"],"linecode":"lc","length":1.0},
         "l2":{"bus_from":"b1","bus_to":"b2","terminal_map_from":["1"],
     "terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
"""; from_string=true)
measurements = [
    Measurement(kind=:vr, bus="b1", terminal="1", reference=nothing,
                value=230.0, sigma=1.0),
    Measurement(kind=:vi, bus="b1", terminal="1", reference=nothing,
                value=0.0, sigma=1.0),
    Measurement(kind=:vmag, bus="b2", terminal="1", reference=nothing,
                value=230.0, sigma=1.0),
]

structure = compile_state_estimator(net, measurements;
                                    neutral="n",
                                    zero_injection=[("b2", "1")])
parameters = SEParameters(structure, measurements)
nf = length(structure.free_state_map)
x0 = vcat(fill(230.0, nf), zeros(nf)) # nonzero magnitude start

# Dense: transparent small-system reference solver.
dense = solve_compiled_state_estimator(structure, parameters, x0)

# Sparse: augmented Hachtel/SuiteSparse-QR step solver.
sparse = solve_sparse_state_estimator(structure, parameters, x0)
```

Successful first-order termination is reported as `:converged_unique` or
`:converged_underobserved`.  The latter has locally unobserved directions. Neither status certifies a global
solution or a strict local minimum.  Other
statuses distinguish, among other conditions, constraint restoration failure,
`invalid_initial_domain`, trust-region stall, `:undefined_derivative` (a `:vmag`
or `:imag` row sitting exactly at zero, where the magnitude derivative does not
exist), and numerical failure. A numerical reduced-Lagrangian curvature check rejects
negative curvature with `:stationary_not_minimum`; `:curvature_check_failed`
means a nearby domain/derivative evaluation prevented the check. These are
non-publishable statuses. The check is performed at nonlinear stationary points
with non-negligible residuals; it is not a formal second-order certificate.

!!! warning "`:converged_unique` is a LOCAL statement"
    It reports that the reduced Jacobian ``HZ`` has full rank **at the returned
    point**.  It is not a global uniqueness certificate.  Two distinct states
    can each earn it while fitting the same data exactly — see
    [current-magnitude measurements](../estimation/current_magnitude.md).

### Step acceptance and convergence

`penalty` is the initial merit weight. Both solvers increase it when a step's
predicted objective cost would overwhelm its predicted feasibility improvement.
The history records the changing penalty and separates measurement from prior
objectives, so merit values can be interpreted with their associated weight.
The dense solver uses SVD least-squares steps for rank-deficient systems; the
sparse system uses a negative primal diagonal so eliminating the residual block
adds positive damping to the reduced least-squares operator. These safeguards
do not constitute a general global-convergence or infeasibility certificate.

### Convergence tolerances are scale relative

`norm(c)` is a difference of ``YV`` products, so the smallest value it can
reach in double precision scales with the network's SI current magnitude — near
``5\times10^3`` A on a 230 V feeder with 20 S conductors, near ``3\times10^9``
A on a 132 kV feeder with short 20 kS segments.  Feasibility is therefore tested
against `max(constraint_tolerance, constraint_rtol * current_scale)`: the
absolute request acts as a floor, so small networks behave exactly as before,
while a high-voltage network is no longer asked for a precision that double
precision cannot deliver.  Stationarity is divided by `norm(H, Inf)`, so
`optimality_tolerance` reads in units of the whitened residual rather than in
whatever units ``\partial h/\partial x`` happens to carry.

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
instrument's meaningful resolution is appropriate. It applies to `:vmag`;
branch `:imag` uses `current_epsilon` in amperes.  Note that a smoothed magnitude row evaluated at
exactly zero is *differentiable but information-free* — its gradient there is
zero — so zero gradient is not evidence of minimisation. If the requested magnitude
is larger than the smoothed prediction, the point can be a local maximum and
the solver reports `:stationary_not_minimum`.  See [current-magnitude measurements](../estimation/current_magnitude.md).

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

!!! note "Branch power is referenced to earth"
    Unlike a node [`Measurement`](@ref), a `BranchMeasurement` takes no
    `reference`: `:pflow`/`:qflow` use the measured conductor's voltage to
    **earth**.  On a four-wire line with a displaced neutral that is not the
    phase-to-neutral power a real meter reports.  A per-measurement reference
    is [on the roadmap](../estimation/state_of_the_art.md#Roadmap).

## Exact devices and continuation

Use [`ExactDeviceEquation`](@ref) only when a device law is truly exact.
`TerminalConnection` is a general oriented branch: phase-neutral connections
model wye/single-phase devices; phase-phase connections model delta devices.

```julia
load = ExactDeviceEquation(ConstantPowerDevice(
    [TerminalConnection(("b1", "1"), nothing)] # earth return in this one-wire example,
    ComplexF64[4_000 + 1_000im],
))

structure = compile_state_estimator(net, measurements; exact_devices=[load])
parameters = SEParameters(structure, measurements; exact_devices=[load],
                           voltage_min_model=1.0)
result = solve_with_continuation(structure, parameters, x0;
                                 alphas=[0.0, 0.5, 1.0])
```

The constant-current term is a fixed complex phasor; it does not rotate with
the terminal voltage to preserve power factor. `ZIPDevice` implements the
explicit law `conj(S/V) + I + YV`, rather than every conventional ZIP-load
parameterisation.

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
requested quantity depends on an unobservable tangent direction. A locally
identifiable quantity can have finite first-order covariance even when other
state directions remain unobserved. This does not establish global identifiability.

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

A prior the caller configured on an `SEParameters` is **preserved**.  When
`previous_state_sigma` is supplied the previous-state prior is layered on top of
it rather than replacing it. Algebraically, duplicated indices add independent
quadratic rows. That is statistically justified only if the information sources
can be treated as independent; a previous estimate usually reuses earlier meter
data and therefore should not be described as a second independent meter.

This driver is temporal regularisation, not a state-space filter: it does not
propagate the previous covariance, specify process noise, or account for
cross-time correlation. Choose `previous_state_sigma` as an explicit movement
assumption and do not interpret the resulting local covariance as a complete
filtered posterior.

## Current limitations

- Diagonal measurement covariance only; no correlated whitening yet.
- Line telemetry only; transformer, uncertain source-phasor, and angle
  measurement models are pending.
- No robust bad-data test, automated bad-constraint ranking, topology-error
  hypothesis search, or full filter globalisation.
- Jacobians use sparse triplets and cached transposed operators for row access.
  Sparse rank/uncertainty methods and cached symbolic factorisations are not yet
  production-scale implementations.

See [`SEStructure`](@ref), [`SEParameters`](@ref), and the solver/result types
in the API reference for the full callable interface.

## Validation and computational scope

Run the focused suite with
`julia --project=. -e 'using Test, PowerOptLab; include("test/constrained_state_estimation_tests.jl")'`.
It includes loaded four-wire finite-difference checks with complex and mutual
impedance, both line ends, an independent analytic noisy reactive-power optimum,
conflicting exact constraints, square rank deficiency, negative-curvature
rejection, identifiable covariance, and seeded Monte Carlo variance of a linear
estimate. These checks complement the existing power-flow recovery tests;
they do not establish robustness on operational feeders or nonlinear coverage.

Run `julia --project=. scripts/benchmark_compiled_state_estimation.jl` to
reproduce sparse-chain allocation and timing measurements. With Julia 1.12.6,
12 BLAS threads, and compilation excluded, the September 2026 check gave:

| States | Residual evaluation allocated | Residual Jacobian allocated |
|---:|---:|---:|
| 100 | 28 kB | 115 kB |
| 400 | 105 kB | 457 kB |
| 1,600 | 420 kB | 1,969 kB |

The previous implementation allocated about 21 MB and 91 MB respectively at
1,600 states. Sparse assembly removes that quadratic allocation in this case.
The complete sparse solve still performs dense rank/null-space diagnostics: an
already exact 1,600-state warm start allocated about 196 MB. No real-time or
large-feeder performance claim follows from these synthetic measurements.

Remaining validation priorities are independent unbalanced feeder oracles,
nonlinear Monte Carlo coverage, bad-data and model-error experiments, and
end-to-end feeder/time-series benchmarks. Remaining implementation priorities
include correlated whitening, robust and leverage-adjusted bad-data processing,
reference-aware branch power and transformer telemetry, multistart, symbolic
factorisation reuse, reusable numerical workspaces, and sparse rank/covariance.

### Compatibility notes

`residual_jacobian` and `constraint_jacobian` now return `SparseMatrixCSC`.
Use `collect(Float64, H)` when calling a dense-only routine such as `svdvals`.
Covariance calls may now succeed for identifiable quantities on a partially
observed state. Failure to reach feasibility at the iteration limit is
`:constraint_not_satisfied`, not a claim that the mathematical problem is
infeasible. Invalid initial device domains and detected negative curvature have
explicit non-publishable statuses described above.
