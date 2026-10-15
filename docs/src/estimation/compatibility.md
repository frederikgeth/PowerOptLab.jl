# Measurement and electrical compatibility

The compiled estimator has a useful separation between topology, exact equations,
measurement functions and numerical solvers. It remains a research prototype:
its tests establish agreement on specific cases, not complete equipment coverage
or production readiness. Prefer this API for new work.

## One model, three solvers

Compile once with `compile_state_estimator`, create `SEParameters`, and solve
with `solve_compiled_state_estimator` (dense composite step),
`solve_sparse_state_estimator` (Hachtel/sparse QR), or
`solve_ipopt_state_estimator` (classic nonlinear programming).
`solve_state_estimation(structure, parameters, x0)` is also the shared Ipopt path.
All three consume the same residuals, exact constraints, first derivatives,
source phasors, customer-power profiles and priors. Solver choice does not change
what a reading means. The Ipopt backend uses limited-memory Hessians and a
complete constraint-Jacobian pattern; it is currently a compatibility/reference
backend, not a sparse scaling solution.

For floating networks, some KCL equations can be identities of the others.
The shared Ipopt backend selects an independent equation basis using the full
affine network/device coefficients plus independent nonlinear branch-current
incidences. It verifies those coefficient identities near roundoff; it does not
drop a physical equation just because its Jacobian is zero at the starting
point. Original equations remain in `evaluation.constraints` and every final
physical feasibility and rank check. `solve_diagnostics` reports both equation
counts and the physical row indices sent to the optimizer. Feasibility scaling
also includes amplification from the independent rows back to all physical
currents, including transformer turns ratios. Inconsistent constant
terms and weak grounding information are retained.

The original `solve_state_estimation(net, measurements)` remains a separate
native-engine compatibility API. It accepts `Measurement` only, requires paired
customer P/Q coverage or explicit zero injection, and does not accept the newer
measurement types or exact-device objects. Its source neutral must be explicitly
grounded: otherwise the native engine would introduce an implicit earth bond.
Use the compiled API for a source neutral connected through finite electrodes.

## Measurement functions

| Type | Quantities | Meaning |
|---|---|---|
| `Measurement` | `:vr`, `:vi`, `:vmag`, `:pinj`, `:qinj` | Voltage or individual customer power; customer powers also declare their injection paths |
| `NetInjectionMeasurement` | `:ire`, `:iim`, `:imag`, `:pinj`, `:qinj` | Total conductor injection into the passive network; no customer allocation and no inferred path |
| `BranchMeasurement` | `:ire`, `:iim`, `:imag`, `:pflow`, `:qflow` | Current or power entering one end/conductor of a finite line |
| `TransformerMeasurement` | same element quantities | Current or power entering a transformer port/conductor |
| `ShuntMeasurement` | same element quantities | Current or power entering a finite shunt, including an electrode |
| `FunctionalMeasurement` | scalar user function | A function of selected rectangular conductor voltages, with automatic or supplied first derivatives |

All readings and standard deviations are SI. Currents into passive elements are
positive; customer injection power is positive into the network. Exact-device
positive power denotes consumption. Element power uses conductor-to-earth
voltage by default; `reference="n"` selects that bus's neutral instead.
`NetInjectionMeasurement` current has no voltage reference. Unlike customer
power, a net meter observes the sum at a merged electrical node. Individual customer power at a source-fixed conductor (including aliases) is
rejected: the source reaction prevents separating that customer from the net
current. Use `NetInjectionMeasurement` for net source telemetry instead. Declare
`injection_connections` separately to retain unknown loads when telemetry drops
out. A custom function also declares no injection path.

A functional callback receives `[real(V[nodes]); imag(V[nodes])]` in the supplied
node order, including repeated nodes and zeros for ideal earth. Return one finite
real scalar. `ForwardDiff` supplies its gradient unless `gradient` is provided;
that callback must return one entry per input component. Functions must be pure,
with immutable captured constants. Changing the function, node order or voltage
reference requires recompilation; updating a reading or its uncertainty does not.
Callbacks cannot introduce hidden current states, equipment equations or a
non-Gaussian error model.

```julia
neutral_current = NetInjectionMeasurement(kind=:ire, bus="load", terminal="n",
    value=2.0, sigma=0.1)
electrode_current = ShuntMeasurement(kind=:ire, shunt="earth", terminal="n",
    value=0.4, sigma=0.05)
secondary_power = TransformerMeasurement(kind=:pflow, transformer="tx",
    variant="single_phase", port=2, terminal="1", reference="n",
    value=2000.0, sigma=30.0)
squared_voltage = FunctionalMeasurement(name=:squared_voltage,
    nodes=[("load", "1"), ("load", "n")],
    prediction=z -> (z[1]-z[2])^2 + (z[3]-z[4])^2,
    value=230.0^2, sigma=100.0)
```

The identifiers in this example must exist in the supplied network. These types
are exercised by `test/state_estimation_compatibility_tests.jl`, including
finite-difference derivative checks and agreement across all three solvers.

## Transformer coverage

Finite-leakage `single_phase`, `center_tap`, `wye_delta`, `delta_wye`,
`single_phase_autotransformer`, `open_delta_regulator`, and `n_winding` are
supported within the electrical preflight contract. Regressions compare all
three shared solvers and compare compatible node telemetry with native WLS.
Ordinary transformer ports are 1/from and 2/to; center-tap port 2 includes both
legs and the center conductor. An n-winding port is its winding-list index.

The open-delta shared phase has an ideal copper bond **in addition to** the finite
winding primitive. Compilation identifies its two nodes exactly, matching the
native engine; it does not insert a large fictitious admittance. Shared-phase
port current includes an ideal-wire reaction and cannot be obtained from the
winding primitive alone, so its element telemetry is rejected. Meter adjacent
finite branches instead. N-winding telemetry is also rejected where different
ports share the same original conductor and the exporter has already summed
those currents. Other port currents use the original primitive order, including
the distinct autotransformer ordering.

This is not support for every transformer model. Ideal/zero-leakage transformers,
unrepresented vector groups/phase shifts, and automatic tap control remain
outside the voltage-only contract. Fixed supported taps are network parameters;
recompile after changing them. Preflight rejects ignored or unsupported fields
rather than treating them as implemented physics.

## Neutral grounding and telemetry

Floating downstream neutrals, solid earth bonds, complex finite grounding
admittances and multiple electrodes are regression-tested. Model finite electrodes
as BMOPF shunts; their current and absorbed power are distinct from the current
in a neutral line. A source-free case is tested with absolute voltage telemetry
and with a remaining common-mode gauge. Grounding and a phasor reference are
separate pieces of information.

An ideal grounded conductor has zero earth-referenced voltage and power, but
may carry nonzero element current. Its net earth-reaction current is eliminated
from reduced Ybus and cannot be measured using `NetInjectionMeasurement`.
Element current telemetry remains available where a finite primitive determines
it. Constant zero magnitudes have zero derivatives; a *variable* magnitude at
zero is nonsmooth. Supply explicit `magnitude_epsilon` or `current_epsilon` when
smooth approximations are appropriate, and account for that change in the model.

Compiled source phasors fix only the listed conductor-to-earth voltages. A
floating phase-to-neutral voltage source with unknown common-mode voltage is a
different model and is not implemented. Soil potentials, frequency-dependent
electrodes and separate ideal-bond reaction-current states are also absent.

## What the consistency tests establish

`test/state_estimation_consistency_tests.jl` adds references that do not simply
generate an answer from the estimator's own Ybus:

| Test | Independent evidence |
|---|---|
| Four-wire linear WLS | Series Ohm's law and electrode KCL give a closed-form constrained optimum and covariance |
| Representation changes | Renaming buses, reversing a line, permuting conductors/readings, and consistent voltage/current rescaling preserve the statistical problem |
| Noisy nonlinear fits | Native SI and per-unit WLS agree with all three shared solvers on the same readings and tolerances |
| Floating exact devices | Hand-constructed constant-current, constant-power and ZIP equilibria expose redundant KCL; inconsistent constants must remain nonpublishable |
| Gauge freedom | Solvers agree on identified voltage differences, objective and derived covariance; individual conductor voltages need not agree |
| Three-winding transformer | An independent OpenDSS deck supplies voltages and winding-port currents/powers for a grounded single-phase WYE/WYE/WYE case; excitation and artificial anti-float shunts are disabled explicitly |

The existing OpenDSS component/feeder matrix covers the other finite transformer
families. These tests do not establish all configurations of each family or
global uniqueness of a nonlinear fit. A common bug can survive agreement among
solvers sharing one evaluator; independent circuit equations and external
simulations are essential.

Use the same initial state, residual weights, injection paths, priors and
stopping tolerances when comparing algorithms. The nonlinear comparison uses
`optimality_tolerance=1e-6` explicitly: dense/sparse currently default to `1e-8`
and shared Ipopt to `1e-6`. An overstrict request can correctly return a
nonpublishable stall near floating-point resolution even when voltages and
objectives already agree. Do not loosen a test merely to hide a physical
constraint error.

For underobserved problems, compare predicted measurements and identifiable
derived quantities, not arbitrary gauge choices. The tests verify that an
independent stochastic prior resolves the gauge while preserving the uncertainty
of measured voltage differences. An absolute voltage covariance that depends
on an unresolved gauge must be rejected.

Remaining quality work should add more independently specified transformer and
neutral configurations, nonlinear uncertainty-coverage studies, realistic sparse
telemetry and dropout sequences, and multi-start/conditioning benchmarks. A
missing physical state or equipment equation needs a model extension; additional
unit tests cannot manufacture that support.

## Priorities for further work

1. An augmented voltage/current formulation for ideal transformers, bond/source
   reaction telemetry, and floating phase-to-neutral source laws.
2. Correlated measurement covariance, robust losses and statistically calibrated
   bad-data detection/identification. Raw standardized residuals alone do not
   identify bad data.
3. Timestamp alignment, missing-channel management and explicit temporal process
   models. Sequential warm starts and priors are not a Kalman filter.
4. Joint parameter/topology hypotheses, meter bias calibration and reproducible
   multistart studies. Local rank and convergence do not prove global uniqueness.
5. Sparse rank/covariance factorizations, symbolic reuse and large-feeder
   benchmarks, followed by an import/export contract for meter identity, units,
   direction, references and covariance provenance.

Always inspect `structure.preflight`, `solve_status(result).publishable` and
`observability_diagnostics`. A converged but underobserved result is feasible;
it does not identify every voltage direction. A regularized continuation stage
is not a publishable solution of the physical device law.
