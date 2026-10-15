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
`solve_state_estimation(structure, parameters, x0)` defaults to the sparse solver.
The explicit solver names are preferred when comparing algorithms. All three
consume the same residuals, exact constraints, first derivatives, source phasors,
customer identities and priors. They report the same `objective = sum(abs2,r)`
and `dof = residual_count - observable_dimension - nuisance_dimension`.
Internally they may minimize half that sum; the public statistic is always J.

Ipopt receives sparse derivative patterns compiled from dependencies and the
Gauss–Newton curvature `H' * H`. This approximation omits residual and nonlinear
constraint second derivatives; it is not the exact Lagrangian Hessian. Residual,
Jacobian and constraint evaluations are cached separately at each iterate.
Every available primal point is checked for physical feasibility, reduced
stationarity, rank and negative curvature, regardless of the optimizer's
termination label. Its original termination and primal labels remain separate
in diagnostics. Finite local checks do not establish a global optimum.

For floating networks, some KCL equations can be identities of the others.
The shared Ipopt backend selects an independent equation basis using the full
affine network/device coefficients plus independent nonlinear branch-current
incidences. Sparse QR verifies those coefficient identities near roundoff. A
bounded, synchronized cache reuses the proof when only readings, uncertainties
or priors change; physical source/device changes require a fresh proof. It does not
drop a physical equation just because its Jacobian is zero at the starting
point. Original equations remain in `evaluation.constraints` and every final
physical feasibility and rank check. `solve_diagnostics` reports both equation
counts and the physical row indices sent to the optimizer. Feasibility scaling
also includes amplification from the independent rows back to all physical
currents, including transformer turns ratios. Inconsistent constant
terms and weak grounding information are retained.

Compilation also omits return-current identities proved directly from explicit
n-winding coil connections. Summing all phase and neutral terminal currents
is exactly zero for a complete isolated winding, but summing rounded Ybus rows
can manufacture a tiny coefficient. The structural proof applies only when all
coil endpoints have equal component weights and no other touched element or
device breaks the identity. Unknown element groups fail closed. Any finite electrode, however weak, retains its
physical return equation. This conservative reduction leaves other element
families to the numerical coefficient checks above.

The original `solve_state_estimation(net, measurements)` remains a separate
native-engine compatibility API. It accepts `Measurement` only, requires paired
customer P/Q coverage or explicit zero injection, and does not accept the newer
measurement types or exact-device objects. Its source neutral must be explicitly
grounded: otherwise the native engine would introduce an implicit earth bond.
The native current-variable model supports overlapping delta/wye paths and
merged customers with distinct returns. Only conflicting returns for the same
original customer are rejected. The native observability diagnostic uses
ForwardDiff, with analytic current elimination for overlapping paths.

A second contract difference remains explicit: native KCL imposes zero injection
at every non-source conductor without an injection variable, including neutrals.
It warns about undeclared conductors. The compiled model leaves such injections
free and reports `:undeclared_injection` in preflight. Declare zero injections
explicitly when comparing these APIs; missing telemetry alone is not a physical
zero-injection declaration. Use the compiled API for finite source-neutral
electrodes or incomplete conductor coverage.

## Injection states and customer identities

`compile_state_estimator(...; formulation=:auto)` retains voltage elimination
for simple disjoint customer paths. It adds one complex current per original
unknown customer for merged customers, overlapping returns or mixtures of known
and unknown devices on shared conductors. The latter include a known PV/charger
on one phase and an unknown load on another phase sharing a neutral. Distinct
original terminals preserve customer identity even when ideal switches merge
their electrical nodes.

Select `formulation=:current` to use explicit unknown currents throughout, or
`:voltage` to require the optimized elimination contract. Voltage entries remain
`[real(V_free); imag(V_free)]`; current entries follow them. Use
`injection_state_indices(s)` rather than guessing their ordering. Branch current
is positive into the network; KCL is `YV + I_exact - A*I_customer = 0`.
An explicit `injection_connections` path denotes an additional unknown customer,
even when a known device has the same endpoints; its P/Q meters describe that
unknown customer. Without an explicit unknown declaration, matching telemetry
meters the exact device. Individual power is `U*conj(I_customer)`, so its evaluator remains defined at
zero connection voltage with finite currents. In contrast, the optional merged
power profile requires nonzero voltage and retains its explicit domain check.

With explicit currents, missing customer channels may leave individual currents
unobserved while total voltage remains identifiable. Compare identified quantities
and selective covariance, not arbitrary unobserved current coordinates.

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
    port=2, terminal="1", reference="n",
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
Transformer IDs are looked up across subtypes; `variant` is optional and only
needed to disambiguate a duplicated ID. Ordinary transformer ports are 1/from and 2/to; center-tap port 2 includes both
legs and the center conductor. An n-winding port is its winding-list index.

The open-delta shared phase has an ideal copper bond **in addition to** the finite
winding primitive. Compilation identifies its two nodes exactly, matching the
native engine. The bond is assembled upstream in BMOPFTools; PowerOptLab no
longer inserts synthetic switches or maintains a separate admittance patch. Shared-phase
port current includes an ideal-wire reaction and cannot be obtained from the
winding primitive alone, so its element telemetry is rejected. Meter adjacent
finite branches instead. N-winding telemetry is also rejected where different
ports share the same original conductor and the exporter has already summed
those currents. Other port currents use the original primitive order, including
the distinct autotransformer ordering.

This is not support for every transformer model. Ideal/zero-leakage transformers,
unrepresented vector groups/phase shifts, and automatic tap control remain
outside the finite-primitive contract. Fixed supported taps are network parameters;
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
| Mixed winding connections and grounding | Independent hidden-emf/star-leakage equations cover WYE/DELTA/WYE, DELTA/WYE/DELTA and WYE/WYE/DELTA, both DELTA rotations, two winding reference orders, unbalanced loads and floating/solid/complex-electrode grounding; OpenDSS terminal equations independently check every phase and neutral port current, with zero and nonzero excitation |
| Loaded 5/20/40-bus four-wire feeders | Independent series Ohm's law and electrode KCL generate unbalanced operating points; noisy P/Q/magnitude fits compare dense, sparse, shared Ipopt and native WLS states, J and degrees of freedom |
| Nonlinear statistical calibration | 500 seeded independent draws with unequal variances, repeated readings and merged customers check E[J] against dof and empirical voltage covariance against the local calculation; explicit and profiled formulations agree |
| Floating winding rank and covariance | A completely uninformative projected Jacobian must not turn roundoff into observability; all shared solvers and the native rank diagnostic preserve gauges, reject absolute-voltage covariance, and retain weak genuine reference information |

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

Projection roundoff is assessed against the original measurement Jacobian as
well as the reduced map. This matters when exact device laws fix every
differential voltage and only unobserved common-mode directions remain:
`H*Z` is then zero mathematically, but its numerical entries need not be zero.
The same rank convention governs shared solver statuses, diagnostics,
unobservable directions and covariance. A quantity fixed by exact equations
has zero local covariance. A weak real reference or prior remains information
when it exceeds the projection's numerical resolution.

The mixed-winding statistical comparisons use explicit absolute-voltage
readings when comparing native SI/per-unit WLS with the shared solvers. Those
readings supply stochastic gauge information; no implicit neutral bond is
added. Separate floating-circuit tests omit that information and compare only
identified quantities. External component checks specify excitation explicitly, disable
OpenDSS anti-float shunts and evaluate terminal equations at specified voltages,
so an external solver's arbitrary gauge choice does not become ground truth.

`ConstantCurrentDevice` defaults to a fixed earth-frame phasor. Set
`angle_reference=:voltage` for fixed magnitude and power factor relative to the
connection voltage. Its angle is undefined at zero voltage unless the coefficient
is zero. `ZIPDevice` retains a fixed-phasor current term; compose exact models
when a voltage-following current term is required.

Remaining quality work includes nonlinear core/saturation models, further
transformer configurations, broader uncertainty-coverage studies, realistic sparse
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

## Statistical interpretation and numerical limits

J includes measurement and independent prior residuals; diagnostics also report
the two components separately. In the optional voltage profile, fitted customer
powers are nuisance parameters and are subtracted from dof. In the current
formulation their estimated directions already enter `rank(HZ)`, so they must
not be subtracted twice. A chi-square interpretation is local/asymptotic and
requires the stated independent Gaussian errors, correct physical model and
properly counted independent prior information. Reusing correlated past data
as an independent prior invalidates that interpretation.

The 500-draw regression covers one low-noise nonlinear circuit, not general
calibration or bad-data processing. Final rank/covariance diagnostics still use
dense SVD. Sparse derivatives and cached equation reduction do not establish
large-network scalability; the demonstrated feeder range is 5–40 downstream
buses.
