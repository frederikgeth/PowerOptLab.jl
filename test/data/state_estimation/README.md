# Finite transformer and regulator SE validation

Run `julia --project -e 'using Pkg; Pkg.test()'` from the repository root.
The matrix is in `test/state_estimation_network_tests.jl`; fixtures, including
its generated SE-IEEE13-finite OpenDSS deck, are in
`test/state_estimation_feeder_fixtures.jl`. See `license.md` for provenance.

Both Gauss–Newton solvers start from the load-free state, using independently
constructed OpenDSS voltage and load-injection readings. JuMP/Ipopt starts
from the same initial guess and checks the fitted state and objective using
independent differentiation, but shares the compiled admittance matrices.

## Reference run

Julia 1.12.6, OpenDSSDirect 0.9.9, BMOPFTools at
`56e3b2f0261488acd30f1c98ad812c4eff4611f1`, 15 October 2026.
Maximum dense-estimator voltage-component error relative to the oracle phasor
magnitude (100 V floor), in percent. Noise is deterministic with seed 349,
independent Gaussian perturbations at the specified measurement standard
deviations. This is one realization, not a coverage study.

| Case | Clean (%) | Noisy (%) |
|---|---:|---:|
| single_phase | 0.0551539 | 0.122061 |
| single_phase_tap | 0.0544997 | 0.120528 |
| wye_delta | 0.0597075 | 0.313663 |
| delta_wye | 0.0648071 | 0.166337 |
| regulator | 0.0380304 | 0.131529 |
| open_delta | 0.00709892 | 0.0838431 |
| center_tap | 2.49465e-06 | 0.00971615 |
| SE-IEEE13-finite | 0.054711 | 0.142925 |

Acceptance bounds are 0.2% clean and 0.6% noisy. The nonzero clean errors in
several fixtures mean the models agree to an engineering tolerance, not to
machine precision. Both estimator objectives agree with JuMP within absolute
and relative tolerances of 1e-5; states agree within 0.01 V. Both solvers must
report unique convergence, satisfy exact constraints and yield finite selected
covariance. The separate package tests exercise underobserved cases.

SE-IEEE13-finite is a modified benchmark with finite leakage/losses, constant
power loads and fixed unequal phase taps, not the original IEEE case. The
open-delta component includes 1 kOhm phase-earth shunts in both models for a
physical common-mode path. Its shared phase B uses one OpenDSS node on both
sides, preserving the exact copper bond; secondary A/C are src.11/src.13, mapped
to reg.1/reg.3. Component neutrals use their declared grounding.

The separate `state_estimation_compatibility_tests.jl` suite exercises n-winding
and other finite transformer-port telemetry across dense, sparse and shared
Ipopt backends, plus native WLS comparisons. This OpenDSS matrix does not
validate n-winding transformers, automatic regulator
controls, unknown taps, ideal transformers, transformer branch telemetry,
large-feeder scaling, bad-data processing or statistical confidence coverage.

## Three-winding independent port oracle

`three_winding.dss` is a synthetic grounded single-phase WYE/WYE/WYE transformer
with 2400/1200/600 V windings, 100 kVA per winding, and two constant-power loads.
Winding resistances are 1/.25/.0625 ohm. Pair reactances X12/X13/X23 are 2/3/2
ohm on winding 1's 57.6-ohm base; the deck expresses those as percentages.
Grounds use node 0 exactly; excitation losses and OpenDSS artificial anti-float
shunts (`ppm_antifloat`) are disabled to match the stated electrical model.

`state_estimation_consistency_tests.jl` reads conductor voltages and transformer
port currents/powers directly from OpenDSS, rather than constructing those meter
readings from the compiled Ybus. It checks all three shared solvers and native
SI/per-unit WLS. This deck covers one three-winding configuration, not arbitrary
winding counts, connections or ungrounded windings.

## Mixed winding and grounding circuit oracles

`state_estimation_transformer_grounding_tests.jl` builds separate synthetic
three-phase WYE/DELTA/WYE, DELTA/WYE/DELTA and WYE/WYE/DELTA cases. The reference
uses three physical star leakage impedances, referred to the high winding's
coil base: `1+0.8j`, `1+1.2j`, `1+2.2j` ohm. Coil nominal voltages are
2400/1200/600 V. Per core it solves for a hidden common emf from
`sum((u[k]/N[k]-E)/z[k])=0`, then obtains winding and conductor currents by
Ohm's law and the explicit WYE/DELTA connections. It never calls the package's
ZB, Yprim or Ybus code to generate expected readings.

Unbalanced current loads return to a WYE neutral or another DELTA phase.
Cases have floating windings, ideal neutral bonds, or explicit complex
electrodes with a separate phase-earth path. Both DELTA rotations and a change
of reference winding are checked. Tests compare phase/neutral port currents,
voltage differences, phase power, finite-electrode current/power and neutral
earth-power where earth reference information exists. Floating tests retain
the physical common-mode gauges and reject absolute state covariance; exact
affine laws give zero covariance for their fixed differential voltages.

OpenDSS independently checks every transformer terminal-current equation at
specified conductor voltages, avoiding arbitrary external load-flow gauges.
The inline DSS commands disable artificial anti-float shunts and test both zero
excitation and a winding-2 coil admittance of `0.0001-0.0002j` S. The independent
star circuit adds this branch explicitly; OpenDSS supplies a separate check of
its loss and reactive-current base conversion. Unused DELTA fourth conductors
are mapped explicitly. Native SI/per-unit WLS
comparisons use the same P/Q and absolute rectangular-voltage readings as all
three shared solvers, with those absolute readings supplying explicit
stochastic reference information.
