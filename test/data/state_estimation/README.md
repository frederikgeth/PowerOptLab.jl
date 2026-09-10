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
`72e6cec22a66cf376c37ec4d64aef350b9f1100d`, 10 September 2026.
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
| open_delta | 0.0055748 | 0.180888 |
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
physical common-mode path. Component neutrals use their declared grounding.

The matrix does not validate n-winding transformers, automatic regulator
controls, unknown taps, ideal transformers, transformer branch telemetry,
large-feeder scaling, bad-data processing or statistical confidence coverage.
