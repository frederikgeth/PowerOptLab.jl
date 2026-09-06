# Generator circuit and numerical validation

The implementation's equations should be checked through an independent
calculation, and its local numerical properties should be measured explicitly.
This reproducible study does both on a small unbalanced feeder. It does not claim
speed superiority, global optimality, dynamical stability, or reliability on a
population of distribution networks.

Run `scripts/generalized_generators/validation.jl` from the project environment:

```sh
julia --project=. scripts/generalized_generators/validation.jl
```

The script checks 26 component/control cases in three coordinate settings: SI,
per-unit with ``S_b=10^4`` VA, and per-unit with ``S_b=10^8`` VA. These 78 solves
cover all five voltage laws for wye, delta and grounded sources. Eleven physical
circuits have independent complex linear-system references, repeated in all
three coordinate settings for 33 oracle comparisons. The remaining nonlinear
control cases check solver feasibility, power accounting and equality rank.

## 1. A separate circuit calculation

The upstream voltage is prescribed, loads are removed, and the native four-wire
feeder has a diagonal conductor impedance ``L``. Generator injection gives

```math
U=U^{\mathrm{grid}}+LJ.
```

For a generalized generator, substitute ``J=CI`` and ``E=C^{\mathrm T}U+ZI``:

```math
\left(Z+C^{\mathrm T}LC\right)I=E-C^{\mathrm T}U^{\mathrm{grid}}.
```

The reference script solves this complex linear system directly, without JuMP or
the stamped component expressions. Wye and delta use their different incidence
matrices. The tested matrices include zero impedance, unequal diagonals, mutual
coupling, and a singular nonzero primitive. For an exact ideal delta loop, the
reference uses a pseudoinverse and selects zero winding circulation; the OPF
receives the same explicit zero-sequence winding-current restriction. The
comparison therefore does not confuse nonunique ideal-loop currents with error.

For a dedicated source, let ``M=L+Z^c``, keep the four conductor currents ``J``,
and use a zero upstream neutral potential. Its reference equations are

```math
\left(M_{1:3,:}-\mathbf1_3 M_{n,:}\right)J
=E-U^{\mathrm{grid}}_{1:3},
```

```math
\begin{cases}
\mathbf1_4^{\mathrm T}J=0,&\text{open earth},\\
\left(M_{n,:}+z_g\mathbf1_4^{\mathrm T}\right)J=0,&\text{finite or ideal earth}.
\end{cases}
```

The source reference includes phase-neutral mutual impedance and separately
checks open, finite, and ideal grounding. Terminal voltages, winding currents
and conductor currents are compared in physical units. Assertions require
errors below ``10^{-5}`` V/A and active-power accounting error below ``10^{-4}`` W.
These are test tolerances, not physical uncertainty estimates.

## 2. Equality rank and coordinate scaling

At each published solution, the script differentiates every affine/quadratic
JuMP equality, including fixed-variable rows, using its polynomial coefficients.
It scales each nonzero Jacobian row to unit Euclidean norm and computes singular
values. A singular value greater than ``10^{-9}`` counts toward rank. Every
sampled model must have full **row** rank under that convention.

Row rank checks for locally dependent equality constraints; it does not establish
uniqueness. Several controls deliberately leave free dispatch directions, so
there can be more variables than equality rows. The check also omits active
inequalities, second-order optimality, continuation and different operating
points. Those omissions prevent interpreting it as a general LICQ or stability
certificate. In particular, hard capability boundaries need their own active-set
analysis.

The printed smallest singular value is coordinate-dependent even after row
normalization. A very large power base can separate current and voltage column
scales and reduce it. Compare the SI operating point across bases to establish
physical invariance, then choose practical voltage/current bases and inspect
solver convergence; the smallest number alone does not rank physical models.

```@example generator_validation
using PowerOptLab
include(joinpath(pkgdir(PowerOptLab), "scripts", "generalized_generators", "validation.jl"))
checks = generator_validation_study()
oracles = filter(r -> isfinite(r.voltage_error), checks)
(solves=length(checks), oracle_comparisons=length(oracles),
 max_voltage_error_V=maximum(r.voltage_error for r in oracles),
 max_current_error_A=maximum(r.current_error for r in oracles),
 all_equalities_independent=all(r.rank == r.equalities for r in checks))
```

The standalone script prints variable counts, equality ranks, smallest singular
values and build/solve times for each case. It warms one ordinary case; later
specializations can still compile. Timings include process and machine effects
and are not a controlled benchmark against another formulation. Every case also
checks that no internal voltage variables were introduced for its series drop.

## 3. What this adds to the unit tests

The unit suite exercises the capabilities under unbalance, deliberately excludes
known physical points with individual ratings, compares overlapping native
BMOPFTools ownership/cost behavior, and checks data round trips. This script adds
an independently assembled network oracle and explicit sampled equality-rank
checks. The [trade-off tutorial](generalized_generator_tradeoffs.md) adds loaded
network comparisons and demonstrates hidden delta circulation.

A publication-level performance/reliability study remains separate work: select
representative feeders and operating conditions, specify reference algorithms
and tolerances, measure repeated timings after full warmup, and report failure
rates and active-set behavior. No new OpenDSS or other external-software parity
claim is made for this generator extension by these tests.
