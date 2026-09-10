# EV model verification

The production NLP path has no optional dependencies:

```sh
julia --project=. examples/ev_charging.jl
julia --project=. examples/ev_workplace.jl
julia --project=. -e 'using Test, PowerOptLab, JuMP, Ipopt; include("test/fixtures.jl"); include("test/ev_tests.jl"); include("test/multiperiod_tests.jl"); include("test/evse_tests.jl")'
```

`test/runtests.jl` also includes these regressions. They check analytical energy
accounting, disconnection under cancelling phase powers, SI/per-unit current
limits, true replenishment, rating and permission intersections, assignment
conflicts, and acceptance-envelope refinement. The source fixture uses an ideal
230 V supply where an exact energy-price answer is useful; separate regressions
use unbalanced feeders. Source prices are currency/kWh, not currency/Wh.

For the optional complete-AC-network MPCC experiment:

```sh
julia scripts/formulations/setup.jl /tmp/pol-ev-optional
POL_EV_RESULTS=/tmp/pol-ev-results.toml \
  julia --project=/tmp/pol-ev-optional scripts/ev/ccopt.jl
```

The environment path must be new or empty. The setup preserves the production
BMOPFTools revision and installs the tested optional solver versions. The
experiment writes package versions, source SHA-256 fingerprints, the upstream
revision, solver settings, statuses and physical residuals. CI retains this
artifact alongside the other optional formulation results.
`reference_results.toml` records one local run, including unsuccessful raw
candidates; its source fingerprints identify the exact implementation tested.
It is evidence from this machine, not a golden cross-platform status baseline.

The replenished cycle requires a successful solve and physical residual checks.
The nonconcave acceptance case also requires a successful solve, agreement with
its analytical 3600 W optimum, and less than 0.1 W violation of the original
envelope. It uses `tol=1e-7`, recorded in that run's `nlp_tolerance`; the other
cases use the shared `tol=1e-8`. This is case-specific numerical configuration,
not a guarantee that one tolerance works for every network or curve.
The biactive full-battery case is a solver characterization: the exact solution
is zero charge/discharge, but local numerical failure must not be interpreted as
a proof of infeasibility. A rejected raw candidate is recorded separately and
never published as dispatch. No independent MPCC stationarity certificate or
global optimality claim is made.

The bounded zero-dispatch sensitivity study is separate:

```sh
POL_EV_ZERO_RESULTS=/tmp/pol-ev-zero-results.toml \
  julia --project=/tmp/pol-ev-optional scripts/ev/zero_dispatch.jl
```

It compares initialization, network bases, continuation, and bound relaxation
in SI and per-unit coordinates. An analytically reduced Ipopt reference is
labeled separately from the unreduced CCOpt experiments. Local observations
are saved in `zero_reference_results.toml`; source fingerprints identify the
tested code. See `docs/src/ev/ccopt_investigation.md` for the derivation and
limitations. A small raw residual never overrides a rejected solver status.

`test/ev_workplace_tests.jl` verifies the feeder tutorial's equipment-dependent
cost, deadline bound, SI/per-unit agreement, and physical residuals.

Documentation: `docs/src/ev/` contains scientific definitions, references,
executable tutorials, and an explicit roadmap for mobility, shared DC cabinets,
thermal/aging physics and stochastic control. Those future features should not
be inferred from the current fixed-assignment AC API.
