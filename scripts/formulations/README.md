# PWL formulation experiments

See `docs/src/formulations/index.md` for the mathematical contracts and executable
tutorial. The production inverter-control defaults are unchanged.

```sh
julia --project=. scripts/instantiate_pinned.jl
julia --project=. scripts/formulations/comparison.jl
julia scripts/formulations/setup.jl /tmp/pol-pwl-env
POL_FORMULATION_RESULTS=/tmp/pol-pwl-results.toml \
  julia --project=/tmp/pol-pwl-env scripts/formulations/optional_tests.jl
```

The optional destination must be empty. Use a manifest resolved by the Julia
version running the setup: a Julia 1.12 manifest is not necessarily consumable by
Julia 1.10. The optional setup preserves the BMOPFTools commit, uses DiffOpt 0.6.2
for bridge compatibility, and records the resolved solver stack in the results.
It neither changes the production environment nor adds a PowerOptLab solve loop.

The experiment solves `V = Vs + I` and `I = f(V)` on a one-ohm resistive feeder.
The canonical controller has knots `[220,240,250,270]` V and values `[20,20,0,0]` A.
Segment enumeration supplies an independent analytic equilibrium. The core suite
compares softplus and local C2 at equal 0.01 A uniform error budgets. The optional
suite covers 20 exact graphs, 20 complementarity encodings, and a hull witness.

`reference_results.toml` is a measured macOS / Julia 1.12.6 snapshot, with source
fingerprints and package versions. It is evidence, not a golden floating-point
fixture or a platform-independent solver guarantee. First-call elapsed times can
include compilation and should not be interpreted as performance benchmarks.

The exact graph / HiGHS runs all solve strictly and meet microunit physical
residual checks. CCOpt 0.1.0 returns `ALMOST_LOCALLY_SOLVED` in all 20 cases, and
only eight meet the microunit canonical equation checks. None is promoted to
strict success. The normalized-coordinate zero-current breakpoint reaches about
0.0376 A controller error, exceeding the smooth methods' 0.01 A budget.

Consequently, optional CI tests exact-graph correctness, CCOpt integration and the
accurate subset, while retaining all broader CCOpt runs as characterization. A
passing optional job is not a claim that every MPCC case converges or meets the
physical accuracy target. The hull intentionally fails canonical controller
agreement (16 A at 246 V, where the controller commands 8 A).

Follow-up research: physical normalization of complementarity variables, external
solver termination and relaxation controls, stationarity assessment, and realistic
three-phase/fleet studies. Those should be compared through the same physical
oracles without weakening solver outcome classifications.
