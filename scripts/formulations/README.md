# PWL formulation experiments

See [the configurable experiment tutorial](../../docs/src/formulations/experiments.md)
to supply your own cases, configurations, methods, metrics and acceptance criteria.
The scripts provide reproducible scalar and control-lowering examples. The
[mathematical guide](../../docs/src/formulations/theory.md) explains the contracts.

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
Setup resolves and develops the complete optional environment before explicitly
precompiling it once. This avoids incremental precompilation requesting newly
resolved transitive packages by names outside the active project's dependencies.

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

`ComplementarityGraph(scale=...)` exposes physical hinge normalization
independently of model coordinates. The optional script also runs three focused
cases through the configurable API and exports them to a separate
`*-configurable.toml` bundle. CI uploads both bundles. No improvement in MPCC
reliability is inferred from these integration checks.

The controller adapter and configurable runner enable studies of external solver
controls, stationarity, normalization and realistic three-phase/fleet cases.
Researchers choose their own physical oracles and acceptance policies while raw
solver outcomes remain available.

## Control intent and bounded relations

```sh
julia --project=. scripts/formulations/control_lowering.jl
julia --project=. scripts/formulations/bounded_relations.jl
POL_FORMULATION_RESULTS=/tmp/pol-control-results.toml \
  julia --project=/tmp/pol-pwl-env scripts/formulations/control_lowering_optional.jl
```

The optional control tutorial compares smooth MadNLP, exact HiGHS, hull HiGHS
and CCOpt paths. The bounded example varies voltage domains, equality versus
upper-limit use, and specialization. Unsupported automatic specializations are
recorded explicitly. Optional exports include `*-bounded.toml` and
`*-bounded-mpcc.toml`; graph and relation observations keep their own residual
semantics. See [bounds and relations](../../docs/src/formulations/bounds_and_relations.md)
for the geometry and reproduction instructions for the figures.

## Offline inverter current allocation

```sh
julia --project=. scripts/formulations/inverter_current_projection.jl /tmp/current_projection.json
```

This fixed-boundary physical reference tests a feasible current request, a
rail-infeasible request with a feasible alternative, and energized zero-PV loss
import on the same synthetic plant. It includes source bounds, independent
per-bound checks and source fingerprints. The measured
[`inverter_current_projection_reference.json`](inverter_current_projection_reference.json)
snapshot accepts all three allocations. It is a local nonconvex reference,
not a global nearest-point certificate, autonomous VV+VW law or EMT simulation.
The hardware is a regression fixture and does not replace the Australian
case-study parameters. See the
[operating-constraint guide](../../docs/src/components/inverter_operating_constraints.md).
`python3 scripts/formulations/plot_current_projection.py` uses Matplotlib/NumPy
to regenerate the standalone capability-slice and full-cycle voltage figure
from that saved reference.

## Local reactive-priority allocation

```sh
julia --project=. scripts/formulations/inverter_service_allocation.jl /tmp/service_allocation.json
```

`allocate_inverter_services` uses staged Q tracking, P export and requested I₂
support on the full physical plant. `evaluate_voltage_aware` supplies exact
mean/worst/sequence VV+VW requests from separate sensed and physical voltage
phasors. The [measured reference](inverter_service_allocation_reference.json)
contains six accepted local commands and one explicitly unresolved strict policy
with a feasible unrestricted-plant diagnostic. Four cases use the unchanged
Australian hardware; the 620 V synthetic rail fixture is separate. Every stage
has solver provenance and independently checked physical/policy bounds.
These are local control-reference results, not global lexicographic optima,
network equilibria, measured response times or EMT.
