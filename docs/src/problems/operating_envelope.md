# Dynamic operating envelopes

> **Kind:** Problem specification · **Maturity:** research prototype · **Direction:** forward · **Temporal:** per-interval

[`solve_operating_envelope`](@ref) allocates an active-power import or export
capacity to each participating connection while retaining the nonlinear,
unbalanced four-wire network model and all operational limits declared in the
BMOPFTools case.

The DOE layer does **not** treat ordinary LV reactive power as freely
dispatchable:

- loads retain the known P/Q values in each DSSE-derived snapshot;
- PV and battery connection points can bind to an existing BMOPFTools `ibr`, so
  its mandatory Volt-VAr/Volt-Watt or fixed-power-factor law remains enforced;
- a STATCOM is a separate network device and is dispatchable within its own
  converter limits.

All capacities returned by this API are positive watts. `direction=:export`
means injection into the network and `direction=:import` means withdrawal.

## Security semantics

There are two deliberately distinct modes:

| Mode | AC dispatches included in the optimization | Claim |
|---|---|---|
| `security=:bound_point` | all participants simultaneously at their allocated bound | Feasible upper-bound allocation only |
| `security=:corners` | every zero/full-utilisation corner of the allocated box | Local AC feasibility at every represented corner |

The corner mode creates ``2^N`` network contexts per forecast scenario and is
therefore limited by `max_exact_corners` (10 participants by default). It is a
much stronger test than the bound point, but it is **not** labelled a global
robust certificate: the AC feasible set is non-convex and an interior hole may
exist even when every corner is feasible. The exact scope and solver primal
status are recorded in `result.diagnostics`.

Diagnostics also report worst tested voltage, ampacity and negative-sequence
margins, their network locations, and constraints within reporting tolerance of
binding. For `security=:corners`, these are aggregated across every scenario and
corner rather than only the displayed representative snapshot.

If an interval has no feasible primal point, its capacities and total are `NaN`.
The package never publishes values from an infeasible solver iterate.

## Basic export envelope

The lightweight connection port is useful for teaching cases. It injects active
power at aggregate unity power factor.

```julia
using PowerOptLab

cps = [ConnectionPoint(id="der1", bus="bus1", export_max=10e3),
       ConnectionPoint(id="der2", bus="bus2", export_max=10e3)]

# One network Dict is one interval. A Vector{Dict} is a time series.
r = solve_operating_envelope(net, cps;
    direction=:export,
    fairness=:equal,
    security=:bound_point)

r.envelope["der1"]       # positive W
r.total_capacity
r.diagnostics[1]["security_scope"]
```

Use `security=:corners` when the connection may operate anywhere between zero
and its advertised bound and the participant count is small enough for exact
corner enumeration.

## Forecast and model scenarios

A vector of vectors groups alternative forecasts or network models by interval.
One envelope is shared by all scenarios in an interval.

```julia
# Two time intervals, each with three load/source/network scenarios.
scenario_nets = [
    [net_t1_central, net_t1_low, net_t1_high],
    [net_t2_central, net_t2_low, net_t2_high],
]

r = solve_operating_envelope(scenario_nets, cps; security=:corners)
r.diagnostics[1]["scenario_count"]  # 3
```

Scenarios may represent demand/PV forecast error, source-voltage uncertainty,
topology alternatives, or candidate impedance models. In particular, candidates
or profile intervals produced by [`solve_inverse_carson`](@ref) can be
materialized into alternative network scenarios.

## PV and batteries with mandatory Q-V control

For realistic DER behaviour, place the converter in the network's `ibr` block
and attach the required BMOPFTools `control_profile`. Bind the connection point
to that IBR:

```julia
cp = ConnectionPoint(id="customer_17", bus="lv17", ibr_id="pv17",
                     export_max=10e3, import_max=5e3)

export_doe = solve_operating_envelope(net, [cp]; direction=:export)
import_doe = solve_operating_envelope(net, [cp]; direction=:import)
```

The active-power equality used by the DOE is added to the existing IBR model.
Its per-phase topology, current limit, apparent-power circle, DC coupling and
Volt-VAr/Volt-Watt equality are not replaced. There is intentionally no separate
Q-envelope decision.

## With and without a STATCOM

STATCOMs are normal BMOPFTools network devices, so the comparison is made by
solving two otherwise identical cases:

```julia
using BMOPFTools: add_statcom!, augment_case

without_statcom = deepcopy(net)
with_statcom = deepcopy(net)
add_statcom!(with_statcom, "lv17"; s_max=50e3)
with_statcom, _ = augment_case(with_statcom)

r0 = solve_operating_envelope(without_statcom, cps)
r1 = solve_operating_envelope(with_statcom, cps)

gain = r1.total_capacity[1] - r0.total_capacity[1]
q_statcom = r1.snapshots[1]["ibr"]["statcom_lv17"]
```

The default STATCOM is reactive-only. BMOPFTools also supports a four-wire,
DC-link-coupled D-STATCOM that circulates active power between phases while its
net active exchange remains zero; this can be important in resistive,
unbalanced LV feeders.

## Parameterized fairness

The legacy symbols remain available:

- `:equal` — equal absolute kW;
- `:sum` — maximum aggregate capacity;
- `:proportional` — proportional fairness.

Use [`FairnessPolicy`](@ref) for explicit policy design:

```julia
policy = FairnessPolicy(
    kind=:max_min,
    normalization=:capacity,
    weights=Dict("der1"=>1.0, "der2"=>1.5))

r = solve_operating_envelope(net, cps; fairness=policy)
```

Available allocation objectives are `:equal`, `:max_total`, `:proportional`,
`:alpha`, `:max_min`, and `:equal_curtailment`. Normalization can use:

| Normalization | Reference |
|---|---|
| `:none` | 1 W: absolute allocation/curtailment |
| `:capacity` | directional import/export nameplate |
| `:request` | `ConnectionPoint.requested` forecast/request |
| `:custom` | `ConnectionPoint.normalization` |

For example, equal allocation with `normalization=:capacity` gives equal
fractions of nameplate rather than equal kW. `:max_min` performs a second local
solve to maximize total capacity while retaining the best normalized minimum.
`kind=:alpha` exposes the standard alpha-fair family: alpha 0 is weighted sum,
alpha 1 is proportional fairness, and larger alpha increasingly emphasizes the
least-served participant.

`OperatingEnvelopeResult.fairness_metrics` records normalized allocations,
curtailment fractions, total capacity and Jain's index for every published
interval. To compare an explicit policy frontier, solve the same study under
multiple policies:

```julia
results = compare_operating_envelope_policies(net, cps,
    ["equal" => :equal,
     "efficient" => :sum,
     "proportional" => FairnessPolicy(kind=:proportional,
                                       normalization=:capacity)])
```

## Rolling fairness and publication

For operational issuance, `temporal_fairness=:cumulative_max_min` carries a
normalized service history between intervals and, at each new interval,
prioritises the least-served participant before maximizing the remaining
weighted allocation. It is a rolling policy—not a horizon-wide optimiser—and
therefore does not assume perfect future forecasts.

```julia
using Dates

r = solve_operating_envelope(nets, cps;
    fairness=FairnessPolicy(kind=:max_min, normalization=:capacity),
    temporal_fairness=:cumulative_max_min,
    fairness_history=Dict("der1" => 2.4), # prior normalized service
    temporal_dt_h=5 / 60,
    issued_at=DateTime(2026, 7, 16, 9),
    interval_seconds=300,
    validity_seconds=600)
```

`r.schedule` contains issue and validity times plus the publication source.
On an infeasible solve, `fallback=:missing` (default), `:zero`, or
`:last_feasible` controls what is published. A fallback is explicitly marked as
not freshly network-safe; use verification before relying on it.

## Verify an issued envelope

Verification fixes an issued allocation rather than re-optimising it. It can
check the simultaneous upper point, every box corner, or custom utilisation
points across all scenarios:

```julia
check = verify_operating_envelope(nets, cps, r; utilizations=:corners)
all(check.feasible) || @warn "issued DOE needs review" check.diagnostics
```

The check retains the network's prescribed Q-V controls and any STATCOM. Its
diagnostics report local nonlinear feasibility and inherited voltage, thermal
and negative-sequence margins at the tested points.

## Network constraints and use cases

The DOE adds no parallel approximation of network limits. It inherits the
constraints already present in the BMOPFTools case, including:

- phase-to-ground, phase-to-neutral and phase-to-phase voltage bounds;
- line, transformer, switch, neutral and converter current limits;
- positive-, negative- and zero-sequence voltage limits such as `vneg_max`;
- prescribed inverter controls and converter ratings;
- controllable network assets already represented by the OPF.

The tests include voltage-limited, thermally limited, negative-sequence-limited,
import, multi-scenario, prescribed-Q-V and STATCOM-assisted examples. They also
cover invalid inputs and infeasible baseline/corner cases.

## DSSE-to-DOE validation runner

The repository includes `scripts/validate_doe_from_dsse.jl` for feeder studies
that are too expensive or data-dependent for the unit-test suite. It keeps three
questions distinct: whether DSSE reconstructs the observed state, whether the
DOE is feasible under its own nonlinear model, and whether a separately fixed
AC power flow reproduces the issued upper-bound point.

Create a case-builder file that defines `doe_validation_case()` and run:

```sh
julia --project=. scripts/validate_doe_from_dsse.jl path/to/my_case.jl
```

The returned named tuple must contain:

- `physics_net`: a passive network suitable for `solve_state_estimation`;
- `operational_net`: the matching DSSE snapshot with known P/Q loads and
  controllable network assets;
- `measurements`: the DSSE `Measurement` vector;
- `connection_points`: DOE connection points bound to existing, single-phase
  IBRs (`ibr_id` is required for the independent power-flow check).

It may also provide `truth_net` (the measurement-generating operational state),
`with_statcom_net` (the otherwise-identical STATCOM case), and `doe_keywords`
(a `NamedTuple` forwarded to `solve_operating_envelope`). The runner reports
DSSE voltage error against the truth power flow, DOE capacity, DOE verification,
the maximum voltage difference between DOE and a separate fixed-setpoint
`solve_pf`, and STATCOM capacity gain where supplied. It deliberately does not
enter the unit suite: use it for real feeder exports, a representative time
series, and solver/runtime logging. For a controllable STATCOM or other flexible
asset, record and replay the controller setpoint selected at issuance; a plain
power flow may otherwise select another feasible reactive dispatch, which is an
operational-policy difference rather than a DOE-model comparison.

## Current limitations

- Ipopt returns local nonlinear solutions; diagnostics never imply global
  optimality or global robust feasibility.
- Exact corner enumeration scales exponentially.
- Rolling fairness is causal rather than globally horizon-optimal; temporal
  storage scheduling and forecast co-optimisation remain later work.
- The independent DSSE validation runner presently fixes single-phase IBR
  active-power setpoints. Multi-phase dispatch replay should preserve the DOE's
  per-phase allocation before being treated as an independent validation.
- Harmonic RMS/THD constraints require a harmonic network model and are outside
  the current fundamental-frequency prototype.
