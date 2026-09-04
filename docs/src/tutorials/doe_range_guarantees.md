# What does a dynamic operating envelope guarantee?

> **Claim:** a feasible simultaneous upper point need not imply that partial
> utilization of the advertised box is feasible.
>
> **Model:** synthetic, unbalanced nonlinear AC, fundamental-frequency RMS
>
> **Expected runtime:** under one minute after Julia compilation
>
> **Data license:** no external data; the complete synthetic case is below

This example isolates a range-safety failure. Three equal single-phase export
limits are connected to a balanced three-phase feeder with a tight
negative-sequence-voltage limit. Simultaneous equal export remains balanced,
but unequal customer utilization creates voltage unbalance.

## Build the complete case

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

net = parse_bmopf(raw"""
{
  "bus": {
    "src": {
      "terminal_names": ["1", "2", "3", "n"],
      "perfectly_grounded_terminals": ["n"]
    },
    "b1": {
      "terminal_names": ["1", "2", "3", "n"],
      "perfectly_grounded_terminals": ["n"],
      "v_min": [200.0, 200.0, 200.0],
      "v_max": [260.0, 260.0, 260.0],
      "vneg_max": 1.0
    }
  },
  "voltage_source": {
    "vs": {
      "bus": "src",
      "terminal_map": ["1", "2", "3"],
      "v_magnitude": [230.0, 230.0, 230.0],
      "v_angle": [0.0, -2.0943951024, 2.0943951024]
    }
  },
  "linecode": {
    "lc": {
      "R_series_1_1": 0.4,
      "R_series_2_2": 0.4,
      "R_series_3_3": 0.4,
      "R_series_4_4": 0.4
    }
  },
  "line": {
    "l1": {
      "bus_from": "src",
      "bus_to": "b1",
      "terminal_map_from": ["1", "2", "3", "n"],
      "terminal_map_to": ["1", "2", "3", "n"],
      "linecode": "lc",
      "length": 1.0
    }
  }
}
"""; from_string=true)

cps = [
    ConnectionPoint(
        id="phase_$phase",
        bus="b1",
        phase_terminals=[string(phase)],
        neutral="n",
        export_max=20e3,
    ) for phase in 1:3
]
```

The case is intentionally small enough to inspect. It is a methodological
fixture, not a calibrated distribution feeder.

## Allocate only at the simultaneous bound

```julia
bound = solve_operating_envelope(
    net, cps;
    security=:bound_point,
    control_policy=PerfectRecourse())

bound.total_capacity
bound.diagnostics[1]["security_scope"]

bound_check = verify_operating_envelope(
    net, cps, bound;
    utilizations=:bound_point,
    control_policy=PerfectRecourse())

@assert bound_check.feasible == [true]
```

At ``u=(1,1,1)``, equal phase exports preserve symmetry. The returned result
supports only the claim “locally feasible at the simultaneous bound point.”
It does not establish box containment.

## Ask the range question

For three participants, exact corner enumeration is still inexpensive:

```julia
corner_check = verify_operating_envelope(
    net, cps, bound;
    utilizations=:corners,
    control_policy=PerfectRecourse())

corner_check.feasible
corner_check.context_results[1]
```

Several asymmetric corners violate the negative-sequence constraint. The same
failure can be found without enumerating every corner by starting from the zero
and upper points and refining around the most stressed point:

```julia
search = search_operating_envelope_adversarial(
    net, cps, bound;
    seed_samples=0,
    refinement_rounds=1,
    restarts=2,
    initial_step=0.5,
    control_policy=PerfectRecourse())

@assert search.outcome == :candidate_counterexample

search.diagnostics["candidate_locations"]
search.worst_interval
search.worst_context.utilization
```

The first refinement evaluates points such as ``(0.5,1,1)``. Their unequal
phase injections expose a failure hidden by the balanced bound point. The
search retains both verification rounds and every normalized-headroom score.

## Interpret the result

This experiment distinguishes four statements:

1. The simultaneous upper point is locally feasible.
2. The advertised independent box is falsified by a partial-utilization point.
3. Exact corner verification also finds unsafe points in this particular case.
4. Neither the local infeasibility status nor the heuristic search is a global
   mathematical certificate for a general nonlinear AC model.

For a publication-quality study, rerun the candidate from multiple starts,
independently recompute the binding network quantity, preserve the complete
study manifest, and repeat with operationally defensible control recourse.

## What this does not prove

- It does not show that all worst cases occur at corners.
- It does not assign a violation probability or uncertainty coverage level.
- It does not certify global AC infeasibility at the candidate point.
- It does not validate the synthetic impedance or limit values against field
  data.
- `PerfectRecourse()` is used because this case has no free network controller;
  studies with taps, STATCOMs, or other controls must state their information
  structure explicitly.
