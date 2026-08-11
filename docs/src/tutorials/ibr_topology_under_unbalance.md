# [Choosing an IBR topology under unbalance](@id ibr-topology-under-unbalance)

This tutorial compares the implemented three-phase converter categories on the
same balanced and unbalanced four-wire grid. Its purpose is to show when
`:THREE_LEG`, `:FOUR_LEG`, and `:SPLIT_DC` cease to be interchangeable.

The experiment holds the AC rating and filter fixed. Only the physical neutral
path and DC structure change.

## 1. Build two study networks

The source is deliberately stiff so the converter limits, rather than feeder
voltage drop, dominate the comparison.

```julia
using PowerOptLab
using BMOPFTools: parse_bmopf

function four_wire_grid(; magnitudes, angles)
    parse_bmopf("""
    {"bus":{
      "grid":{"terminal_names":["a","b","c","n"],
              "perfectly_grounded_terminals":["n"]},
      "poc":{"terminal_names":["a","b","c","n"],
             "perfectly_grounded_terminals":["n"],
             "v_min":[180.0,180.0,180.0],
             "v_max":[270.0,270.0,270.0]}},
     "voltage_source":{"vs":{"bus":"grid",
       "terminal_map":["a","b","c"],
       "v_magnitude":$(magnitudes),"v_angle":$(angles)}},
     "linecode":{"lc":{"R_series_1_1":0.05,"R_series_2_2":0.05,
       "R_series_3_3":0.05,"R_series_4_4":0.05}},
     "line":{"l1":{"bus_from":"grid","bus_to":"poc",
       "terminal_map_from":["a","b","c","n"],
       "terminal_map_to":["a","b","c","n"],
       "linecode":"lc","length":1.0}}}
    """; from_string=true)
end

balanced = four_wire_grid(
    magnitudes=[230.0,230.0,230.0],
    angles=[0.0,-2.0944,2.0944],
)
unbalanced = four_wire_grid(
    magnitudes=[245.0,215.0,230.0],
    angles=[0.05,-2.15,2.0],
)
```

The unbalanced case contains both magnitude and angle asymmetry. A magnitude-
only test is not enough to exercise every sequence pathway.

## 2. Define the physical alternatives

```julia
common = (
    bus="poc", phase_terminals=["a","b","c"], neutral="n",
    s_max=20e3, i_max=40.0,
    r_filter=0.05, x_filter=0.15, m_max=0.96,
)

devices = [
    AdvancedInverter(; id="three-leg", topology=:THREE_LEG,
        v_dc=700.0, c_dc=1.1e-3, common...),
    AdvancedInverter(; id="four-leg", topology=:FOUR_LEG,
        v_dc=700.0, c_dc=1.1e-3, In_max=40.0, common...),
    AdvancedInverter(; id="split-dc", topology=:SPLIT_DC,
        v_dc=800.0, c_dc=2.8e-3, In_max=21.0, common...),
]
```

The split-link unit uses a higher DC voltage because each phase is constrained
against a half bus. This is a physical utilisation penalty, not a numerical
adjustment intended to equalise the answers.

## 3. Solve and construct an auditable comparison

```julia
function study(network, device)
    r = solve_advanced_inverter(network, device)
    @assert r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    return (
        topology=r.topology,
        status=r.termination_status,
        export_kw=r.p_poc/1e3,
        neutral_a=r.i_neutral,
        zero_a=r.i_zero,
        negative_a=r.i_negative,
        bus_ripple_v=r.dv2,
        switching_headroom_v=r.switching_margin,
    )
end

balanced_results = study.(Ref(balanced), devices)
unbalanced_results = study.(Ref(unbalanced), devices)

foreach(println, balanced_results)
foreach(println, unbalanced_results)
```

Expected qualitative findings:

- the balanced case makes all three topologies look similar and produces almost
  no neutral current or double-frequency DC ripple;
- the 3-leg bridge still serves an unbalanced voltage set, but enforces
  `Ia + Ib + Ic = 0` and therefore has no zero-sequence current;
- the 4-leg bridge uses its semiconductor neutral leg and must respect `In_max`;
- the split link returns neutral current through its two capacitors, coupling
  AC unbalance to midpoint and bank stress.

Do not rank the topologies from `p_poc` alone. Compare current paths, required
DC voltage, capacitor stress, and switching headroom.

## 4. Make a sequence-current policy explicit

Hardware capability and grid-code policy are separate limits. For example, add
a negative-sequence ceiling to the four-leg unit without changing its neutral-
leg rating:

```julia
sequence_limited = AdvancedInverter(; id="four-leg-sequence-limited",
    topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3, In_max=40.0,
    i_negative_max=0.5, common...)

r_sequence = solve_advanced_inverter(unbalanced, sequence_limited)
@assert r_sequence.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
@show r_sequence.p_poc r_sequence.i_zero r_sequence.i_negative
```

`i_negative_max` does not stand in for `In_max`: negative sequence does not
require neutral current, while zero sequence does. Keeping the two constraints
separate is essential for a defensible four-wire study.

## 5. Publication checks

For every reported point, record:

1. `termination_status` and `solve_status(result)`;
2. `switching_margin`, especially when DC utilisation affects the conclusion;
3. all three sequence currents rather than only `i_neutral`;
4. the physical reason each topology has its selected `v_dc`, `In_max`, and
   capacitor ratings; and
5. sensitivity to the unbalance scenarios, not just the balanced base case.

This remains a fundamental-frequency capability comparison. It does not prove
dynamic neutral control, fault ride-through, or small-signal stability.
