# Bilevel PV export and utility tap control

> **Kind:** Problem specification · **Maturity:** proof of concept · **Direction:** hierarchical · **Temporal:** single-snapshot

This page is a proof of concept for a hierarchical distribution-network model:

1. the utility chooses a continuous transformer tap;
2. the PV systems respond through native BMOPFTools Volt-var and Volt-watt
   curves; and
3. the utility evaluates the lower-level response and minimises a smooth voltage
   stress metric plus tap movement. The upper-level search is safeguarded and
   one-dimensional; each trial uses a fresh lower-level solve.

The important implementation detail is that the lower model remains a live
staged BMOPFTools model. Its utility tap is a `JuMP.Parameter` bound to the
native transformer tap with `bind_opf_parameter!`. After the lower model is
solved, `DiffOpt.forward_differentiate!` gives the local derivative of monitored
voltages with respect to the tap. The outer loop uses that derivative in a
safeguarded derivative search.

```julia
using PowerOptLab

net = bilevel_demo_network()
bilevel = solve_bilevel_pv_tap(net;
    transformer_id="reg",
    pv_ids=["pv1", "pv2"],
    monitored_buses=["lv1", "lv2"],
    max_iterations=12)

centralized = solve_single_level_pv_tap(net;
    transformer_id="reg",
    pv_ids=["pv1", "pv2"],
    monitored_buses=["lv1", "lv2"],
    export_weight=0.25)

(bilevel.tap, bilevel.exported_power_W, bilevel.voltages_V)
(centralized.tap, centralized.exported_power_W, centralized.voltages_V)
```

For a fixed tap, the lower-level response and its local voltage sensitivity can
be inspected directly. This is useful for reproducing the DiffOpt-versus-finite
difference checks in the accompanying case-study repository:

```julia
response = solve_bilevel_pv_response(net;
    transformer_id="reg", pv_ids=["pv1", "pv2"],
    monitored_buses=["lv1", "lv2"], tap=1.0,
    lower_level=:local_controller)

response.voltage_sensitivity_V_per_tap
```

Voltage outputs and sensitivities are keyed consistently by `(bus, phase)`
tuples. The default monitored quantity is phase-to-neutral (`voltage_measurement=:pn`);
use `:pg` for phase-to-ground or `:pp` for cyclic phase-to-phase monitoring on
three-wire or delta configurations.

There are two lower-level interpretations:

```julia
aggregate = solve_bilevel_pv_tap(net;
    transformer_id="reg", pv_ids=["pv1", "pv2"],
    monitored_buses=["lv1", "lv2"], lower_level=:aggregate)

local_controller = solve_bilevel_pv_tap(net;
    transformer_id="reg", pv_ids=["pv1", "pv2"],
    monitored_buses=["lv1", "lv2"], lower_level=:local_controller)
```

`:aggregate` is the shared greedy `Max sum(p_pv)` formulation. It is useful
when the lower-level abstraction is a coordinated export maximiser. The
`:local_controller` option removes that shared objective and adds one
deterministic response equation per PV: native Volt-var supplies reactive
power, while active power follows Volt-watt and is smoothly capped by the
remaining apparent-power capability. This is the closer model of independent
PV controllers responding to their observed local voltage. DiffOpt is
appropriate here because the controller equations and network equilibrium are
still one smooth differentiated lower-level solve; it is not being used to
differentiate through a discrete best-response or a nonsmooth kink.

The centralized solve can optionally include an export reward with
`export_weight`; the default is zero so the voltage/tap objective is not
silently changed. It chooses the tap and PV decisions together and answers a
different question from the hierarchy.

## Modelling interpretation

The IBR `control_profile` carries the native smooth Volt-var/Volt-watt equations.
For the current BMOPFTools DiffOpt wrapper, use `softplus=:builtin`; the stable
user-defined softplus is not accepted by that wrapper. The demo uses
`q_ref="VAR_MAX"` for Volt-var and `p_ref="S_MAX"` for Volt-watt.
For Volt-Watt, BMOPFTools interprets `p_limits` as `[p_low, p_high]`, so the
demo's `[0.0, 1.0]` setting retains full output at the low-voltage breakpoint
and curtails toward the high-voltage breakpoint.

Monitored voltages are constructed using the actual bus terminals. Result keys
use the consistent `(bus, phase)` tuple form for both single- and multi-phase
measurements. The default is phase-to-neutral; `voltage_measurement=:pg` and
`:pp` support phase-to-ground and cyclic phase-to-phase measurements.

The aggregate lower-level objective is `Max sum(p_pv)`. This is an aggregate
greedy model: all consumers value export positively and there is no export
price or penalty. It is not a full non-cooperative Nash model. The local
controller option is different: it prescribes each consumer's response rather
than solving a shared dispatch problem. A true consumer-by-consumer game would
require a best-response loop (or a generalized Nash/MPEC formulation) because
consumers share network voltage constraints.

On the small nominal demo, the aggregate and local-controller responses may
coincide. That is expected when the aggregate objective selects the same native
Volt-Watt/apparent-power-limited export point as the deterministic controller;
it is not evidence that the two formulations are interchangeable in general.
Use a stressed feeder or the centralized free-dispatch benchmark to expose
authority differences without manufacturing a gap between mathematically
aligned lower-level objectives.

The outer objective is

```math
\sum_i \left((V_i - V_\mathrm{ref})/\Delta V\right)^4
  + \rho (t-1)^2,
```

subject to the lower-level network bounds. The fourth-power term penalises
voltage excursions smoothly while remaining differentiable for DiffOpt. It is
not a probabilistic robustness certificate; for uncertainty, solve several
scenarios and take the worst-case margin.

## Limitations to keep explicit

- BMOPFTools' IVR-EN OPF is nonconvex, so both lower-level solutions and
  sensitivities are local.
- DiffOpt derivatives are valid only while the lower-level active set remains
  stable. Volt-watt corners and voltage-bound transitions need perturbation
  checks.
- Continuous tap control is used here; discrete regulator positions need a
  mixed-integer outer layer or enumeration.
- `BilevelPVResult.converged` and `termination_reason` report whether the
  upper-level derivative search reached a stationary point, a tap bound, or a
  numerical/iteration limit. A lower-level `LOCALLY_SOLVED` status alone is
  not an upper-level convergence certificate.
- The consumer objective does not include negative export prices, curtailment
  compensation, or fairness. Those are policy choices that can be added as
  lower-level weights or a game layer later.
