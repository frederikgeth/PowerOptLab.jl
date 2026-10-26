# Constraint-based volt-var and volt-watt operation

[`ConstrainedInverter`](@ref) adds operating constraints to an
[`AdvancedInverter`](@ref). It is an alternative to the forward
[`SequenceController`](@ref), which retains its existing current limiters and
defaults. Choose according to the question being studied:

| Formulation | Selects the operating point through | When a current limit binds |
|:--|:--|:--|
| `ControlledDevice(plant, SequenceController(...))` | Prescribed current-command equations | The selected limiter backs off the command |
| `ConstrainedInverter(plant, VoltVarWattConstraints(...))` | An explicit dispatch objective, or a specified dispatch | The current inequality restricts feasible P while the Q equality remains |

The second formulation is an OPF feasible set. Maximising export over a network
does not establish that an autonomous local controller will select that point,
or that it will be dynamically stable. Neither option is a certification model.

## Physical and operating relationships

With positive injection, the new formulation imposes total POC reactive power

```math
Q_{POC}=q_{scale} f_{VV}(V),\qquad P_{DC}\le P_{available}.
```

The voltage input and reference come from the declared positive policy and
`sensing_reference`. For `sensing=:average_voltage`, the input is the arithmetic
mean of the three RMS phase magnitudes. Curve encodings use the same smooth,
flat-extended PWL functions as the existing controller. A missing Q curve imposes
zero total POC Q. A missing watt curve gives a fraction of one.

For `volt_watt_basis=:rated`,

```math
P_{POC}\le P_{rated} f_{VW}(V).
```

For `:available`, define the solved AC availability
``A=P_{available}-(P_{DC}-P_{POC})`` and impose
``P_{POC}\le f_{VW}(V) A``. Together with the independent DC bound this is
equivalent to ``P_{POC}\le\min(A,f_{VW}A)``. When ``A<0``, the full loss import
is retained. DC availability is not incorrectly treated as lossless AC output.

The plant retains per-conductor RMS current, converter apparent-power,
switching-voltage, capacitor and DC ripple limits. There is no radial P/Q
allocation, common current scale, or automatic Q slack. An impossible reactive
request gives an unsuccessful solve. Physical limits are not relaxed to make a
voltage curve feasible.

This path normalises squared current and apparent-power capability rows by the
device rating. This keeps the solver's constraint-bound relaxation small in
physical units when a small inverter is embedded in a large network power base.
Existing forward-controller stamping retains its original numerical defaults.

The initial implementation supports three-leg, grid-following plants with
balanced **converter** current: the wrapper imposes ``I_{2,conv}=0`` and the
plant already supplies ``I_{0,conv}=0``. This allocation is explicit; three-leg
topology by itself does not require balanced current. LCL shunt current can make
grid-side current unbalanced. Q tracking includes the filter and optional POC
shunt through the plant's total POC exchange. `pwm_strategy=:NONE` is required;
fixed calibrated RMS current reserves are supported, but automatic carrier
reserve iteration is not part of this solve.

## Compare saturation at a fixed 207 V POC

This ideal, lossless example fixes POC voltage to isolate capability allocation
from network feedback. Zero filter impedance is a deliberate idealisation.

```@example constrained_inverter
using PowerOptLab
using BMOPFTools: parse_bmopf

net = parse_bmopf("""
{"bus":{"poc":{"terminal_names":["a","b","c","n"],
 "perfectly_grounded_terminals":["n"],"v_min":[180,180,180],"v_max":[280,280,280]}},
 "voltage_source":{"vs":{"bus":"poc","terminal_map":["a","b","c"],
 "v_magnitude":[207,207,207],
 "v_angle":[0,-2.0943951023931953,2.0943951023931953]}}}
"""; from_string=true)

plant = AdvancedInverter(id="pv", bus="poc", topology=:THREE_LEG,
    phase_terminals=["a","b","c"], neutral=nothing,
    s_max=6000., i_max=1.1*6000/(3*230),
    v_dc=800., c_dc=2e-3, r_filter=0., x_filter=0.)
intent = VoltVarWattIntent(
    volt_var=PWLFunction([207.,220.,240.,258.], [.44,0.,0.,-.60]),
    volt_watt=PWLFunction([253.,260.], [1.,.2]),
    sensing=:average_voltage, volt_watt_basis=:rated)
encoding = VoltVarWattEncoding(
    volt_var=LocalC2Formulation(.01), volt_watt=LocalC2Formulation(.01))
request = InverterControlRequest(p_available=6000., p_rated=6000., q_scale=6000.)

policy = VoltVarWattConstraints(intent, encoding)
device = ConstrainedInverter(plant, policy)
constrained = solve_constrained_inverter(net, device, request;
    objective=:max_export, solver_options=("tol" => 1e-9,))
forward = solve_controlled_inverter(net,
    ControlledDevice(plant, SequenceController(lower_positive_policy(intent, encoding);
        limiter=CommonScaleLimiter(pq_priority=:var))), request;
    solver_options=("tol" => 1e-9,))
@assert constrained.solve.publishable && forward.solve.publishable
@assert abs(constrained.operation.q_constraint_residual_var) < .05

[(model="forward controller", P_W=forward.plant.p_poc, Q_var=forward.plant.q_poc),
 (model="operating constraints", P_W=constrained.plant.p_poc, Q_var=constrained.plant.q_poc)]
```

Ignoring the small curve-smoothing effect at the endpoint, the forward model
delivers approximately 5334 W and 2614 var. The constrained model delivers
approximately 5321 W and 2640 var. In this balanced idealisation,
``S_{effective}=\min(S_{rated},3VI_{max})=5940`` VA, so the independent
capability oracle is ``P=\sqrt{5940^2-2640^2}``. The constrained model curtails
about 13 W more active power and retains the requested reactive power.

## Fixed dispatch and composition

```@example constrained_inverter
fixed = solve_constrained_inverter(net, device, request;
    objective=:feasibility, p_set=4000., solver_options=("tol" => 1e-9,))
@assert fixed.solve.publishable
(P_W=fixed.plant.p_poc, Q_var=fixed.plant.q_poc)
```

For an external OPF, call `stamp_device!(ctx, device; request)` inside the model
hook. The returned handles expose physical plant quantities through
`inverter_handles(device, handles)`. The stamp adds no objective, so the caller
can select cost, export, a fixed dispatch, or another criterion. Use
`extract_device(device, handles, solve_status)` after solving. For a fleet,
stamp each device in the shared model; a sum-of-export objective represents
coordinated fleet optimisation. The existing `ControlledInverterFleetSpec`
continues to describe the forward-controller formulation.

## Read the diagnostics

`result.plant` reports physical SI quantities. `result.operation` reports:

- `q_target_var` and `q_constraint_residual_var`: the encoded curve equality;
- `q_exact_target_var` and `q_tracking_residual_var`: the exact PWL policy at the
  solved voltage, separating curve approximation from electrical feasibility;
- `volt_watt_margin_W` and `exact_volt_watt_margin_W`: encoded and exact cap slack;
- `dc_availability_margin_W` and `ac_available_W`: the physical power budget;
- current, converter apparent-power, switching-voltage and DC ripple margins.

Margins are positive inside their limits. Undeclared optional limits return
`nothing`. Failed solves retain `SolveStatus` and return NaN numerical outputs.
There are no invented limiter scales or exact-current replay fields for this
feasible-set formulation. As with the other nonlinear solvers, a successful
solve provides a local solution and does not prove global optimality or global
infeasibility.
