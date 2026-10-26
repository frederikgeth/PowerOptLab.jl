# Constraint-based inverter operation

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

## Simulation and optimisation formulations

Here, **simulation** means the algebraic forward controller at a fundamental-frequency
network equilibrium. It does not mean time-domain or switched EMT simulation.
The physical plant is shared; operating relationships and the selection objective
determine which feasible point is selected.

| Concept | Forward simulation (`SequenceController`) | Optimisation (`ConstrainedInverter`) | Support / equivalence condition |
|:--|:--|:--|:--|
| Network, reduced filter, primitive matrices, LCL, POC shunt | Plant KCL/KVL and prescribed current | Same plant KCL/KVL and dispatchable current | Implemented; compare powers at the same location |
| Three-leg zero sum | Plant requires converter `I0=0` | Same topology equality, including `:free` allocation | Implemented; topology does not require `I2=0` |
| Current allocation | Positive/negative sequence commands at converter or grid target | `:balanced_converter`, `:balanced_grid`, or `:free` | Implemented; LCL capacitor current makes the two balance choices distinct |
| Volt-var | Curve produces a Q command before limiting | `VoltVarWattConstraints` imposes total POC Q equality | Implemented; saturation can intentionally produce different Q |
| No reactive curve | Zero Q command | `VoltVarWattConstraints`: Q=0; `InverterOperatingConstraints`: Q free | Explicitly different defaults |
| Fixed Q / signed PF / free Q | No dedicated fixed-PF policy; a constant curve can prescribe Q | `FixedReactivePower`, `FixedPowerFactor`, `FreeReactivePower` | Implemented on total POC powers; PF does not select current sequences |
| Volt-watt | Active command, rated or AC-available basis | Upper bound on POC P with the same curve/basis | Implemented; maximum export can reproduce an unsaturated command |
| P/Q dispatch bounds | Request and controller prescribe operation | Optional `p_min`, `p_max`, `q_min`, `q_max` | Implemented; inconsistent bounds and service remain infeasible |
| DC source availability and charging | Request availability with loss-aware command | Request upper bound plus signed `DCSourceBounds` | Implemented; `p_min=0` declares non-charging PV; default has no lower bound |
| Converter/grid conductor current and neutral RMS | Command backoff plus hard plant inequalities | Hard plant inequalities | Implemented, including fixed switching RMS reserves |
| Converter apparent power | P/Q allocation and plant-aware backoff | Hard converter S circle | Implemented; total POC Q can differ from converter Q |
| Proportional current limiter | Scales both current sequences / P and Q together | No implicit limiter; optional explicit current-ray restriction plus maximise scale | Reproduction example provided; the ray is an extra operating assumption |
| Sequence priority | Reserve one sequence, fill headroom with the other | Fix admissible reserved sequence and maximise other-sequence scale | Fixed-measurement reproduction provided; arbitrary network-equilibrium equivalence is not established |
| Negative-sequence admittance / ripple cancellation | Prescribed negative-sequence current law | `:free` exposes sequence choices; physical sequence/ripple caps still apply | Law-specific optimal objectives or exact command graphs remain future work |
| Internal voltage and switching utilisation | Same sampled switching hull and voltage boxes | Same hard bounds | Implemented; sampled hull is an outer approximation, audit `switching_margin` |
| DC 2ω ripple and finite line-frequency source | Plant ripple circuit and backoff | Same hard circuit/capability constraints | Implemented; source impedance here differs from signed DC power bounds |
| Four-leg / split-link neutral, midpoint, capacitor heating, ESR | Physical plant; controller wrapper remains three-leg | THREE_LEG, FOUR_LEG, SPLIT_DC plants; common/individual bank bounds compose | Implemented for three-phase optimisation; split midpoint balancing remains bounded by configured charge authority |
| Semiconductor and capacitor losses | Same physical loss equations | Same physical loss equations under either power direction | Implemented; no objective-dependent inflatable loss epigraph |
| Carrier PWM reserve closure | Plant solver has post-solve carrier iteration/audits | Fixed reserves supported; wrapper requires `pwm_strategy=:NONE` | Automatic optimisation closure remains future work |
| Single phase / grid forming | Plant supports legacy single phase and fixed balanced internal EMF; not `SequenceController` | Wrapper rejects both | Dedicated operating laws remain future work; no GFM droop or dynamic-stability claim |
| Curve encoding | Smooth local/global encodings and separate exact numeric evaluator | Same smooth voltage encodings, exact residual audit | Exact PWL graphs in complete inverter models remain future work |
| Objective / fleet / time coupling | Objective should not change prescribed operation; forward fleet/study APIs exist | Caller chooses objective and can stamp devices in a shared model | Composition implemented; dedicated constrained-fleet/DOE/time APIs remain future work |
| Hardware sizing / variable DC rail | Existing controller parameter/sizing studies | Fixed physical hardware in this wrapper | Continuous design variables and DC energy/rail dynamics remain future work |

## Select an operating mode

```julia
# Free Q, bounded active dispatch, PV source, balanced grid-arm current.
operating = InverterOperatingConstraints(p_min=0., p_max=4000.,
    q_min=-2000., q_max=2000.)
device = ConstrainedInverter(plant, operating;
    current_allocation=:balanced_grid, dc_source=DCSourceBounds(p_min=0.))

# Alternatives: FixedReactivePower(600.) or FixedPowerFactor(.95).
operating = InverterOperatingConstraints(reactive=FixedPowerFactor(.95))
```

For positive P, positive PF absorbs Q and negative PF supplies Q. The relation
is ``|pf|Q+sign(pf)\sqrt{1-pf^2}P=0`` for either direction of active power.
It imposes neither a direction nor a balanced-current assumption; PF at zero
complex power is undefined. Zero PF is rejected. Optional P/Q bounds compose
with any reactive mode. Use `VoltVarWattConstraints` for a voltage-dependent Q
equality; the general policy rejects a second, ambiguous Q law.

`DCSourceBounds(p_min=0)` permits loss imports at the AC port, while preventing
net DC charging. Signed bounds can model a bidirectional source. Omitted bounds
retain the previous availability-only behaviour. Source bounds are physical DC
power bounds, not bounds on POC active power or a battery energy model.

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
``A=P_{upper}-(P_{DC}-P_{POC})``, where ``P_{upper}`` is the tighter of request
availability and an optional DC source upper bound, and impose
``P_{POC}\le f_{VW}(V) A``. Together with the independent DC bound this is
equivalent to ``P_{POC}\le\min(A,f_{VW}A)``. When ``A<0``, the full loss import
is retained. DC availability is not incorrectly treated as lossless AC output.

The plant retains per-conductor RMS current, converter apparent-power,
switching-voltage, capacitor and DC ripple limits. There is no radial P/Q
allocation, common current scale, or automatic Q slack. An impossible reactive
request gives an unsuccessful solve. Physical limits are not relaxed to make a
voltage curve feasible.

This path normalises capability rows and operating power rows by the
device rating. This keeps the solver's constraint-bound relaxation small in
physical units when a small inverter is embedded in a large network power base.
Existing forward-controller stamping retains its original numerical defaults.

The wrapper supports three-phase THREE_LEG, FOUR_LEG and SPLIT_DC plants with
`grid_forming=false`. The default `:balanced_converter` imposes zero and negative
sequence converter current (without duplicating the three-leg zero-sum row).
`:balanced_grid` imposes balance on the grid arm; `:free` leaves current allocation
to the optimiser within topology and capability bounds. Balanced grid current
and three-leg converter zero sum can be incompatible with a grounded LCL shunt
under zero-sequence voltage: infeasibility is a valid physical result.
Q tracking includes the filter and optional POC
shunt through the plant's total POC exchange. `pwm_strategy=:NONE` is required;
fixed calibrated RMS current reserves are supported, but automatic carrier
reserve iteration is not part of this solve.

## [Compare saturation at a fixed 207 V POC](@id inverter-saturation-comparison)

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

General policies also report `reactive_mode`, `current_allocation`, source
bounds and source/P/Q margins. Fixed PF has `fixed_pf_residual_VA`; free Q and PF
have `nothing` for Q targets and tracking residuals. A missing watt ceiling in
`InverterOperatingConstraints` similarly returns `nothing` for cap margins.
Four-leg/split-link results include neutral RMS, midpoint and capacitor margins.

Margins are positive inside their limits. Undeclared optional limits return
`nothing`. Failed solves retain `SolveStatus` and return NaN numerical outputs.
There are no invented limiter scales or exact-current replay fields for this
feasible-set formulation. As with the other nonlinear solvers, a successful
solve provides a local solution and does not prove global optimality or global
infeasibility.

## Reproduce forward operation through optimisation

The executable [reproduction cases](../ibr/formulation_reproductions.md) compare
independently solved feeder equilibria, saturated proportional limiting and
sequence-priority allocation. They distinguish agreement under matched
assumptions from the intentional saturation difference above.
