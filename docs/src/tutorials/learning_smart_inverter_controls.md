# Learning and validating local smart-inverter controls

This tutorial gives a reproducible workflow for learning a local inverter
control law from data, then evaluating it with PowerOptLab. The target laws are
constant power factor, Volt-VAr, and Volt-Watt.

The central rule is simple: a deployed controller is a **local policy**, not a
free OPF decision. It observes signals at its point of connection (POC), applies
a prescribed response, and thereby changes the voltages it will observe next.
An OPF can be a useful benchmark or teacher, but its setpoints often depend on
non-local information unavailable to a physical inverter.

## What PowerOptLab evaluates

PowerOptLab does not yet train a learning model. It provides the physically
consistent evaluation stage: an IBR can reference a BMOPFTools `control_profile`,
and that prescribed control remains active in power flow, OPF, and dynamic
operating-envelope (DOE) calculations. Fit policy parameters externally, then
attach the fitted profile to the IBR.

Keep these three tasks separate:

| Task | Target | Deployable result? |
| --- | --- | --- |
| Behavioural identification | Measured inverter P/Q from local voltage | Yes, if the data describe the installed controller |
| Policy design | A parameterised local policy optimised over scenarios | Potentially, after safety validation |
| OPF imitation | Centralized optimal setpoints | Usually no; targets may use non-local information |

## Choose a deployable policy class first

A constant-power-factor law imposes, per phase,

```math
\operatorname{sign}(pf) Q + \tan(\arccos |pf|) P = 0.
```

Positive `pf` is lagging (absorbing VAr) and negative `pf` is leading
(injecting VAr). It is appropriate for a deliberately fixed operating mode, but
it does not provide voltage feedback.

Volt-VAr uses four voltage breakpoints and two reactive limits: it injects VAr
at low voltage, has a zero-VAr deadband, and absorbs VAr at high voltage.
Volt-Watt uses two high-voltage breakpoints to cap active export. Together they
make a closed loop:

```text
local voltage → Volt-VAr / Volt-Watt curve → P/Q injection → feeder voltage → local voltage
```

A profile normally contains one law. If `power_factor` and a droop are both
present, fixed power factor takes precedence with a warning. Do not fit both and
mistake the result for a hybrid controller.

## Construct an identification data set

For each inverter and time index, retain local voltage, measured P and Q,
available active power, controller mode, rating, terminal/phase, timestamps,
curtailment commands, and quality flags. The voltage must be the inverter's POC
measurement—not a substation or another bus. Without available PV power,
Volt-Watt curtailment cannot be distinguished from low irradiance.

Split training and test data by day, weather event, feeder state, or customer
cluster. Randomly splitting adjacent timestamps leaks the same voltage trajectory
into both sets and greatly overstates generalisation.

## Fit constrained, interpretable curves

For Volt-VAr, fit ordered breakpoints and VAr limits subject to apparent-power
capability and grid-code constraints. For Volt-Watt, fit ordered breakpoints and
a non-increasing cap in the curtailment region. These constraints are essential:
a flexible regressor can fit noise with a non-monotone response that creates
positive feedback once it is placed in the network loop.

Use only features that the installed device can measure. If the learner consumes
remote voltage, peer injections, or a state estimate, it is a centralized policy
and must be evaluated with its communication and delay assumptions.

## Replay the learned profile

Translate fitted parameters into a control profile, then bind it to the IBR.

```julia
net["control_profile"] = Dict(
    "learned_pv" => Dict(
        "volt_var" => Dict(
            "voltage_reference" => "PN_PER_PHASE",
            "breakpoints" => [207.0, 220.0, 240.0, 258.0],
            "q_limits" => [-0.60, 0.44],
            "q_unit" => "VA_FRACTION",
            "q_ref" => "VAR_MAX",
        ),
        "volt_watt" => Dict(
            "voltage_reference" => "PN_PER_PHASE",
            "breakpoints" => [253.0, 260.0],
            "p_limits" => [0.20, 1.00],
            "p_unit" => "VA_FRACTION",
            "p_ref" => "P_AVAILABLE",
        ),
    ),
)
net["ibr"]["pv17"]["control_profile"] = "learned_pv"
```

The policy is an equality/control constraint, not a suggestion that the OPF may
override. To test it in an active-power DOE, bind the point to its IBR:

```julia
cp = ConnectionPoint(id="customer_17", bus="lv17", ibr_id="pv17",
                     export_max=10e3)
r = solve_operating_envelope(net, [cp]; direction=:export)
```

The DOE changes active capacity only. It retains the IBR topology,
apparent/current limits, and Volt-VAr/Volt-Watt equality; it does not invent a
separately dispatchable Q source.

### Pitfall: units and signs

`VA_FRACTION` scales a limit by rating; `VAR` and `W` are absolute. Volt-VAr
limits are ordered `[q_absorb ≤ 0, q_inject ≥ 0]`. Check the fitted curve at a
low-voltage point, a deadband point, and a high-voltage point before feeder runs.

## Validate the closed loop

One-step P/Q prediction error is not sufficient. Compare the fitted policy with
a mandated/default curve and fixed-PF baseline on held-out feeder scenarios:

1. high-PV/low-demand and low-PV/high-demand cases;
2. source-voltage, measurement, and parameter uncertainty;
3. simultaneous nearby-inverter response, not one-at-a-time replay; and
4. voltage, thermal, neutral, apparent-power, curtailment, and VAr-throughput
   metrics.

Use voltage-violation frequency and magnitude, total curtailed energy, and
constraint margins as primary metrics. A small error near a steep knee can move
the closed-loop equilibrium substantially.

## Multi-phase pitfalls

The voltage reference is physical. `PN_PER_PHASE` is phase-to-neutral,
`PG_PER_PHASE` phase-to-ground, and `PP_PER_PHASE` phase-to-phase; each can be
averaged across phases. In a four-wire feeder with neutral displacement,
phase-to-ground and phase-to-neutral are not interchangeable, while averaging
can hide the phase that needs support.

Volt-VAr/Volt-Watt droop is supported for `SINGLE_PHASE` and `FOUR_LEG` IBRs. A
`THREE_LEG` IBR lacks the degrees of freedom for the per-phase droop model and
falls back to box bounds with a warning. Do not validate a per-phase learned
controller with that topology.

## A defensible study protocol

1. State whether the task is identification, local-policy design, or OPF
   imitation.
2. Declare deployment-time measurements and actuation signals.
3. Fit a constrained, grid-code-compliant policy with event-wise data splits.
4. Replay the exact profile in the nonlinear feeder model.
5. Validate simultaneous closed-loop behaviour on held-out, stressed, and
   unbalanced cases where relevant.
6. Report the centralized OPF only as a benchmark unless its information and
   communications are actually deployable.

This separation between learned data behaviour, network feasibility, and local
device information is what turns an attractive fitted curve into a credible
smart-inverter control study.
