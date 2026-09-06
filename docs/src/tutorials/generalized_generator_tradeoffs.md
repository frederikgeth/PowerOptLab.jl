# Comparing generator trade-offs under unbalance

This tutorial runs actual network solves using the new components. Start with the
[phasor and grounding examples](generalized_generator_models.md) if the distinction
between an internal EMF, PCC voltage, and an angle reference is unfamiliar. The
[API guide](../components/generalized_generator_api.md) documents the implementation.

## The experiment

One four-wire feeder supplies unequal phase loads of 1.5, 3.5, and 0.8 kW, with
reactive demands 0.3, 0.6, and 0.1 kvar. Its upstream phase magnitudes are
240, 225, and 232 V, with slightly noncanonical phase angles. The PCC neutral is
not perfectly grounded; the upstream neutral is. Phase and neutral line
impedances are retained, so voltage unbalance and neutral displacement arise
from the same network equations used for every case.

The tested local component attaches at the PCC. Cases a-e compare voltage laws;
the terminal-PV case regulates positive-sequence PCC voltage; f leaves phase
allocation to a loss-minimizing OPF subject to sequence bounds. Dedicated-source
cases compare open, finite, and ideal internal-star grounding. The complete
data and computation are in `scripts/generalized_generators/tradeoffs.jl`.

Run it independently:

```sh
julia --project=. scripts/generalized_generators/tradeoffs.jl
```

The following results are generated during the documentation build, so the
displayed comparison cannot silently drift from the implementation:

```@example generator_tradeoffs
using PowerOptLab
include(joinpath(pkgdir(PowerOptLab), "scripts", "generalized_generators", "tradeoffs.jl"))
rows = generator_tradeoff_study()
print_generator_tradeoffs(rows)
```

All phasors are fundamental-frequency RMS values. `Iw` is maximum winding
current, `IL` maximum external conductor current, `In` neutral-lead current,
`Ig` earth current, and `Ic` delta circulating current. VUF is
``100|V_2|/|V_1|`` at the PCC with the 1/3 Fortescue transform. The delta section
below explains its line-line sequence quantities and the `V1*` conversion. The loss column contains the local
component's physical series and grounding losses; it excludes feeder losses.

## What is held fixed?

The network, loads, and upstream source remain identical. The regular generator
cases use a common full port impedance, including a reduced neutral contribution.
Their source laws/control constraints differ deliberately:

| Case | Specified | Allowed to respond |
|:--|:--|:--|
| a | Each phase P and Q | Internal EMF shape |
| b | Entire internal phasor | P, Q and phase currents |
| c | Internal magnitude template and total P | Common source rotation and Q |
| d | Total P/Q, balanced internal shape | Common magnitude/rotation |
| e | Total P and each phase Q, fixed relative angles | Separate phase magnitudes and common rotation |
| Terminal PV | Total P and PCC ``|V_1|=230`` V | Common internal magnitude, rotation and Q |
| f | Total P/Q and sequence current limits | Phase allocation and internal phasor, selected by custom-loss minimization |

These rows are different physical/control feasible regions, not alternative
formulations of one identical feasible region. In particular, b has no active
dispatch target. Its power exchange must not be described as a fair economic
improvement over a P-constrained case. The source variants use explicit conductor
primitives rather than the regular generator's reduced one-return impedance.

## Read the trade-offs

First compare c with terminal PV. Holding the internal magnitude is an excitation
constraint; regulating ``|V_1|`` measures and controls a terminal quantity. Series
voltage drop separates those quantities. Check Q and the maximum phase current
alongside VUF: voltage support uses reactive/current capability, and a perfectly
balanced internal source does not enforce balanced PCC voltage.

Then compare d with e. Both share a relative-angle pattern, but only e can change
the individual magnitudes. Its extra flexibility requires a physical actuator
interpretation; a conventional machine's common field winding is not automatically
three independent voltage actuators. Comparing phase current and loss explains
the consequence of that modeling assumption without claiming e is universally
preferable.

In f, limits bound the sequence currents but do not prescribe their angles or a
controller priority. The OPF chooses a feasible allocation minimizing local losses.
It should be interpreted as dispatch capability. It does not prove an actual
converter would settle there; compare it with the advanced IBR and an explicit
controller if that is the scientific question.

For the grounded sources, compare `In` and `Ig` separately. The component's
external conductor-current sum is balanced by earth current. A ground connection
can redistribute neutral and phase-to-neutral voltage stress. It is not merely a
numerical reference, and a lower local loss is not sufficient to assess grounding
design, protection, or electrode safety.

Finally compare the two a cases. They specify identical PCC phase P/Q and differ
only in series impedance. With unconstrained internal EMF, the terminal solutions
should agree while internal voltage and losses differ. This is a useful invariant:
adding an unconstrained internal node cannot make a terminal PQ injection behave
like a voltage regulator.

```@example generator_tradeoffs
finite = only(filter(r -> r.case == "a: phase PQ", rows))
ideal = only(filter(r -> r.case == "a: ideal series", rows))
@assert isapprox(finite.result.port_voltage, ideal.result.port_voltage; atol=1e-5)
@assert finite.loss_w > ideal.loss_w
@assert maximum(r.power_error_w for r in rows) < 1e-4
println("PQ terminal equivalence and all device power balances verified.")
```

## Change a capability rather than hiding a failed target

To study a limit transition, reduce a PV generator's `i_max` or internal
`magnitude_max` while retaining its P and terminal voltage target. A hard target
can become infeasible. This implementation reports that outcome and masks
unpublished measurements; it does not silently relax the voltage or switch to PQ.
An AVR/current-limiter experiment requires an explicitly modeled limit policy.

For a controlled comparison, vary only one physical parameter and preserve the
measurement location, dispatch target, grounding, and other ratings. For example,
sweep neutral impedance with fixed grounding, or sweep a negative-sequence current
limit with the same objective and total P/Q. Record solver status, physical
residuals, active constraints and alternative starts, as well as the decision
metric. Local NLP results do not establish global optimality or dynamic stability.


## Delta: the terminal meter cannot see every winding current

The final three cases use ordered ab,bc,ca windings. `delta: fixed EMF` uses
line-line EMFs obtained from the same balanced 230 V phase template. The
`delta: circulation` case adds the same ``10+j5`` V to all three internal EMFs.
With equal winding impedance ``z=0.8+j1.2`` Ω, this adds

```math
\Delta I_{\mathrm{circ}}=\frac{10+j5}{0.8+j1.2}\ \mathrm A
```

to every winding current. Their differences, and therefore external line currents
and voltages, are unchanged. The executable tutorial checks that equality and the
increase in loss. This is why a line-current limit cannot replace winding-current
limits, even if the externally delivered P and Q are identical. Neither delta
case has a neutral or earth terminal.

`Iw` is maximum winding current, `IL` maximum external conductor current,
`Ic` the magnitude of delta circulating current, and `In`/`Ig` the neutral/earth
currents. For a common-return component, `IL` also includes its neutral.
`V1*` reports phase-neutral positive-sequence magnitude for wye and line-line
positive-sequence magnitude divided by ``\sqrt3`` for delta, allowing comparison
on the same voltage scale. This conversion does not reconstruct a delta's
unobservable phase-to-earth zero-sequence voltage. VUF uses the same magnitude
ratio for either convention.

`delta: ideal PQ` instead fixes three unequal winding P/Q injections with zero
series impedance. Winding voltages close around the triangle automatically;
internal nodes and extra cycle equations are unnecessary. These winding P/Q
setpoints differ from wye phase-to-neutral setpoints, so that comparison changes
the connection and the physical location of the prescribed phase powers.
