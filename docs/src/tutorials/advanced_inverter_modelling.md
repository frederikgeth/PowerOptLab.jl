# Advanced inverter modelling: when the POC IBR is not enough

This tutorial helps choose between PowerOptLab's ordinary inverter-based-resource
(IBR) model and `AdvancedInverter`. The distinction is important: an ordinary IBR
is a controlled, bounded current injection at the point of connection (POC), while
an advanced inverter makes the converter's internal AC node and DC-side constraints
explicit.

For many distribution studies, the ordinary IBR is the more defensible model. Use
the advanced model only when the omitted internal physics changes the conclusion.

## Model-selection guide

| Study question | Recommended model | Why |
| --- | --- | --- |
| Feasible POC injections, fixed power factor, or mandatory Volt-VAr / Volt-Watt behaviour | Ordinary IBR | The control law and limits are imposed at the POC; an internal converter circuit is not needed. |
| Hosting capacity or a DOE where DERs are represented by POC setpoints | Ordinary IBR | It is compact, transparent, and matches the available control interface. |
| Does an LCL/L output filter change the POC voltage, reactive power, or export? | `AdvancedInverter` | The filter separates the POC from converter voltage and current. |
| Can a grid-forming converter sustain its requested internal EMF? | `AdvancedInverter` | Internal voltage, modulation, and DC limits must be represented. |
| Is a four-wire unit limited by neutral current or DC-link ripple? | `AdvancedInverter` | These are topology- and DC-capacitor-dependent constraints. |

The ordinary IBR should not be treated as a low-fidelity version of every physical
converter. It is a different abstraction: a POC-level resource model. In particular,
it is often the right choice for a PV inverter following a prescribed Q-V curve. The
advanced model exposes electrical capability; it does not automatically reproduce a
manufacturer's supervisory or mandatory grid-support control law.

## The internal-node model

`AdvancedInverter` inserts an internal AC node behind an optional output filter:

```text
POC bus ──[ r + jx, optional grid-side shunt ]── internal node ── converter ── DC link
```

Converter current and apparent-power limits apply at the internal node. Therefore,
the converter's rating is not necessarily the POC export rating. With a non-zero
filter impedance, part of the converter power supplies filter losses and reactive
exchange; the POC voltage can also differ materially from internal voltage.

Here is the basic pattern. The exact network construction is deliberately omitted;
the object can be passed to `solve_advanced_inverter` with the network used elsewhere
in PowerOptLab.

```julia
inv = AdvancedInverter(
    id = "pv-1",
    bus = "load_bus",
    s_max = 5_000.0,
    i_max = 25.0,
    r_filter = 0.08,
    x_filter = 0.15,
    p_loss_fixed = 15.0,
    a_loss = 0.2,
    c_loss = 0.03,
)

result = solve_advanced_inverter(network, inv;
    objective = :min_loss,
    p_set = 4_500.0,
)
```

Inspect `result.p_poc`, `result.q_poc`, `result.p_conv`, `result.p_loss`, and
`result.v_int_mag` together. A result can be feasible at the converter while falling
short of a requested POC export, because the filter consumes real power. Losses are
modelled without a complementarity branch: `p_dc = p_conv + p_loss`, so the same
convention is retained in charging and exporting regimes.

### Pitfall: applying nameplate limits at the wrong terminal

If a datasheet's current limit is a converter-side limit, applying it directly to the
POC ignores filter current and reactive flow. Conversely, a contractual POC export
limit still needs a POC constraint in addition to internal converter limits. Report
both terminals rather than calling either one simply "inverter power".

## Grid-forming operation is an internal-voltage statement

Set `grid_forming = true` when the converter controls a balanced internal EMF behind
its filter. The model constrains the three internal phase voltages to a balanced,
120-degree set and chooses their magnitude within its specified bounds.

```julia
gfm = AdvancedInverter(
    id = "bess-gfm",
    bus = "pcc",
    topology = :THREE_LEG,
    phase_terminals = ["a", "b", "c"],
    s_max = 50_000.0,
    i_max = 100.0,
    r_filter = 0.02,
    x_filter = 0.08,
    grid_forming = true,
    v_int_min = 220.0,
    v_int_max = 260.0,
    v_dc = 800.0,
    c_dc = 10e-3,
)
```

This is not equivalent to declaring the POC a slack bus. The surrounding network
still needs a voltage reference. In an unbalanced grid, a balanced internal EMF can
produce unbalanced POC voltages and currents through the filter. Treating a
grid-forming flag as a network reference silently hides that distinction.

## DC modulation is a capability limit, not a post-processing check

At a fixed DC voltage, an inverter cannot create an arbitrary internal AC voltage.
For a single-phase unit, the internal-voltage magnitude is limited by the modulation
index and DC-link voltage. For three-phase units, PowerOptLab can use an exact,
sampled switching polytope to enforce the topology-specific modulation limit.

```julia
single_phase = AdvancedInverter(
    id = "single-phase-pv",
    bus = "service",
    phase_terminals = ["a"],
    topology = :SINGLE_PHASE,
    s_max = 7_000.0,
    i_max = 32.0,
    v_dc = 400.0,
    modulation_max = 0.95,
)
```

For a three-phase device, choose the physical bridge topology explicitly:

```julia
four_wire = AdvancedInverter(
    id = "four-wire-bess",
    bus = "lv_bus",
    phase_terminals = ["a", "b", "c"],
    neutral = "n",
    topology = :FOUR_LEG,
    s_max = 30_000.0,
    i_max = 60.0,
    In_max = 35.0,
    v_dc = 750.0,
    c_dc = 8e-3,
    m_max = 0.95,
    n_samples = 36,
)
```

`n_samples` controls the sampled switching polytope resolution. The default 36
samples gives a close but outer approximation to the true switching boundary. It is
appropriate for most planning studies, but do not interpret a point very near its
edge as a hardware guarantee. Increase the resolution and compare the result when
the conclusion depends on modulation headroom.

### Three-leg, four-leg, and split-DC are not interchangeable

`:THREE_LEG` has no independent neutral leg and cannot carry zero-sequence neutral
current. `:FOUR_LEG` provides an explicit fourth leg and therefore permits bounded
neutral current. `:SPLIT_DC` uses the midpoint of two DC-link capacitors; it can serve
four-wire loads, but midpoint utilisation reduces modulation headroom and couples
unbalance to capacitor stress.

The modelling choice is consequential. A balanced test network may make all three
topologies appear equally capable. Under phase-voltage or load unbalance, a three-leg
bridge rejects neutral current, a four-leg bridge may bind `In_max`, and a
split-DC bridge may require a larger DC voltage to deliver the same AC operating
point. Avoid using a balanced feeder as evidence that topology does not matter.

## Neutral current and double-frequency ripple

Single-phase conversion transfers pulsating power to the DC link. The model can
limit this with `p_ripple_max`; in three-phase four-wire and split-DC topologies it
also reports a two-times-line-frequency DC-bus ripple `result.dv2`. Specify `c_dc`
and, where relevant, `dv2_max` to make capacitor sizing a feasibility condition
rather than an after-the-fact calculation.

```julia
ripple_limited = AdvancedInverter(
    id = "unbalanced-bess",
    bus = "lv_bus",
    phase_terminals = ["a", "b", "c"],
    neutral = "n",
    topology = :SPLIT_DC,
    s_max = 30_000.0,
    i_max = 60.0,
    In_max = 30.0,
    v_dc = 800.0,
    c_dc = 5e-3,
    dv2_max = 20.0,
)

result = solve_advanced_inverter(network, ripple_limited;
    objective = :max_export,
)

@show result.i_neutral result.ripple result.dv2
```

With balanced voltage and current, a three-phase bridge has little low-frequency
DC-link pulsation. In contrast, phase imbalance can create neutral current and a
substantial 2ω ripple. A smaller DC capacitor increases `dv2`; a binding `dv2_max`
can reduce feasible export even when RMS current and apparent-power limits look
comfortable.

### Pitfall: checking only RMS quantities

`s_max` and `i_max` do not protect the DC capacitor or neutral conductor. A feasible
RMS operating point can still violate modulation headroom, neutral-current rating,
or ripple tolerance. Conversely, imposing all of these limits on a POC-only IBR
without evidence can create a falsely conservative model. Use the advanced model
when the hardware information is available and the omitted constraint is material to
the study decision.

## A practical modelling workflow

1. Begin with an ordinary IBR whose POC setpoint and control law match the asset's
   operational specification. This is usually the correct baseline for DOE and
   volt-var studies.
2. Introduce `AdvancedInverter` for a small set of representative locations or
   operating points where filters, grid-forming behaviour, or unbalance are expected
   to bind. Keep POC-level contractual limits explicit.
3. Compare POC and internal quantities: `p_poc`/`q_poc`, `p_conv`, `p_loss`,
   `v_int_mag`, phase currents, `i_neutral`, `ripple`, and `dv2`.
4. Stress the model with voltage unbalance, weak-grid voltage excursions, and DC-link
   variations. A balanced nominal case is a useful sanity check, not a validation of
   four-wire capability.
5. Repeat solutions from more than one sensible initial point for critical cases.
   This is a smooth nonlinear prototype, so a solver status alone is not proof of a
   globally best operating point.

The solver accepts SI or per-unit networks, but reports `InverterResult` quantities
in SI in either case. Make the chosen network scaling explicit when comparing cases,
especially alongside DC-link capacitance and voltage-ripple limits.

## What to report in a research study

At minimum, state the POC/control abstraction, filter parameters, which terminal
each rating applies to, bridge topology, DC voltage and capacitance, modulation
assumption, neutral-current rating, and ripple criterion. Then report which
constraint binds in each claimed capability result. This makes it possible to tell a
genuine converter limitation from a network-voltage, POC-contract, or control-law
limitation.

For the complete API and the underlying equations, see the
[advanced inverter component reference](@ref AdvancedInverter).
