# [Advanced inverter modelling: when the POC IBR is not enough](@id advanced-inverter-modelling)

This tutorial helps choose between PowerOptLab's ordinary inverter-based-resource
(IBR) model and `AdvancedInverter`. The distinction is important: an ordinary IBR
is a controlled, bounded current injection at the point of connection (POC), while
an advanced inverter makes the converter's internal AC node and DC-side constraints
explicit.

For many distribution studies, the ordinary IBR is the more defensible model. Use
the advanced model only when the omitted internal physics changes the conclusion.

After this overview, continue with the focused studies on [topology under
unbalance](@ref ibr-topology-under-unbalance), [finite DC-source and split-link
carrier stress](@ref ibr-dc-source-and-split-link), and [carrier harmonics through
L/LCL filters](@ref ibr-ac-harmonics-lcl).

## Model-selection guide

| Study question | Recommended model | Why |
| --- | --- | --- |
| Feasible POC injections, fixed power factor, or mandatory Volt-VAr / Volt-Watt behaviour | Ordinary IBR | The control law and limits are imposed at the POC; an internal converter circuit is not needed. |
| Hosting capacity or a DOE where DERs are represented by POC setpoints | Ordinary IBR | It is compact, transparent, and matches the available control interface. |
| Do a series filter or explicit fundamental-frequency LCL network change converter voltage, arm currents, reactive exchange, or loss? | `AdvancedInverter` | The filter separates converter, midpoint, and POC quantities. Dynamic resonance and control interaction remain outside this model. |
| Can a grid-forming converter sustain its requested internal EMF? | `AdvancedInverter` | Internal voltage, modulation, and DC limits must be represented. |
| Is a four-wire unit limited by neutral current or DC-link ripple? | `AdvancedInverter` | These are topology- and DC-capacitor-dependent constraints. |

The ordinary IBR should not be treated as a low-fidelity version of every physical
converter. It is a different abstraction: a POC-level resource model. In particular,
it is often the right choice for a PV inverter following a prescribed Q-V curve. The
advanced model exposes electrical capability; it does not automatically reproduce a
manufacturer's supervisory or mandatory grid-support control law.

## The internal-node model

`AdvancedInverter` inserts an internal AC node behind an optional output filter.
The reduced circuit is

```text
POC bus ──[ r + jx, optional grid-side shunt ]── internal node ── converter ── DC link
```

Converter current and apparent-power limits apply at the internal node. Therefore,
the converter's rating is not necessarily the POC export rating. With a non-zero
filter impedance, part of the converter power supplies filter losses and reactive
exchange; the POC voltage can also differ materially from internal voltage.

Use scalar `r_filter`/`x_filter` and optional neutral values for identical,
uncoupled conductors. When the cable or filter data are supplied as a primitive,
pass the complete matrix instead. For terminal order `a,b,c,n`:

```julia
R = [0.05 0.004 0.003 0.006;
     0.004 0.05 0.004 0.006;
     0.003 0.004 0.05 0.006;
     0.006 0.006 0.006 0.12]
X = [0.15 0.020 0.018 0.025;
     0.020 0.15 0.020 0.025;
     0.018 0.020 0.15 0.025;
     0.025 0.025 0.025 0.22]

inv = AdvancedInverter(
    id = "four-wire",
    bus = "lv_bus",
    phase_terminals = ["a", "b", "c"],
    neutral = "n",
    topology = :FOUR_LEG,
    s_max = 30_000.0,
    In_max = 35.0,
    v_dc = 750.0,
    c_dc = 8e-3,
    r_filter_matrix = R,
    x_filter_matrix = X,
)
```

The matrix conductor order must exactly match the terminal order followed by the
neutral. Resistance must be symmetric positive semidefinite; reactance must be
symmetric. Scalar and matrix parameterisations cannot be mixed, which prevents
an easy-to-miss double count of the diagonal terms.

### When to use the explicit LCL midpoint

Use the explicit midpoint when capacitor current, separate converter/grid
inductor ratings, damping loss, or fundamental reactive exchange can affect the
operating point:

```julia
lcl = AdvancedInverter(
    id = "four-wire-lcl",
    bus = "lv_bus",
    phase_terminals = ["a", "b", "c"],
    neutral = "n",
    topology = :FOUR_LEG,
    s_max = 30_000.0,
    i_max = 60.0,          # converter-side arm
    i_grid_max = 55.0,     # grid-side arm
    In_max = 35.0,
    v_dc = 750.0,
    c_dc = 8e-3,
    r_filter = 0.02,
    x_filter = 0.06,
    r_filter_grid = 0.03,
    x_filter_grid = 0.09,
    c_filter_mid = 30e-6,
    r_filter_damping = 0.5,
)

r = solve_advanced_inverter(network, lcl)
@show r.v_int_mag r.v_filter_mag
@show r.i_mag r.i_grid_mag r.i_filter_shunt_mag
@show r.p_filter_loss r.filter_resonance_hz
```

With `c_filter_mid = 0`, converter and grid currents are identical and the two
series arms reduce to their sum. With a capacitor, `i_mag` is the semiconductor/
converter-inductor current while `i_grid_mag` reaches the POC. For a three-leg
converter, a grounded-wye filter capacitor can exchange zero-sequence current
with the grid neutral even though the converter still enforces zero-sequence
current at its own terminals.

Do not interpret `filter_resonance_hz` as a stability result. It is the undamped
scalar LCL estimate. Matrix-valued arms require modal analysis, and controller
delay, PLL/GFM control, active damping, grid impedance, ESR/ESL, and switching
effects require a frequency-sweep or dynamic model.

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
    r_filter_neutral = 0.01,
    x_filter_neutral = 0.04,
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
For `:SINGLE_PHASE`, PowerOptLab currently retains a scalar `/√3` convention;
it is not a bridge-specific full- or half-bridge model, so treat it as a user-supplied
engineering cap. For three-phase units, PowerOptLab uses a sampled outer
approximation of the topology-specific continuous-time switching hull.

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
    i_zero_max = 10.0,
    i_negative_max = 8.0,
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
the conclusion depends on modulation headroom. Also inspect
`result.switching_margin`: a negative value means the independent dense audit
found a rail violation between optimisation samples. A positive value is useful
numerical headroom, but is not an analytic continuous-time certificate.

### Three-leg, four-leg, and split-DC are not interchangeable

`:THREE_LEG` is a three-leg, three-wire bridge with no zero-sequence neutral-current
path. `:FOUR_LEG` is a four-leg bridge with a monolithic DC link and bounded
fourth-leg current. `:SPLIT_DC` is specifically a **three-leg, four-wire
split-capacitor** bridge: it can serve four-wire loads, but midpoint utilisation
reduces modulation headroom and couples unbalance to capacitor stress. A hybrid
fourth-leg-plus-split-link converter is a different topology and is not represented.

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
    r_filter = 0.05,
    x_filter = 0.15,
    v_dc = 800.0,
    c_dc = 5e-3,
    c_dc_upper = 4.5e-3,
    c_dc_lower = 5.5e-3,
    dv2_max = 20.0,
    dv_mid_max = 5.0,
    v_mid_mean_max = 5.0,
    q_mid_balance_max = 0.25,
    i_cap_upper_max = 12.0,
    i_cap_lower_max = 14.0,
    cap_thermal_weights = (1.4, 1.0, 0.7), # neutral, 2ω, switching ESR ratios
    esr_dc_upper = 0.030,
    esr_dc_lower = 0.025,
    i_sw = 2.0,                  # independent residual allowance
    f_sw = 12e3,
    pwm_strategy = :SPWM,        # split midpoint excludes centered injection
    pwm_dc_source_r = 0.08,      # optional upstream carrier-frequency R–L
    pwm_dc_source_l = 80e-6,
    pwm_dc_harmonics = 64,
    pwm_ac_ripple = true,
)

result = solve_advanced_inverter(network, ripple_limited;
    objective = :max_export,
)

@show result.i_neutral result.i_zero result.i_negative
@show result.ripple result.dv2 result.dv_mid result.v_mid_mean
@show result.i_cap_upper result.i_cap_lower
@show result.i_cap_thermal_upper result.i_cap_thermal_lower
@show result.q_mid_balance result.p_cap_loss result.switching_margin
@show result.i_cap_switching result.i_cap_switching_reserved
@show result.i_dc_bridge_switching_rms result.i_dc_source_switching_rms
@show result.p_dc_source_switching_loss result.pwm_dc_network_margin
@show result.dv_switching_rms result.dv_switching_pp
@show result.dv_switching_upper_rms result.dv_switching_lower_rms
@show result.pwm_reserve_margin result.pwm_modulation_margin result.pwm_iterations
@show result.i_ac_switching_rms result.i_neutral_switching_rms
@show result.i_ac_total_rms result.i_neutral_total_rms
```

With balanced voltage and current, a three-phase bridge has little low-frequency
DC-link pulsation. In contrast, phase imbalance can create neutral current and a
substantial 2ω ripple. A smaller DC capacitor increases `dv2`; a binding `dv2_max`
can reduce feasible export even when RMS current and apparent-power limits look
comfortable.

`result.ripple` is the sinusoidal 2ω power amplitude. The physical upper/lower
RMS currents include capacitance-proportional neutral shares, the common 2ω
component `result.ripple/(√2*v_dc)`, and `i_sw`. Their `i_cap_thermal_*`
counterparts apply the squared-current ESR-ratio weights and are the quantities
checked against the ratings. `p_cap_loss` converts those weighted currents back
to watts using the supplied reference ESRs.

With PWM enabled, `i_cap_switching` is reconstructed from ideal leg switch
states using one shared triangular carrier, while `i_cap_switching_reserved` is
the conservative current-norm allowance closed around the optimisation. The
manual `i_sw` above remains an independent residual term and combines in
quadrature. Treat the result as publishable only when `pwm_reserve_margin ≥ 0`,
`pwm_modulation_margin ≥ 0`, and a tighter carrier/fundamental sampling study
does not materially change the diagnostics. `dv_switching_rms` and
`dv_switching_pp` report the associated high-frequency DC-link voltage ripple.
With a finite `pwm_dc_source_*` branch, bridge ripple divides between the source
and capacitor at each retained harmonic. Inspect the two current diagnostics,
the separate source-resistor loss, and `pwm_dc_network_margin`; values near zero
indicate an inadequately damped source-capacitor antiresonance. For a split link,
the upper/lower voltage diagnostics expose the larger ripple across the smaller
half-bank. Increase `pwm_dc_harmonics` and `pwm_carrier_samples` together before
using these quantities for component selection.

With `pwm_ac_ripple=true`, the pole-voltage carrier harmonics are also passed
through the phase and neutral conductor inductances. `i_ac_switching_rms` and
`i_neutral_switching_rms` are the predicted high-frequency components;
`i_ac_total_rms` and `i_neutral_total_rms` combine them in quadrature with the
fundamental solution and are the quantities that should be compared with
physical conductor/device ratings. For an explicit LCL filter, also inspect
`i_grid_switching_rms` and `i_filter_shunt_switching_rms` rather than assuming
the converter-side ripple reaches the grid unchanged.

Unequal capacitances create a natural mean midpoint shift even at zero neutral
current. `q_mid_balance_max` supplies bounded quasi-static charge authority to
correct it, while `v_mid_mean_max` states the admissible residual. This is a
steady-state feasibility abstraction: it does not demonstrate that a physical
balancing controller has adequate bandwidth or remains within its transient
voltage limits.

The RMS distinction matters because some papers tabulate Fourier/phasor
amplitudes as "ripple current" without the `√2` conversion.

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
   `v_int_mag`, `v_filter_mag`, converter/grid/shunt currents, symmetrical-
   component currents, `i_neutral`, `p_filter_loss`, `ripple`, `dv2`, `dv_mid`,
   `v_mid_mean`, physical and thermal capacitor currents, `p_cap_loss`,
   `filter_resonance_hz`, and `switching_margin`.
4. Stress the model with voltage unbalance, weak-grid voltage excursions, and DC-link
   variations. A balanced nominal case is a useful sanity check, not a validation of
   four-wire capability.
5. Repeat solutions from more than one sensible initial point for critical cases.
   This is a smooth nonlinear model, so a solver status alone is not proof of a
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
[advanced inverter component reference](@ref AdvancedInverter) and
[scientific foundations](@ref ibr-scientific-foundations).
