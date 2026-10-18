# [Choosing inverter and controller parameters](@id ibr-parameter-selection)

Every number passed to [`AdvancedInverter`](@ref) or a
[`SequenceController`](@ref) is a modelling claim about hardware or firmware.
This page explains how to choose each one so that a reader can audit it. Each
parameter has one of three kinds of justification:

1. **derived** from the model's own equations and a stated requirement, as
   below (preferred);
2. **sourced** from a cited design reference or standard; or
3. **assumed**, which then needs a sensitivity sweep in the study.

Bibliography entries marked † in [IBR references](@ref ibr-references) have not
been checked against live sources, and standard clauses have not been read
from the licensed text. Treat every value attributed to them here as
*to be verified* before publication.

## Work in per unit of each device's rating

For a three-phase device with rating ``S`` and nominal phase-to-neutral voltage
``V_n``,

```math
I_r=\frac{S}{3V_n},\qquad Z_b=\frac{V_n}{I_r}=\frac{3V_n^2}{S},\qquad
Y_b=\frac1{Z_b},\qquad C_b=\frac1{\omega Z_b}.
```

A 6 kVA, 230 V inverter has ``I_r=8.70`` A, ``Z_b=26.5`` Ω and
``C_b=120`` µF. Specify filters, gains and current headroom as fractions of
these bases and convert per device. A fixed ohmic filter or a fixed A/V gain
applied across a fleet of mixed ratings silently gives small and large
inverters different physics and different control effort.

## Output filter

The filter sets the internal-voltage drop, the converter-side carrier ripple and
the LCL resonance. Liserre, Blaabjerg and Hansen give the usual design
constraints (reference 14):

- total series inductance at most about 0.1 pu, which bounds the fundamental
  voltage drop and hence the DC voltage needed;
- midpoint capacitance at most about 5% of ``C_b``, which bounds the capacitor's
  reactive power at fundamental frequency;
- resonance between roughly ``10f`` and ``f_{sw}/2``; and
- converter-side inductance large enough to keep the carrier ripple current
  within the semiconductor and inductor ratings.

Check the last constraint with the model itself: enable `pwm_ac_ripple=true`
and read `i_ac_switching_rms` against the leg rating rather than relying on a
closed-form ripple estimate. `filter_resonance_hz` reports the undamped LCL
estimate.

For a reduced single-arm study, ``X_f=0.05``–``0.10`` pu with
``X/R\approx10``–``30`` is a defensible **assumption**: inductor resistance at
line frequency is small, but not zero. For the 6 kVA example, 5% is
``X_f=1.32`` Ω (4.2 mH). A value an order of magnitude smaller makes the
internal voltage essentially equal to the POC voltage and removes the filter
from the study; state that if it is intended. [`AdvancedInverter`](@ref) warns
when a declared filter is below ``10^{-3}`` pu.

## DC-link voltage

From the switching hull in [Scientific foundations](@ref
ibr-scientific-foundations), a balanced internal RMS voltage ``|U|`` needs

```math
V_{dc}\ge\frac{\sqrt6\,|U|}{m}\quad(\text{three/four-leg}),\qquad
V_{dc}\ge\frac{2\sqrt2\,|U|}{m}\quad(\text{split link}).
```

Size ``|U|`` at the highest supply voltage the inverter must operate at plus the
filter drop at rated current. At 253 V (230 V +10%), a 5% filter, full reactive
absorption and ``m=0.96``, ``|U|\approx266`` V. That requires ``V_{dc}\ge679`` V for a
three-leg bridge and ``784`` V for a split link (**derived**). Unbalance increases
the largest line-to-line or phase voltage; `switching_margin` reports the
remaining headroom. Use the bus voltage of the actual product when it is known,
and treat the model's limit as a check rather than a design.

## DC-link capacitance and ripple-current rating

Capacitance and ripple-current rating are separate constraints. For a
three-wire converter the 2ω power is

```math
\widetilde S=3(U_1I_2+U_2I_1),\qquad
C\ge\frac{|\widetilde S|}{2\omega V_{dc}\,\Delta V_{2}},\qquad
I_{2\omega,rms}=\frac{|\widetilde S|}{\sqrt2V_{dc}},
```

for an open line-frequency DC source (**derived**; the source-branch
generalisation is in the foundations page). The controller choice dominates
this requirement. Take the 6 kVA example at 1% voltage unbalance with
``V_{dc}=750`` V and a 1% (7.5 V peak) ripple allowance:

| operating point | ``\lvert\widetilde S\rvert`` | required ``C`` | ``I_{2\omega,rms}`` |
| --- | ---: | ---: | ---: |
| balanced current (``I_2=0``), rated ``I_1`` | ≈ 60 VA | ≈ 17 µF | 0.06 A |
| plus ``I_2=3`` A of negative-sequence support | up to ≈ 2.1 kVA | up to ≈ 600 µF | up to 2.0 A |

A three-phase product sized for balanced operation therefore cannot be assumed
to tolerate a negative-sequence law. Studies should either
(a) use a **sourced** capacitance and ripple rating for a named product class and
declare `dv2_max` and `i_cap_max`, or (b) treat capacitance as the hardware
counterfactual and sweep it (see [Fleet and hardware studies](network_control_studies.md)). A
capacitance chosen so large that `dv2` never binds removes the trade-off this
model is built to expose; if that is deliberate, say so. Switching ripple,
hold-up energy and control requirements add further lower bounds that the 2ω
formula does not cover.

## Current rating and headroom

`i_max` is the per-leg semiconductor (or inductor) rating, distinct from
``I_r=S/(3V_n)``. To deliver rated apparent power at the lowest operating
voltage ``V_{min}``, it must satisfy ``i_{max}\ge S/(3V_{min})``. At 207 V this is
``1.11I_r`` (**derived**). A headroom of 10–20% is a common **assumption**;
anything above the supply-voltage requirement is current available for
negative-sequence support. With `SequencePriorityLimiter(priority=:power)` a
phase can always carry at least ``I_{max}-|I_1|`` of negative-sequence current.

## Losses

Fit `p_loss_fixed`, `a_loss` and `c_loss` to a datasheet efficiency curve at
several load points (for example 10, 20, 30, 50, 75 and 100% of rating) by least
squares on ``P_{loss}=P_0+a\sum|I|+c\sum|I|^2``. Report the fit residual.
Losses now enter the DC-availability limit, so a lossless plant overstates
delivered energy by the full converter loss. Capacitor ESR (`esr_dc*`) and its
frequency weights come from the capacitor datasheet.

## Volt-var and Volt-watt curves

Use the curve and the *basis* of the jurisdiction being studied. The AS/NZS
4777.2:2020 Australia A settings are commonly reported as follows; the
remaining regional settings differ (**sourced, to be verified** against the
licensed text, reference 47):

| curve | knots (V) | ordinates | basis |
| --- | --- | --- | --- |
| Volt-var | 207, 220, 240, 258 | +44%, 0, 0, −60% | rated VA (`q_scale = S`) |
| Volt-watt | 253, 260 | 100%, 20% | rated power (`volt_watt_basis=:rated`) |

A rated basis caps the curve ordinate at availability. A Volt-watt curve
applied to *available* power, or with an earlier or lower knee, curtails more
than the standard at partial irradiance. Such curves are legitimate
sensitivity cases but should not be labelled as the standard.

## Negative-sequence gain and angle

The static attenuation result in [Phase-aware local control laws](@ref
ibr-phase-aware-control-laws) gives the design rule:

```math
\frac{|U_2|}{|E_2|}=\frac1{|1+g\,e^{j\delta}|},\qquad
g=\kappa_{tot}|Z_2|,\qquad \delta=\arg Z_2-\phi_2.
```

1. **Measure ``Z_2``.** Use the feeder's negative-sequence Thévenin impedance at
   the connection point, for example from
   [`inverter_control_network_voltage_sensitivity`](@ref) in sequence
   coordinates. Do not assume a typical LV angle: a short feeder behind a
   transformer is much more inductive than the cable alone.
2. **Choose the angle.** ``\phi_2=\arg Z_2`` aligns the response (``\delta=0``); any
   ``|\delta|>90°`` can amplify unbalance. Report ``\delta`` for every device.
3. **Choose the attenuation.** For a target ratio ``a<1``, ``g=1/a-1`` when
   aligned; for example ``a=0.8`` needs ``g=0.25``. With a common per-inverter
   gain ``\kappa``, bus ``i`` sees ``g_i=\kappa\left|\sum_jZ_{2,ij}\right|``, where
   the sum over the fleet's measured mutual impedances already accounts for the
   number of inverters. Do not divide by ``N`` again: for ``N`` inverters behind
   one shared impedance this is ``g=N\kappa|Z_2|``.
4. **Check the current.** ``|I_2|=\kappa|U_2|`` must fit the headroom above.
   Express ``\kappa`` in per unit of ``Y_b`` so that heterogeneous fleets receive
   comparable effort.

A gain chosen without these steps can be a no-op (``g\ll0.1``), saturate every
inverter, or point in the wrong direction. Each of these is a different study
conclusion.

## Sensing reference

Declare what the product measures: `sensing_reference=:device` (the plant's
neutral, or ground for a three-wire plant), `:ground`, or a terminal name such
as `"n"` for a three-wire power stage with a neutral sensing lead. Worst-phase
laws respond to the zero-sequence voltage to first order and mean-phase laws to
second order, so this choice changes results on feeders with neutral
displacement.

## Smoothing widths

Smoothing widths are approximation budgets, not tuning knobs. The defaults
(extrema 0.05 V, curve corners 0.05 V, current widths ``2.5\times10^{-5}`` of the
rating, power widths ``5\times10^{-8}`` of the rating) keep the exact-versus-smooth
current residual far below 0.02 A in the regression suite. Use identical widths
for every variant in a matched comparison, record them, and show that halving
them does not change the conclusions.

## What to record

For each study, give a table with: parameter, value, kind (derived / sourced /
assumed), source or derivation, and the sensitivity run that covers it. The
[Australian LV study](https://github.com/frederikgeth/PowerOptLab.jl/tree/main/studies/australian_pv)
README follows this format.
