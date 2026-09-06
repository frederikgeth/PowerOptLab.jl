# Understanding generalized generator models

These worked examples accompany the [scientific model](../components/generalized_generator.md).
They isolate elementary RMS phasor circuits. Follow them with the
[implemented network tutorial](generalized_generator_tradeoffs.md) to compare
the actual component modes under load and network unbalance.
Run their reproducible calculations from the repository root:

```sh
julia --project=. scripts/generalized_generators/analytic_examples.jl
```

The script uses only Julia standard libraries. Its assertions check the numerical
identities below; the separate component unit tests cover BMOPFTools integration.
Neither is a dynamic simulation.
All impedances are ohms, voltages volts RMS, and currents amperes RMS.

## 1. A balanced internal source can have unbalanced terminal voltage

Imagine three balanced internal EMFs feeding a network through a series
impedance. The following sketch shows each phase of a simple wye case; use the
full circuit when a neutral or earth impedance is present.

```text
internal EMF e_k ---- series impedance ---- PCC voltage u_k ---- network
                         current i_k -->
```

With ``a=\exp(\mathrm{j}2\pi/3)``, choose internal source
``e=230[1,a^2,a]``. Use sequence impedances in order zero, positive, negative:

```math
(Z_0,Z_1,Z_2)=(0.4+0.3\mathrm{j},\ 0.2+0.4\mathrm{j},\ 0.1+0.2\mathrm{j}),
\qquad (I_0,I_1,I_2)=(0,10,2).
```

From ``u=e-Zi`` the terminal sequences are
``U_1=228-4\mathrm{j}`` and ``U_2=-0.2-0.4\mathrm{j}``. Thus
``|U_1|\approx228.035`` V and ``|U_2|\approx0.447`` V. The balanced source's
internal negative-sequence voltage is zero, but the PCC negative-sequence voltage
is nonzero. Passive negative-sequence current is not the same as a controlled
negative-sequence EMF.

The same calculation gives a source-to-terminal active-power difference of
61.2 W, equal to ``i^\mathsf{H}Hi``. This is an independent power-ledger check.

**Interpretation:** a fixed internal voltage magnitude does not regulate terminal
voltage. For this imposed-current example, restoring ``|U_1|=230`` V with real
positive internal ``E_1`` requires

```math
E_1=\Re(Z_1I_1)+\sqrt{230^2-[\Im(Z_1I_1)]^2}
\approx231.965\ \mathrm{V}.
```

This isolates the distinction between fixed excitation and terminal regulation.
It is not a complete PV power-flow solution: in a network solve, current and
source angle must also satisfy KCL and the selected active-power law.

## 2. A source's four conductor currents need not sum to zero

A grounded source needs a fifth current path in its accounting:

```text
                         phases a, b, c --> network
                        /
        EMFs referenced to internal neutral *
                                          | \
                            neutral lead  |  \ grounding impedance
                                          |   \
                                     PCC neutral   earth
```

Use a real-valued phasor example. Let the sum of outward phase currents be
10 A, the neutral lead resistance be 0.2 ohm, the grounding resistance 0.8 ohm,
and both remote earth and the PCC neutral potential be zero. The two return
paths share current. Internal neutral KCL gives

```math
10+v_*/0.2+v_*/0.8=0,
\quad v_*=-1.6\ \mathrm{V},
\quad I_n=-8\ \mathrm{A},\quad I_g=-2\ \mathrm{A}.
```

The four-conductor sum is **2 A**, while the sum including earth is zero.
The source is receiving 8 A through neutral and 2 A through earth. For example,
take ``I_a=10`` A and ``I_b=I_c=0``. Then ``I_0=10/3`` A and the identity is
``3I_0+I_n+I_g=0``. Inferring ``I_n=-3I_0`` would give the wrong neutral current.

With phase-a terminal voltage 230 V and zero phase-series drop, its internal
EMF relative to the shifted star is 231.6 V. Source power is 2316 W, terminal
power is 2300 W, and the difference is
``0.2(8)^2+0.8(2)^2=16`` W of neutral and grounding losses.

Opening the earth connection sets ``I_g=0`` and recovers zero four-conductor
sum. An ideal ground fixes star voltage and lets KCL determine earth current.
Neither operation defines the system's phasor angle reference.

## 3. Independent sequence limits do not protect every phase

Set ``I_0=0``, ``I_1=6`` A, and ``I_2=6`` A, with the two nonzero sequences
aligned in phase a. The inverse Fortescue transform gives

```math
(I_a,I_b,I_c)=(12,-6,-6)\ \mathrm{A}.
```

Both nonzero sequences fit an 8 A sequence limit. Phase a violates a 10 A
conductor limit. Keep the phase limits even when sequence limits are independent
configuration fields. Neutral and grounding ratings need their own constraints.

Similarly, fixed relative angles do not guarantee balanced magnitude: the
template ``[230,220a^2,240a]`` has both zero- and negative-sequence magnitude
``10/\sqrt3\approx5.774`` V despite exact 120-degree phase separation.

## 4. Relative rotation is different from an absolute reference

For a fixed-magnitude rotating template, the allowed EMFs are
``e(\delta)=\exp(\mathrm{j}\delta)\bar e``. That is one real source-voltage
degree of freedom. A common free magnitude adds one more; independent magnitudes
add one per independent phase port.

Rotating **every** voltage and current in an otherwise rotation-invariant
connected model by the same angle changes no voltage magnitude, power, or
impedance relation. Rotating only one generator changes its angle relative to
the network and generally changes its power transfer. Fix a global reference
without fixing every generator's physically meaningful relative angle.

The script verifies invariance of the circuit drop and complex power under a
common rotation. This is why a numerical angle reference and a stiff physical
source should be separate concepts.

## 5. Split phase is a connection, not three phases with one missing

For terminals ordered ``(L_1,L_2,N)``, the two leg-to-neutral port map is

```math
C=\begin{bmatrix}1&0\\0&1\\-1&-1\end{bmatrix}.
```

Choose ``v=(120,-120,0)`` V and port injections ``i=(10,-6)`` A. Then
``j=Ci=(10,-6,-4)`` A and the total power is 1920 W. The legs are 120 V RMS,
but the line-to-line voltage is 240 V RMS. The script verifies both port and
conductor power accounting. If the source also has an earth connection, extend
this circuit as in example 2; do not force all return current through N.

The phase convention is ``(0,\pi)``, and ordinary three-phase positive/negative/
zero-sequence labels do not apply. A different modal transform can be supported
later with explicit definitions, normalization, and a physical interpretation.

For delta, also check winding-cycle consistency. The unequal 120-degree template
in example 3 does not sum to zero, so it cannot be three ideal voltage differences
around a closed delta. The implementation validates ideal-loop closure. Finite loop impedance can
support that template by carrying circulating winding current. The
[network tutorial](generalized_generator_tradeoffs.md) includes delta cases.

## From examples to a scientific study

Start with one network and change one assumption at a time: ideal versus finite
source, fixed internal EMF versus terminal positive-sequence regulation, a single
neutral return versus a grounded source, then sequence-only versus full conductor
capability. Preserve the same demands and physical ratings where the comparison
permits it. Distinguish a change of physical feasible region from a change of
algebraic formulation of the same region.

Measure PCC phase/sequence voltage, neutral/earth current, source/terminal power,
active limits, physical residuals, and the decision of interest. When studying
controllers, freeze their specified laws before comparing responses. The
[implementation acceptance plan](../components/generalized_generator.md#Implementation-increments-and-acceptance-evidence)
lists the further tests needed before these examples become runtime claims.
