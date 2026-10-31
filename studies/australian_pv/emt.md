# Representative time-domain circuit validation

> **Historical campaign.** These results describe a pinned earlier model revision. See the [current benchmark](engineering_benchmark.md) and [evidence/replay index](evidence/README.md).

This extension integrates a physical bridge, L-filter currents and DC-capacitor
voltage. It complements the frozen-phasor carrier audits in the package. Three
fixed-boundary cases use the unchanged Australian phase/DC sizing:

1. Balanced 230 V, maximum available export.
2. A declared synthetic 2% negative-sequence boundary, −1800 var Q request,
   70% rated VW ceiling and a 2 A negative-sequence request. The staged allocator
   reaches a binding DC-link ripple boundary.
3. The saved Australian four-leg 2% VUF point at the customer with greatest
   neutral current. Its explicit 10 A neutral and equal phase/neutral filter
   impedances are the same assumptions as the four-leg study.

POC voltage is fixed at each accepted phasor point. The network and the VV+VW
service controller are not simulated dynamically. The independent circuit file
imports JSON/SHA only; it does not call PowerOptLab's optimizer, filter stamps,
rail oracle or carrier replay functions.

For the four-wire bridge, eliminating the physical neutral return by KCL gives

```math
(L I + L_n 11^T)\dot i
 = v_{dc}(s_{abc}-s_n1)-v_{g,PN}-(R I+R_n11^T)i.
```

The three-wire version projects pole and grid voltages onto the sum-zero
subspace. Both integrate the DC capacitor through the actual bridge current:

```math
C\dot v_{dc}=(P_{source}-P_{loss})/v_{dc}-s^Ti_{bridge}.
```

The phase/neutral filter structure is consistent with the circuit KVL and
neutral-inductor treatment in [Mandrioli et al., Energies 2023](https://doi.org/10.3390/en16041710).
This experiment does not reproduce their complete analytical or experimental
matrix. RMS/peak and fundamental/switching conventions are explicit in the files.

The source supplies constant input power at the regulated DC-link boundary:
there is no 2ω input **power**, while source **current** can ripple. This is an
explicit temporal assumption; it does not validate the optional small-signal
R–L source branch. Loss is a fixed lumped mean from the fundamental RMS point,
not a semiconductor switching-loss or instantaneous conduction-loss model.

Centered PWM runs at 10 kHz with regular sampled sinusoidal feedforward,
carrier-midpoint sampling and measured DC compensation. The neutral is another
physical switched pole. Duty saturation is an error, never silent clipping.
The initial capacitor energy is chosen so the ideal averaged orbit has the
declared mean DC voltage; initializing voltage by a first-order sinusoid alone
would shift its mean at second order. No current feedback, PLL or service-loop
dynamics are claimed.

RK4 steps end exactly at each switch event. Each of 32/64/128 maximum substep
resolutions is run for ten line cycles, discarding six and using the final four
for Fourier and RMS quantities. Statistics use RK4 quadrature too. Stored
electrical energy is checked against source power, grid power, filter resistance
loss and the declared converter loss. Finite-window DC drift is retained in
the records; this is not a claim of periodic steady state or stability.

The input fixture declares comparison bands before execution: 1% plus 0.01 A
for fundamental current, 2% plus 0.05 V for the complex 2ω bus-voltage coefficient,
the corresponding capacitor-current band, 0.01 A / 0.01 V mesh differences and
1 ppm relative energy closure. Engineering bounds are checked separately with
unchanged ratings and explicit failure flags. A comparison can pass while a
hardware rating fails. Peak currents are reported without inventing a peak
rating from the RMS rating.

```sh
julia --project=. scripts/australian_pv/test_emt_circuit.jl
julia --project=. scripts/australian_pv/prepare_emt.jl /tmp/au_emt_inputs.json
julia --project=. scripts/australian_pv/emt_circuit.jl /tmp/au_emt_inputs.json /tmp/au_emt
python scripts/australian_pv/audit_emt.py /tmp/au_emt_inputs.json /tmp/au_emt
python scripts/australian_pv/test_emt_records.py /tmp/au_emt_inputs.json /tmp/au_emt
julia --project=. scripts/australian_pv/crosscheck_emt.jl /tmp/au_emt_inputs.json /tmp/au_emt
python scripts/australian_pv/plot_emt.py /tmp/au_emt_inputs.json /tmp/au_emt
```

Reference-source fingerprints are separate from the executed circuit digest;
the result pins both its exact input and circuit source. Waveform digests,
complete mesh/mode coverage, independent energy reconstruction from endpoint
states and RMS losses, physical flags and comparison bands are rechecked by a
separate Python audit. Corruption tests exercise these publication gates.

The measured non-2ω capacitor RMS component is mapped to the existing `i_sw`
engineering bound as a **point-specific diagnostic**. It is not a safe reserve
over changing dispatch or a substitute for a PWM-aware allocation model.
LCL/split-link dynamics, ESR/ESL, device loss/thermal dynamics, dead time,
measurement/control delays, protection and feeder stability require further
validation. [Generated results and figures](emt_results.md) retain failed
engineering checks without changing the sizing to make them pass.
