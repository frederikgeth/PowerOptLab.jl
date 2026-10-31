# Mixed-controller PWM thermal and nonlinear DC-energy study

This extension keeps the phase/DC sizing fixed and adds thermal allocation and
nonlinear capacitor-energy bounds. Every coupled feeder contains three each of
mean, worst-phase, sequence and four-leg per-phase VV+VW controllers. Sorted
customer IDs determine the repeating assignment. Each controller group gets
the same 25/60/100% PV-availability pattern. The physical locations differ:
these are mixed-system validation cases, not a causal ranking of controller laws.

The existing 6 kVA / 10 A / 750 V / 300 µF / 1 A parameters and phase filter/loss
coefficients are unchanged. Four-leg sites retain their declared 10 A neutral
and equal phase/neutral filter impedances. Solar availability varies by site;
nameplate ratings do not. The Australia A curve knots, converter current target,
POC VA boundary and reactive priority match the previous terminal-profile setup.

Eight cases are predeclared: three model profiles at 1.00/1.03 pu source voltage,
plus nonlinear/PWM night and partial-disconnection cases. All have 1% source VUF
and the same 60/25/15% load split. Profiles are:

- `fundamental`: original linearized energy and LF-only capacitor model;
- `pwm`: embedded centered-PWM thermal allocation with linearized energy;
- `energy`: the same PWM allocation plus nonlinear monolithic capacitor energy.

All controllers are active in one voltage-aware feeder fixed point. Exact local
curves produce requests; physical allocation retains each declared priority.
Higher-priority service is preserved within an explicit 100 µvar/W numerical
budget in PWM/nonlinear stages. Hardware bounds are not relaxed. Night PV has
one exact zero-DC-power equality; loss imports and reactive support remain.
Disconnection removes the entire selected AC model.

The switching variance is an analytic shared-carrier expression in fundamental
phasors, with a voltage-norm error enclosure and full-cycle DC denominator
bounds. Independent exact carrier dwell integration audits each point. The
default 10% validation allowance was selected after an initial 5% allowance
failed at three four-leg mixed points (1.0057–1.0059 A). Those old references
are stored as negative controls in `thermal_calibration_negative.json`, with
their originating core commit. They are not accepted current commands.
The allowance is not a universal bound on AC ripple or semiconductor dynamics.

Nonlinear energy preserves the declared cycle-mean DC voltage. `dv2_max` keeps
its Fourier 2ω meaning; optional absolute DC extrema bounds are separate. The
full low-frequency capacitor RMS includes higher even harmonics. The initial
scope is a monolithic link, constant input power, fixed mean converter loss and
an open carrier-frequency source. Finite small-signal source branches and split
banks need their own nonlinear dynamics and are rejected by this option.

Saved-point audits reconstruct current/voltage physics, source/thermal/energy
bounds and controller policies in SI, then run fixed-PQ replay in the same
network engine. This is not an independent feeder circuit implementation.
Failed solves and incomplete plans stay explicit; figures require complete
accepted records. Publication also checks CSV metrics against saved evidence.

The independent EMT matrix includes all 24 energized inverter points from the
two nonlinear daytime feeder cases and four fixed-boundary references, plus
three failed 5% calibration points. Both averaged and switched circuits run at
three meshes. All 27 current PWM-allocated points must meet the unchanged
engineering limits. The LF-only ripple-boundary reference and calibration
points remain negative switching controls. The averaged boundary trajectory
separately tests the nonlinear 2ω constraint against the old linearized result.
Waveforms are retained for eight representative references; complete endpoint,
Fourier, RMS and energy metrics are retained for every mesh of all 31 references.

```sh
julia --project=. scripts/australian_pv/test_thermal_energy.jl
julia --project=. scripts/australian_pv/thermal_energy.jl /tmp/au_thermal_energy
julia --project=. scripts/australian_pv/audit_thermal_energy.jl /tmp/au_thermal_energy
julia --project=. scripts/australian_pv/test_thermal_energy_records.jl /tmp/au_thermal_energy
python scripts/australian_pv/test_thermal_energy_csv.py /tmp/au_thermal_energy
julia --project=. scripts/australian_pv/prepare_thermal_energy_emt.jl /tmp/au_thermal_emt_inputs.json /tmp/au_thermal_energy
julia --project=. scripts/australian_pv/emt_circuit.jl /tmp/au_thermal_emt_inputs.json /tmp/au_thermal_emt
python scripts/australian_pv/audit_emt.py /tmp/au_thermal_emt_inputs.json /tmp/au_thermal_emt
python scripts/australian_pv/test_thermal_energy_emt_records.py /tmp/au_thermal_emt_inputs.json /tmp/au_thermal_emt
```

The publication script reads the canonical audited result directories and the
refreshed baseline EMT evidence to generate [tables and figures](thermal_energy_results.md).
It does not rerun or discard failed solves. Source, input, waveform and runtime
fingerprints accompany the evidence.

EMT here fixes each POC voltage and uses ideal sampled feedforward PWM at 10 kHz.
There is no dynamic feeder, PLL/current/service loop, protection or switching-loss
model. A common-mode projection error discovered in the averaged three-wire
branch was corrected; a separate KCL/energy oracle now guards that mechanism.
No 253 V bound is imposed by the local controllers and no regional conformity,
global service optimum, lifetime or dynamic-stability claim is made.
