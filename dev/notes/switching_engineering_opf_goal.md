# Switching-aware engineering and matched Australian OPF

Active objective: complete AC phase/neutral RMS and declared peak-current
limits, combine carrier ripple with nonlinear DC voltage bounds, then benchmark
mixed local controllers against coordinated OPF under identical physical and
network constraints. Hardware sizing stays fixed. The completed thermal/energy
campaign at `cdc19d2` is the baseline, not completion of this objective.

## Requirements and evidence to deliver

1. Carrier-event AC/DC diagnostics without a truncated switching spectrum.
   Independent analytic and time-domain checks must cover P/Q, sequence
   magnitudes/angles and modulation over a declared operating domain.
2. Enforced phase/neutral total RMS and optional declared peak limits, and
   combined LF/carrier DC voltage extrema and modulation. Apply the same
   constraints in local service allocation and shared-network OPF. Separate
   mathematical bounds within the carrier approximation from calibrated
   allowances and errors of that approximation.
3. Predeclared matched Australian cases with mean, worst-phase, sequence and
   four-leg local control and coordinated dispatch. Use identical sites,
   ratings, availability, loads, source, thermal/energy models and network
   limits. Save failed results as well as successes. Do not claim a global
   nonconvex OPF bound.
4. Independent saved-phasor/network audits, representative averaged/switched
   EMT including current and DC-voltage extrema, corruption gates and source
   fingerprints. Run affected tests, the complete package suite and docs.
5. Reproducible tables and pedagogical figures showing binding constraints,
   achieved services, losses, curtailment and the benefit of coordination.
   Commit implementation and validation checkpoints on the existing branch.

## Current implementation checkpoint

The initial carrier-event calculation integrates exact gate-event intervals
and piecewise-affine primitives for lossless inductive AC ripple and
frozen-current DC charge. It retains physical neutral coupling. Separate
all-duty current bounds cover positive inductance and proportional passive
R/L. They use a bounded-forcing primitive and variance inequality, not the
existing empirical 10% DC-thermal allowance.

These calculations freeze the slow bus and reference over each carrier. They
do not prove a bound for a coupled time-varying converter or an LCL resonance.
The diagnostics are internal foundations; enforcement, declared peak fields,
combined DC bounds, domain sweeps, the matched OPF campaign and final validation
remain required. This file is progress tracking, not an acceptance certificate.

Initial verification passes 895 assertions: a closed-form square/triangle
orbit, dimensional scaling, exact-event versus independent 16,384-step carrier
integration over 108 operating points, all-duty bounds for phase/neutral RMS
and peak, independent periodic damped-RL orbits, and unsupported-domain checks.
The existing PWM allocation and nonlinear energy regressions also pass all
242 assertions. The full package suite is reserved for the integrated milestone;
no previously published campaign is relabeled with the new source fingerprints.
