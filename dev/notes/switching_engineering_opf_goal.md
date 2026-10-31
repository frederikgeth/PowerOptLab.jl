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

## Current-limit enforcement checkpoint

`pwm_current_envelope=true` now enforces the analytical AC RMS reserve in the
shared plant used by local allocators and constrained OPF. A declared
`v_dc_max` supplies the maximum frozen carrier bus; the existing DC API requires
nonlinear energy for that absolute voltage declaration. Optional `i_peak_max`,
`i_grid_peak_max` and `In_peak_max` declare peak limits explicitly. No peak
rating is inferred from the unchanged RMS ratings. Phase/neutral limits,
independent SI physical checks, extraction and operation margins agree.

Manual RMS reserves add by Minkowski, since their spectral correlation is not
known. They cannot certify a declared peak waveform. Unavailable AC switching
bounds are `nothing` in extracted results. The new option does not yet add the
carrier contribution to the existing low-frequency DC extrema constraints;
that combined guard remains the next required implementation step.

The current-limit tests add 56 assertions covering binding converter/grid peaks
across coordinate bases, binding neutral peaks under per-phase service,
three/four-leg RMS and peak allocation, invalid configurations and unavailable
bounds. The affected operating-constraint, service, phase, feasibility and new
envelope tests pass **1,684/1,684** assertions. The 242 PWM/energy assertions
also passed the first current-limit integration run. The full documentation
build passes. No full-suite or new campaign completion is claimed yet.

Benchmark audit finding: `_voltage_aware_network_checks` reconstructs native
zero-equality rows; it is not a network-inequality certificate. The matched OPF
work must explicitly reconstruct the same voltage/conductor inequalities in
SI for local and coordinated records, in addition to network balance and replay.

## Combined DC guard checkpoint

`pwm_dc_voltage_envelope=true` adds the carrier charge reserve to nonlinear
voltage extrema, full-cycle modulation and the capacitor thermal denominators.
It rejects finite carrier-source sharing and unresolved RMS-only allowances.
The reserve separates frozen fundamental leg current from the variation of
gate/AC-ripple-current products. Its inequalities are conditional on the
periodic carrier approximation; no transient or cumulative-drift claim is made.

The first 243 DC-guard assertions pass, with 2,055 affected assertions passing
in a fresh Julia process. Tight four-leg per-phase service exposed an adaptive
barrier cycle at the 100 kVA coordinate base. Positive voltage lifts replace
nearly dependent squared energy gaps. Combined-voltage service defaults use
monotone continuation; caller overrides remain effective. One-sided voltage and
current-norm padding and the existing 100 µW fitted-loss error budget remain
explicit. Engineering limits, source acceptance and raw solver gates are not
relaxed. Tests cover 10/100/1000 kVA coordinates.

Three development references complete 18 averaged/switched integrations and
742 independent record conditions. The LF-only 747–753 V boundary reaches
753.170120 V in the switched circuit and fails its upper limit. The combined
three-leg point stays between 748.697633 and 751.504819 V; the combined four-leg
point between 749.354610 and 750.869469 V. Their capacitor RMS currents stay
below the unchanged 1 A rating. These are development probes in `/tmp`, not the
final committed campaign or the full goal's acceptance evidence.

The circuit now retains current/neutral peaks, DC extrema and their sampled
state/time witnesses. Its record audit reconstructs optional peak/absolute-
voltage flags and checks mesh convergence of extrema as well as RMS/Fourier
quantities. The original 32 circuit assertions pass. Formal reproduction,
expanded EMT/domain sweeps and the matched OPF benchmark still remain.

The integrated package suite passes **22,100/22,100** assertions (8m14.5s) and
the documentation build passes. The independent circuit now passes 53 analytic,
engineering-flag and extremum-witness assertions. The new SI network oracle
passes eight corruption/coverage tests and reconstructs all seven fully
energized historical mixed points (826 conditions each), including line and
switch limits, ground-shunt losses, delta-wye KVL, source voltage and total
active-power balance. These archived phasors are input fixtures, not relabeled
new-model studies. Maximum active-power mismatch is below 0.7 mW.

`switching_engineering.jl` centralizes fixed-hardware templates. New declared
operating ceilings are common 735–765 V DC and 14 A phase/neutral peak limits,
with a narrower 747–753 V band for boundary validation. These are explicit
study operating assumptions, not manufacturer ratings or sizing choices. Its
constructor checks every other physical field remains unchanged. The boundary
preparation driver records all stages, requests, source hashes and runtime
fingerprints. The full-domain and coordinated campaigns remain outstanding.

## Audited boundary publication checkpoint

The clean-source `693b2cb` boundary campaign is now preserved under
`studies/australian_pv/results/engineering_boundary_emt*`. It contains all 18
averaged/switched integrations and passes **747 independent record checks**.
The publication script reaudits source/input/waveform fingerprints before
regenerating its table and three-panel voltage/capacitor/peak figure. It retains
the LF-only upper-voltage failure and explicitly separates finite-window DC
drift from the periodic carrier approximation. The four-leg reference has a
different declared service request and is not labeled as a matched performance
comparison. The wider operating-domain campaign remains required.

The coordinated prototype still rejects pure-export attempts. A regularized
initialization can achieve raw `LOCALLY_SOLVED`; it is not substituted for the
intended final export objective. Failed/ALMOST/iteration-limit candidates remain
rejected. Investigation now includes primal/dual continuation and independent
candidate-voltage retention (native failed extraction masks those voltages).
