# Inverter-control implementation workflow

> Historical developer workflow record, excluded from the documentation build.
> Extracted from [control_study_methodology.md](../../docs/src/ibr/control_study_methodology.md).
> Check the published reference guides for current capabilities and evidence.

## Implementation sequence

Steps 1--5 are implemented by the public BMOPFTools smooth-PWL API and the
PowerOptLab phase-aware control component. The exact and smooth evaluators,
three voltage-reference policies, converter/grid current targets, dominant
conflict policy, converter-terminal power extraction, three-leg plant
composition, supported-unit regression, and extracted exact-versus-smooth current
residual are covered by unit tests. Fleet construction, matched-case execution
and initial hardware counterfactuals are also implemented, as recorded in steps
6--9 below. Dataset-specific campaigns and the independent exact-law network
oracle remain separate work; batteries are outside this programme.

1. **Complete:** submit the BMOPFTools public smooth-PWL API PR.
2. **Complete:** add PowerOptLab exact sequence/curve/limiter numeric helpers and tests.
3. **Complete:** add `SequenceController`, `ControlledDevice{AdvancedInverter}`, and semantic
   control handles.
4. **Complete:** stamp the smoothed worst-phase and negative-sequence admittance laws into one
   advanced inverter; verify exact/smooth agreement.
5. **Complete:** add continuous normalised conflict handling, rated/available
   Volt-watt bases, P/Q priority, plant-side converter/grid current and
   apparent-power backoff, optional `dv2_max` backoff, and binding-limit tests.
6. **Complete:** add multi-device single-snapshot construction using
   per-component `OpfDeviceBuilder` replacement for dataset IBRs.
7. **Complete:** add deterministic tidy device/phase extraction, ownership
   provenance, physical-balance tests, and supported-unit checks. Large-network
   scaling benchmarks remain part of the case-study campaign rather than the
   unit suite.
8. **Infrastructure complete:** add deterministic matched-case execution,
   failure retention, duration/weight-aware summaries, and paired variant
   differences. The dataset-specific fixed-hardware penetration campaign is the
   next research run.
9. **Initial infrastructure complete:** expand matched cases over explicit
   converter-current, grid-current, DC-capacitance, and capacitor-current
   counterfactuals; report achieved-current and closed-form 2ω requirements.
   Dataset-specific service predicates and complete capability frontiers are
   the next research run.
10. **Open:** implement the exact-law fixed-point oracle described in
    [Verification and benchmark cases](../../docs/src/ibr/verification.md) and report
    ``\|\Delta V\|_\infty`` between the firmware-law and smooth-law equilibria on
    at least one unbalanced feeder per control law. Until this exists, the
    reported exact/smooth residual bounds controller-algebra smoothing only, and
    equilibrium agreement is an assumption of every comparison above.
11. **Open:** audit the starts and scales of every lifted magnitude, decide
    whether the controller selectors become rating-relative, and re-measure the
    `per_unit=true` restriction against the shifted-expression representation.
    The evidence and the three specific defects are in the numerical-policy
    section of [Phase-aware local control laws](../../docs/src/ibr/phase_aware_control_laws.md). This gates any study on devices whose ratings
    differ by more than about a decade, and any MV study.
12. **Open:** read the licensed AS/NZS 4777.2 and IEEE 1547-2018 clauses that
    specify the monitored voltage for a multi-phase DER, and settle the novelty
    framing of the worst-phase envelope before writing it up.

## Diagram development notes

Extracted from [the IBR overview](../../docs/src/ibr/index.md).

The next most useful figures are a constraint-dependency map from parameters to
`InverterResult` diagnostics, and a validation pyramid linking exact identities,
published PLECS cases, EMT sweeps, and hardware measurements. A later
switching-state figure should focus on dead time, discontinuous/interleaved PWM,
and common-mode earth paths only when those mechanisms are present in the equations.
