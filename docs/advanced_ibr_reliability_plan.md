# Advanced IBR reliability implementation

Branch: `codex/advanced-ibr-reliability`. Baseline: `556d547`.

All five findings from the scientific review remain present at the baseline.
The existing 513 advanced-inverter assertions pass, but do not cover the review
counterexamples. This work keeps the fundamental-frequency model and its
declared small-ripple/frozen-carrier approximations.

1. **Publication contract:** distinguish inner NLP termination from PWM audit
   validity and reserve convergence. Failed audits and exhausted outer iterations
   must be non-publishable, retain diagnostic reasons, and mask operating results.
   Test both an insufficient iteration budget and an unrealizable SPWM reference.
2. **Rating normalization:** normalize apparent-power, conductor, sequence,
   AC voltage, single-phase power-ripple, and capacitor-current limits by their own
   physical ratings. Keep the 2ω DC-voltage row in SI for controller conditioning. Exercise
   multiple power bases and independently check extracted physical residuals.
3. **DC spectral tail:** bound the omitted sampled-waveform energy with the
   maximum source/capacitor transfer gain over omitted integer carrier harmonics.
   Include resonances beyond the cutoff; retain the open/resistive-source limits.
   Compare low-order bounds against refined spectra and test singular tails.
4. **Split-link physics:** include the first-order 2ω midpoint displacement
   `(Cu-Cl)/(2(Cu+Cl))*d(t)` in the hull and both PWM audits. Verify against
   independent capacitor charge balance, including a limiting rail example.
5. **Three-leg projection:** support the declared conductor dimension in
   validation and the AC audit. Preserve grounded-midpoint LCL paths rather than
   silently deleting a neutral connection. Test default-neutral and floating
   configurations, with scalar and coupled filter primitives.
6. **Documentation and integration:** update the equations, publication checks,
   spectral-tail guarantees, and numerical guidance. Run advanced-inverter,
   controller, and device integration regressions and build the documentation
   where the installed environment permits it. Record results below.

Acceptance requires the original suite plus independent counterexamples to pass;
NLP success alone is not evidence of PWM convergence or physical feasibility.
Finite AC harmonic resolution, network solver tolerances, local optimality, and
small-ripple validity remain explicit limitations.

## Execution record

- Created the implementation branch and rechecked all five baseline defects.
- Implemented all five fixes. `InverterResult.inner_solve` retains the NLP
  status; `pwm_status` and `solve_status(result)` describe full PWM validity.
  Failed closure masks operating quantities while retaining carrier diagnostics.
- Added regressions for iteration exhaustion, SPWM overmodulation, an omitted
  singular DC harmonic, and a between-sample peak that passes both coarse grids.
  Added physical rating base sweeps, spectral-tail bounds, unequal-bank charge
  checks, and three-leg scalar/coupled conductor cases.
- Kept the 2ω voltage-ripple row in its original SI scaling. Normalizing that
  already base-independent row caused a tightly saturated controller case to
  reach 500 iterations; restoring SI scaling converged in 16 iterations with
  an independently measured residual below `1e-12`. No physical test tolerance
  or solver acceptance criterion was relaxed.
- Replaced an SI/per-unit comparison of an unpriced balancing-charge decision
  with checks of its charge-balance identity and physical bounds. That decision
  is not uniquely determined when the rail constraints have headroom.
- Validation on Julia 1.12.6:
  - advanced-inverter suite plus reliability regressions: **849 passed**;
  - controller, controller numerics, fleet/study/sizing, and multiperiod device
    integration: **1,004 passed**, one optional OpenDSS test skipped because
    OpenDSSDirect is not installed in this environment;
  - all Julia code blocks in the three focused IBR tutorials executed successfully;
  - Documenter build passed (deployment omitted), with only HTML/search-index
    size warnings; `git diff --check` passed.
- Tests used `--compiled-modules=existing` to avoid writing to the read-only
  user cache. Documentation used a temporary writable depot for package logs.
- Remaining model limits: the DC tail bound applies to the sampled frozen-carrier
  spectrum, not unmodelled sidebands; AC spectra still require refinement;
  the dense rail audit is not a continuous-time certificate; and a rejected
  modulation reference is reported as an audit failure rather than repaired
  by a strategy-constrained reoptimization.

## Dual controller encodings (follow-up implementation)

Implemented on the same branch:

1. `SequenceController(...; encoding=:smooth)` preserves the default graph;
   `encoding=ComplementarityGraph()` selects exact flat-tail curves, extrema,
   clipping, norms, P/Q/current limits and plant capability allocation. Both
   paths share the physical controller graph and retain physical voltage floors,
   priority headroom and firmware conflict blending.
2. Added shared complementarity selectors with simplex weights. Exact norms
   materialize their components before squaring to avoid cancellation around
   balanced negative-sequence voltage. Zero-state degeneracies are documented.
3. Added `stamp_control!` dispatch and a pre-attachment `configure!` callback for
   optional solver bridges. Added finite complementarity and exact-current-replay
   publication checks; the inner solve cannot override a failed encoding audit.
4. Added independent exact-graph witnesses across sensing policies, flat tails,
   unbalance and model bases; P/Q priority and zero-headroom allocator checks; and
   a deliberately misreported-success regression for the publication gate.
5. Added a matched optional-backend comparison and CI artifacts. Per the user's
   request, CCOpt v0.1.0 failures are characterized for the PR rather than hidden,
   promoted, or used to reject the requested complementarity functionality.

Validation:

- Controller/controller-study regressions: **8,088 passed**, one optional
  OpenDSS check skipped because OpenDSSDirect is unavailable.
- Functional formulation regressions: **3,382 passed**.
- Public snapshot API checks: **10 passed**; both successful CCOpt cases pass
  the full publication audit and return finite physical operating quantities.
- Matched backend comparison: **20 runs / 23 assertions passed**. Smooth Ipopt
  strictly solved 5/5; smooth MadNLP 3/5; default CCOpt 0/5; the explicitly
  tighter CCOpt relaxation 2/5. These counts characterize this environment and
  are not cross-platform convergence guarantees.
- Documenter built successfully with only the existing size warnings; no
  deployment. `git diff --check` passed.

The [PR-ready failure list](diagnostics/ccopt-inverter-controller-failures.md)
separates CCOpt-only failures, shared MadNLP difficulties, and resolved cases.
Its [compact result bundle](diagnostics/ccopt-inverter-controller-results.json)
contains raw statuses, candidate residuals, versions and source hashes. The
comparison script and CI retain the full TOML evidence. The solver failure
report is prepared for the PR discussion.

### Zero bound-relaxation follow-up

At the user's request, extended the same comparison with three paired methods
setting only `bound_relax_factor=0.0`: smooth MadNLP, default-relaxation CCOpt,
and tighter-relaxation CCOpt. The original 20 runs reproduced their outcomes;
the complete **35-run comparison passed 38 record/baseline assertions**.
Strict success counts changed from 3/5 to 2/5 for smooth MadNLP, 0/5 to 1/5
for default-relaxation CCOpt, and 2/5 to 0/5 for tighter-relaxation CCOpt.
The report and compact evidence bundle now include every pair and its residuals.
Production solver defaults and publication criteria were not changed.
