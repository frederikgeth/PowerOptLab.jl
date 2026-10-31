# Advanced IBR branch: final release review

Reviewed 31 October 2026. The four release gates are complete. No unresolved
blocking correctness, API, numerical, documentation or reproducibility finding
remains within the declared model scope. No PR has been opened or published.

The tested implementation is `35d5cac`; subsequent commits package its evidence
and documentation. `src`, package tests and `Project.toml` remain byte-identical
to that tested revision. Hardware sizing and validation tolerances were not
changed during release preparation.

## Validation

| Gate | Final evidence | Result |
|---|---|---|
| Supported Julia package suite | Julia 1.10.12 and 1.13.3, independently resolved compatible environments | 22,372/22,372 package assertions pass on each, plus package-quality checks |
| Generation-runtime regressions | Julia 1.12.6 PCC review and terminal-profile fleet cases | All pass |
| Documentation | Final `docs/make.jl` on Julia 1.10.12 and 1.13.3 | Both builds pass; existing size warnings only; deployment skipped |
| Executable tutorials | Eight study tutorials on both supported versions, including mixed local PF/coordinated OPF; eight DOE tutorials on latest Julia | All pass |
| Tutorial figures | Three generated study figures checked against fresh solved data on both supported versions | Pass |
| Mixed-controller campaign | Four conditions × local/coordinated modes, 96 inverter rows | Eight accepted points; 375 record, 3,960 SI inverter and 6,608 independent SI feeder checks |
| Coordinated OPF regression | Raw solver, final objective, neutral release, coordinate invariance and rejected-record gates | 79 assertions pass on both supported versions |
| Saved-point corruption gates | Reconstruct saved points and fixed-PQ network replay; reject altered currents, source conditions, hardware, raw outcomes, objectives and dependency records | 27 assertions pass on both supported versions |
| Independent circuit regression | Averaged/switched energy, three-wire common mode, peaks, voltage extrema and witnesses | 53 assertions pass on both supported versions |
| Matched EMT | 32 maximum-modeled-stress representatives; three meshes, averaged and switched | 192 integrations, 391 parent checks and 7,276 circuit checks pass |
| Thermal enclosure | Current all-duty AC-induced DC-current term plus passive-damping residual | All 32 predictions enclose switched measurements; largest selected RMS 0.880210 A against 1 A; minimum prediction margin 0.119437 A |
| Boundary experiment | Three references and 18 integrations | 749 checks pass; LF-only switched voltage failure remains rejected |
| Declared carrier domain | 972 predeclared P/Q, sequence, modulation and angle points | 36,033 independent checks pass; all 855 modeled-constraint rejections retained |
| Python integrity/publication gates | Archive, domain, EMT selection and publication checks | 35 tests pass |
| Compact artifacts | Independent audit replay and exact semantic comparison in `summarize_evidence.py --check` | Pass; original receipts remain unchanged |
| Clean exported checkout | No `.git` and no evidence cache; restore actual distributed bundles | All 182 current and 546 historical files verify; 27 saved-record and 35 Python gates plus compact-artifact check pass |
| Publication | Regenerated figures, report/receipt hashes and visual inspection | Pass |

The full-suite logs are `/tmp/ibr-release-julia110-v2-full.log` and
`/tmp/ibr-release-julia113-v2-full.log`. The repository CI workflows encode the
repeatable commands. GitHub-hosted CI has not been run for a PR yet; the local
checks above exercise its release paths.

## Review findings resolved

- Support-floor numerical failures: implicit controllers preserve explicit
  solver options and retain every raw stage. A bounded cold fallback can try
  alternate default scaling; a zero-objective fallback transfers only an
  accepted loss-objective seed into the final zero objective. Failed raw values
  are not promoted. Zero-iteration and time-limit outcomes remain final.
- Stale evidence: the authoritative PF/OPF, carrier-domain, boundary and matched
  EMT campaigns were regenerated from committed implementation `35d5cac`.
  All four source manifests match current files (67, 68, 64 and 66 entries).
  The generation dependency lock is preserved; runtime audit environments are
  recorded separately.
- Incomplete archival packaging: `current.json` and immutable
  `current-35d5cac.json` identify the complete current bundle. Distributed
  compressed bytes allow restoration after squash merges or branch deletion.
  The actual current and historical bundles were tested without Git history or
  an evidence cache. Every raw payload removed from the review tree matched one
  of the immutable manifests before compaction.
- Inconsistent scientific wording: historical frozen-current thermal results
  are labeled explicitly. Current coupled thermal bounds are distinguished from
  the 10% frozen-current allowance. Rejection by a sufficient engineering bound
  is not described as proof of physical hardware infeasibility.

## Compatibility and interpretation

The coupled embedded PWM thermal allocator requires a centered monolithic
three/four-leg L-filter bridge, proportional passive R/L and an open
carrier-frequency source. LCL and finite-source PWM thermal allocation remain
outside that bound; other documented fundamental/harmonic paths retain their
scope. A missing phase inductance is rejected. The retained four-leg fixture
with zero neutral inductance cannot fit its conservative all-duty reserve into
1 A and remains unresolved; its hardware was not resized.

The implementation and case-study review checked current/power reference
locations, sequence/neutral physics, complete final OPF capability, source and
loss balances, thermal/DC limits, strict raw/physical/replay publication gates,
source freshness, archive integrity and API support documentation.

The coordinated benchmark releases local VV/VW reactive-service obligations.
It is a feasible local nonconvex OPF comparison, not an equal-service dominance
or global-optimality claim. Feeder limits are slack in these cases. Analytical
bounds use the periodic frozen-slow-quantity carrier approximation; the slow
energy model retains fitted mean losses. Fixed-POC EMT includes finite-window
DC drift and does not validate feeder/controller dynamics, periodic steady
state, regional conformity or all-time transient bounds.

## Review package

See `studies/australian_pv/evidence/README.md`, `current_summary.json` and the
mixed PF/OPF report. The review tree retains 49 compact result files (about
0.2 MB), two physical reference fixtures, summaries and teaching figures.
The complete current and historical result bundles contain 182 and 546 files,
compressed to 6,824,341 and 11,967,257 bytes. Historical source hashes, failed
predictions and failed hardware/solver outcomes have not been relabeled.

The next step is to open the PR after reviewing this report. No additional
modeling scope or inverter-sizing work is needed to close this package.
