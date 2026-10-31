# Advanced IBR / Australian PV study: review roadmap

Offline tracker for the PR #75 review (1 November 2026). Status keys:
`[ ]` open, `[~]` in progress, `[x]` done, `[U]` needs a user decision.
Update this file in the same commit as the work it tracks.

## Why

The PR's code and verification infrastructure are careful, but the headline
Australian result was set by one unsourced fixture value: a 1 A RMS DC-link
capacitor budget on a 6 kVA, 10 A, 750 V inverter. Every daytime inverter row
was capacitor-limited, PV size had no effect, and the "coordinated" OPF simply
switched volt-var off (Q ≈ 0) on a feeder with slack voltage and thermal limits.
The local fleet could not deliver its own AS/NZS 4777.2 volt-var requests, so it
was not representative of certifiable hardware.

## P0: required before merge

- [ ] **R1 Fixture capacitor rating.** Derive `i_cap_max` from a declared
  design rule (nameplate deliverable at the AS/NZS 4777.2 volt-var corners,
  checked against Kolar & Round 2006), not from the old 2ω-only model. Keep
  the 1 A value only as a labelled sensitivity.
- [ ] **R2 Rated-capability fixture gate.** A study test that every fixture
  plant delivers rated apparent power at unity PF and the Australia A volt-var
  extremes (−60 % / +44 % of rated S) across 216–253 V, with phase current,
  not an internal budget, as the limiting quantity. Run it before any campaign.
- [ ] **R3 Rerun the matched PF/OPF benchmark** with the corrected fixture and
  regenerate tables/figures.
- [ ] **R4 Equal-service OPF.** Add a coordinated mode that enforces each
  site's local volt-var request (same reactive service) and maximises export,
  so "coordination benefit" and "cost of releasing volt-var" are separated.
- [ ] **R5 PF/OPF agreement test.** OPF with the local laws imposed must
  reproduce the local fixed point (voltages, P, Q) within tolerance.
- [ ] **R6 Feeder-limited case.** At least one predeclared case where the
  253 V planning limit or a native thermal limit binds, so Q has network value.
- [ ] **R7 Report claims generated from data.** Binding-constraint counts and
  limiting-quantity statements in the report come from the CSVs, with a check
  that fails on text/data drift. Fix the "DC source ceiling" contradiction.
- [ ] **R8 Wording.** "Recovery" → "export gained by releasing volt-var";
  "independent EMT" → "separately coded"; drop assertion-count headlines in
  favour of what was tested and the worst margin; state that the thermal
  enclosure is not an out-of-sample test; short glossary for project jargon.
- [ ] **R9 Remove temporary/bloat files** (list below) and consolidate the
  study pages into one current page plus one archive page.

## P1: should follow in this PR if time allows

- [ ] **R10 AS/NZS 4777.2 grounding.** Report the 10-minute average
  sustained-overvoltage quantity (258 V) beside the 253 V planning limit;
  label var priority and multi-phase sensing choices as study settings that
  must be checked against the standard's clauses; add a watt-priority
  sensitivity.
- [ ] **R11 Scalability record.** Save solve times, iterations and NLP size
  for every PF/OPF solve; add a site-count sweep on the LV feeder.
- [ ] **R12 Superseded sensitivity results.** Label or archive the
  `sensitivity_*` folders computed without the carrier thermal term.

## P2: follow-up work (separate PRs)

- [ ] **R13 Held-out thermal-bound validation.** Predeclared random
  operating points and different filter parameters, run through the switched
  EMT; the current 32/32 enclosure was checked after the bound was revised.
- [ ] **R14 External cross-checks.** Fundamental PF against OpenDSS on the
  same feeder with fixed injections; one or two EMT points against PSCAD/PLECS.
- [ ] **R15 Scalable OPF.** Offline per-device capability region in
  (P, Q, V, unbalance) from the rigorous carrier bound, a convex inner
  approximation inside the network OPF, then detailed post-check/EMT.
  Scaling to MV21 with LV feeders.
- [ ] **R16 Australian realism.** CSIRO LV feeder taxonomy feeders, a
  day-long time series with public profiles (e.g. Ausgrid solar home data),
  DOE / CSIP-AUS flexible exports as the coordinated counterpart, and
  optimised versus default Australia A/B/C curves (`bilevel.jl`).
- [ ] **R17 Fixed-point diagnostics.** Report local-mode iteration history
  and a uniqueness/contraction indicator for the volt-var equilibrium.

## User decisions

- [U] **D1 Evidence bundles** (`studies/australian_pv/evidence/bundles`,
  18 MB binary). Move to a GitHub Release or Zenodo and keep only manifests
  with checksums and a download step. Creating a release is an external action.
- [U] **D2 Split PR #75** into core models + tests, and the case study.
- [U] **D3 LCL filters.** The embedded PWM thermal bound supports only L
  filters; LCL is the common commercial topology.

## Temporary/bloat files (R9)

| Item | Action |
|---|---|
| `dev/notes/switching_engineering_opf_acceptance.md` | delete (stale release review) |
| `dev/notes/switching_engineering_opf_goal.md` | delete (agent goal record) |
| `dev/notes/ibr_review_resolution.md` | delete; this roadmap replaces it |
| `dev/notes/inverter_pr_review.md`, `inverter_full_opf_validation.md`, `ibr_formulation_design.md`, `equilibrium_workbook_design.md`, `inverter_control_implementation.md` | delete (pre-implementation plans / old reviews) |
| `dev/notes/documentation_branch_review.md`, `documentation_tutorial_audit.md`, `doe_documentation_workflow.md`, `dsse_review_resolution.md` | delete (historical records, not from this PR) |
| `evidence/current-35d5cac.json` | byte-identical to `current.json`; check references first |
| `results/experimental_*` | delete (debugging artefacts) |
| `studies/australian_pv/*_results.md` twins and superseded study pages | merge into one current page + archive |

## Log

- 2026-11-01: roadmap created from the PR #75 review.
