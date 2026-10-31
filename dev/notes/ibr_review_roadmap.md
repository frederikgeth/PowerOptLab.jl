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

- [x] **R1 Fixture capacitor rating.** Derive `i_cap_max` from a declared
  design rule (nameplate deliverable at the AS/NZS 4777.2 volt-var corners,
  checked against Kolar & Round 2006), not from the old 2ω-only model. Keep
  the 1 A value only as a labelled sensitivity.
- [x] **R2 Rated-capability fixture gate.** A study test that every fixture
  plant delivers rated apparent power at unity PF and the Australia A volt-var
  extremes (−60 % / +44 % of rated S) across 216–253 V, with phase current,
  not an internal budget, as the limiting quantity. Run it before any campaign.
- [x] **R3 Rerun the matched PF/OPF benchmark** with the corrected fixture and
  regenerate tables/figures.
- [x] **R4 Equal-service OPF.** Add a coordinated mode that enforces each
  site's local volt-var request (same reactive service) and maximises export,
  so "coordination benefit" and "cost of releasing volt-var" are separated.
- [x] **R5 PF/OPF agreement test.** OPF with the local laws imposed must
  reproduce the local fixed point (voltages, P, Q) within tolerance.
- [x] **R6 Feeder-limited case.** At least one predeclared case where the
  253 V planning limit or a native thermal limit binds, so Q has network value.
- [x] **R7 Report claims generated from data.** Binding-constraint counts and
  limiting-quantity statements in the report come from the CSVs, with a check
  that fails on text/data drift. Fix the "DC source ceiling" contradiction.
- [x] **R8 Wording.** "Recovery" → "export gained by releasing volt-var";
  "independent EMT" → "separately coded"; drop assertion-count headlines in
  favour of what was tested and the worst margin; state that the thermal
  enclosure is not an out-of-sample test; short glossary for project jargon.
- [x] **R9 Remove temporary/bloat files** (list below) and consolidate the
  study pages into one current page plus one archive page.

## P1: should follow in this PR if time allows

- [~] **R10 AS/NZS 4777.2 grounding.** Report the 10-minute average
  sustained-overvoltage quantity (258 V) beside the 253 V planning limit;
  label var priority and multi-phase sensing choices as study settings that
  must be checked against the standard's clauses; add a watt-priority
  sensitivity.
  *Done:* 258 V sustained quantity reported per case; var priority and
  sensing labelled as study settings in the README. *Open:* watt-priority
  sensitivity needs a watt-priority stage in `allocate_inverter_services`
  (the voltage-aware allocator currently requires var priority).
- [~] **R11 Scalability record.** Save solve times, iterations and NLP size
  for every PF/OPF solve; add a site-count sweep on the LV feeder.
  *Done:* elapsed time, local iterations, equal-service rounds and NLP size
  are recorded in `cases.csv` and summarised in the report. *Open:* site-count
  sweep (needs a second feeder or synthetic replication; see R15).
- [x] **R12 Superseded sensitivity results.** Label or archive the
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
  27 MB binary after dropping the superseded first-benchmark bundle on 2026-11-01). Move to a GitHub Release or Zenodo and keep only manifests
  with checksums and a download step. Creating a release is an external action.
- [x] **D2 Split PR #75** — decided 2026-11-01: keep as one PR.
- [U] **D3 LCL filters.** The embedded PWM thermal bound supports only L
  filters; LCL is the common commercial topology.

## Temporary/bloat files (R9) — done

Deleted: 12 temporary dev notes (agent goal/acceptance/review logs, old plans);
29 superseded study scripts; 10 superseded study pages and 14 figures; all
tracked historical result folders (`experimental_*`, `sensitivity_*`,
`coordinated`, `four_leg`, `thermal_energy`, `voltage_aware_*`, ...). All
historical records remain in `evidence/historical-d2c9a97.json`; pages and
scripts remain in Git history. `HISTORY.md` summarises them. Kept: the DOE and
operability roadmaps (forward-looking, not from this PR).

## New items found during the work

- [ ] **R18 Ripple-limit margin.** Switched EMT exceeds the 7.5 V 2ω DC-link
  ripple limit by 4 mV at one released-OPF four-leg site placed exactly on
  that limit (averaged model 7.500 V). Options: a small OPF margin on active
  ripple limits, or a switching correction in the 2ω model. Reported, not fixed.
- [ ] **R19 Evidence brittleness.** Any edit to a hashed shared script makes
  all saved records stale (full rerun ~15 min). Consider hashing only the
  functions a campaign uses, or separating study helpers from campaign runners.
- [ ] **R20 Equal-service on weak hardware.** Same-service OPF is infeasible
  when hardware cannot deliver its own volt-var request (legacy fixture). A
  "same service fraction" variant (match the local λ) would make the
  comparison defined there too.

## Log

- 2026-11-01: roadmap created from the PR #75 review.
- 2026-11-01: R1–R9, R12 done. Declared fixture 7 A capacitor / 23 A peak
  (model requirement at the Australia A corners + 10 %); legacy 1 A / 14 A
  delivers 0.84–1.25 kVA of 6 kVA and fails the gate. Three modes (local,
  same-service OPF, released OPF) over five conditions, including the
  feeder-limited `oversized_dc`. PF/OPF agreement holds to ~1e-8 V in every
  condition. Releasing volt-var gains only losses unless 253 V binds.
  Superseded scripts/pages/results removed (archived in the historical
  bundle); temporary dev notes deleted.
- 2026-11-01: campaign regenerated from clean commits and sealed as
  `evidence/current-05db2db.json`. 60 EMT references / 360 integrations: every
  thermal bound encloses the switched value (max 4.65 A vs 7 A); one 2ω ripple
  finding retained (R18). Carrier-domain (972 points) and DC-boundary
  validations regenerated under the legacy limits. All study CI checks pass
  locally; package `src/` unchanged; GitHub CI not yet run.
