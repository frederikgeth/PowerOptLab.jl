# Results

Tracked here (compact, reviewable):

- `engineering_benchmark/` — current matched PF/OPF benchmark: `cases.csv`,
  `devices.csv`, `agreement.json`, `rated_capability.csv`, `environment.json`,
  `summary.json`. Report: [engineering_benchmark.md](../engineering_benchmark.md).
- `engineering_benchmark_legacy_fixture/cases.csv` — the same benchmark with the
  legacy 1 A / 14 A fixture, kept as a sensitivity. That fixture fails the
  rated-capability gate.
- `engineering_domain_summary.csv`, `engineering_boundary_emt/comparison.csv` —
  carrier-model validation summaries.
- `historical_night_thermal_misses.csv` — the four night-time thermal-prediction
  misses of the first PWM thermal bound (negative result; see
  [study history](../HISTORY.md)).
- `validation.json`, `modified_validation.json` — network conversion checks for
  the CSIRO feeder import.

Full raw records (saved phasors, EMT waveforms, solver stages) are in
checksum-verified bundles. Restore them from the repository root before running
saved-record tests or the report generator:

```sh
python scripts/australian_pv/evidence.py restore studies/australian_pv/evidence/current.json
```

Earlier campaigns: [study history](../HISTORY.md).
