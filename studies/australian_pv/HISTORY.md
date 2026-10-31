# Study history

Earlier campaigns on the same feeder, in order. They are superseded by the
[current benchmark](engineering_benchmark.md). Their records are kept unchanged
in the historical evidence bundle (`evidence/historical-d2c9a97.json`, 546 files)
and their pages and scripts are in Git history. Restore them into a separate
folder; they must not overwrite current results:

```sh
python scripts/australian_pv/evidence.py restore studies/australian_pv/evidence/historical-d2c9a97.json --output /tmp/australian-pv-history
```

Reproducing a historical campaign needs the source bytes matching its recorded
hashes (see each campaign's `environment.json`), not only the archive commit.

## Why the earlier export results are not used

From the PWM thermal campaign onward, every campaign used a 1 A RMS DC-link
capacitor budget on a 6 kVA / 10 A / 750 V inverter. That value was chosen when
the model had only a line-frequency (2ω) capacitor term. Adding the carrier
ripple term raised the modelled capacitor current about twentyfold, and the
budget then limited each inverter to about 1.2 kVA. The earlier matched benchmark
reported 4–13 kW of export "recovered" by coordination; that came from releasing
volt-var on hardware that could not deliver both P and Q. With a nameplate-capable
fixture the difference is losses only (see the current report's legacy-fixture
sensitivity table).

## Campaigns

| Campaign | What it did | Status of its conclusions |
|---|---|---|
| Fundamental controller comparison (`results/cases.csv`, `summary.csv`) | Mean, worst-phase, sequence, droop and virtual-impedance laws on the 12-site feeder at 0/50/100 % PV, 1.00/1.03 pu source; single-phase counterfactual | Fundamental-frequency only; no PWM thermal model. Valid as a sensing-law comparison. |
| Sensitivities (`results/sensitivity_*`) | DC-link size, 0.5/3 A capacitor budget, smoothing, balanced load/source, pedagogical VW curve | Computed without the carrier thermal term. The capacitor-budget sensitivities do not describe the current model. |
| Terminal-profile and PCC fleet experiments (`results/experimental_*`, `terminal_profile_validation`) | Solver scaling and failure diagnostics | Numerical diagnostics only. |
| Voltage-aware allocation (`voltage_aware_*`) | Staged Q→P→I₂ local allocation; unbalance angle and load sweeps; synthetic 24-period day | Superseded allocator settings; no carrier thermal term. |
| Coordinated capability (`coordinated/`) | First free-Q OPF with balanced and free-I₂ current | Superseded by the matched benchmark. |
| Four-leg counterfactual (`four_leg/`) | Per-phase VV/VW with a neutral leg | Research counterfactual; STNW1170 §4.3.3 requires balanced output. Used the 1 A budget. |
| First switched-circuit EMT (`emt/`) | Three fixed-POC references | Retained capacitor thermal failures under the 1 A budget. |
| PWM thermal and nonlinear DC energy (`thermal_energy/`, `thermal_energy_emt/`) | Embedded PWM thermal bound with fixed sizing; four night prediction misses | The misses led to the all-duty thermal term. 1 A budget. |
| Matched PF/OPF benchmark, first version (`engineering_benchmark` at `35d5cac`) | Local versus free-Q OPF, 4 conditions, 32 EMT references | Superseded: 1 A budget made every site capacitor-limited; no same-service OPF. |

The carrier-domain and DC-boundary validations ([domain](engineering_domain.md),
[boundary](engineering_boundary.md)) check the carrier approximations, not the
fixture, and were regenerated with the current source under the legacy limits.
