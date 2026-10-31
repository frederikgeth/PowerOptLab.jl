# Voltage-aware allocation experiments

These experiments couple the local Q→P→I₂ physical allocator to the Australian
LV3 feeder. They use the existing 6 kVA / 10 A / 750 V / 300 µF plant unchanged,
the sourced Australia A VV/VW knots, converter-arm current limiting, POC apparent
power limiting, neutral-referenced sensing and reactive priority. They are a
distinct allocation policy from the earlier fixed-structure scaling controller;
results must retain their campaign labels.

The three compared sensing laws are mean, worst-phase and positive sequence
VV+VW. Their hardware, limits, settings, available PV and numerical budgets are
matched. No negative-sequence compensation is requested in this comparison.
Within-request Q service is maximized first, achieved Q is preserved while P is
maximized, then any requested I₂ service would be allocated without changing
those powers. An unsuccessful local NLP remains unresolved, not a proof of
physical or controller-policy infeasibility.

## Reproduction and acceptance

From the repository root, with the project environment installed:

```sh
julia --project=. scripts/australian_pv/test_voltage_aware.jl
julia --project=. scripts/australian_pv/voltage_aware.jl snapshots
julia --project=. scripts/australian_pv/voltage_aware.jl chronology
python3 scripts/australian_pv/test_voltage_aware_records.py
python3 scripts/australian_pv/plot_voltage_aware.py
```

An optional second argument selects an empty output directory. Existing nonempty
outputs are rejected to prevent mixing campaigns. Each case gets one 100 kVA
base and a common 1,000-iteration Ipopt profile per NLP, with at most 30 feeder
iterations. There are no per-case retries. Chronological warm starts are numerical
current seeds from the last accepted interval, not physical state transitions.

The plotting command needs Matplotlib and NumPy. It audits both saved campaigns
and their source hashes before generating the [result tables](voltage_aware_results.md)
and standalone PNG/SVG figures. Failed or incomplete campaigns retain their audit
report but withhold the matched figures. Four saved-record regressions cover
complete, missing, nonfinite failed and falsely successful physical records.

Every candidate is audited at its actual feeder voltage. Publication requires
local allocation success, continuous full-cycle rails, independently recomputed
inverter/source/policy limits, converged voltage/current maps, native SI equality
residual checks and a fixed-PQ network replay within 1 mV. That replay changes the
injection formulation in the same network engine; it is not an external or
independent circuit implementation. Native inequality constraints retain their
native enforcement; the SI residual audit does not independently certify them.

`environment.json` records the predeclared cases, complete executed source-file
hashes, revision, environment and fixed hardware/settings. `cases.csv` retains
every completed case, including failures; `devices.csv`, `physical_checks.csv`,
`iterations.csv` and per-case bus phasors supply detailed evidence.
`summary.json` distinguishes completed from planned cases. Any failed or missing
chronological interval withholds the corresponding energy totals.

## Controlled scenarios

The 41 snapshots comprise the 18 matched source/PV/controller combinations,
six source-VUF magnitude cases (0% and 2%), nine U₂/U₁ relative-angle cases
(90°, 180°, 270° at 1% VUF), six load-phase-allocation cases (balanced and a
cyclic redistribution of the original 60/25/15% shares), and two disconnected
baselines. U₁ is held fixed while U₂ changes, and no upstream U₀ is introduced.
The transformer and full four-wire feeder determine the customer phasors; an
upstream sequence perturbation is not a prescribed customer voltage.

Disconnection means an open grid-side contactor: the entire converter/filter/
POC-shunt AC model is absent. Energized zero-PV operation instead retains Q
control and modeled loss imports. Disconnected DC standby is outside the AC
model and is recorded as undefined, rather than assigned zero consumption.
The disconnected baselines are pure power flows, not different controllers.

The synthetic day has six piecewise-constant intervals:

| Start (h) | Duration (h) | Available PV fraction | Load multiplier |
| --- | --- | --- | --- |
| 0 | 6 | 0 | 0.6 |
| 6 | 3 | 0.15 | 0.8 |
| 9 | 3 | 0.65 | 0.9 |
| 12 | 3 | 1 | 1 |
| 15 | 3 | 0.45 | 1.1 |
| 18 | 6 | 0 | 1.2 |

All intervals use a 1.03 pu positive-sequence source with 1% aligned upstream
VUF. Each controller is compared with its inverters energized throughout the
day and with their grid-side contactors open during the zero-PV intervals.
This is an illustrative quasi-static profile, not measured Australian yield,
weather, switching behavior or transient response. Available fleet PV energy is
405 kWh. Export and import are integrated per inverter before aggregation;
signed reactive exchange and converter/filter loss energy are also retained.
The modeled identity is available PV minus net POC energy equals unused DC
availability plus active converter/filter losses. Excluded DC standby must not
be inferred from that identity.

Iteration histories describe numerical convergence, with voltage residuals in
volts and current-map residuals in amperes. Intermediate virtual feeder probes
may violate inverter constraints and are never treated as transient dispatch
or included in energy totals. These experiments do not establish dynamic
stability, global uniqueness, regional connection conformity or EMT validity.
