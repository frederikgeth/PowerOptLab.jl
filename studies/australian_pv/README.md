# Australian LV photovoltaic control study

This is the Australian successor to the ENWL showcase. It starts with a small
real feeder, preserves a converted reference, and makes the teaching experiment's
changes explicit. The first campaign compares **12 three-phase PV systems at
12 three-phase customer connections, with unbalanced loading**. All variants use
the same network, loads, hardware and available solar power.

## Reproduce

From the PowerOptLab repository root, with its pinned Julia environment installed:

```sh
julia --project=. scripts/australian_pv/prepare.jl
julia --project=. scripts/australian_pv/run.jl
julia --project=. scripts/australian_pv/audit.jl
python3 scripts/australian_pv/plot.py
```

`LV3_three_phase_pv.bmopf.json` contains native ownership placeholders for the
advanced PVs; use `run.jl` for their controller requests and ratings. It is not a
standalone fixed-PQ dispatch case. `LV3_three_phase_no_pv.bmopf.json` is the
modified-load power-flow baseline.

The plot script needs Matplotlib and NumPy (`python3 -m pip install matplotlib`).
PNG and SVG figures are generated locally and ignored by Git. Run the plotting
command above to create the figures and customer key linked below. Julia
`--compiled-modules=existing` is useful in environments with a read-only depot.
`run.jl` also exposes `run(modes=..., levels=..., sources=...)` for smaller
experiments; the publication plot script intentionally requires the complete
36-case default campaign. Outputs are overwritten on rerun.

- [Feeder schematic](figures/feeder.png): customer locations along each route;
  [customer key](figures/customer_key.csv) maps presentation labels to source IDs.
- [Phase voltage comparison](figures/phase_voltages.png): identical axes for all
  controllers, with the upper supply limit visible.
- [Benefits and effort](figures/tradeoffs.png): voltage, VUF, neutral displacement,
  active power reduction, reactive absorption and converter current together.
- [Numerical results](results/summary.csv), [case acceptance](results/cases.csv),
  [baseline validation](results/validation.json), and
  [modified-network validation](results/modified_validation.json).

## First results

At full PV, source positive-sequence voltage 1.03 pu and source VUF 1%, the
accepted solutions give the following maxima across the 12 customer buses:

| Control | Maximum PN voltage (V) | Maximum VUF (%) | PV delivered (kW) | Peak converter current (A) |
|---|---:|---:|---:|---:|
| Unity PF | 253.547 | 0.985 | 59.997 | 6.623 |
| Mean Volt-var | 253.295 | 0.987 | 59.997 | 7.320 |
| Mean VV+VW | 253.113 | 0.986 | 50.602 | 6.416 |
| Worst-phase VV+VW | 252.940 | 0.986 | 43.113 | 5.884 |
| Sequence VV+VW | 252.968 | 0.986 | 43.015 | 5.682 |
| Sequence + I₂ droop | 252.565 | 0.803 | 44.645 | 10.000 |

This gives a useful engineering example: mean-voltage control reduces voltage,
but the highest phase remains above 253 V in this snapshot. The worst-phase watt
guards bring the customer maxima below that line while delivering less active
power. Adding negative-sequence droop to the otherwise identical sequence policy
reduces maximum VUF from about 0.986% to 0.803% and increases delivered power by
about 1.63 kW, but brings one converter phase essentially to its 10 A limit.
Neutral displacement does not improve in that matched pair (maximum about
0.255 V in both). These are outcomes of this case and tuning, not universal
rankings or a claim of regulatory compliance.

## Source and deliberate changes

The source is CSIRO's *Realistic Australian Medium Voltage Feeder with Associated
Low Voltage Feeders*, Geth et al. (2025), v1,
[DOI 10.25919/ghnz-bk28](https://doi.org/10.25919/ghnz-bk28), copied from BMOPFTools
revision `8ca84ab12c0c91aaa8ad4c9986d6adbeb969ea0b`. The source and derived data are
**CC BY-NC-SA 4.0**; see the copied `DATA_LICENSE.md` files. The Julia/Python
study code follows the repository's code license. `source/provenance.json`
records copied-file SHA-256 hashes and source paths.

`LV3_55bus` has **12 customers**, four originally on each phase. The imported
model has 56 buses, including the upstream source bus. The original line lengths,
impedance matrices, switches, explicit neutral, grounding and transformer are
retained. The transformer is 500 kVA, 11 kV/433 V with a 1.025 primary tap;
433 V is a transformer nameplate value, not the statutory nominal supply voltage.
The feeder is electrically short and relatively stiff: this limits how dramatic
local voltage-control benefits can be.

The original Master uses `BatchEdit load..* kw=1`. The pinned import leaves loads
at the default 10 kW, so preparation explicitly sets each original customer to
1 kW and power factor 0.88 (0.539743 kvar). This corrects the intended baseline
setpoints rather than accepting the import defaults silently. Conversion warnings
are retained: duplicate linecode definitions, ignored switch fields and retained
source-only fields must not be mistaken for numerical validation. BMOPFTools
restores transformer tap and source angle in its parsed representation. This
study validates BMOPFTools power flows, **not independent OpenDSS equivalence**.

The teaching network then makes these changes:

| Quantity | Teaching setting |
|---|---|
| Customer connections | All 12 upgraded to three-phase; existing a/b/c/n bus terminals retained |
| Demand | 1.5 kW/customer, constant PQ, PF 0.95 lagging |
| Demand phase shares | A/B/C = 60/25/15% at each customer; 18 kW total |
| PV | One three-leg inverter per customer; 5 kW available at full irradiance, 6 kVA rating |
| Inverter | 10 A phase current limit, 750 V DC link, 2 mF capacitance, modulation limit 0.96, filter R=0.05 Ω and X=0.15 Ω |
| Source positive sequence | 1.00 or 1.03 times the original 11 kV source |
| Source negative sequence | 1% of positive sequence, aligned with phase-A positive-sequence phasor |
| Source zero sequence | Zero |
| PV availability | 0, 0.5, 1.0 of 5 kW per customer |

These are transparent hypothetical connection/load changes, not measured household
profiles or a proposed Queensland connection approval. Aligning all demand shares
is a deliberately systematic imbalance. The 1% upstream negative-sequence term
makes the sequence-control mechanism visible on a short feeder; it is not an
assertion about the source dataset. At 1.03 positive-sequence voltage, the largest
source phase magnitude is 1.0403 pu. Three-phase upgrades are a modeling choice,
not a claim that all Australian houses already have three-phase service.

`MV21_328bus` is also copied, converted and run through BMOPFTools. A new standalone
Master includes its MV files only. Its load and transformer files are empty:
loads and distribution transformers belong to the associated LV feeders. Thus
this is an **unloaded MV connectivity/energization check**, not a validated loaded
MV operating point. The later integrated study must attach the LV feeders or
explicit, validated aggregate equivalents. It must not treat this unloaded result
as an MV voltage-impact study.

## Matched controller experiment

The six modes distinguish sensing, positive-sequence response and negative-sequence
actuation. Comparing `sequence_vv_vw` with `sequence_droop` changes only the
negative-sequence policy; comparing `mean_vv` with `mean_vv_vw` adds only Volt-watt.

| Mode | Reactive-power voltage input | Active-power voltage input | Negative-sequence current |
|---|---|---|---|
| `unity` | None; Q command zero | Available PV | None |
| `mean_vv` | Mean of three phase magnitudes | Available PV | None |
| `mean_vv_vw` | Mean phase magnitude | Mean phase magnitude | None |
| `worst_vv_vw` | Worst-phase guards (default conflict policy) | Largest phase magnitude | None |
| `sequence_vv_vw` | Positive-sequence magnitude | Largest phase magnitude | None |
| `sequence_droop` | Same as preceding row | Same as preceding row | Voltage-oriented I₂ droop |

Illustrative Volt-var curve: voltage `[207,220,240,258]` V maps to
`[+0.44,0,0,-0.60]` times 6 kVA. Positive Q is injection; negative Q is
absorption. Volt-watt: `[250,260]` V maps to `[1,0]` times **available power**.
These curves are pedagogical settings, **not certified AS/NZS 4777 settings**.
At zero PV availability the inverter remains energized and may supply reactive
and negative-sequence current, drawing active power for losses. This is a
STATCOM-capable assumption, not a model of an inverter that switches off at night.
A supply limit and a controller knee are different quantities. The negative-sequence
admittance curve maps VUF-like ratio `[0,0.0005,0.003,0.01]` to `[0,0,1.5,3]` A/V,
with configured impedance angle `atan(0.3)` rad and ripple blend zero. It is a
fixed illustrative tuning, not feeder-optimal tuning or a stability guarantee.

All controls target grid current. Common capability limiting can reduce the power
command. Three-leg currents sum to zero, so these converters cannot directly
inject zero-sequence current. **They sense phase-to-ground voltage in the current
network API** (`neutral=nothing`). Customer voltage plots use phase-to-neutral
voltage, calculated explicitly as `Vphase − Vneutral`. “Worst phase” here means
worst controller-sensed phase; it does not guarantee the worst customer
phase-to-neutral voltage is held below a limit.

VUF is `100 |V₂|/|V₁|`, with a/b/c sequence order. Subtracting a common neutral
voltage changes V₀, not V₁ or V₂. Neutral displacement is shown separately.
An unbalanced three-wire PV can change phase power while maintaining zero total
phase current. The study must not equate negative-sequence suppression with
neutral compensation or claim either follows automatically from Volt-var.

## Engineering story

1. **Start at the connection.** Show the feeder schematic, original neutral and
   grounding, all three phase loads and the three-wire inverter. Establish which
   voltage a customer experiences and which voltage the controller measures.
2. **Give every controller the same sunlight and network.** Show the three phase
   voltages at full PV and the same upstream condition. A mean or positive-sequence
   magnitude can hide a high phase. Explain why a worst-phase watt guard can
   curtail earlier, without calling every earlier curtailment beneficial.
3. **Separate maximum voltage from voltage quality.** Compare maximum customer
   voltage and VUF. Balanced-current reactive absorption can lower the maximum
   voltage while changing VUF in a different direction. Compare matched sequence
   policies to isolate the effect of I₂ actuation.
4. **Show what the response costs.** Read delivered power, reactive absorption and
   peak current beside voltage benefits. A voltage reduction on a stiff feeder
   can cost substantial curtailment. Available-minus-delivered P includes any
   converter losses/limiting; it is not an energy figure. These are snapshots,
   not a daily yield calculation.
5. **End with the sensing boundary.** Use neutral displacement to explain why
   PG sensing, PN service voltage and VUF answer different questions. Local
   controls do not certify feeder-wide compliance or solve all unbalance.

Use fixed phase colors, identical axes across controller panels and route-distance
scatter points rather than lines joining unrelated feeder branches. For live
teaching, a later interactive view should expose source VUF, phase-load allocation,
PV availability and droop tuning with a linked phasor diagram and limit/current
indicators. It should show precomputed accepted equilibria or rerun the solver;
interpolated animation must not be presented as physical transient simulation.

## Limits, numerical checks and next steps

Use **207–253 V phase-to-neutral** and **360–440 V phase-to-phase**, corresponding
to the requested 230/400 V ±10% supply band. Energex's
[2024 DAPR, power quality section](https://www.energex.com.au/__data/assets/pdf_file/0011/1492382/Energex_Distribution_Annual_Planning_Report_2024.pdf)
reports this range. Its
[2025 DAPR](https://www.energex.com.au/__data/assets/pdf_file/0010/1846702/Energex_Distribution_Annual_Planning_Report_2025.pdf)
distinguishes the regulatory 207 V lower limit from its preferred lower operating
limit near 216 V. Do not substitute that operating preference for the requested
study band. These snapshot plots do not constitute a complete compliance
assessment with prescribed time aggregation or all power-quality criteria.

Preparation validates serialized/reparsed BMOPF JSON using determined power flow,
with no elastic KCL slack, in SI and per-unit scaling. It requires solved status,
finite voltage phasors and less than 0.01 V cross-scaling disagreement. The modified
three-phase demand case receives the same check. A second check fixes solved
phase P/Q injections from three representative controller cases and replays them
through BMOPFTools power flow; it requires less than 0.01 V disagreement and
verifies each three-wire inverter's phase-current sum is below 1 µA.

The controlled sweep uses Ipopt, tolerance 1e-8, `bound_relax_factor=0`, and a
100 kVA numerical power base. Accepted cases must have a publishable inner solve,
finite customer metrics, exact-versus-smooth current residual below 0.02 A and
converter current no more than 10 A plus 10 µA tolerance. Rejected cases remain
in `cases.csv` and are excluded from customer/device result tables. These are
local nonlinear equilibria: convergence does not prove uniqueness or dynamic
stability. No control mode is presumed universally better.

This branch starts from main before PR #49. The network controller API here uses
smooth encoding and `pwm_strategy=:NONE`. This campaign addresses fundamental
steady-state voltage/current behavior. It does not validate switching ripple,
harmonic emissions, thermal lifetime, certified firmware, or the PWM outer-loop
closure fixed in that PR. Inverter hardware values are illustrative.

Next experiments should add a balanced-load/source control, phase permutation,
load and source-VUF sweeps, objective-invariance and smoothing sensitivity checks,
and independently calibrated negative-sequence tuning. Preserve matched hardware
and input scenarios. Then add mixed single-/three-phase connections, larger
Australian LV feeders, and finally MV+LV coupling with validated aggregation and
transformer/source settings. Dynamic controller stability requires a separate
model and study.
