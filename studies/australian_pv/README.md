# Australian LV photovoltaic control study

This is the Australian successor to the ENWL showcase. It starts with a small
real feeder, preserves a converted reference, and makes the teaching experiment's
changes explicit. The campaign compares **local control laws for 12 PV systems at
12 three-phase customer connections, with unbalanced loading**. Every variant
uses the same network, loads, hardware and available solar power. The status-quo
case keeps a single-phase PV inverter at each customer instead.

## Reproduce

From the PowerOptLab repository root, with its pinned Julia environment installed:

```sh
julia --project=. scripts/australian_pv/prepare.jl
julia --project=. scripts/australian_pv/run.jl        # main campaign + sensitivities
julia --project=. scripts/australian_pv/summarize.jl  # summary.csv and the tables below
julia --project=. scripts/australian_pv/audit.jl
python3 scripts/australian_pv/plot.py                 # optional figures
```

`summarize.jl` regenerates every results table in this README from the committed
CSVs; do not edit the tables between the `<!-- ... -->` markers by hand.
`LV3_three_phase_pv.bmopf.json` contains native ownership placeholders for the
advanced PVs; use `run.jl` for their controller requests and ratings. It is not a
standalone fixed-PQ dispatch case. `LV3_three_phase_no_pv.bmopf.json` is the
modified-load power-flow baseline.

The plot script needs Matplotlib and NumPy (`python3 -m pip install matplotlib`).
PNG and SVG figures are generated locally and ignored by Git. `run.jl` also
exposes `run(modes=..., levels=..., sources=..., curves=..., hw=...)` for smaller
experiments; the plot script requires the complete 54-case main campaign.
Outputs are overwritten on rerun.

- [Feeder schematic](figures/feeder.png), [phase voltages](figures/phase_voltages.png)
  and [benefits and effort](figures/tradeoffs.png) (generated locally).
- [Numerical summary](results/summary.csv), [case acceptance](results/cases.csv),
  [negative-sequence design](results/negative_sequence_design.json),
  [baseline validation](results/validation.json),
  [modified-network validation](results/modified_validation.json), and the
  `results/sensitivity_*` directories.

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
| Source positive sequence | 1.00 or 1.03 times the original 11 kV source |
| Source negative sequence | 1% of positive sequence, aligned with phase-A positive-sequence phasor |
| Source zero sequence | Zero |
| PV availability | 0, 0.5, 1.0 of 5 kW per customer |

These are transparent hypothetical connection and load changes, not measured
household profiles or a proposed Queensland connection approval. Aligning all
demand shares is a deliberately systematic imbalance. The 1% upstream
negative-sequence term makes the sequence-control mechanism visible on a short
feeder; it is not an assertion about the source dataset. At 1.03
positive-sequence voltage, the largest source phase magnitude is 1.0403 pu.
Three-phase upgrades are a modelling choice, not a claim that all Australian
houses already have three-phase service.

`MV21_328bus` is also copied, converted and run through BMOPFTools. A new standalone
Master includes its MV files only. Its load and transformer files are empty:
loads and distribution transformers belong to the associated LV feeders. Thus
this is an **unloaded MV connectivity/energization check**, not a validated loaded
MV operating point. The later integrated study must attach the LV feeders or
explicit, validated aggregate equivalents. It must not treat this unloaded result
as an MV voltage-impact study.

## Parameter assumptions

Each value is *derived* from the model and a stated requirement, *sourced*, or
*assumed*; assumed values that can change a conclusion have a sensitivity run.
The derivations are in the [parameter-selection guide](../../docs/src/ibr/parameter_selection.md).

| Parameter | Value | Kind | Basis | Sensitivity |
|---|---|---|---|---|
| Rating | 6 kVA, I_r = 8.70 A | assumed | residential three-phase inverter class | — |
| Current limit `i_max` | 10 A | derived | ≥ 6 kVA/(3·207 V) = 9.66 A to keep rated VA at 0.9 pu | — |
| DC-link voltage | 750 V | derived | ≥ √6·266 V/0.96 = 679 V for the three-leg hull | — |
| Modulation limit | 0.96 | assumed | utilisation margin | — |
| Series filter | X = 1.32 Ω (5%), R = 0.132 Ω | assumed | Liserre et al.: total L ≤ 0.1 pu; X/R = 10 | — |
| DC-link capacitance | 300 µF (50 µF/kVA) | assumed | film link sized for near-balanced operation | 150 µF and 2 mF |
| 2ω ripple limit `dv2_max` | 7.5 V peak (1%) | assumed | ripple allowance | via capacitance |
| Loss curve | 15 W + 2.0 W/A·Σ\|I\| + 0.3 W/A²·Σ\|I\|² | assumed | ≈ 2.2% at rated, ≈ 3.5% at 10% load | — |
| Volt-var | 207/220/240/258 V → +44/0/0/−60% of VA | sourced (verify) | AS/NZS 4777.2:2020 Australia A; BMOPFTools `Aus_A` preset | — |
| Volt-watt | 253/260 V → 100/20% of rating | sourced (verify) | same; rated basis | earlier 250/260 V → 100/0% of available |
| Sensing reference | phase-to-neutral (`"n"`) | assumed | three-wire power stage with neutral sensing lead | — |
| Negative-sequence angle φ₂ | median measured fleet Z₂ angle | derived | `results/negative_sequence_design.json` | — |
| Negative-sequence gain κ | for 20% aligned fleet attenuation | derived | same file; per-unit value recorded | — |
| Grid-forming virtual impedance | equal to the series filter | assumed | passive comparator | — |
| Status-quo single-phase inverter | 6 kVA, 28.7 A, native Australia A droops, phase-to-neutral | assumed/sourced | same curves; original customer phase | — |

## Matched controller experiment

The modes distinguish sensing, positive-sequence response, negative-sequence
actuation and current priority. Comparing `sequence_vv_vw` with
`sequence_droop` changes only the negative-sequence policy, and comparing
`sequence_droop` with `sequence_droop_power` changes only the limiter.

| Mode | Reactive-power voltage input | Active-power voltage input | Negative-sequence current | Current limiting |
|---|---|---|---|---|
| `unity` | none; Q command zero | available PV | none | common scale |
| `mean_vv` | mean phase magnitude | available PV | none | common scale |
| `mean_vv_vw` | mean phase magnitude | mean phase magnitude | none | common scale |
| `worst_vv_vw` | worst-phase guards | largest phase magnitude | none | common scale |
| `sequence_vv_vw` | positive-sequence magnitude | largest phase magnitude | none | common scale |
| `sequence_droop` | as above | as above | designed admittance droop | common scale |
| `sequence_droop_power` | as above | as above | designed admittance droop | power priority |
| `gfm_virtual_impedance` | as above | as above | passive −U₂/Z_f | power priority |
| `single_phase_native` | per-phase magnitude (native) | per-phase magnitude (native) | n/a (single-phase) | native current limit |

All three-phase controls target grid current and sense phase-to-neutral voltage.
Three-leg currents sum to zero, so these converters cannot inject zero-sequence
current. At zero PV availability the three-phase inverters remain energised and
may supply reactive and negative-sequence current. Their DC source then supplies
nothing, so they import their losses from the grid (negative POC active power).
This is a STATCOM-capable assumption, not a model of an inverter that switches off
at night. A supply limit and a controller knee are different quantities.

The negative-sequence design follows the parameter guide. The fleet's
negative-sequence voltage sensitivity is measured at zero current; φ₂ is the
median aggregate impedance angle; and κ gives a 20% aligned attenuation for the
median aggregate impedance shared by the 12 inverters. The single-phase status
quo uses native BMOPFTools droops. Because BMOPFTools `solve_pf` currently strips ([BMOPFTools.jl#447](https://github.com/frederikgeth/BMOPFTools.jl/issues/447))
IBR `s_max`, which zeroes S_MAX-based droops, that case is solved with the OPF
builder, bus voltage bounds removed and maximum PV export. P is then the tighter
of availability and the Volt-watt cap.

VUF is `100 |V₂|/|V₁|`, with a/b/c sequence order. Subtracting a common neutral
voltage changes V₀, not V₁ or V₂. Neutral displacement is shown separately.

## Results

At full PV, source positive-sequence voltage 1.03 pu and source VUF 1%, maxima
across the 12 customer buses and totals across the 12 inverters:

<!-- results:start -->
| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Peak converter current (A) | Peak I₂ (A) | Peak 2ω ripple (V) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.70 | 0.00 | 1.099 | 6.48 | 0.00 | 0.34 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 |
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.46 | 31.43 | 1.287 | 7.34 | 0.00 | 0.39 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.04 | 1.249 | 7.17 | 0.00 | 0.38 |
| Sequence + I₂ droop (common scale) | 252.58 | 0.925 | 0.287 | 23.57 | 10.77 | 0.510 | 4.18 | 1.44 | 7.50 |
| Sequence + I₂ droop (power priority) | 253.11 | 0.933 | 0.232 | 58.47 | 28.10 | 1.275 | 8.23 | 1.36 | 6.84 |
| Grid-forming virtual impedance | 253.37 | 0.967 | 0.217 | 58.47 | 28.17 | 1.275 | 8.52 | 1.41 | 7.45 |
| Single-phase PV, native AS/NZS 4777.2 | 254.54 | 1.099 | 1.310 | 59.93 | 30.58 | — | 22.56 | — | — |
<!-- results:end -->

<!-- interpretation:start -->
Read the table with the sensitivity blocks below. Six observations hold for
this case:

1. **Standard Volt-watt hardly curtails on this feeder.** With the Australia A
   curve (253 V knee, rated basis), every three-phase law delivers about
   58.5 kW of the 60 kW available. Most of the remainder is the converter loss
   that the plant now draws from the DC side. The pedagogical curve (250 V knee,
   available basis) curtails to about 42 kW. The earlier headline that
   worst-phase guards curtail heavily was a property of that curve, not of
   worst-phase sensing.
2. **Positive-sequence policies barely move the maximum voltage.** The feeder is
   short and stiff: about 28 kvar of absorption lowers the highest customer
   voltage by about 0.25 V, and all three-phase laws end within 0.3 V of 253 V.
   With phase-to-neutral sensing and a small phase spread, mean, worst-phase and
   sequence sensing give almost the same response.
3. **Negative-sequence support is limited by the DC link, not by the
   controller.** The measured design asks for 20% attenuation (κ = 1.63 A/V,
   φ₂ = 19°). With the 300 µF link the 2ω ripple limit binds and maximum VUF
   falls from 0.987% to 0.933%; with 150 µF to 0.961%; with 2 mF it is current
   limited at 0.873%. Capacitance buys unbalance attenuation.
4. **The limiter decides who pays.** Under the same ripple limit, common scaling
   delivers 23.6 kW because it cuts active power to reduce a ripple that comes
   mostly from I₂. Power priority, which now extends to the plant-aware ripple
   backoff, delivers 58.5 kW for almost the same VUF benefit.
5. **The passive grid-forming law is a weaker unbalance actuator here.** A
   balanced EMF behind the series filter gives about 0.75 A/V at 84°, far from
   the feeder's 19° impedance angle, so it attenuates less (0.967%) for a
   similar I₂.
6. **The single-phase status quo is worse on every unbalance measure.** Single-
   phase PV on the original phases raises maximum VUF to 1.10%, neutral
   displacement to 1.31 V (more than five times any three-phase case) and the
   highest customer voltage to 254.5 V. Its delivered power is not comparable:
   the native single-phase model carries no converter loss.
<!-- interpretation:end -->

### Sensitivity to assumptions

Each block changes one assumption at the same operating point.

<!-- sensitivity:start -->
**large_dc_link_2mF**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Peak converter current (A) | Peak I₂ (A) | Peak 2ω ripple (V) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.70 | -0.00 | 1.099 | 6.48 | 0.00 | 0.05 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.06 |
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.06 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.46 | 31.43 | 1.287 | 7.34 | 0.00 | 0.06 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.04 | 1.249 | 7.17 | 0.00 | 0.06 |
| Sequence + I₂ droop (common scale) | 252.89 | 0.855 | 0.232 | 58.30 | 28.16 | 1.387 | 10.00 | 3.49 | 2.70 |
| Sequence + I₂ droop (power priority) | 252.93 | 0.873 | 0.232 | 58.36 | 28.15 | 1.356 | 9.46 | 2.84 | 2.19 |
| Grid-forming virtual impedance | 253.39 | 0.962 | 0.213 | 58.45 | 28.20 | 1.291 | 8.91 | 1.82 | 1.44 |
| Single-phase PV, native AS/NZS 4777.2 | 254.54 | 1.099 | 1.310 | 59.93 | 30.58 | — | 22.56 | — | — |

**pedagogical_volt_watt**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Peak converter current (A) | Peak I₂ (A) | Peak 2ω ripple (V) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.70 | 0.00 | 1.099 | 6.48 | 0.00 | 0.34 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 |
| Mean VV+VW | 253.09 | 0.986 | 0.245 | 49.64 | 27.72 | 1.058 | 6.32 | 0.00 | 0.33 |
| Worst-phase VV+VW | 252.92 | 0.986 | 0.254 | 42.14 | 30.83 | 0.955 | 5.81 | 0.00 | 0.31 |
| Sequence VV+VW | 252.95 | 0.985 | 0.256 | 42.06 | 27.42 | 0.912 | 5.60 | 0.00 | 0.29 |
| Sequence + I₂ droop (common scale) | 252.47 | 0.925 | 0.295 | 17.73 | 10.67 | 0.438 | 3.52 | 1.44 | 7.50 |
| Sequence + I₂ droop (power priority) | 252.80 | 0.930 | 0.254 | 42.64 | 27.51 | 0.950 | 6.61 | 1.38 | 6.98 |
| Grid-forming virtual impedance | 253.04 | 0.966 | 0.241 | 41.55 | 27.54 | 0.932 | 6.97 | 1.41 | 7.47 |
| Single-phase PV, native AS/NZS 4777.2 | 253.66 | 1.053 | 1.284 | 45.98 | 29.52 | — | 20.46 | — | — |

**small_dc_link_150uF**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Peak converter current (A) | Peak I₂ (A) | Peak 2ω ripple (V) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.70 | 0.00 | 1.099 | 6.48 | 0.00 | 0.68 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.75 |
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.75 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.46 | 31.43 | 1.287 | 7.34 | 0.00 | 0.77 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.04 | 1.249 | 7.17 | 0.00 | 0.75 |
| Sequence + I₂ droop (common scale) | 252.47 | 0.953 | 0.307 | 11.27 | 5.06 | 0.311 | 2.01 | 0.72 | 7.50 |
| Sequence + I₂ droop (power priority) | 253.19 | 0.961 | 0.232 | 58.50 | 28.07 | 1.256 | 7.67 | 0.65 | 6.25 |
| Grid-forming virtual impedance | 253.32 | 0.977 | 0.225 | 58.50 | 28.10 | 1.255 | 7.81 | 0.67 | 7.19 |
| Single-phase PV, native AS/NZS 4777.2 | 254.54 | 1.099 | 1.310 | 59.93 | 30.58 | — | 22.56 | — | — |

<!-- sensitivity:end -->

These are outcomes of this case and tuning, not universal rankings or a claim of
regulatory compliance.

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
   policies to isolate the effect of I₂ actuation, and the two limiters to see
   whether unbalance support is paid for in active power.
4. **Show what the response costs.** Read delivered power, reactive absorption and
   peak current beside voltage benefits. A voltage reduction on a stiff feeder
   can cost substantial curtailment. Available-minus-delivered P includes any
   converter losses/limiting; it is not an energy figure. These are snapshots,
   not a daily yield calculation.
5. **Compare with the status quo and with hardware.** The single-phase case
   shows how much of the unbalance problem is created by single-phase PV in
   the first place. The capacitance sensitivities show when the DC link, not the
   controller, sets the limit.
6. **End with the sensing boundary.** Use neutral displacement to explain why
   phase-to-ground sensing, phase-to-neutral service voltage and VUF answer
   different questions. Local controls do not certify feeder-wide compliance or
   solve all unbalance.

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
100 kVA numerical power base. A solve that stops at a relaxed (acceptable) tolerance
is retried at 1 MVA and then 10 kVA bases; the base that produced each accepted
equilibrium is recorded in `cases.csv`. Accepted cases must have a publishable inner solve,
finite customer metrics, exact-versus-smooth current residual below 0.02 A,
converter current no more than 10 A plus 10 µA, 2ω ripple within `dv2_max`, and
DC draw no more than the available PV power (10 mW tolerance). Rejected cases remain
in `cases.csv` and are excluded from customer/device result tables. These are
local nonlinear equilibria: convergence does not prove uniqueness or dynamic
stability. No control mode is presumed universally better.

The network controller uses smooth encoding and `pwm_strategy=:NONE`. This
campaign addresses fundamental steady-state voltage/current behaviour. It does
not validate switching ripple, harmonic emissions, thermal lifetime, certified
firmware, or PWM reserve closure. Hardware values follow the assumptions table
above; they are not a specific product.

Next experiments should add a balanced-load/source control, phase permutation,
load and source-VUF sweeps, objective-invariance and smoothing sensitivity checks,
a gain/angle sweep around the measured design, and a four-leg per-phase
comparator. Preserve matched hardware and input scenarios. Then add mixed
single-/three-phase fleets, larger
Australian LV feeders, and finally MV+LV coupling with validated aggregation and
transformer/source settings. Dynamic controller stability requires a separate
model and study.
