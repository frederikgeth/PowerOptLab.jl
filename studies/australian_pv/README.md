# Australian LV photovoltaic control study

The [evidence index](evidence/README.md) identifies the current coupled thermal
campaign and gives checksum-verified restoration commands. Earlier result
tables below are pinned historical experiments; their source hashes have not
been changed to match the updated model.

The [matched switching-aware PF/OPF benchmark](engineering_benchmark.md)
compares all four local controller families with coordinated dispatch under
unchanged hardware and network limits. Its figures show per-inverter limits,
export recovery, the reactive-service tradeoff and thermal prediction coverage
from 32 representative EMT references. The former low-current prediction
failures remain in the historical archive.

The [combined switching/DC-voltage boundary experiment](engineering_boundary.md)
adds phase/neutral peak witnesses and compares averaged and switched voltage
extrema with fixed sizing. It retains the LF-only failure and documents the
carrier approximation and finite-window drift.

The [declared carrier validation domain](engineering_domain.md) covers 972
P/Q, sequence and modulation points with independent current/charge integration
and passive R/L orbit checks. Engineering-limit failures remain in the records.

## Historical extensions and fundamental comparison

The [mixed-controller thermal/energy extension](thermal_energy.md) adds
dispatch-dependent PWM capacitor limits and nonlinear DC-energy bounds with
fixed sizing. Mean, worst-phase, sequence and four-leg controls operate together
in each feeder case, with independent saved-point and EMT validation.

The [voltage-aware allocation extension](voltage_aware.md) adds controlled
unbalance magnitude/angle/load experiments, explicit energized/disconnected
operation and a synthetic 24-hour comparison using unchanged hardware. Its
[generated results and pedagogical figures](voltage_aware_results.md) cover 77
accepted cases. It uses staged Q→P→I₂ allocation and retains its own campaign
labels, distinct from the fixed-structure controller results below.

The [coordinated capability benchmark](coordinated.md) adds matched balanced-
current and free-I₂ dispatch references with fixed hardware and an explicit
customer phase-to-neutral voltage criterion. Its [generated comparison and
figure](coordinated_results.md) distinguish feasible local OPF results from
global upper bounds and distinguish extra sequence actuation from coordination.

The [four-leg research-control counterfactual](four_leg.md) adds an explicit
neutral path and per-phase VV+VW with one shared DC source. Its
[eleven accepted cases and phase/neutral teaching figure](four_leg_results.md)
keep phase/DC sizing fixed and state the additional neutral assumptions.

The [representative switched-circuit extension](emt.md) integrates AC currents
and DC voltage. Its [results and figures](emt_results.md) show close fundamental
agreement while retaining capacitor thermal failures under the fixed 1 A budget.
The fundamental study does not provide switching-thermal certification.

This is the Australian successor to the ENWL showcase. It starts with a small
real feeder, preserves a converted reference, and makes the teaching experiment's
changes explicit. The campaign compares **local control laws for 12 PV systems at
12 three-phase customer connections, with unbalanced loading**. Every three-phase variant
uses the same network, loads, hardware and available solar power. The single-phase
counterfactual keeps PV on each customer's original phase, with the same modified
three-phase loading.

## Reproduce historical experiments

These commands use the checked-out implementation. For replay of an archived
claim, use its recorded revision, environment and source hashes as described in
the [evidence index](evidence/README.md).

From the PowerOptLab repository root, with its pinned Julia environment installed:

```sh
julia --project=. scripts/australian_pv/prepare.jl
julia --project=. scripts/australian_pv/run.jl        # main, sensitivities, terminal validation
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
Outputs are overwritten on rerun. New campaigns also write
`physical_reports.csv` and `physical_checks.csv`: one report per advanced
inverter/base attempt and one row per independently recomputed bound/equality.
Rejected attempts remain in these diagnostics. Native single-phase devices are
outside the advanced-inverter certificate and have no entries in those files.
Acceptance requires each advanced inverter's physical report as well as the
existing network and exact/smooth checks.

Zero PV availability keeps the inverters energized: the voltage controls and
I₂ service remain active, and AC loss import supports them. This is a declared
operating-mode assumption, not a disconnected-inverter reference or a claim
about every product's night-time behaviour.

The terminal comparison and its 18-case validation sweep use one fixed 100 kVA
base and solver profile, with no per-case retries. Earlier rejected runs are
retained separately under `results/experimental_*`. Any future rejected campaign
still raises an assertion after saving its results; `summarize.jl` prints its
status table and withholds a partial matched metric comparison.

The refreshed checkpoint uses source commit `e19cbae`: all 152 canonical
main/sensitivity/terminal snapshots pass. Each environment records 52 source-file
hashes, and the advanced-device reports independently check declared physical
bounds. The modified-network replay audit passes 380 checks. These checks cover
the declared averaged model; this campaign disables PWM and is not EMT validation.
The separate 300-iteration 1 MVA scaling diagnostic retains a mean-control
iteration limit and is excluded from matched scientific comparisons. A
[separate 1 MVA repeat](results/experimental_terminal_common1000/cases.csv) with
one common 1,000-iteration budget accepts all three laws. Its
[comparison with the 100 kVA reference](results/experimental_terminal_common1000/base_comparison.json)
retains the numerical differences; this additional repeat changes no hardware.

A separate [local voltage-aware allocation reference](../../scripts/formulations/inverter_service_allocation_reference.json)
uses these unchanged three-leg hardware parameters at a fixed unbalanced terminal
voltage. Mean, worst-phase and sequence VV+VW all issue accepted Q-priority
commands, and energized zero-PV operation imports its losses. The allocator uses
staged physical solves and differs from the main campaign's algebraic limiters.
These local measurements are teaching inputs, not additional feeder equilibria;
source/phase perturbation campaigns and chronological integration remain pending.

- [Feeder schematic](figures/feeder.png), [phase voltages](figures/phase_voltages.png)
  and [benefits and effort](figures/tradeoffs.png) (generated locally).
- [Sensing and requests](figures/sensing_and_requests.png) explains a deliberately
  unbalanced local-voltage input using the public exact laws. It is a teaching
  fixture, not a feeder equilibrium. Its inputs and requests are saved in
  [teaching_fixture.json](results/teaching_fixture.json).
- [Power budget and headroom](figures/power_budget_and_headroom.png) separates
  delivered power, converter/filter losses and unused PV power, then shows
  physical limit utilization. Nearness to a limit alone does not identify a
  binding controller constraint.
- [Matched terminal controls](figures/terminal_control_comparison.png) shows the
  accepted mean/worst/sequence cases with reactive priority and a common PCC
  power boundary. This figure is separate from the research-profile comparison.
- [Numerical summary](results/summary.csv), [case acceptance](results/cases.csv),
  [negative-sequence design](results/negative_sequence_design.json),
  [baseline validation](results/validation.json),
  [modified-network validation](results/modified_validation.json), and the
  `results/sensitivity_*` directories.
- [Terminal-profile validation sweep](results/terminal_profile_validation/cases.csv)
  and its [solver/environment settings](results/terminal_profile_validation/environment.json).

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
*assumed*. Sensitivities cover selected influential assumptions; the remaining
unswept assumptions are limitations, not calibrated product values.
The derivations are in the [parameter-selection guide](../../docs/src/ibr/parameter_selection.md).

| Parameter | Value | Kind | Basis | Sensitivity |
|---|---|---|---|---|
| Rating | 6 kVA, I_r = 8.70 A | assumed | residential three-phase inverter class | — |
| Current limit `i_max` | 10 A | derived | ≥ 6 kVA/(3·207 V) = 9.66 A to keep rated VA at 0.9 pu | — |
| DC-link voltage | 750 V | derived | ≥ √6·266 V/0.96 = 679 V for the three-leg hull | — |
| Modulation limit | 0.96 | assumed | utilisation margin | — |
| Series filter | X = 1.32 Ω (5%), R = 0.132 Ω | assumed | Liserre et al.: total L ≤ 0.1 pu; X/R = 10 | — |
| DC-link capacitance | 300 µF (50 µF/kVA) | assumed | illustrative link | 150 µF and 2 mF, with independent thermal ratings |
| Capacitor thermal rating `i_cap_max` | 1 A RMS | assumed | independent line-frequency thermal budget; PWM disabled | 0.5 A; 2 mF at 1 A and 3 A |
| 2ω ripple limit `dv2_max` | 7.5 V peak (1%) | assumed | ripple allowance | via capacitance |
| Loss curve | 15 W + 2.0 W/A·Σ\|I\| + 0.3 W/A²·Σ\|I\|² | assumed | ≈ 2.2% at rated, ≈ 3.5% at 10% load | — |
| Volt-var | 207/220/240/258 V → +44/0/0/−60% of VA | sourced | Energex STNW1170 v10 (2025), Table 8 | — |
| Volt-watt | 253/260 V → 100/20% of rating | sourced | same, Table 9; rated basis | earlier 250/260 V → 100/0% of available |
| Sensing reference | phase-to-neutral (`"n"`) | assumed | three-wire power stage with neutral sensing lead | — |
| Negative-sequence angle φ₂ | median measured fleet Z₂ angle | derived | `results/negative_sequence_design.json` | — |
| Negative-sequence gain κ | for 20% aligned fleet attenuation | derived | same file; per-unit value recorded | — |
| Passive negative-sequence impedance | equal to the series filter | assumed | specific untuned comparator | matched droop gain and angle |
| Single-phase hardware counterfactual | 6 kVA, 28.7 A, native Australia A droops, phase-to-neutral | assumed/sourced | same curves; original customer phase | — |

Australia A knots and bases are confirmed in
[Energex STNW1170 v10 (2025), Tables 8–9](https://www.energex.com.au/__data/assets/pdf_file/0010/1072549/STNW1170-Standard-for-Small-IES-Connections.pdf).
The product standard is AS/NZS 4777.2:2020 including Amendment 2:2024
(effective August 2025). This source does not establish full conformity of
our multiphase sensing, priority or negative-sequence laws.

Checked 31 October 2026: Energex's
[installer standards index](https://www.energex.com.au/our-services/connections/residential-and-commercial-connections/solar-connections-and-other-technologies/connection-information-for-installers)
still lists the February 2025 STNW1170 release. Its Tables 8–9 retain the knots
used here. Section 4.3.3 also requires balanced AC output for multiphase IES;
the three-phase-inverter definition addresses simultaneous connection and
disconnection. Negative-/zero-sequence support and future four-leg comparisons
are research capability experiments. Using the sourced voltage curves does not
establish that those operating modes meet the connection requirements.

The two 2 mF sensitivities intentionally separate stored energy from thermal
capability. At the baseline 300 µF/750 V, 7.5 V peak 2ω ripple corresponds to
about 1.06 kW oscillating power and 1 A capacitor RMS current for the open source.
Increasing C while retaining a 1 A rating reduces voltage ripple but cannot
remove the thermal constraint. These are hardware counterfactuals, not a
capacitor product recommendation. `pwm_strategy=:NONE`, zero capacitor ESR and
an open line-frequency DC source remain explicit limitations.
The controller's protection geometry matters too: the research profile uses a
conservative triangle bound for voltage ripple and pure 2ω thermal current.
Joint carrier/ripple budgets use an RMS-sphere intersection to preserve their
orthogonal channels. A coordinated capability optimum remains a useful future
bound on utilisation.

## Matched controller experiment

The full-controller research baseline explicitly retains
`pq_priority=:proportional`, `current_target=:grid` and
`power_limit_location=:converter`. This preserves the matched common-scale versus
power-priority sequence comparison. It is an experimental controller profile,
not a claim about certified residential firmware.

The proposed equivalence-demo profile is named `TERMINAL_CONTROL_SETTINGS`:
`pq_priority=:var`, `current_target=:converter`, `power_limit_location=:poc`.
It allocates feasible Q before P and compares mean, worst-phase and
sequence-voltage VV+VW in `sensitivity_terminal_var_priority`. All three cases
pass, and `terminal_profile_validation` additionally covers both source levels
and zero, half and full PV with the same profile. The mean and
worst-phase regression demos use the same profile and check separately binding
current and apparent-power limits against independent oracles.
This profile also changes the interpretation of the 6 kVA boundary; it is not
a pure sensing-policy comparison. Equivalence to a hard Q-curve equality requires
that the requested reactive service is feasible.

Changing `current_target` changes the commanded arm. `i_max` always protects
converter-arm current; `i_grid_max` would separately protect grid-arm current.
The series-only filter here makes those currents equal, while converter and PCC
powers remain different. In PCC/var-priority mode, the common and power-priority
constructors both allocate Q→P→I₂, so constructor names alone must not be used
as evidence of different service priorities.

| Mode | Reactive-power voltage input | Active-power voltage input | Negative-sequence current | Research-baseline sequence limiting |
|---|---|---|---|---|
| `unity` | none; Q command zero | available PV | none | common scale |
| `mean_vv` | mean phase magnitude | available PV | none | common scale |
| `mean_vv_vw` | mean phase magnitude | mean phase magnitude | none | common scale |
| `worst_vv_vw` | worst-phase guards | largest phase magnitude | none | common scale |
| `sequence_vv_vw` | positive-sequence magnitude | largest phase magnitude | none | common scale |
| `sequence_droop` | as above | as above | designed admittance droop | common scale |
| `sequence_droop_power` | as above | as above | designed admittance droop | power priority |
| `gfm_virtual_impedance` | as above | as above | passive −U₂/Z_f | power priority |
| `single_phase_native` | per-phase magnitude (native) | per-phase magnitude (native) | n/a (single-phase) | native capability |

Baseline three-phase controls target grid current and sense phase-to-neutral
voltage. Their curves use target-current power semantics; the terminal profile
instead uses net PCC P/Q. Read cross-profile differences with that distinction.
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
quo uses native BMOPFTools droops. That case is solved with the OPF
builder, bus voltage bounds removed and maximum PV export. The earlier
`solve_pf` rating-base bug was fixed upstream in BMOPFTools #448; the explicit
export objective also resolves dispatch beneath the Volt-watt ceiling. P is the tighter
of availability and the Volt-watt cap.

VUF is `100 |V₂|/|V₁|`, with a/b/c sequence order. Subtracting a common neutral
voltage changes V₀, not V₁ or V₂. Neutral displacement is shown separately.

## Results

At full PV, source positive-sequence voltage 1.03 pu and source VUF 1%, maxima
across the 12 customer buses and totals across the 12 inverters:

<!-- results:start -->
| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.70 | -0.00 | 1.099 | 6.48 | 0.00 | 0.34 | 0.200 | 0.045 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.46 | 31.43 | 1.287 | 7.34 | 0.00 | 0.39 | 0.255 | 0.051 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.04 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Sequence + I₂ droop (common scale) | 252.58 | 0.925 | 0.287 | 23.57 | 10.77 | 0.510 | 4.18 | 1.44 | 7.50 | 0.049 | 1.000 |
| Sequence + I₂ droop (power priority) | 253.11 | 0.933 | 0.232 | 58.47 | 28.10 | 1.275 | 8.23 | 1.36 | 6.84 | 0.253 | 0.911 |
| Passive I₂ admittance (filter Z) | 253.37 | 0.967 | 0.217 | 58.47 | 28.17 | 1.275 | 8.52 | 1.41 | 7.45 | 0.253 | 0.993 |
| Single-phase PV, native AS/NZS 4777.2 | 254.54 | 1.099 | 1.310 | 59.93 | 30.58 | — | 22.56 | — | — | — | — |
<!-- results:end -->

<!-- interpretation:start -->
Read these as conditional results for this feeder, tuning and hardware:

1. **The ordinary voltage-control benefit is small.** Worst-phase VV+VW changes
   maximum voltage from 253.52 to 253.23 V while absorbing
   31.43 kvar. The short, stiff feeder limits the voltage benefit.
2. **Sequence support has a hardware cost.** With power-priority sequence allocation,
   maximum VUF changes from 0.987% to
   0.933%, with 58.47 kW delivered
   and 0.911 A maximum capacitor thermal RMS current.
3. **Capacitance and thermal rating are independent.** Compare `large_dc_link_2mF`
   (1 A) with `large_dc_link_3A` (3 A). Increasing capacitance alone cannot remove
   the thermal constraint; the unchanged 7.5 V peak allowance is a separate limit.
4. **The limiter decides how support is paid for.** In the matched research
   baseline, common scaling delivers 23.57 kW versus
   58.47 kW with power priority. Those cases retain the old
   converter-power/grid-current/proportional-PQ semantics. They do not describe
   the separate reactive-priority PCC voltage-control comparison.
5. **A passive law is only as well tuned as its impedance.** The filter-Z case
   differs in both gain and angle from the designed droop. The matched virtual
   impedance sensitivity removes those differences at this operating point;
   neither comparison establishes dynamic grid-forming performance.
6. **The single-phase case is a connection/hardware counterfactual.** Its phase allocation,
   native control implementation and absent converter losses differ from the
   three-phase experiment. Read voltage quality separately from delivered power.

Removing the imposed source negative-sequence voltage reduces maximum VUF for
sequence VV+VW to 0.059%;
balancing the loads while retaining that source gives
1.022%.
The headline unbalance result therefore depends strongly on the assumed source;
site representativeness requires measured source conditions.
The half-smoothing case changes curve, extrema and conflict widths; it
does not halve every numerical width in the controller. Available-minus-delivered
power includes converter/filter losses and limiting; these snapshots do not
measure daily energy or certify feeder compliance.

`terminal_var_priority` compares mean, worst-phase and sequence-voltage VV+VW
with reactive priority, converter-current actuation and a PCC VA circle. It records
3/3 accepted cases using one 100 kVA base and the same solver
profile. The separate `terminal_profile_validation` sweep covers the same three
laws at zero, half and full PV and both source levels
(18/18 accepted snapshots).
The regression also checks that a zero selection objective preserves dispatch.

The scalar capability root is normalized by service magnitude, preserving the
quadratic and smoothing law as service tends to zero. Earlier failed runs are
archived in `experimental_terminal_var_priority_before_fix` and
`experimental_pcc_fleet`. Numerical conditioning remains sensitive to other
solver/base configurations. The current `experimental_terminal_scaling_current`
diagnostic rejects mean VV+VW at the 300-iteration limit at 1 MVA, while worst-phase
and sequence VV+VW pass. No partial matched metric comparison is published for
that group. `experimental_terminal_numerical_scaling` retains an older sequence-case
failure; its source checkpoint differs. These diagnostics establish no physical
infeasibility or hardware-sizing limit.

A separate `experimental_terminal_common1000` repeat uses one common
1,000-iteration budget for all three laws at 1 MVA. All three pass; the maximum
device P/Q difference from the 100 kVA reference is below 0.00003 W/var and the
maximum phase-neutral voltage difference below 0.3 µV. The 300-iteration failure
remains recorded. This is a separately declared numerical profile, not a per-case
retry or a change to the inverter hardware.
<!-- interpretation:end -->

### Sensitivity to assumptions

Each block uses the same stressed operating point with the named parameter or
controller-comparison changes. The 2 mF/3 A case changes both hardware values
relative to the main case; compare it with 2 mF/1 A to isolate the rating effect.

<!-- sensitivity:start -->
**balanced_load**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.37 | 1.015 | 0.066 | 58.70 | -0.00 | 1.099 | 6.48 | 0.00 | 0.35 | 0.200 | 0.047 |
| Mean VV | 253.17 | 1.022 | 0.079 | 58.50 | 28.01 | 1.249 | 7.17 | 0.00 | 0.39 | 0.244 | 0.052 |
| Mean VV+VW | 253.17 | 1.022 | 0.079 | 58.50 | 28.01 | 1.249 | 7.17 | 0.00 | 0.39 | 0.244 | 0.052 |
| Worst-phase VV+VW | 253.15 | 1.023 | 0.082 | 58.46 | 31.27 | 1.285 | 7.33 | 0.00 | 0.40 | 0.255 | 0.053 |
| Sequence VV+VW | 253.17 | 1.022 | 0.079 | 58.50 | 27.99 | 1.249 | 7.17 | 0.00 | 0.39 | 0.244 | 0.052 |
| Sequence + I₂ droop (common scale) | 252.43 | 0.944 | 0.040 | 22.75 | 10.37 | 0.497 | 4.01 | 1.44 | 7.50 | 0.047 | 1.000 |
| Sequence + I₂ droop (power priority) | 253.07 | 0.953 | 0.088 | 58.47 | 28.06 | 1.275 | 8.19 | 1.36 | 6.81 | 0.253 | 0.908 |
| Passive I₂ admittance (filter Z) | 253.09 | 0.986 | 0.097 | 58.47 | 28.13 | 1.275 | 8.53 | 1.40 | 7.44 | 0.253 | 0.992 |
| Single-phase PV, native AS/NZS 4777.2 | 254.46 | 1.138 | 1.380 | 59.98 | 30.51 | — | 22.64 | — | — | — | — |

**balanced_source**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 252.31 | 0.051 | 0.254 | 58.70 | -0.00 | 1.099 | 6.48 | 0.00 | 0.02 | 0.200 | 0.002 |
| Mean VV | 252.11 | 0.059 | 0.239 | 58.50 | 28.03 | 1.249 | 7.17 | 0.00 | 0.02 | 0.244 | 0.003 |
| Mean VV+VW | 252.11 | 0.059 | 0.239 | 58.50 | 28.03 | 1.249 | 7.17 | 0.00 | 0.02 | 0.244 | 0.003 |
| Worst-phase VV+VW | 252.10 | 0.059 | 0.238 | 58.50 | 28.69 | 1.256 | 7.21 | 0.00 | 0.02 | 0.246 | 0.003 |
| Sequence VV+VW | 252.11 | 0.059 | 0.239 | 58.50 | 28.03 | 1.249 | 7.17 | 0.00 | 0.02 | 0.244 | 0.003 |
| Sequence + I₂ droop (common scale) | 252.11 | 0.059 | 0.239 | 58.50 | 28.03 | 1.249 | 7.18 | 0.01 | 0.06 | 0.244 | 0.008 |
| Sequence + I₂ droop (power priority) | 252.11 | 0.059 | 0.239 | 58.50 | 28.03 | 1.249 | 7.18 | 0.01 | 0.06 | 0.244 | 0.008 |
| Passive I₂ admittance (filter Z) | 252.12 | 0.056 | 0.239 | 58.50 | 28.03 | 1.249 | 7.27 | 0.11 | 0.56 | 0.244 | 0.074 |
| Single-phase PV, native AS/NZS 4777.2 | 253.45 | 0.212 | 1.332 | 60.00 | 30.55 | — | 22.40 | — | — | — | — |

**half_smoothing**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.46 | 31.42 | 1.287 | 7.34 | 0.00 | 0.39 | 0.255 | 0.051 |
| Sequence + I₂ droop (power priority) | 253.11 | 0.933 | 0.232 | 58.47 | 28.10 | 1.275 | 8.23 | 1.36 | 6.84 | 0.253 | 0.911 |

**large_dc_link_2mF**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.69 | -0.00 | 1.099 | 6.48 | 0.00 | 0.05 | 0.200 | 0.045 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.06 | 0.244 | 0.050 |
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.06 | 0.244 | 0.050 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.45 | 31.43 | 1.287 | 7.34 | 0.00 | 0.06 | 0.255 | 0.051 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.03 | 1.249 | 7.17 | 0.00 | 0.06 | 0.244 | 0.050 |
| Sequence + I₂ droop (common scale) | 252.58 | 0.925 | 0.287 | 23.58 | 10.77 | 0.510 | 4.18 | 1.44 | 1.13 | 0.049 | 1.000 |
| Sequence + I₂ droop (power priority) | 253.11 | 0.933 | 0.232 | 58.46 | 28.09 | 1.275 | 8.23 | 1.36 | 1.03 | 0.253 | 0.912 |
| Passive I₂ admittance (filter Z) | 253.37 | 0.967 | 0.217 | 58.46 | 28.17 | 1.275 | 8.52 | 1.41 | 1.12 | 0.253 | 0.993 |
| Single-phase PV, native AS/NZS 4777.2 | 254.54 | 1.099 | 1.310 | 59.93 | 30.58 | — | 22.56 | — | — | — | — |

**large_dc_link_3A**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.70 | -0.00 | 1.099 | 6.48 | 0.00 | 0.05 | 0.200 | 0.045 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.06 | 0.244 | 0.050 |
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.06 | 0.244 | 0.050 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.45 | 31.43 | 1.287 | 7.34 | 0.00 | 0.06 | 0.255 | 0.051 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.03 | 1.249 | 7.17 | 0.00 | 0.06 | 0.244 | 0.050 |
| Sequence + I₂ droop (common scale) | 252.89 | 0.855 | 0.232 | 58.30 | 28.16 | 1.387 | 10.00 | 3.49 | 2.70 | 0.294 | 2.396 |
| Sequence + I₂ droop (power priority) | 252.93 | 0.873 | 0.232 | 58.36 | 28.15 | 1.355 | 9.46 | 2.84 | 2.19 | 0.283 | 1.943 |
| Passive I₂ admittance (filter Z) | 253.39 | 0.962 | 0.213 | 58.45 | 28.20 | 1.291 | 8.91 | 1.82 | 1.44 | 0.259 | 1.280 |
| Single-phase PV, native AS/NZS 4777.2 | 254.54 | 1.099 | 1.310 | 59.93 | 30.58 | — | 22.56 | — | — | — | — |

**legacy_allocation**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.70 | -0.00 | 1.099 | 6.48 | 0.00 | 0.34 | 0.200 | 0.045 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.46 | 31.43 | 1.287 | 7.34 | 0.00 | 0.39 | 0.255 | 0.051 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.04 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Sequence + I₂ droop (common scale) | 252.58 | 0.925 | 0.287 | 23.57 | 10.77 | 0.510 | 4.18 | 1.44 | 7.50 | 0.049 | 1.000 |
| Sequence + I₂ droop (power priority) | 253.11 | 0.933 | 0.232 | 58.47 | 28.10 | 1.275 | 8.23 | 1.36 | 6.84 | 0.253 | 0.911 |
| Passive I₂ admittance (filter Z) | 253.37 | 0.967 | 0.217 | 58.47 | 28.17 | 1.275 | 8.52 | 1.41 | 7.45 | 0.253 | 0.993 |
| Single-phase PV, native AS/NZS 4777.2 | 254.54 | 1.099 | 1.310 | 59.93 | 30.58 | — | 22.56 | — | — | — | — |

**matched_virtual_impedance**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Sequence + I₂ droop (power priority) | 253.11 | 0.933 | 0.232 | 58.47 | 28.10 | 1.275 | 8.23 | 1.36 | 6.84 | 0.253 | 0.911 |
| Passive I₂ admittance (filter Z) | 253.37 | 0.967 | 0.217 | 58.47 | 28.17 | 1.275 | 8.52 | 1.41 | 7.45 | 0.253 | 0.993 |
| Passive I₂ admittance (matched gain/angle) | 253.11 | 0.933 | 0.232 | 58.47 | 28.10 | 1.275 | 8.23 | 1.36 | 6.84 | 0.253 | 0.911 |

**pedagogical_volt_watt**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.70 | -0.00 | 1.099 | 6.48 | 0.00 | 0.34 | 0.200 | 0.045 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Mean VV+VW | 253.09 | 0.986 | 0.245 | 49.64 | 27.72 | 1.058 | 6.32 | 0.00 | 0.33 | 0.188 | 0.044 |
| Worst-phase VV+VW | 252.92 | 0.986 | 0.254 | 42.14 | 30.83 | 0.955 | 5.81 | 0.00 | 0.31 | 0.158 | 0.041 |
| Sequence VV+VW | 252.95 | 0.985 | 0.256 | 42.06 | 27.42 | 0.912 | 5.60 | 0.00 | 0.29 | 0.147 | 0.039 |
| Sequence + I₂ droop (common scale) | 252.47 | 0.925 | 0.295 | 17.73 | 10.67 | 0.438 | 3.52 | 1.44 | 7.50 | 0.035 | 1.000 |
| Sequence + I₂ droop (power priority) | 252.80 | 0.930 | 0.254 | 42.64 | 27.51 | 0.950 | 6.61 | 1.38 | 6.98 | 0.159 | 0.931 |
| Passive I₂ admittance (filter Z) | 253.04 | 0.966 | 0.241 | 41.55 | 27.54 | 0.932 | 6.97 | 1.41 | 7.47 | 0.154 | 0.996 |
| Single-phase PV, native AS/NZS 4777.2 | 253.66 | 1.053 | 1.284 | 45.98 | 29.52 | — | 20.46 | — | — | — | — |

**small_dc_link_150uF**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.70 | -0.00 | 1.099 | 6.48 | 0.00 | 0.68 | 0.200 | 0.045 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.75 | 0.244 | 0.050 |
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.05 | 1.249 | 7.17 | 0.00 | 0.75 | 0.244 | 0.050 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.46 | 31.43 | 1.287 | 7.34 | 0.00 | 0.77 | 0.255 | 0.051 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.50 | 28.04 | 1.249 | 7.17 | 0.00 | 0.75 | 0.244 | 0.050 |
| Sequence + I₂ droop (common scale) | 252.47 | 0.953 | 0.307 | 11.27 | 5.06 | 0.311 | 2.01 | 0.72 | 7.50 | 0.011 | 0.500 |
| Sequence + I₂ droop (power priority) | 253.19 | 0.961 | 0.232 | 58.50 | 28.07 | 1.256 | 7.67 | 0.65 | 6.25 | 0.246 | 0.417 |
| Passive I₂ admittance (filter Z) | 253.32 | 0.977 | 0.225 | 58.50 | 28.10 | 1.255 | 7.81 | 0.67 | 7.19 | 0.246 | 0.479 |
| Single-phase PV, native AS/NZS 4777.2 | 254.54 | 1.099 | 1.310 | 59.93 | 30.58 | — | 22.56 | — | — | — | — |

**terminal_var_priority**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.51 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.46 | 31.43 | 1.287 | 7.34 | 0.00 | 0.39 | 0.256 | 0.051 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.51 | 28.04 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |

**thermal_0p5A**

| Control | Max PN voltage (V) | Max VUF (%) | Max neutral (V) | PV delivered (kW) | Q absorbed (kvar) | Converter loss (kW) | Max phase RMS current (A) | Max I₂ RMS (A) | Peak 2ω ripple (V) | Filter loss (kW) | Max capacitor thermal RMS (A) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unity PF | 253.52 | 0.984 | 0.249 | 58.69 | -0.00 | 1.099 | 6.48 | 0.00 | 0.34 | 0.199 | 0.045 |
| Mean VV | 253.27 | 0.987 | 0.233 | 58.49 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Mean VV+VW | 253.27 | 0.987 | 0.233 | 58.49 | 28.05 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Worst-phase VV+VW | 253.23 | 0.988 | 0.231 | 58.44 | 31.42 | 1.287 | 7.34 | 0.00 | 0.39 | 0.255 | 0.051 |
| Sequence VV+VW | 253.27 | 0.987 | 0.233 | 58.49 | 28.03 | 1.249 | 7.17 | 0.00 | 0.38 | 0.244 | 0.050 |
| Sequence + I₂ droop (common scale) | 252.47 | 0.953 | 0.307 | 11.27 | 5.06 | 0.311 | 2.01 | 0.72 | 3.75 | 0.011 | 0.500 |
| Sequence + I₂ droop (power priority) | 253.19 | 0.961 | 0.232 | 58.48 | 28.06 | 1.255 | 7.66 | 0.65 | 3.13 | 0.246 | 0.417 |
| Passive I₂ admittance (filter Z) | 253.32 | 0.977 | 0.225 | 58.48 | 28.09 | 1.255 | 7.81 | 0.67 | 3.60 | 0.246 | 0.479 |
| Single-phase PV, native AS/NZS 4777.2 | 254.54 | 1.099 | 1.310 | 59.93 | 30.58 | — | 22.56 | — | — | — | — |

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
5. **Compare connection hardware.** The single-phase case
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
phase P/Q injections from four representative controller cases and replays them
through BMOPFTools power flow; it requires less than 0.01 V disagreement and
verifies each three-wire inverter's phase-current sum is below 1 µA.
The campaign records Git revision/dirty state, tracked-source diff hash,
Project/Manifest hashes and all control settings in `environment.json`.
The final replay audit records per-file hashes of the model and study sources in
`modified_validation.json`, so the validated equations can be checked independently
of the dirty-checkout snapshot recorded at campaign start.

The terminal-profile audit independently replays each of its three stressed
equilibria with fixed phase P/Q injections. The regression covers all 18 source/
availability/policy combinations and checks selection-objective invariance.

The terminal runs fix a 100 kVA base, adaptive Ipopt barrier, automatic NLP
scaling disabled, tolerance 1e-8, `bound_relax_factor=0` and a 300-iteration cap.
All three policies use this same profile. Other numerical bases can still
fail to converge: the 1 MVA sequence-case diagnostic is retained separately in
`results/experimental_terminal_numerical_scaling`. Physical unit equivalence
does not imply identical nonlinear-solver convergence.

The baseline controlled sweep uses Ipopt, tolerance 1e-8,
`bound_relax_factor=0`, and a 100 kVA numerical power base. A non-publishable solve
or failed physical gate is retried at 1 MVA and then 10 kVA bases; the base that produced each accepted
equilibrium is recorded in `cases.csv`. Accepted cases must have a publishable inner solve,
finite customer metrics, exact-versus-smooth current residual below 0.02 A,
converter current no more than 10 A plus 10 µA, 2ω ripple within `dv2_max`,
capacitor thermal current within its independent rating (10 µA tolerance),
selected-location apparent power within 0.01 VA and dense switching margin
nonnegative to 10 µV, and
DC draw no more than the available PV power (10 mW tolerance). Rejected cases remain
in `cases.csv` and are excluded from customer/device result tables. These are
local nonlinear equilibria: convergence does not prove uniqueness or dynamic
stability. No control mode is presumed universally better.

The network controller uses smooth encoding and `pwm_strategy=:NONE`. This
campaign addresses fundamental steady-state voltage/current behaviour. It does
not validate switching ripple, harmonic emissions, thermal lifetime, certified
firmware, or PWM reserve closure. Hardware values follow the assumptions table
above; they are not a specific product.

The committed sensitivities now include separately balanced source and load,
half-width curve/extrema/conflict smoothing, the unrated legacy reference and terminal-var-priority comparison, independent
capacitor thermal ratings and a matched-gain/angle passive admittance.
The incomplete `results/experimental_pcc_fleet` diagnostic retains accepted and
rejected attempts from promoting the terminal profile across the full fleet
campaign. Its metadata declares the completed subset. This historical run
predates the scalar-root normalization fix. The earlier three-policy failures
are archived in `results/experimental_terminal_var_priority_before_fix`.
The present validation covers the three positive voltage-control policies;
promoting the terminal profile across the full negative-sequence-controller
benchmark still requires its own matched campaign.

Next experiments should add phase permutation, load and source-VUF sweeps,
broader objective-invariance checks,
a gain/angle sweep around the measured design, and a four-leg per-phase
comparator. Preserve matched hardware and input scenarios. Then add mixed
single-/three-phase fleets, larger
Australian LV feeders, and finally MV+LV coupling with validated aggregation and
transformer/source settings. Dynamic controller stability requires a separate
model and study.
