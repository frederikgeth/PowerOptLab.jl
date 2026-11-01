# Australian LV photovoltaic study: local inverter control versus OPF

**Current results:** [matched PF/OPF benchmark](engineering_benchmark.md) (generated
from audited records). **Earlier campaigns:** [study history](HISTORY.md).
**Evidence and restore commands:** [evidence index](evidence/README.md).

## Question

Twelve 6 kVA three-phase PV inverters sit on a real Australian LV feeder. Each
runs a local volt-var/volt-watt (VV/VW) controller with Australia A settings
from AS/NZS 4777.2. Four controller variants are mixed in one feeder:
mean-voltage, worst-phase and positive-sequence sensing on three-leg bridges,
and per-phase control on four-leg bridges. The inverter model includes converter
current, filter, DC link, PWM ripple and nonlinear DC-link energy.

The study asks:

1. Do the simulation (local fixed point) and the optimisation (OPF) give the
   same answer when they are given the same reactive service?
2. With the same reactive service, does central dispatch export more?
3. What does releasing the volt-var obligation buy, and when does the feeder
   voltage limit bind?

## Network

CSIRO's *Realistic Australian Medium Voltage Feeder with Associated Low Voltage
Feeders*, Geth et al. (2025), v1, [DOI 10.25919/ghnz-bk28](https://doi.org/10.25919/ghnz-bk28),
feeder `LV3_55bus` (CC BY-NC-SA 4.0; see `source/provenance.json`). It has 12
customers on a 500 kVA, 11 kV/433 V transformer with explicit neutral and
grounding. The feeder is short and stiff, so voltage differences between
customers are small.

Deliberate study changes: every customer gets a three-phase connection and
one inverter; demand is 1.5 kW per customer (PF 0.95) with a systematic
60/25/15 % A/B/C phase split; the 11 kV source has a 1 % negative-sequence
component. These are transparent hypothetical settings, not measured profiles.
Power flows are solved by BMOPFTools; equivalence with OpenDSS is not claimed.

## Inverter fixture

Each value is *sourced*, *derived* from a stated requirement, or *assumed*.

| Parameter | Value | Kind | Basis |
|---|---|---|---|
| Rating | 6 kVA (8.70 A at 230 V) | assumed | residential three-phase class |
| Phase current limit | 10 A RMS | derived | 6 kVA at 207 V needs 9.66 A |
| DC link | 750 V, 300 µF, band 735–765 V | derived / assumed | modulation hull; illustrative link |
| Filter | L only, X = 5 % (4.2 mH), X/R = 10 | assumed | Liserre et al.: total L ≤ 0.1 pu |
| Switching | 10 kHz centred PWM | assumed | |
| **Capacitor RMS rating** | **7 A** | derived | nameplate at the AS/NZS 4777.2 corners needs 6.22 A under the model's conservative carrier bound; +10 % |
| **Phase/neutral peak rating** | **23 A** | derived | nameplate needs 20.1 A (fundamental plus conservative ripple envelope); +10 % |
| Volt-var | 207/220/240/258 V → +44/0/0/−60 % of rated VA | sourced | Energex STNW1170 v10 (2025), Table 8 (Australia A) |
| Volt-watt | 253/260 V → 100/20 % of rating | sourced | same, Table 9 |
| Priority | reactive (var) priority | assumed | study setting; check against the AS/NZS 4777.2 clause before relying on it |

**Rated-capability gate.** Before any feeder run,
[`rated_capability.jl`](../../scripts/australian_pv/rated_capability.jl) checks
that each plant delivers 6 kVA at 207, 230 and 258 V with the Australia A
volt-var extremes. The limiting quantity must be apparent power or phase
current, never an internal budget. The capacitor current at rated unity PF
(5.7 A modelled) is close to the Kolar & Round (2006) estimate of about
0.59 × 9.2 A = 5.4 A; the difference is filter ripple.

The earlier fixture used a 1 A capacitor budget, chosen when the model had
only a 2ω term, and a 14 A peak (≈ √2 × 10 A, no ripple headroom). Under the
current carrier model that fixture delivers only 0.84–1.25 kVA of its 6 kVA
nameplate and fails every gate corner. It is kept as a labelled sensitivity
(`run_engineering_benchmark.jl ... legacy`), not as a result.

## Conditions and modes

| Condition | Source | PV per site (kW) | Purpose |
|---|---|---|---|
| `nominal` | 1.00 pu, 1 % U₂ | 1.25 / 3 / 5 | baseline |
| `high_source` | 1.03 pu | 1.25 / 3 / 5 | high upstream voltage |
| `rotated_unbalance` | 1.00 pu, 2 % U₂ rotated 90° | 1.25 / 3 / 5 | unbalance sensing |
| `night` | 1.03 pu | 0 | volt-var with no PV |
| `oversized_dc` | 1.03 pu | 2 / 4.8 / 8 | feeder-limited: arrays up to 133 % DC:AC |

| Mode | What is solved |
|---|---|
| `local` | Fixed point of each site's local VV/VW allocation and the network power flow. Local controllers do not enforce 253 V. |
| `opf_equal_service` | OPF maximising export with each site's exact local VV request and VW ceiling imposed at the OPF's own voltages, plus the 253 V planning limit. |
| `opf_vv_released` | OPF maximising export with Q free (volt-var released), plus the 253 V limit. Under AS/NZS 4777.2 this needs network-operator agreement. |

**PF/OPF agreement.** For every condition, the equal-service OPF *without*
the 253 V row must reproduce the local fixed point (voltages within 1 mV, site
P and Q within 0.5 W/var). Two separately built formulations give one
equilibrium; the audit fails if they do not.

## Findings

The generated [benchmark report](engineering_benchmark.md) has the numbers.

- With nameplate-capable hardware and slack feeder limits, volt-var costs
  almost no export: every site is limited by its PV availability. Releasing
  volt-var gains only converter/filter losses.
- Equal-service OPF matches local control wherever local control meets
  253 V. Central dispatch has no export advantage at equal service on this
  feeder.
- In the feeder-limited `oversized_dc` case, local control slightly exceeds
  the 253 V planning limit. Equal-service OPF must curtail to meet it;
  released OPF meets it with less absorption and more export.
- The earlier headline (4–13 kW "recovered" by coordination) was an artefact of
  the 1 A capacitor budget. See [history](HISTORY.md).

## What this does not show

- Conformity with AS/NZS 4777.2 or Energex STNW1170. STNW1170 §4.3.3 requires
  balanced AC output from multiphase inverters, so per-phase four-leg control is
  a research counterfactual.
- Dynamic behaviour: controller response times, PLL, feeder dynamics and the
  10-minute averaging of the 258 V sustained-overvoltage setting. Snapshots are
  treated as their own 10-minute average.
- Global optimality: OPF results are local optima of a nonconvex problem.
- Out-of-sample validation of the thermal bound. The EMT enclosure check uses
  the same selection rule after the bound was revised (roadmap R13).
- Scalability: one 55-bus feeder with 12 sites (roadmap R11/R15).

## Validation

- **Saved-point audit** (`audit_engineering_benchmark.jl`): rebuilds every
  saved point in SI, checks hardware, policy and the PF/OPF agreement, and
  rejects corrupted records (`test_engineering_benchmark.jl`).
- **Network check** (`network_si_audit.py`): separately coded Python
  reconstruction of feeder currents and limits at the fundamental frequency.
- **EMT** (`emt_circuit.jl`): separately coded averaged and switched circuit at a
  fixed POC voltage, one highest-capacitor-stress site per condition, mode and
  controller group. It checks that the modelled capacitor current encloses the
  switched one. It is not a PSCAD/PLECS comparison or a hardware test.
- **Carrier-model validation** ([domain](engineering_domain.md),
  [DC boundary](engineering_boundary.md)): checks the carrier approximations
  across P/Q, sequence and modulation, under the legacy fixture limits. Those
  limits only decide which points are labelled infeasible.

## Reproduce

```sh
julia --project=. scripts/australian_pv/rated_capability.jl
julia --project=. scripts/australian_pv/run_engineering_benchmark.jl /tmp/bench python3
julia --project=. scripts/australian_pv/audit_engineering_benchmark.jl /tmp/bench python3
julia --project=. scripts/australian_pv/prepare_engineering_benchmark_emt.jl /tmp/bench /tmp/emt_inputs.json
julia --project=. scripts/australian_pv/emt_circuit.jl /tmp/emt_inputs.json /tmp/emt
python3 scripts/australian_pv/audit_engineering_benchmark_emt.py /tmp/emt_inputs.json /tmp/emt /tmp/bench
```

The benchmark takes a few minutes on a laptop. Restore the committed records
from the [evidence index](evidence/README.md) before running saved-record tests.

## Glossary

- **Fixed point (local mode):** the operating point where every inverter's
  local controller output is consistent with the voltages it causes.
- **Equal service:** the same reactive power at every site as its local
  volt-var curve requests at the solved voltages.
- **Carrier bound:** an upper bound on PWM ripple currents over a switching
  period, used for capacitor and peak limits; conservative by construction.
- **Rated-capability gate:** the nameplate check described above.
- **Planning limit:** 253 V phase-neutral at every LV bus (AS 61000.3.100 upper
  steady-state value), used as an illustrative study criterion.
