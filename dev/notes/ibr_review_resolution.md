# Advanced IBR review: resolution plan

Branch: `ibr-review-fixes` (from `main` @ `2c49436`). This note tracks every
finding of the October 2026 critical review of the advanced IBR model, the
phase-aware controllers and the two network studies. Each item stays here until
it has code, a regression test and documentation, or an explicit decision to
defer it. Update the status column in the same commit as the change.

Status: `todo`, `doing`, `done`, `deferred` (with reason).

## A. Correctness

| ID | Finding | Plan | Status |
|---|---|---|---|
| A1 | DC power availability is never enforced. `P_dc` is unbounded, so a controlled PV supplies its own losses from nothing: +57 W at 5 kW available with typical loss coefficients, and 20 W at zero PV with POC P = 0. The docs list DC availability as a protection priority. | The controller limits the positive-sequence active-power command to `p_available - L`, where `L = P_dc - 3Re(U1 conj(I1_target))` is the solved plant offset (all losses plus negative-sequence active power). The exact evaluator receives the same offset from the solved plant, as the existing capability backoff does. Report `dc_availability_margin_W = p_available - p_dc` in the device rows. At zero PV the inverter then imports its losses from the grid. Regression: lossy plant at full and zero PV. | todo |
| A2 | The common-scale limiter scales down I1 when I2 saturates a phase, so negative-sequence support curtails active power. This confounds the sequence-droop comparison. | Add `SequencePriorityLimiter(priority=:power | :balance)`. The reserved sequence is limited alone; the other sequence gets the largest exact per-phase scale `t` with `|a + t b| <= I_max`, using the rationalised quadratic root. The smooth model overestimates the root, so it stays conservative. Plant-aware backoff remains the final common protection layer. Exact and smooth tests, plus a network test showing power priority keeps P where common scaling curtails it. | todo |
| A3 | The split-DC model is internally inconsistent: neutral current divides by `C_u/C_Σ` (stiff source at ω) while the 2ω and carrier models assume an open source. Only unequal banks are affected. | Use one DC-source model at every frequency. Derivation: `i_s = ½J(C_l-C_u)/(C_Σ + jωC_uC_l Z_s)`, `i_cu = i_s - ½J`, `i_cl = i_s + ½J`. This also gives a fundamental total-bus ripple that enters the half-rail hull. Equal banks are unchanged for every `Z_s`. Decide on parameter naming against the existing `pwm_dc_source_*` fields during implementation. | todo |
| A4 | Sensing reference is tied to wire count: `neutral=nothing` forces phase-to-ground sensing, although many three-wire products sense phase-to-neutral. | Add `sensing_neutral` to `SequenceController`: `:device` (current behaviour), `:ground`, or a bus terminal name. Route the smooth stamp, the exact evaluation, the fixed-point oracles and the sensitivity helpers through one helper. Test the common-mode metamorphic contract with PN sensing on a three-wire plant. | todo |
| A5 | The `:dominant` conflict rule changes the Volt-var curve without a conflict: near the deadband edge q ∝ s², up to about 0.15·ε. | Scale the transition width by a conflict measure `ρ = min(l,h)/(min(l,h)+ε)`. The firmware law then reproduces the standard curve exactly when only one branch is active, while keeping continuity and zero at the tie. The smooth model uses the zero-preserving `nonnegative_min`. | todo |
| A6 | Docs say the positive-sequence policy is invariant to U0. Only its Volt-var is; the default worst-phase Volt-watt guard is not. | Correct the wording in the control-law page, the `SequenceController` docstring and the Australian README. | todo |
| A7 | `screen_fixed_point_gain` reports only Picard criteria (spectral radius, ∞-norm). Those are stricter than needed for lagged/incremental droop. | Also report `max_real_eigenvalue`. The criterion for a small-step relaxed iteration is `Re λ < 1`; document both. | todo |

## B. Parameters and studies

| ID | Finding | Plan | Status |
|---|---|---|---|
| B1 | The Australian Volt-watt curve (250→1, 260→0 of available) is stricter than AS/NZS 4777.2:2020 Australia A (253→1, 260→0.2 of rated). It drives the headline curtailment result. | Default to the Australia A curve on a rated basis; keep the pedagogical curve as a named sensitivity variant. Re-run and update the README tables. | todo |
| B2 | 2 mF for 6 kVA with no `dv2_max` means DC-link stress never binds. | Use a representative film link (about 50–80 µF/kVA, stated as an assumption) and declare `dv2_max` and `i_cap_max`. Report 2ω ripple and capacitor current. | todo |
| B3 | The filter is fixed at 0.05 + j0.15 Ω (about 0.6 % per unit for 6 kVA, roughly 10× below typical) and is not scaled with rating in ENWL. | Specify filters in per unit of each device's rating (stated assumption, e.g. X = 5 %, X/R = 10) in both scripts. | todo |
| B4 | The negative-sequence gain κ is in A/V and not rating-scaled. ENWL's κ is a no-op (g ≈ 0.02). | Both scripts declare κ in per unit of `I_rated/V_nom` and convert per device. Calibrate the ENWL gain to the Australian per-unit value. | todo |
| B5 | `impedance_angle = atan(0.3)` is not linked to the feeder. | Measure the local negative-sequence Thévenin angle with `inverter_control_network_voltage_sensitivity` (or a direct Z2 estimate). Record it per customer and choose `phi2` from it. Keep the fixed angle only as a sensitivity. | todo |
| B6 | ENWL variants use different smoothing widths, which confounds a matched comparison. | Use one smoothing configuration across all variants. | todo |

## C. Comparison set

| ID | Missing comparator | Plan | Status |
|---|---|---|---|
| C1 | Power-priority sequence allocation | Covered by A2; add as an Australian study mode. | todo |
| C2 | Grid-forming passive negative-sequence admittance (balanced EMF behind `Z_v` gives `I2 = -U2/Z_v`) | Constructor `VirtualImpedanceNegativeSequence(r, x)` for a flat admittance at `arg Z_v`, documented as the steady-state negative-sequence equivalent of a balanced-EMF voltage source. Add a study mode. | todo |
| C3 | Status quo: single-phase PV with per-phase Volt-var | Australian study variant keeping native single-phase BMOPF PV (PN_PER_PHASE) at the original customer phases. | todo |
| C4 | Four-leg per-phase Volt-var | Needs a FOUR_LEG controller with I0 enabled (new policy and stamp). Larger API change; planned as a follow-up PR unless time allows. | todo |

## D. Literature

| ID | Gap | Plan | Status |
|---|---|---|---|
| D1 | No BPSC/PNSC/AARC/IARC taxonomy | Cite Wang, Duarte & Hendrix (TPEL 2011) and map every implemented and candidate law onto it. | todo |
| D2 | Peak-current-constrained voltage support in R–X grids | Cite Camacho et al. (TIE 2013; TPEL 2018) and Afshari et al. (TIE 2017); position A2 against them. | todo |
| D3 | LV three-phase inverter unbalance predecessor | Cite Meersman et al. (EPSR 2011) as the predecessor of the conductance (`phi2 = 0`) law. | todo |
| D4 | Grid-code negative-sequence injection | Add k-factor context (VDE-AR-N 4110/4120, NC RfG), marked unverified like the other standards. | todo |
| D5 | Grid-forming natural unbalance response | Add a paragraph in the foundations and control-law pages connecting `grid_forming=true` and C2. | todo |

All new bibliography entries carry the existing † "not checked against live
sources" marker until verified.

## Order of work

1. A1, A2 (largest effect on study conclusions), then A5–A7.
2. A4, A3.
3. C2, then D1–D5.
4. B1–B6 plus C1–C3 in the study scripts; re-run the Australian study, and the
   ENWL dry run (its data checkout is external).
5. C4 as a follow-up unless scope allows.
