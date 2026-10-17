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
| A1 | DC power availability is never enforced. `P_dc` is unbounded, so a controlled PV supplies its own losses from nothing: +57 W at 5 kW available with typical loss coefficients, and 20 W at zero PV with POC P = 0. The docs list DC availability as a protection priority. | The controller limits the positive-sequence active-power command to `p_available - L`, where `L = P_dc - 3Re(U1 conj(I1_target))` is the solved plant offset (all losses plus negative-sequence active power). The exact evaluator receives the same offset from the solved plant, as the existing capability backoff does. Report `dc_availability_margin_W = p_available - p_dc` in the device rows. At zero PV the inverter then imports its losses from the grid. Regression: lossy plant at full and zero PV. Done: `dc_power_offset`, `p_dc_W` and `dc_availability_margin_W` rows; the oracles refresh the offset from each plant re-solve. Follow-up: the study gate caught a grid-forming night case overdrawing 4–8 W because the backoff scaled the loss import. The import is now an unscaled offset through every limiter (1e-3-of-rating selector, one-sided bias), and the `dv2_max` backoff uses the line-frequency source factor. Regression added. | done |
| A2 | The common-scale limiter scales down I1 when I2 saturates a phase, so negative-sequence support curtails active power. This confounds the sequence-droop comparison. | Add `SequencePriorityLimiter(priority=:power Done: `current_scale` reports the I1 scale; rows add `negative_sequence_scale`. \| :balance)`. The reserved sequence is limited alone; the other sequence gets the largest exact per-phase scale `t` with `\|a + t b\| <= I_max`, using the rationalised quadratic root. The smooth model overestimates the root, so it stays conservative. Plant-aware backoff remains the final common protection layer. Exact and smooth tests, plus a network test showing power priority keeps P where common scaling curtails it. Follow-up: the corrected study showed the binding 2ω ripple limit in the plant-aware backoff still common-scaled both sequences (power priority kept only 25.5 of 58.5 kW). The backoff now honours the priority in the stamp and the exact oracle; regression added. | done |
| A3 | The split-DC model is internally inconsistent: neutral current divides by `C_u/C_Σ` (stiff source at ω) while the 2ω and carrier models assume an open source. Only unequal banks are affected. | Use one DC-source model at every frequency. Derivation: `i_s = ½J(C_l-C_u)/(C_Σ + jωC_uC_l Z_s)`, `i_cu = i_s - ½J`, `i_cl = i_s + ½J`. This also gives a fundamental total-bus ripple that enters the half-rail hull. Equal banks are unchanged for every `Z_s`. Decide on parameter naming against the existing `pwm_dc_source_*` fields during implementation. Done: new `dc_source_r`/`dc_source_l` (line-frequency branch, default open) govern both the 2ω ripple and the split fundamental; the carrier branch stays `pwm_dc_source_*`. Added `dv_bus_fundamental` and the 2ω midpoint term to the hull and all audits. Behaviour change: unequal banks now default to equal neutral shares (open source). The old stiff-ω/open-2ω pairing is removed as inconsistent; `dc_source_r=0` gives the fully stiff limit. A nodal-circuit oracle covers four source models. | done |
| A4 | Sensing reference is tied to wire count: `neutral=nothing` forces phase-to-ground sensing, although many three-wire products sense phase-to-neutral. | Add `sensing_neutral` to `SequenceController`: `:device` (current behaviour), `:ground`, or a bus terminal name. Route the smooth stamp, the exact evaluation, the fixed-point oracles and the sensitivity helpers through one helper. Test the common-mode metamorphic contract with PN sensing on a three-wire plant. Done as `sensing_reference` (`:device`, `:ground` or a terminal name) plus the exported `sensing_terminal`. Also documented: the mean-phase law responds to U0 only at second order. | done |
| A5 | The `:dominant` conflict rule changes the Volt-var curve without a conflict: near the deadband edge q ∝ s², up to about 0.15·ε. | Scale the transition width by a conflict measure `ρ = min(l,h)/(min(l,h)+ε)`. The firmware law then reproduces the standard curve exactly when only one branch is active, while keeping continuity and zero at the tie. The smooth model uses the zero-preserving `nonnegative_min`. | done |
| A6 | Docs say the positive-sequence policy is invariant to U0. Only its Volt-var is; the default worst-phase Volt-watt guard is not. | Correct the wording in the control-law page, the `SequenceController` docstring and the Australian README. | done |
| A7 | `screen_fixed_point_gain` reports only Picard criteria (spectral radius, ∞-norm). Those are stricter than needed for lagged/incremental droop. | Also report `max_real_eigenvalue`. The criterion for a small-step relaxed iteration is `Re λ < 1`; document both. Review finding withdrawn: `maximum_real_eigenvalue`, `continuous_time_margin` and `alpha_max` already existed. Re-reading found `alpha_max = 2/\|1-λ\|`, which is wrong for complex eigenvalues; fixed to `2Re(1-λ)/\|1-λ\|²` with a rotation regression. | done |
| A8 | Found while verifying C2: `grid_forming=true` on `:THREE_LEG` imposed a ground-referred balanced EMF. With I0 = 0 that forces the POC zero-sequence voltage to zero, so any network with zero-sequence unbalance was infeasible. | For a three-leg bridge, constrain only the internal negative-sequence EMF to zero and bound the positive-sequence magnitude; U0 follows the network. Regression: plant I2 = -U2/Z_f on a source with V0 ≠ 0. Docstring and foundations §10 updated. | done |
| A9 | Found while preparing B5: network voltage sensitivity used 6×6 phase coordinates for a three-leg grid-current target. Independent phase perturbations violate ΣI = 0, so the perturbed plant solves cannot succeed. | Use sequence coordinates for every three-leg target. Regression: Z2 self-impedance of the resistive fixture at both targets. | done |
| A10 | Upstream (BMOPFTools `solve_pf`): `_strip_operational_limits!` deletes IBR `s_max`, which zeroes S_MAX/VAR_MAX droop bases. Native Volt-var/Volt-watt then gives zero output in a power flow, contrary to the docstring's claim that profile-governed IBRs remain determined. | Not fixed here (other repository). Report upstream. The status-quo study case (C3) uses the OPF builder with bus voltage bounds removed and a maximum-export objective. | deferred (upstream) |

## B. Parameters and studies

| ID | Finding | Plan | Status |
|---|---|---|---|
| B1 | The Australian Volt-watt curve (250→1, 260→0 of available) is stricter than AS/NZS 4777.2:2020 Australia A (253→1, 260→0.2 of rated). It drives the headline curtailment result. | Default to the Australia A curve on a rated basis; keep the pedagogical curve as a named sensitivity variant. Re-run and update the README tables. Australia A (BMOPFTools `Aus_A` preset) on a rated basis is the default; the pedagogical curve is the `sensitivity_pedagogical_volt_watt` run. | done |
| B2 | 2 mF for 6 kVA with no `dv2_max` means DC-link stress never binds. | Use a representative film link (about 50–80 µF/kVA, stated as an assumption) and declare `dv2_max` and `i_cap_max`. Report 2ω ripple and capacitor current. 300 µF (50 µF/kVA) with `dv2_max` = 7.5 V; 150 µF and 2 mF sensitivity runs. | done |
| B3 | The filter is fixed at 0.05 + j0.15 Ω (about 0.6 % per unit for 6 kVA, roughly 10× below typical) and is not scaled with rating in ENWL. | Specify filters in per unit of each device's rating (stated assumption, e.g. X = 5 %, X/R = 10) in both scripts. 5% series reactance at X/R = 10 in per unit of each rating (both scripts); per-unit loss curve added. | done |
| B4 | The negative-sequence gain κ is in A/V and not rating-scaled. ENWL's κ is a no-op (g ≈ 0.02). | Both scripts declare κ in per unit of `I_rated/V_nom` and convert per device. Calibrate the ENWL gain to the Australian per-unit value. Australian κ is derived for 20% aligned fleet attenuation and recorded in per unit. The first design divided the fleet-aggregate impedance by N twice, making κ 12× too small; this was caught in the results and corrected to κ = 1.63 A/V (43 pu). ENWL uses the same per-unit value, marked as an assumption. | done |
| B5 | `impedance_angle = atan(0.3)` is not linked to the feeder. | Measure the local negative-sequence Thévenin angle with `inverter_control_network_voltage_sensitivity` (or a direct Z2 estimate). Record it per customer and choose `phi2` from it. Keep the fixed angle only as a sensitivity. Measured with the fleet sequence sensitivity. Aggregate Z2 angles are 10–21° (median 19.4°): the feeder is resistive, so the earlier 16.7° angle was reasonable. The review was wrong to expect a transformer-dominated angle; what it lacked was the measurement. | done |
| B6 | ENWL variants use different smoothing widths, which confounds a matched comparison. | Use one smoothing configuration across all variants. ENWL variants now use default smoothing throughout. | done |

## C. Comparison set

| ID | Missing comparator | Plan | Status |
|---|---|---|---|
| C1 | Power-priority sequence allocation | Covered by A2; add as an Australian study mode. Mode `sequence_droop_power`. | done |
| C2 | Grid-forming passive negative-sequence admittance (balanced EMF behind `Z_v` gives `I2 = -U2/Z_v`) | Constructor `VirtualImpedanceNegativeSequence(r, x)` for a flat admittance at `arg Z_v`, documented as the steady-state negative-sequence equivalent of a balanced-EMF voltage source. Add a study mode. Done as `negative_sequence_virtual_impedance(r, x)`; verified against a grid-forming plant solve. | done |
| C3 | Status quo: single-phase PV with per-phase Volt-var | Australian study variant keeping native single-phase BMOPF PV (PN_PER_PHASE) at the original customer phases. Mode `single_phase_native`: native BMOPFTools single-phase PV on the original phases with the same loads, solved around upstream bug A10. | done |
| C4 | Four-leg per-phase Volt-var | Needs a FOUR_LEG controller with I0 enabled (new policy and stamp). Larger API change; planned as a follow-up PR unless time allows. Deferred to a follow-up PR: needs a FOUR_LEG controller with I0 enabled. | deferred (follow-up PR) |

## D. Literature

| ID | Gap | Plan | Status |
|---|---|---|---|
| D1 | No BPSC/PNSC/AARC/IARC taxonomy | Cite Wang, Duarte & Hendrix (TPEL 2011) and map every implemented and candidate law onto it. | done |
| D2 | Peak-current-constrained voltage support in R–X grids | Cite Camacho et al. (TIE 2013; TPEL 2018) and Afshari et al. (TIE 2017); position A2 against them. New references are numbered 54–59 in a separate section so existing item numbers stay stable. | done |
| D3 | LV three-phase inverter unbalance predecessor | Cite Meersman et al. (EPSR 2011) as the predecessor of the conductance (`phi2 = 0`) law. | done |
| D4 | Grid-code negative-sequence injection | Add k-factor context (VDE-AR-N 4110/4120, NC RfG), marked unverified like the other standards. | done |
| D5 | Grid-forming natural unbalance response | Add a paragraph in the foundations and control-law pages connecting `grid_forming=true` and C2. | done |

All new bibliography entries carry the existing † "not checked against live
sources" marker until verified.

## E. Documentation: parameter justification and code sync

No page justifies the hardware or control parameters. The studies label values
"illustrative" or "pedagogical" without sources, ranges or sensitivities, and
the docs-to-code check is manual.

| ID | Gap | Plan | Status |
|---|---|---|---|
| E1 | No parameter-selection guidance | New page `docs/src/ibr/parameter_selection.md`. Cover: filter impedance in per unit (L and LCL ranges, X/R); DC-link voltage by topology and the modulation limit; DC-link capacitance (film vs electrolytic, µF/kVA, ripple-current rating); current headroom over rated current; loss-curve coefficients from a datasheet efficiency curve; AS/NZS 4777.2 Volt-var/Volt-watt settings by region and basis; κ from the target `g = κ\|Z2\|` and rating; `phi2` from a measured Z2 angle; smoothing widths from stated error budgets. Each value gets a source (marked † until verified) or is explicitly a stated assumption. Done as `docs/src/ibr/parameter_selection.md`; the Australia A values are cross-checked against the BMOPFTools `Aus_A` preset. | done |
| E2 | Study parameters are not justified per study | Both study READMEs get an assumptions table: value, source or rationale, and the sensitivity run that covers it (links B1–B6). Assumptions table in the Australian README; per-unit assumptions listed in the ENWL README. | done |
| E3 | Docs made stale by this branch | Update as each item lands: A1 changes the `:rated` basis text (`inverter_controls.md` ≈141), the study boundary (≈596), the protection-priority text in `phase_aware_control_laws.md` (≈192, ≈492) and the Australian README "drawing active power for losses" (now true and stated precisely). Done alongside each item; stale statements updated in the control, foundations, verification and study pages. | done |
| E4 | No systematic docs/code sync check | Audit every IBR page against the code: API names, defaults, keyword lists and "implemented/future" status tables. Hard-coded numeric result tables (study READMEs) must be regenerated from the committed result CSVs, not hand-copied. Record the audit result here. Done: a script checked every backticked identifier on the IBR pages against package names, fields and keywords (no drift found); verification rows were updated for the new oracles. Five table cells whose inline pipes broke rendering were fixed and confirmed in the docs build. Study README tables are generated by `summarize.jl`. | done |

## Order of work

1. A1, A2 (largest effect on study conclusions), then A5–A7.
2. A4, A3.
3. C2, then D1–D5.
4. E1 (parameter guide), then B1–B6 plus C1–C3 in the study scripts with E2;
   re-run the Australian study, and the
   ENWL dry run (its data checkout is external).
5. E3 continuously with each item; E4 audit before opening the PR.
6. C4 as a follow-up unless scope allows.


## Outcome

Study results after all fixes (Australian LV3, full PV, 1.03 pu source, 1%
source VUF) are generated in `studies/australian_pv/README.md`. All 54 main
and 27 sensitivity cases are accepted and the replay audit passes 80/80. The
review's main concerns are borne out:

- the headline curtailment came from the pedagogical Volt-watt curve;
- unbalance attenuation is limited by the DC link (capacitance sensitivity);
- common scaling gives up most of the active power for that attenuation,
  while power priority (extended to the plant backoff) keeps it;
- the single-phase status quo has the worst unbalance and neutral displacement.

One review claim was withdrawn: the measured negative-sequence impedance angle
is resistive (about 19°), so the original fixed angle was reasonable.
