# [Verification and benchmark cases](@id ibr-verification)

A circuit-aware optimisation model should be tested at four levels: algebraic
identities, closed-form topology limits, published cases, and a higher-fidelity
oracle. A solver returning `LOCALLY_SOLVED` is necessary, but not sufficient.

## Tests implemented in the package

| Mechanism | Test oracle | Expected invariant |
| --- | --- | --- |
| Phase and neutral filter KVL | conductor power balance | `Pconv - Ppoc = Rp sum(Ix^2) + Rn In^2` |
| Primitive filter matrix | direct complex matrix product and diagonal equivalence | `ΔVxn = (Zf J)x - (Zf J)n`; a diagonal primitive reproduces scalar solves |
| LCL zero-capacitance limit | paired reduced/explicit solves | with `Cmid=0`, the two arms collapse to their summed impedance and `Iconv=Igrid` |
| LCL midpoint KCL | damped-capacitor admittance | `|Ish|=|Vmid|/|Rd+1/(jωCmid)|` and converter/grid currents differ |
| LCL real-power balance | terminal powers and branch currents | `p_filter_loss=Pconv-Ppoc=IcᴴRcIc+IgᴴRgIg+Rd Σ|Ish|²` |
| Scalar resonance estimate | closed form | reported `fres = sqrt((Lc+Lg)/(Lc Lg Cmid))/(2π)`; matrix mode reports no scalar estimate |
| Fortescue transform | prescribed sequence phasors | `i_neutral = 3 i_zero`; pure negative sequence has zero neutral current |
| Sequence limits | unbalanced voltage sources | each of `i_zero_max`, `i_positive_max`, `i_negative_max` bounds its reported component |
| 3-wire topology | KCL | the three phase-current phasors sum to zero under balanced and unbalanced voltage |
| Balanced switching hull | closed form | `Vdc_min = sqrt(6) U/m` for 3-leg and `2sqrt(2) U/m` for split DC |
| Sample resolution and audit | nested sample grids plus independent dense evaluation | increasing `n_samples` tightens the outer approximation; `switching_margin < 0` reveals a missed between-sample peak |
| 2ω power | Deakin et al. Table V | six prescribed current patterns produce `[0,1,1.5,3,0,1] Vph Iref` |
| 2ω PLECS comparison | Deakin et al. Table VII | all six rows are retained; non-cancelling analytical powers agree with PLECS within 0.3% and reported currents within 1% |
| Bus capacitance | capacitor energy balance | `dv2 = ripple/(2ω Ceq Vdc)`, with `Ceq=Cu Cl/(Cu+Cl)` for a split link |
| Split midpoint | charge balance | `dv_mid = i_neutral/[ω(Cu+Cl)]`; symmetric banks recover `i_neutral/(2ω Cdc)` |
| Split-bank mismatch | series charge and midpoint KCL | natural mean offset is `Vdc(Cu-Cl)/[2(Cu+Cl)]`; neutral shares are `Cu/(Cu+Cl)` and `Cl/(Cu+Cl)` |
| Bounded midpoint balancing | differential-charge identity | `v_mid_mean - v_mid_natural = 2q_mid_balance/(Cu+Cl)` and insufficient charge authority is infeasible |
| Capacitor RMS budget | orthogonal components | reported physical upper/lower currents combine their neutral share, common 2ω current, and reserved switching term |
| Capacitor thermal and ESR loss | weighted squared-current sum | rating constraints use `i_cap_thermal_*`; `p_cap_loss = ESRu Iu,th² + ESRl Il,th²` and enters `p_dc` exactly |
| Balanced carrier PWM | Mandrioli et al. Eqs. (40), (42) | numerical shared-carrier integration reproduces the published SPWM and centered-PWM DC-link ripple RMS at modulation indices 0.1, 0.3, and 0.5 within 0.2% |
| PWM scaling | carrier-period charge balance | switching-current RMS scales with AC current and is independent of `Cdc` and `f_sw` under the frozen-current assumption; switching-voltage ripple scales as `1/(Cdc f_sw)` |
| PWM reserve closure | post-solve carrier oracle | allocated switching current covers the predicted value, the capacitor bound remains satisfied, and a binding bank rating curtails export |
| PWM topology | paired 4-leg/split-link per-unit solves | the fourth leg participates in shared-carrier DC current; split-link SPWM uses the series-equivalent bus capacitance |
| DC harmonic KCL | finite source R–L in parallel with `Ceq` | `Ibridge,h + Icap,h + Isource,h = 0` at every retained carrier harmonic; the open-source and high-impedance limits converge |
| DC source loss | Parseval plus a per-unit closure solve | `Psource,sw = Rsource Isource,rms²`; bridge/source/capacitor diagnostics and the normalized network margin remain in physical units |
| Split rail switching voltage | series-capacitor charge | `Vupper=(Ceq/Cu)Vbus`, `Vlower=(Ceq/Cl)Vbus`, and the two rail ripples sum to total-bus ripple |
| DC antiresonance screen | parallel-admittance cancellation | a lossless `Lsource` satisfying `Ω²LsourceCeq=1` drives `pwm_dc_network_margin` to zero and returns non-finite voltage rather than a misleading bounded result |
| Split-link AC ripple | Mandrioli et al. Eqs. (15), (27) | carrier-harmonic phase and neutral RMS values agree across modulation indices 0.1, 0.3, and 0.5 within 0.3% |
| AC-ripple scaling | inductive harmonic circuit | phase and neutral ripple scale as `1/(L f_sw)` |
| Harmonic filter topology | reduced-L/LCL limiting cases | two series arms with `Cf=0` reproduce their summed impedance; a physical midpoint branch separates converter, shunt, and grid ripple |
| Neutral inductance and 3-wire projection | paired topology audits | increasing `Ln/L` suppresses four-leg neutral ripple; the 3-leg sum-zero subspace has exactly zero neutral ripple |
| AC-ripple reserve closure | physical total-RMS ratings in a per-unit solve | allocated carrier RMS covers every predicted phase/neutral value and `hypot(Ifund, Isw)` respects `i_max`/`In_max` |
| Split rating composition | paired bound cases | common and individual half-bank ratings compose; enlarging any rating cannot shrink the feasible set |
| Fourth-leg loss | fitted loss equation | phase currents plus `i_neutral` enter the per-leg loss sum |
| Unit system | representative paired solves | per-unit and raw-SI formulations agree on well-scaled cases, including LCL midpoint quantities |
| Controller Fortescue algebra | prescribed `U1/U2` phasors | the transform recovers both sequences and reconstructed three-leg current sums to zero |
| Controller closed forms | direct formulas | common power/current scales, sequence powers, voltage-oriented `I2`, and regularized ripple residual agree near machine precision |
| Controller metamorphic tests | rotated and cyclically relabelled phasors | scalar decisions are invariant and current phasors transform with the reference frame |
| Controller exact/smooth agreement | independent numeric evaluators plus solved local voltage | maximum phase-current residual remains below the declared tolerance |
| Controller conflict continuity | dense minimum-voltage sweep through equal branch severity | deployed `:dominant` command changes continuously and reverses sign without a winner-take-all jump |
| PWL smoothing refinement | widths 0.1, 0.05, and 0.025 V at a knot | bias is positive, approximately first order, and below 0.5% at 0.05 V for the documented curve |
| P/Q and Volt-watt bases | exact/smooth numeric policies plus partial-irradiance OpenDSS point | all priorities satisfy the capability circle; rated and available bases separate as specified |
| PV oversizing with P/Q priority | exact/smooth local-law sweep at DC/AC ratios 0.9, 1.0, 1.1, 1.25, and 1.4; stamped ratio 1.1 for all priorities | the full local capability sweep respects rating and priority identities; each representative stamped solve is publishable |
| Plant-aware saturation | binding converter `s_max` and binding `dv2_max` cases | the controller backs off to a publishable solution at the physical location rather than creating an infeasible equality |
| Current-target phasors | converter- and grid-target LCL cases | full complex current phasors—not only magnitudes—equal the stamped target |
| Selection objective | loss and zero objectives on the same controlled power flow | control requests and current phasors are invariant within solver tolerance |
| Controller formulation units | raw-SI fail-fast check and SI-valued numeric evaluators | unsupported raw-SI stamping cannot silently yield a scientific result; controller parameters and extracted results retain SI semantics |
| Balanced AC topology reduction | one four-leg model versus three independent 1φ models | aggregate `P,Q` and per-phase current/internal voltage agree; neutral and zero sequence vanish |
| Volt-watt external oracle | OpenDSS `PVSystem`/`InvControl` | balanced slope/saturation and a partial-irradiance rated-basis point agree in POC voltage and active power |

These tests live in `test/advanced_inverter_tests.jl` and
`test/inverter_control_tests.jl`. They intentionally include balanced cases,
pure sequence excitations, strong mixed unbalance, and binding
limits. A balanced feeder alone cannot validate a four-wire topology because it
removes the very neutral and ripple mechanisms being tested.

Numerically demanding PWM-reserve and DC-link closure regressions use the
per-unit formulation. Raw SI is retained for representative, well-scaled
`AdvancedInverter` unit-conversion checks. The coupled phase-aware controller
requires per-unit stamping and rejects raw SI because its saturated convergence
is not sufficiently reliable.

The controller derivations and calling conventions are documented in the
[phase-aware inverter-control API](@ref inverter-control-api). This page is the
canonical claim-to-test map and publication gate. In particular, neither a balanced
OpenDSS comparison nor the four-leg/three-single-phase reduction validates
negative-sequence control or shared-DC-link dynamics; those are mechanism tests,
not universal model equivalences.

## Infeasibility triage for fleet studies

Never discard a failed or non-publishable snapshot. Retain scenario identifiers,
termination status, starts, and all available residuals. Re-run it in the
following diagnostic order: zero controller objective, wider smoothing via
continuation, relaxed experimental service commands, and finally individually
relaxed physical limits. Classify the first relaxation that restores a solution
as controller-equilibrium, numerical-conditioning, or physical-capability
failure. Report each class as a fraction of every penetration cohort. Only a
physical-limit relaxation may support a hardware-sizing conclusion, and only
after multiple starts and a higher-fidelity representative case.

## Published cases worth lifting next

### Heidari and Geth (2024)

Their four-conductor algebraic model provides equation-level cases for:

1. a non-zero neutral series impedance under zero-sequence current;
2. equality between internal and external powers plus full conductor losses;
3. 3-leg versus 4-leg feasible regions under the same unbalanced POC voltage;
4. independent positive-, negative-, and zero-sequence current caps; and
5. balanced internal voltage under a GFM constraint.

The present suite covers the first, second, fourth, and fifth as invariants. The
first two now exercise both diagonal and mutually coupled primitive conductor
matrices. A direct digitisation of their plotted three-leg/four-leg capability
boundaries would add a valuable paper-to-code regression.

### Deakin, Heidari, and Deng (2025)

Table V is especially useful because the voltage and six current phasor patterns
are fully prescribed. It tests the unconjugated-product convention without
depending on an OPF optimum. Tables VI-VII then separate formula error from
PLECS controller/filter effects. All six Table VII rows are now literal
regression data. For the four non-cancelling cases, analytical and PLECS 2ω
powers agree within 0.3%, and the associated reported currents within 1%; the
balanced active/reactive rows retain their small numerical residuals. The
package separately tests its **time-domain RMS** capacitor current, which
divides a sinusoidal Fourier magnitude by ``\sqrt2``. The table's
``P_{2\omega}/V_{dc}`` current convention must not be used as an RMS bank rating
without that conversion. A machine-readable PLECS waveform archive would still
be stronger than copied summary values.

### Deakin, Heidari, and Deng (2026)

The OPF parameter sweeps in Figs. 4-6 are useful **monotonic and ordering tests**:

- increasing capacitor ripple rating cannot shrink the feasible set;
- a monolithic 4-leg link does not route neutral current through its capacitor;
- a split link trades neutral capability against capacitor ripple headroom; and
- reconfigurable four-leg-plus-split hardware dominates a fixed allocation when
  all other ratings are equal.

The first three properties are covered, including unequal split half-banks and
frequency-weighted bank ratings. The final one must wait for an explicit
hybrid topology; using the current `:FOUR_LEG` result as its surrogate would be a
model error.

### Mandrioli et al. (2021) and Hammami et al. (2020)

Mandrioli et al. give balanced-load RMS closed forms for sinusoidal PWM and
centered continuous PWM. These are unusually clean unit-test oracles because
they isolate shared-carrier correlation from the OPF optimum. The package now
evaluates both formulas at three modulation indices, then separately checks
current, capacitance, and switching-frequency scaling. Hammami et al. provide
the complementary three-phase four-wire split-capacitor derivation. The present
split-link test exercises asymmetric half-banks and the series-equivalent bus
capacitance, but digitising several of their modulation/load-angle maps remains
a useful independent regression target.

### Mandrioli et al. (2021), Viatkin et al. (2021), and the 2023 generalization

The split-capacitor paper provides closed-form phase and neutral switching RMS
currents under SPWM. Equations (15) and (27) are now literal tests at three
modulation indices, including the predicted `1/(L f_sw)` scaling. Viatkin et
al. add an explicit neutral inductor to the four-leg circuit; the package tests
their central monotonic result that increasing the neutral-to-phase inductance
ratio suppresses neutral ripple. The 2023 generalized treatment covers many
common-mode injections and hundreds of simulation/experimental operating
conditions. Its published numeric archive would be the best next fixture for
DPWM, THIPWM, and strategy-dependent neutral-inductor validation.

## Recommended higher-fidelity test harness

For each topology, build one parameter-identical averaged or switched model in
PLECS, PSCAD, or an open EMT tool and export machine-readable sweeps. A compact
validation matrix should vary:

- positive/negative/zero-sequence voltage independently;
- active/reactive current angle and current magnitude;
- DC voltage, capacitance, and half-bank mismatch;
- midpoint-balancer charge authority, bandwidth, and loss;
- phase and neutral inductance/resistance;
- converter-side/grid-side LCL split, capacitance, damping, and grid impedance;
- modulation index, sampling density, and switching strategy; and
- capacitor ESR versus frequency and temperature.

Compare converter/midpoint/POC fundamental phasors, both arm currents, filter-
capacitor current and damping loss, mean `P/Q`, ``|I_0|,|I_1|,|I_2|``, neutral RMS,
``|\widetilde S|``, DC-bus 2ω amplitude and phase, midpoint fundamental RMS and
mean offset, each half-bank's current spectrum, temperature, and loss. Store
inputs and expected outputs in a versioned
CSV/JSON fixture rather than recreating a plotted curve by eye.

Suggested acceptance bands are exact/near-machine tolerance for algebraic
identities, below 0.5% for closed-form ideal cases, and study-specific bands for
EMT comparisons. For switched simulations, report the Fourier window, transient
discard interval, RMS/peak convention, switching frequency, controller bandwidth,
and whether the DC source absorbs 2ω power.

## Model-risk checks before publication

For every claimed boundary point:

1. require non-negative `switching_margin`, increase `n_samples`, and demonstrate convergence;
2. for PWM-enabled points, require non-negative `pwm_reserve_margin` and
   `pwm_modulation_margin`, require every AC `*_reserved` value to cover its
   predicted ripple, then refine both sampling grids and `pwm_ac_harmonics`;
3. resolve from more than one physically sensible initial point;
4. check all reported conductor, sequence, rail, and capacitor residuals;
5. perturb the operating point inside and outside the claimed boundary;
6. state whether each value is RMS, peak amplitude, or Fourier magnitude; and
7. compare at least one representative point with a higher-fidelity oracle.

For an explicit LCL case, also verify that the fundamental frequency is well
separated from every passive modal resonance and that the resonance remains
acceptably damped after controller delay and grid impedance are included. The
reported scalar `filter_resonance_hz` is only an initial screening number.

The maintained sources behind these tests are listed in [IBR references](@ref
ibr-references).
