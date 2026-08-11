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
| Split rating composition | paired bound cases | common and individual half-bank ratings compose; enlarging any rating cannot shrink the feasible set |
| Fourth-leg loss | fitted loss equation | phase currents plus `i_neutral` enter the per-leg loss sum |
| Unit system | paired solves | SI and per-unit solves return the same SI diagnostics, including LCL midpoint quantities |

These tests live in `test/advanced_inverter_tests.jl`. They intentionally include
balanced cases, pure sequence excitations, strong mixed unbalance, and binding
limits. A balanced feeder alone cannot validate a four-wire topology because it
removes the very neutral and ripple mechanisms being tested.

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
2. resolve from more than one physically sensible initial point;
3. check all reported conductor, sequence, rail, and capacitor residuals;
4. perturb the operating point inside and outside the claimed boundary;
5. state whether each value is RMS, peak amplitude, or Fourier magnitude; and
6. compare at least one representative point with a higher-fidelity oracle.

For an explicit LCL case, also verify that the fundamental frequency is well
separated from every passive modal resonance and that the resonance remains
acceptably damped after controller delay and grid impedance are included. The
reported scalar `filter_resonance_hz` is only an initial screening number.

The maintained sources behind these tests are listed in [IBR references](@ref
ibr-references).
