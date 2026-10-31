# [Verification and benchmark cases](@id ibr-verification)

A circuit-aware optimisation model should be tested at four levels: algebraic
identities, closed-form topology limits, published cases, and a higher-fidelity
oracle. A solver returning `LOCALLY_SOLVED` is necessary, but not sufficient.

## Independent physical acceptance

Every standalone advanced, constrained or controlled inverter solve, including
fleet devices, now produces `physical_feasibility(result)`. The report
recomputes SI quantities from RMS voltage/current phasors, rather than trusting
squared-current, sequence-power or ripple auxiliaries. It checks phase and
neutral current reserves; converter sequence ratings; primitive/LCL KVL and
KCL; internal/GFM voltage bounds; selected-boundary VA; source-dependent DC
ripple; split midpoint charge/voltage; capacitor thermal budgets; and DC-source
availability supplied by the wrapper. Reported power and ripple quantities are
also checked against the independent calculation. Controlled-device reports
add actual-versus-commanded current residuals in amperes; constrained-device
reports add declared POC P/Q bounds. The exact fitted linear-loss
norm is used in the audit, making the known one-sided smoothing deficit
visible rather than hiding it in a reported auxiliary.

`InverterFeasibilityTolerances` declares physical floors of 10 µA, 10 µV,
0.01 VA and 1 nC, plus `1e-8` times the relevant rating. These do not depend on
the network power base or silently grow with solver tolerances. The rail
minimum is enclosed to 1 µV by default. All tolerances are inspectable in the
individual `InverterLimitCheck` records.

The full-cycle rail check bounds each signed inequality as a degree-two
trigonometric polynomial. On each interval, Taylor's theorem supplies a lower
bound using `abs(f′′) ≤ hypot(c₁,s₁) + 4 hypot(c₂,s₂)`; evaluated points supply
upper bounds. Adaptive subdivision encloses the global minimum over the full
electrical cycle, including split-link fundamental and 2ω displacement. Float64
roundoff is padded. This is a numerical enclosure of the declared averaged
waveform, not an interval-arithmetic proof of hardware behavior.

The standalone solvers refine the sampled hull by adding violating angles and
re-solving, for at most 12 refinements. `rail_refinement=false` retains the
sampled research formulation; the physical acceptance gate still rejects a
missed rail violation. Direct stamping remains composable: call
`inverter_physical_feasibility(inv, handles, solve_status(model))` after an
external solve and add/refine constraints in the host problem as appropriate.

`solve_status(result).publishable` requires solver success **and** physical
acceptance. MOI termination/primal statuses and the solver's `feasible` and
`optimal` flags remain unchanged for diagnosis. Candidate quantities may remain present
when a physically rejected solver-success point is returned; a termination
string alone is insufficient to select study results. Candidate reports can
also diagnose unsuccessful solves without promoting them to accepted results.
No-primal outcomes have an unchecked, non-passing report.

The fixed-current single-device/fleet replays used by network sensitivities and
fixed-point comparisons use the same publication gate, full-cycle rail refinement
and independently checked command residuals. A coarse six-angle zero-current
replay once returned publishable `LOCALLY_SOLVED` with roughly -35 V rail
headroom; that candidate now fails physical acceptance even when refinement is
explicitly disabled for diagnosis. One rejected device rejects the fleet.

Sensitivity perturbations and intermediate fixed-point iterates are virtual
imposed-current capability probes. They enforce the declared plant limits but
do not assert a PV availability bound. The returned fixed-point plant additionally
audits request availability using an upper enclosure of independently recomputed
DC power. An overdrawn final point cannot publish or report convergence. This
distinction permits perturbed-voltage numerical starts without calling their
virtual injections transient PV dispatch. Explicit source bounds can also be
applied to an individual current replay. `physical_feasibility` exposes the final
single-device or per-device reports; equilibrium use requires both `converged`
and a publishable final solve.

The report checks the declared inverter fundamental model and calibrated RMS
reserves. Standalone PWM-enabled solves additionally require finite carrier
predictions, converged reserves and actual predicted capacitor-thermal and
AC total-RMS ratings. A numerical closure tolerance cannot relax those ratings;
under-allocation at a binding rating triggers another update. Nonconvergence,
nonfinite harmonic networks and modulation violations prevent publication while
retaining the final candidate for diagnosis. These are model audits, not EMT.
Network-wide feasibility, regional compliance and transient/dynamic stability
remain outside this report; network replay and controller-policy checks still
provide separate evidence.

The balanced three-phase fixture uses exact ±120° angles. Rounded decimal
angles introduce a small U₀ and can make a strict balanced-grid LCL controller
incompatible with its passive shunt, even when a converter-target plant is
feasible. A separate common-mode case retains and tests that distinction.

Constrained-device acceptance additionally reconstructs POC P/Q from terminal
phasors, evaluates the declared smooth VV/VW relationship numerically, checks
fixed Q/PF and P/Q bounds, and checks zero/negative sequence on the selected
balanced-current arm. Standalone fixed-dispatch solves also check the requested
POC watts. These checks do not rely on the constrained model's lifted power or
sequence auxiliaries. Exact PWL tracking errors remain separate diagnostics:
the encoded smooth curve is the declared optimization policy, not an assertion
of exact firmware or regional-standard behavior.

A deliberately coarse Ipopt reproduction returns raw `LOCALLY_SOLVED` at
60 W / 0 var for a 2,500 W / −600 var request. All hardware checks pass, but
the independently reconstructed Q and fixed-dispatch checks reject publication.
Separate raw-success reproductions similarly reject violated PF and volt-watt
relationships. Solver tolerances do not change the physical acceptance floors.

Omitting a conductor rating does not create a hidden current bound. Numerical
current characteristics use declared ratings when present, otherwise physical
terminal-voltage starts (with declared voltage-rating fallback), rather than
the numerical coordinate voltage base. Three-leg KCL, LCL KCL and balanced-arm
equations retain physical scaling in SI and at different per-unit bases.
The shifted current norm's one-sided loss deficit is capped at 0.0001 W in
total across conducting legs. Phase paths retain any smaller rating-relative
width. The fourth leg instead uses its share of the absolute watt budget,
capped at 1 mA: its current often vanishes while the three phase currents are
nonzero. This avoids needless neutral-channel curvature without changing its
exact current bounds or the independent physical acceptance floors.
Composite current norms evaluate conductor sums before squaring them. Expanding
the square into a quadratic polynomial can otherwise produce a small negative
number through floating-point cancellation at balanced neutral current. Tests
exercise this cancellation in plain and staged models at different bases.
The unrated 1 mA fallback alone once understated demand by about 0.30 W for
`a_loss=100 W/A`; the unchanged exact-norm physical gate correctly rejected
that raw solver success. Regression solves now accept the corrected point in
SI and at 10 kVA, 1 MVA and 100 MVA bases without widening physical tolerances.

`dc_source_line_frequency_diagnostics` separately reports ω/2ω source currents,
their RMS/peak conventions, source resistance dissipation and a normalized
antiresonance cancellation screen. Independent capacitor-node KCL and period
integration check these quantities. The regulated DC-link input is the mean
power boundary for `p_dc`, request availability and `DCSourceBounds`; converter
and modeled capacitor ESR losses are included. Upstream source/DC–DC mean-path
loss, source-impedance harmonic dissipation and DC standby are outside that
boundary. The source diagnostic does not add an upstream power bound, a
component thermal rating or a dynamic-stability guarantee.

## Voltage-aware feeder coupling

`solve_voltage_aware_inverter_fleet(net, spec)` couples the local staged
Q→P→I₂ allocator to the native network equations. Every local allocation uses
separate sensed and physical terminal voltages. The feeder map imposes grid-arm
currents and the declared POC shunt exactly once. Intermediate feeder solves
are virtual numerical probes without inverter hardware constraints; they are
not transient dispatch samples. Returned phasors are reconstructed through the
full primitive/LCL filter and independently audited against hardware, source,
volt-watt and current-command limits. Publication also requires convergence of
the voltage/current map and successful local allocation stages.

Native zero-equality rows are checked in SI with physical tolerances independent
of the power base. This residual audit evaluates registered model rows; it is
not an independent circuit implementation or an independent audit of all native
operational inequalities. Small-feeder regressions separately calculate Ohm's
law and KCL from returned phasors. A failed network solve, local allocation,
controller voltage-domain check or iteration budget remains unpublishable and
retains its explicit diagnostic reason. Failure of a local NLP does not prove
plant or policy infeasibility.

`operating_modes=:energized` retains reactive control and AC loss imports at
zero PV. `:disconnected` represents an open grid-side contactor, removing the
entire converter/filter/POC-shunt AC model. A complete per-device dictionary
supports mixed fleets. Disconnected devices have zero AC injections and no
inverter physical report or modeled DC standby power. Unselected native IBRs
keep their native physics. `voltage_aware_inverter_rows(result)` records modes,
acceptance, powers, losses, current/ripple/rail quantities and sensed voltages.
Only accepted rows support equilibrium comparisons.

Three-leg reactive-priority control uses Q→P→I₂ ordering and optional
voltage-admittance negative-sequence requests. Four-leg
[`PerPhaseVoltVarWatt`](@ref) instead preserves the direction of its three Q
requests, maximizes their common service fraction, then maximizes equal POC
phase active powers. The shared DC source, neutral path and full filter remain
coupled physical constraints. Independent regressions cover a closed-form
binding neutral rating, mixed three/four-leg fleets, neutral KCL, full LCL losses,
night operation and an extreme zero-sequence-only voltage boundary. This is
solver-based quasi-static reference allocation; implicit controller stamping is
not provided for this family. Current-dependent ripple blending, dynamic
stability, global uniqueness and real-time guarantees remain outside this capability.
Warm starts and damping are numerical choices, not a transient inverter model.

## Tests implemented in the package

| Mechanism | Test oracle | Expected invariant |
| --- | --- | --- |
| Phase and neutral filter KVL | conductor power balance | `Pconv - Ppoc = Rp sum(Ix^2) + Rn In^2` |
| Primitive filter matrix | direct complex matrix product and diagonal equivalence | `ΔVxn = (Zf J)x - (Zf J)n`; a diagonal primitive reproduces scalar solves |
| LCL zero-capacitance limit | paired reduced/explicit solves | with `Cmid=0`, the two arms collapse to their summed impedance and `Iconv=Igrid` |
| LCL midpoint KCL | damped-capacitor admittance | `abs(Ish)=abs(Vmid)/abs(Rd+1/(jωCmid))` and converter/grid currents differ |
| LCL real-power balance | terminal powers and branch currents | `p_filter_loss=Pconv-Ppoc=IcᴴRcIc+IgᴴRgIg+Rd Σabs(Ish)²` |
| Scalar resonance estimate | closed form | reported `fres = sqrt((Lc+Lg)/(Lc Lg Cmid))/(2π)`; matrix mode reports no scalar estimate |
| Fortescue transform | prescribed sequence phasors | `i_neutral = 3 i_zero`; pure negative sequence has zero neutral current |
| Sequence limits | unbalanced voltage sources | each of `i_zero_max`, `i_positive_max`, `i_negative_max` bounds its reported component |
| 3-wire topology | KCL | the three phase-current phasors sum to zero under balanced and unbalanced voltage |
| Balanced switching hull | closed form | `Vdc_min = sqrt(6) U/m` for 3-leg and `2sqrt(2) U/m` for split DC |
| Sample resolution and audit | nested grids, full-cycle Taylor enclosure and adaptive separating cuts | a coarse-grid false optimum cannot publish; feasible dispatch remains available after refinement |
| Offline current allocation | analytic VA-circle projection across bases; primitive LCL reconstruction; rail-conflict witness | accepted alternatives meet all declared bounds; a failed solve never establishes plant infeasibility; inconclusive rail enclosures remain unresolved |
| Local Q→P→I₂ allocation | analytic Q-priority circle; fixed-boundary Australian mean/worst/sequence laws; primitive LCL/grid commands | every stage passes physical and policy checks; lower-priority stages preserve achieved powers; explicit override is distinct from a failed strict request |
| 2ω power | Deakin et al. Table V | six prescribed current patterns produce `[0,1,1.5,3,0,1] Vph Iref` |
| 2ω PLECS comparison | Deakin et al. Table VII | all six rows are retained; non-cancelling analytical powers agree with PLECS within 0.3% and reported currents within 1% |
| Bus capacitance | capacitor energy balance | open line-frequency source: `dv2 = ripple/(2ω Ceq Vdc)`, with `Ceq=Cu Cl/(Cu+Cl)` for a split link; ideal stiff source: `dv2 = 0` |
| Split midpoint | charge balance | equal banks: `dv_mid = i_neutral/[ω(Cu+Cl)]` for every source model, recovering `i_neutral/(2ω Cdc)` |
| Split-bank mismatch | series charge, midpoint KCL, and an independent three-node nodal solve over open, stiff, resistive and R–L sources | natural mean offset is `Vdc(Cu-Cl)/[2(Cu+Cl)]`; open source: equal neutral shares, `dv_mid = In(1/Cu+1/Cl)/(4ω)` and `dv_bus_fundamental = In abs(1/Cl-1/Cu)/(2ω)`; stiff source: shares `Cu/(Cu+Cl)`, `Cl/(Cu+Cl)` and no bus ripple |
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
| DC availability | lossy controlled plant at full and zero irradiance, both current targets | `p_dc <= p_available`; at zero irradiance the plant imports exactly its losses; a filter-only offset equals `3 R abs(I)^2` |
| Sequence-priority limiting | exact quadratic fill, smooth/exact comparison and a network solve | `abs(a + t b) <= Imax` with equality when binding; power priority keeps `I1` while common scaling curtails it |
| Grid-forming negative sequence | three-leg grid-forming plant with a zero-sequence source unbalance | the plant solves (no forced `V0 = 0`) and `abs(I2) = abs(U2)/abs(Zf)`, matching `negative_sequence_virtual_impedance` |
| Three-leg network sensitivity | converter- and grid-current targets on a resistive fixture | sequence coordinates at both targets; the negative-sequence self-impedance equals the line resistance |
| Sensing reference | three-wire plant sensing `"n"`, `:ground` and `:device` on a floating neutral | sensed phasors equal the declared terminal differences; `U1`/`U2` laws are invariant, worst-phase laws move at first order and mean-phase laws only at second order |
| Relaxed fixed-point step | rotation map with complex eigenvalues | `alpha_max = 2Re(1-λ)/abs(1-λ)^2` is the exact contraction boundary of the relaxed iteration |
| Split rating composition | paired bound cases | common and individual half-bank ratings compose; enlarging any rating cannot shrink the feasible set |
| Fourth-leg loss | fitted loss equation | phase currents plus `i_neutral` enter the per-leg loss sum |
| Unit system | representative paired solves | per-unit and raw-SI formulations agree on well-scaled cases, including LCL midpoint quantities |
| Controller Fortescue algebra | prescribed `U1/U2` phasors | the transform recovers both sequences and reconstructed three-leg current sums to zero |
| Controller closed forms | direct formulas | common power/current scales, sequence powers, voltage-oriented `I2`, and regularized ripple residual agree near machine precision |
| Controller metamorphic tests | rotated and cyclically relabelled phasors | scalar decisions are invariant and current phasors transform with the reference frame |
| Controller zero-sequence reference | common-mode offset added to all three phasors | `U1`/`U2`, both curve outputs of the sequence comparator, and the whole command are invariant; the physical extrema and every phase-magnitude comparator move |
| Capability-backoff sign safety | offsets at, below, and above the per-leg limit | the smooth allocation factor stays in `[0,1]` and matches the exact allocator within the declared selector width, so an exhausted limit cannot reverse the command |
| Controller exact/smooth agreement | independent numeric evaluators plus solved local voltage | maximum phase-current residual remains below the declared tolerance; this bounds controller-algebra smoothing at the solved point, **not** the distance between the two network equilibria |
| Controller conflict continuity | dense minimum-voltage sweep through equal branch severity | deployed `:dominant` command changes continuously and reverses sign without a winner-take-all jump |
| PWL smoothing refinement | widths 0.1, 0.05, and 0.025 V at a knot | bias is positive, approximately first order, and below 0.5% at 0.05 V for the documented curve |
| P/Q and Volt-watt bases | exact/smooth numeric policies plus partial-irradiance OpenDSS point | all priorities satisfy the capability circle; rated and available bases separate as specified |
| PV oversizing with P/Q priority | exact/smooth local-law sweep at DC/AC ratios 0.9, 1.0, 1.1, 1.25, and 1.4; stamped ratio 1.1 for all priorities | the full local capability sweep respects rating and priority identities; each representative stamped solve is publishable |
| Plant-aware saturation | binding selected-location `s_max` and binding `dv2_max` cases | the controller backs off to a publishable solution at the physical location rather than creating an infeasible equality |
| Current-target phasors | converter- and grid-target LCL cases | full complex current phasors—not only magnitudes—equal the stamped target |
| Selection objective | loss and zero objectives on the same controlled power flow | control requests and current phasors are invariant within solver tolerance |
| Controller formulation units | raw-SI fail-fast check and SI-valued numeric evaluators | unsupported raw-SI stamping cannot silently yield a scientific result; controller parameters and extracted results retain SI semantics |
| Balanced AC topology reduction | one four-leg model versus three independent 1φ models | aggregate `P,Q` and per-phase current/internal voltage agree; neutral and zero sequence vanish |
| Volt-watt external oracle | OpenDSS `PVSystem`/`InvControl` | balanced slope/saturation and a partial-irradiance rated-basis point agree in POC voltage and active power |
| Strict-Q feasible-set operation | independent `sqrt(min(Smax,3VImax)^2-Q^2)` oracle | current and S saturation curtail P while Q remains fixed; proportional forward backoff deliberately differs |
| Operating modes and source direction | signed PF equation, fixed Q, external free-Q objective and charging/loss budget | PF signs work for import/export/zero; explicit PV DC bound forbids charging but permits AC loss import |
| Operating-row scaling | contradictory dispatch at 6 kVA, 1 MVA and 100 MVA bases | incompatible fixed P / DC bounds cannot publish a false feasible result |
| Current allocation | independent Fortescue transform and LCL arm phasors | `:free` permits I2; converter/grid balance differ with midpoint current; three-leg zero sum remains physical |
| Four-leg / split capability optimisation | neutral-current, split-bank RMS and midpoint charge closed forms | tight neutral, individual bank and midpoint limits bind at predicted currents |
| Independent feeder reproduction | separately solved forward/maximum-export equilibria at 215, 230, 259 V | matched unsaturated laws agree within 0.2 W/var, 1 mA and 1 mV |
| Proportional limiter reproduction | current-ray maximum scale plus closed-form S/current oracle | optimiser scale agrees within `1e-7`; forward conservative backoff within 0.4 W and 0.2 var |
| Sequence-priority reproduction | fixed reserved sequence, ray fill and independent disk intersection | both priorities preserve the reserved sequence and agree with exact local allocation within 1 µA; fixed measurement only |
| P–Q teaching geometry | independently computed circle/ray intersections at 230/207 V | proportional, strict Q, signed fixed PF and free-Q solves match closed forms |
| Unbalanced feeder/LCL reproduction | separate equilibria, feeder KVL and capacitor admittance | matching Q/current locations agrees within 0.2 W/var, 1 mA and 1 mV; unmatched locations deliberately differ |
| Strict-service versus allocation infeasibility | apparent-power lower bound and zero-sequence KCL contradiction | impossible Q and incompatible balanced-grid allocation cannot publish; releasing grid balance yields the analytically predicted grid I0 |
| Sequential-to-equilibrium workbook | guarded affine-branch enumeration, independent fixed-point formula and sequential energy replay | cycling numerical iteration can coexist with a unique equilibrium; ceilings change dispatch freedom; state and hysteresis examples preserve declared history |
| Full multi-inverter network OPF | seven-bus four-wire feeder, three matched controls, native/hook voltage limits and independent SI radial replay | control agreement within 0.2 W / 0.1 var / 1 mV; coordinated 244 V and 12 A limits bind; KCL, KVL, power balance, base invariance and failed-output masking pass |

Physical certificate and full-cycle enclosure tests live in
`test/inverter_feasibility_tests.jl`. The other tests live in
`test/advanced_inverter_tests.jl` and
`test/inverter_control_tests.jl`, with operating constraints in
`test/inverter_operating_constraints_tests.jl`,
`test/inverter_formulation_variants_tests.jl` and
`test/inverter_formulation_reproduction_tests.jl`. The executable
[reproduction cases](formulation_reproductions.md) state each equivalence
condition and tolerance. They intentionally include balanced cases,
pure sequence excitations, strong mixed unbalance, and binding
limits. A balanced feeder alone cannot validate a four-wire topology because it
removes the very neutral and ripple mechanisms being tested.

Numerically demanding PWM-reserve and DC-link closure regressions use the
per-unit formulation. Raw SI is retained for representative, well-scaled
`AdvancedInverter` unit-conversion checks. The coupled phase-aware controller
requires per-unit stamping and rejects raw SI because its saturated convergence
is not sufficiently reliable.

The controller derivations and calling conventions are documented in the
[phase-aware control modeling guide](@ref inverter-control-api). This page is the
canonical claim-to-test map and publication gate. In particular, neither a balanced
OpenDSS comparison nor the four-leg/three-single-phase reduction validates
negative-sequence control or shared-DC-link dynamics; those are mechanism tests,
not universal model equivalences.

### Numerical regression settings

The DC-ripple saturation and LCL current-target regressions use `tol=1e-8`
with `nlp_scaling_method="none"`; the ripple fixture additionally uses
`mu_strategy="adaptive"` to avoid monotone-barrier stagnation on Julia 1.10.
Their equations already include per-unit and
device normalization. In the grid-target fixture, Ipopt's additional gradient
scaling stopped at acceptable tolerance; solving with the declared row scales
reaches the requested convergence criterion and satisfies the physical
current/power comparisons. These tests require publishable solutions and keep
their physical tolerances; the ripple assertions are active rather than skipped.
This is a measured configuration for these fixtures, not a universal solver
setting. Further tightening `tol` can encounter the formulation's numerical
accuracy floor instead of improving the physical answer.

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

## Comparing exact-law and smooth-law equilibria

The `exact_smooth_current_residual` evaluates both control laws at the smooth
model's solved operating point. It measures controller-algebra smoothing error
at that point, rather than the difference between two network equilibria.

A separate calculation is available through
[`solve_inverter_control_network_fixed_point`](@ref) and
[`solve_controlled_inverter_fleet_network_fixed_point`](@ref). Each iteration
evaluates the exact local law, fixes the resulting converter- or grid-side
current targets, and re-solves the physical plant and network. The fleet routine
updates all selected targets in one network solve, retaining coupling between
inverters. The default initial voltage comes from the smooth controlled solve;
other initial voltages can be supplied.

The routines return voltage and current residuals, convergence and cycling
information, and the status of the final plant solve. The
[closed-loop comparison guide](closed_loop_evidence.md) describes these results,
network sensitivities and comparisons across initial states. Regression tests
in `test/closed_loop_evidence_tests.jl` cover single-inverter and fleet cases,
including curved characteristics, a floating neutral and failed plant solves.

These calculations compare quasi-static equilibria under the specified model.
Convergence of the iteration does not establish dynamic stability or global
uniqueness. Conversely, cycling or a failed plant solve does not establish that
no equilibrium exists. Hardware and dynamic-model validation remain separate
requirements.

## Australian study checks

The Australian LV study (`studies/australian_pv/README.md`) combines four kinds
of check. Numbers are in the generated report
(`studies/australian_pv/engineering_benchmark.md`); earlier campaigns are listed
in `studies/australian_pv/HISTORY.md`.

- **Rated-capability gate.** Before any feeder run, each fixture plant must
  deliver its nameplate at the AS/NZS 4777.2 Australia A volt-var corners with
  apparent power or phase current, not an internal budget, as the limit. The
  earlier 1 A capacitor budget failed this gate (about 1.2 kVA of 6 kVA) and
  made every earlier export comparison capacitor-limited.
- **PF/OPF agreement.** The OPF with each site's local volt-var request and
  volt-watt ceiling imposed at its own voltages, and without the 253 V row, must
  reproduce the local fixed point. The two formulations are built separately.
- **Saved-point and network audits.** Every saved point is rebuilt in SI, and a
  separately coded Python script reconstructs feeder currents and limits.
- **Fixed-POC EMT.** A separately coded averaged/switched circuit checks that
  the modelled capacitor-current bound encloses the switched value at the
  highest-stress site per condition, mode and controller group. The bound was
  revised after earlier misses and is checked with the same selection rule, so
  this is not out-of-sample validation. It omits PLL, controller and feeder
  dynamics, and is not a PSCAD/PLECS or hardware comparison.

The carrier-domain and DC-boundary experiments
(`studies/australian_pv/engineering_domain.md`, `engineering_boundary.md`) check
the carrier approximations across P/Q, sequence and modulation under the legacy
fixture limits; those limits only decide which points are labelled infeasible.

## Further higher-fidelity validation

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

Before comparing controllers or dispatch on a feeder, check that the device
fixture can deliver its declared nameplate (the rated-capability gate above).
An internal budget that binds below nameplate dominates every comparison.

For every claimed boundary point:

1. check the per-unit base against the device ratings. Ratings are stamped as
   squared per-unit inequalities, and Ipopt's `bound_relax_factor` floors its
   relaxation at `max(1, |bound|)`, so the admissible physical violation grows
   as `s_base²`: a 20 kVA rating is honoured to ~0.25 VA at `s_base=1e6` but
   exceeded by ~2.4 kVA (12 %) at `s_base=1e8`. Keep `s_base` within a couple of
   decades of the ratings or pass `bound_relax_factor=0`. The mechanism and the
   durable fix are documented on [`AdvancedInverter`](@ref). No hardware-sizing
   claim should be published from a base that has not been checked this way;
2. require non-negative `switching_margin`, increase `n_samples`, and
   demonstrate convergence;
3. for PWM-enabled points, require non-negative `pwm_reserve_margin` and
   `pwm_modulation_margin`, require every AC `*_reserved` value to cover its
   predicted ripple, then refine both sampling grids and `pwm_ac_harmonics`;
4. resolve from more than one physically sensible initial point;
5. check all reported conductor, sequence, rail, and capacitor residuals;
6. perturb the operating point inside and outside the claimed boundary;
7. state whether each value is RMS, peak amplitude, or Fourier magnitude; and
8. compare at least one representative point with a higher-fidelity oracle.

For an explicit LCL case, also verify that the fundamental frequency is well
separated from every passive modal resonance and that the resonance remains
acceptably damped after controller delay and grid impedance are included. The
reported scalar `filter_resonance_hz` is only an initial screening number.

The maintained sources behind these tests are listed in [IBR references](@ref
ibr-references).

PCC rating regressions in `test/inverter_power_location_tests.jl` compare the
legacy converter circle with net PCC capability through reduced and LCL filters,
including the PCC shunt. They check SI/per-unit behavior, optional current
limits, exact/smooth controller replay at both current targets, and location-aware
operating, fleet and hardware-sizing reports.

### PCC terminal saturation regressions

`inverter_sequence_limit_tests.jl` separately forces known converter-sequence
currents and compares contradictory and feasible ratings at 10 kVA, 1 MVA and
100 MVA numerical bases. This checks the plant's hard bounds independently of
controller backoff. Controller cases cover I₁/I₂ saturation, both targets and
power boundaries, common/power/balance allocation, LCL shunt offsets, equal
positive/negative ratings and fleet margin reporting. Physical Fortescue sums,
DC availability and exact/smooth command agreement are checked independently.

The Australian terminal-profile fleet regression covers mean, worst-phase and
sequence VV+VW at both source levels and zero/half/full PV. It uses a single
100 kVA base, adaptive Ipopt profile and 300-iteration cap without retries,
then checks that a zero selection objective preserves the stressed dispatch.
The network audit replays all three stressed equilibria with fixed phase P/Q.
The larger PCC graph retains squared-norm coefficients as expressions and
normalizes its scalar root by `sqrt(d + epsilon^2)`. The quadratic, smoothing
width and service order are unchanged; the root stays scaled to headroom for
vanishing service. Simple constant-service graphs retain their prior root
coordinates. The smooth maximum is anchored at a positive argument to prevent
cancellation in its denominator on trial points. Convergence at other network
bases remains a separate conditioning question, with failed probes retained.

The Australia A review fixture also checks mean and worst-phase sensing with
the same reactive-priority allocator, converter-arm current limit and PCC
apparent-power circle, under separately binding current and power limits.
Capacitor thermal regressions independently reconstruct weighted RMS stress
from physical 2ω and declared carrier currents, cover both power boundaries,
both current targets and common/power-priority sequence allocation, and compare
the exact and smooth controller commands. These tests protect allocation
semantics; they do not imply that the sensing policies agree under unbalance.
They also do not establish equivalence to an optimisation that retains a hard
Q-curve equality when that Q request itself exceeds the available capability:
the forward reactive-priority allocator clips infeasible Q before allocating P.
The PCC stamp lifts the curve outputs and phase-current service directions
once before reusing them in capability intersections. This preserves the
equations while avoiding repeated copies of the worst-phase conflict graph
when exporting a fleet model to the nonlinear solver.

PCC mode is checked against independently calculated terminal-circle and
arm-current intersections for watt, var and proportional priority. Both current
targets must deliver the same net P/Q with filter compensation; positive and
negative PCC shunts are included once. Exact replay includes the physical plant
measurements and compares the quadratic segment law with its smooth interior
approximation. A plant-backed network replay refreshes those measurements at
each iteration. These checks validate the idealised policy, not manufacturer
firmware or AS/NZS 4777.2 certification.

The review-derived fixture in `test/inverter_pcc_review_tests.jl` specifies a
neutral-connected LCL filter on an ideal balanced four-wire source. It checks
the declared Australia A curves at knots, slopes and the analytic rating/volt-watt
crossover, with DC availability at 1.2, 0.5 and zero times the rating. Independent
softplus formulas separate curve bias from terminal saturation and component
losses. Width-halving tests check the signed corner error, and separate arm
ratings check feasible-Q priority and the absence of an implied grid-arm cap.
These reproducible cases do not establish general unbalanced-grid behaviour,
charging control, dynamic performance or certification.
