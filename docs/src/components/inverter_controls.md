# [Phase-aware inverter-control API](@id inverter-control-api)

PowerOptLab implements the first manufacturer-facing law from the
[phase-aware control design](@ref ibr-phase-aware-control-laws) as a closed-form
controller around the existing [`AdvancedInverter`](@ref) plant. The
voltage-curve law uses only the local three-phase RMS voltage phasor; the final
plant-aware capability backoff additionally uses converter-terminal voltage and
both filter-arm currents, as described below. Neither stage solves an online
optimization problem.

> **Kind:** Component/control composition · **Maturity:** research prototype ·
> **Direction:** forward controlled power flow · **Temporal:** single snapshot

!!! warning "Validity boundary"
    This is a quasi-steady-state, fundamental-frequency equilibrium model with
    ideal instantaneous phasor and sequence extraction. It contains no PLL,
    measurement filtering or delay, ramp-rate dynamics, hysteresis,
    anti-windup, ride-through/cessation state machine, or dynamic-stability
    model. A solved equilibrium is not evidence that the physical closed loop
    is stable or meets a standardised response time.

This initial API deliberately supports three-leg, three-wire inverters only.
The Fortescue transform ensures `I0 = 0`, and a common algebraic scale enforces
the conductor-current limit after positive- and negative-sequence requests are
combined. Four-leg and split-link policies need an explicit zero-sequence and
neutral allocation and are not silently approximated by this implementation.

## Construct a controller

Curve breakpoints retain SI units. Volt-watt values are fractions of available
or rated active power according to `volt_watt_basis`, Volt-var values are
fractions of `q_scale`, and the
negative-sequence curve returns admittance in A/V.

```@example phase_control
using PowerOptLab

volt_watt = PiecewiseLinearLaw(
    [230.0, 240.0, 250.0], [1.0, 1.0, 0.2];
    smoothing_epsilon=0.05)

volt_var = PiecewiseLinearLaw(
    [210.0, 220.0, 240.0, 250.0], [0.3, 0.0, 0.0, -0.3];
    smoothing_epsilon=0.05)

negative_gain = PiecewiseLinearLaw(
    [0.0, 0.01, 0.10], [0.0, 0.0, 0.08];
    smoothing_epsilon=1e-4)

controller = SequenceController(
    WorstPhaseVoltVarWatt(
        volt_watt=volt_watt,
        volt_var=volt_var,
        conflict_policy=:dominant),
    NegativeSequenceAdmittanceDroop(
        negative_gain;
        impedance_angle=0.0,
        ripple_blend=0.5,
        voltage_floor=1.0),
    CommonScaleLimiter(pq_priority=:var);
    current_target=:converter,
    power_voltage_floor=1.0)

typeof(controller)
```

Worst-phase Volt-watt observes the largest phase magnitude. Volt-var evaluates
the low-voltage request at the smallest magnitude and the high-voltage request
at the largest. This prevents an average phase voltage from hiding a phase that
is already in a different curve regime.

Three positive-sequence comparators are available:

| Policy | Volt-var input | Volt-watt input | Purpose |
| --- | --- | --- | --- |
| `AverageVoltageVoltVarWatt` | arithmetic mean of the three magnitudes | arithmetic mean | legacy comparator; can hide a phase in a different regime |
| `PositiveSequenceVoltVarWatt` | ``|U_1|`` | maximum phase magnitude by default | inexpensive sequence-coordinate comparator with a safety guard |
| `WorstPhaseVoltVarWatt` | minimum and maximum envelopes | maximum phase magnitude | direction-aware candidate for unbalanced LV systems |

`PositiveSequenceVoltVarWatt(worst_phase_watt_guard=false)` makes both curves
observe ``|U_1|``. That mode is useful as a scientific comparator, but the guard
should normally remain enabled because positive-sequence voltage can conceal one
high phase.

For every policy, result fields `voltage_min` and `voltage_max` mean the physical
minimum and maximum phase magnitudes: hard extrema in `evaluate_exact`, and the
same pairwise smooth extrema in `evaluate_smooth` and the stamped graph. They do
not change meaning to match a policy's curve input.

!!! note "Which reference the phase magnitudes are measured against"
    The measurement is the POC voltage referred to the plant's declared
    `neutral` terminal, and to the network's ground when the composed inverter
    is three-wire (`neutral=nothing`). A common-mode offset therefore leaves
    `PositiveSequenceVoltVarWatt` Volt-var and the negative-sequence droop
    unchanged — they see only ``U_1`` and ``U_2`` — while moving every
    phase-magnitude comparator, including the recommended
    `WorstPhaseVoltVarWatt`. A three-leg bridge cannot control the ``U_0`` it is
    then reacting to. Report `voltage_sequence[1]` with any study that uses a
    phase-magnitude comparator, and see
    [the design note](@ref ibr-phase-aware-control-laws) for why this is a
    sensor specification rather than a modelling detail.

### Simultaneous low- and high-voltage Volt-var requests

Under unbalance, the minimum phase can request positive reactive injection while
the maximum phase requests negative reactive absorption. A three-leg inverter
cannot independently satisfy those two scalar requests with balanced current, so
the conflict rule is part of the control law rather than a numerical detail:

| `conflict_policy` | Result | Intended interpretation |
| --- | --- | --- |
| `:dominant` | continuously blend using branch-normalised severity | recommended default; avoids cancellation and a discontinuous feedback switch |
| `:net` | add positive low-voltage response and negative high-voltage response | legacy/symmetric comparator; may cancel to zero even while two phases violate opposite limits |
| `:low_voltage` | retain the minimum-phase response | explicit undervoltage/service-continuity priority |
| `:high_voltage` | retain the maximum-phase response | explicit overvoltage/export-compliance priority |

For `:dominant`, let ``q_\ell\geq0`` and ``q_h\leq0`` be the two branch
ordinates and let ``q_+`` and ``q_-`` be their respective maximum magnitudes.
The deployed law and network model use

```math
d={q_\ell\over q_+}+{q_h\over q_-},\qquad
w={1\over2}\left(1+{d\over\sqrt{d^2+\epsilon_c^2}}\right),\qquad
q=wq_\ell+(1-w)q_h.
```

Thus a normalised-severity tie returns zero continuously rather than choosing
one side. At that single tie, `:dominant` and `:net` coincide; away from it,
`:dominant` suppresses cancellation by moving rapidly toward the more severe
normalised branch. `conflict_epsilon` is dimensionless and sets a real firmware
transition width; it is not model-only smoothing. This deliberately rejects the
old winner-take-all rule because a discontinuous closed-loop map may have no
equilibrium near the tie surface. The priority rule never overrides physical limits:
converter-leg current, grid-side current, switching, and DC-link constraints
remain hard plant constraints.

`volt_watt_basis=:available` multiplies the curve ordinate by
`request.p_available`. `:rated` multiplies by `request.p_rated` and then caps at
availability. The second convention is needed for rated-power curves at partial
irradiance; the first is retained as an explicit scientific comparator.

`CommonScaleLimiter(pq_priority=...)` declares the positive-sequence
capability allocation: `:watt` preserves active power, `:var` preserves reactive
power, and `:proportional` preserves the requested power factor. These names are
behavioural, not claims of compliance with a particular jurisdictional profile.
Watt and var priority retain the declared `priority_headroom_fraction` (default
``10^{-3}``) before the corresponding capability-circle axis. This 0.1% margin
keeps the remaining-component capacity well conditioned under sustained
saturation. It is an intentional protection margin used by both firmware and
network laws—not an epsilon inside a magnitude or square root.

The negative-sequence policy evaluates

```math
\eta=\frac{|U_2|}{\sqrt{|U_1|^2+U_{floor}^2}},\qquad
I_2^v=-\kappa(\eta)e^{-j\phi_2}U_2.
```

`ripple_blend=0` uses this voltage-oriented request. `ripple_blend=1` uses the
regularized ripple-cancelling target; intermediate values expose the tradeoff
between voltage-unbalance actuation and twice-frequency DC-link stress.

With the RMS Fortescue convention

```math
U_0=(U_a+U_b+U_c)/3,\quad
U_1=(U_a+aU_b+a^2U_c)/3,\quad
U_2=(U_a+a^2U_b+aU_c)/3,
```

where ``a=e^{j2\pi/3}``, the reconstructed three-leg currents are

```math
(I_a,I_b,I_c)=(I_1+I_2,\ a^2I_1+aI_2,\ aI_1+a^2I_2).
```

Consequently ``I_a+I_b+I_c=0`` identically. The requested positive-sequence
current uses

```math
I_1^{req}=\frac{(P-jQ)U_1}{3(|U_1|^2+U_{floor}^2)}.
```

`power_voltage_floor` is a declared low-voltage regularization, not a
curve-smoothing parameter. It makes the division well-defined without replacing
any physical magnitude by an epsilon-perturbed square root. Before other
limiters, its relative positive-sequence power residual is exactly

```math
\frac{|S^{req}-S_1|}{|S^{req}|}
=\frac{U_{floor}^2}{|U_1|^2+U_{floor}^2}.
```

For the one-volt default this is approximately `0.0019%` at 230 V, `0.0076%` at
115 V, and `0.99%` at 10 V. Thus 1 V is acceptable for ordinary 230 V LV
steady-state studies, but it is not a physical ride-through law: results near
voltage collapse depend materially on the declared floor and should be excluded
or paired with an explicit undervoltage cessation/ride-through policy. The test
suite verifies the closed-form residual at 0.1 V, 1 V, and 10 V floors.

`NegativeSequenceAdmittanceDroop.voltage_floor` is separate. It regularizes the
unbalance index and ripple-cancelling ratio; it does not control the ``P,Q`` to
``I_1`` conversion. For `ripple_blend=1`,

```math
I_2^{ripple}=-\frac{U_2U_1^*}{|U_1|^2+U_{floor}^2}I_1,
```

so the remaining ripple has the closed form

```math
\widetilde S=
3U_2I_1\frac{U_{floor}^2}{|U_1|^2+U_{floor}^2}.
```

This residual, rather than an assertion of exact zero, is the unit-test oracle.

## Exact and smooth evaluators

Use [`evaluate_exact`](@ref) as the firmware oracle:

```@example phase_control
measurement = InverterControlMeasurement(ComplexF64[
    245cis(0.05), 215cis(-2.15), 230cis(2.0)])
request = InverterControlRequest(
    p_available=12e3, p_rated=12e3, q_scale=8e3)
ratings = InverterControlRatings(s_max=20e3, i_max=40.0)

exact = evaluate_exact(controller, measurement, request, ratings)
smooth = evaluate_smooth(controller, measurement, request, ratings)

(p_exact=exact.p_request, q_exact=exact.q_request,
 max_exact_smooth_current_error=
     maximum(abs, exact.phase_current .- smooth.phase_current))
```

The exact path uses exact PWL corners, extrema, and saturation. The smooth path
uses the same fixed computation graph as the JuMP model, including the public
BMOPFTools smooth-PWL primitive. Both paths return [`InverterControlResult`](@ref)
with identical SI semantics. Comparing `phase_current` gives a direct
exact-versus-smooth residual without building an OPF model.

The smoothing and regularisation parameters are:

| Parameter | Units | Default/example | Role and sensitivity test |
| --- | --- | ---: | --- |
| `smoothing_epsilon` | curve-input units | user supplied; 0.05 V in example | PWL corner bias/curvature; halve repeatedly and report observable convergence |
| `extrema_epsilon`, `guard_epsilon` | V | 0.05 | smooth phase min/max; sweep relative to voltage-data precision |
| `conflict_epsilon` | normalised severity | 0.01 | width of deployed dominant transition; sweep the equal-severity manifold |
| `current_epsilon` | A | 0.001 | current-selector width; refine around binding legs |
| `power_epsilon` | VA | 0.001 | power-selector width; refine around the capability circle |
| `priority_headroom_fraction` | 1 | ``10^{-3}`` | declared watt/var-priority axis reserve; verify priority results across the intended DC/AC-ratio range |
| `power_voltage_floor` | V | 1 | regularises ``P,Q\mapsto I_1``; report low-voltage sensitivity |
| `voltage_floor` | V | 1 | regularises unbalance/ripple ratios; report balanced and low-voltage sensitivity |
| ``\kappa`` | A/V | study specific | size from expected negative-sequence impedance and required attenuation |
| ``\phi_2`` | rad | study specific | align with the negative-sequence network impedance; sweep angle uncertainty |
| ``\lambda`` (`ripple_blend`) | 1 | 0–1 | trades voltage actuation against 2ω ripple; publish a Pareto sweep |

The near-exact limiter defaults minimise controller bias but may be too sharp
for a difficult fleet solve. There is no universal setting: use continuation
from wider selectors and retain the refinement evidence for reported cases.

`current_epsilon` and `power_epsilon` are **absolute** SI widths. A 1 mVA
selector width is a different fraction of a 5 kVA inverter than of a 500 kVA
one, and once divided by a megavolt-ampere base it is ``10^{-9}`` per unit —
small enough that the selector is effectively a hard `max` with a curvature of
``1/\epsilon`` at the tie. A heterogeneous fleet therefore does not share one
relative controller bias. Declare these widths explicitly per study rather than
inheriting the defaults across devices of different ratings; making them
rating-relative, as [`AdvancedInverter`](@ref) already does for its own
smoothing, is an open design decision recorded in
[the design note](@ref ibr-phase-aware-control-laws).

The documented Volt-watt knot has first-order smoothing bias. This executable
table makes its magnitude visible rather than treating smoothing as exact:

```@example phase_control
refinement = map((1.0, 0.5, 0.1, 0.05, 0.01)) do epsilon
    law = PiecewiseLinearLaw(
        [230.0, 240.0, 250.0], [1.0, 1.0, 0.2];
        smoothing_epsilon=epsilon)
    c = SequenceController(AverageVoltageVoltVarWatt(volt_watt=law))
    m = InverterControlMeasurement([
        240 + 0im, 240cis(-2pi/3), 240cis(2pi/3)])
    e = evaluate_exact(c, m, request, ratings)
    s = evaluate_smooth(c, m, request, ratings)
    (epsilon_V=epsilon, p_exact_W=e.p_request,
     p_smooth_W=s.p_request,
     relative_error=(e.p_request-s.p_request)/e.p_request)
end
refinement
```

The equal-severity conflict is also an intentional smooth control transition:

```@example phase_control
conflict_only = SequenceController(WorstPhaseVoltVarWatt(
    volt_var=volt_var, conflict_policy=:dominant))
[(vmin_V=v, q_var=evaluate_exact(
    conflict_only,
    InverterControlMeasurement([
        245cis(0.05), v*cis(-2.15), 230cis(2.0)]),
    request, ratings).q_request) for v in (214.8, 215.0, 215.2)]
```

PowerOptLab converts these widths with the model bases inside the supported
per-unit formulation. The controller parameters and extracted results retain
SI semantics.

`solve_controlled_inverter` requires `per_unit=true`. Raw-SI stamping of this
coupled nonlinear controller is deliberately rejected because saturated cases
do not have sufficiently reliable Ipopt convergence for scientific use. This
restriction does not change the SI units of controller parameters or returned
results, and it does not restrict the uncoupled exact and smooth numeric
evaluators. The underlying `AdvancedInverter` retains its separate raw-SI
support.

### Differentiable magnitude representation

The JuMP formulation does not replace a magnitude by
``\sqrt{x+\epsilon^2}``. Every required magnitude introduces a nonnegative
auxiliary variable

```math
y\geq0,\qquad (y/y_b)^2=x/y_b^2,
```

where ``y_b`` is a fixed scaling quantity used only to condition the equality.
Thus ``y=\sqrt{x}`` exactly at every feasible point — to the solver's constraint
tolerance, which is why the choice of ``y_b`` matters. This representation is
used for phase-voltage, sequence-voltage, apparent-power, and phase-current
magnitudes, and for the square roots inside smooth min/max selectors. Epsilon is
retained only where it defines an intended smooth approximation: BMOPFTools PWL
corners and algebraic selection between two values.

There is one deliberate exception. The zero-clamp on remaining capability
headroom inside the plant-aware backoff uses the shifted expression
``(a+\sqrt{a^2+\epsilon^2})/2`` rather than another lifted variable. Lifting it
was measured to push the `s_max` and `dv2_max` saturation regressions from
`LOCALLY_SOLVED` to non-publishable — the same failure mode that moved
[`AdvancedInverter`](@ref)'s current-magnitude loss term to the expression form.
The clamp exceeds ``\max(a,0)`` by at most ``\epsilon/2``, so it relaxes the
conservative triangle bound by less than the selector width already declared.
See the numerical-policy discussion in
[the design note](@ref ibr-phase-aware-control-laws).

Implicit magnitude equalities are degenerate at an exact zero vector: the
equality pins ``y=0`` while its gradient vanishes there, so LICQ fails at the
model's own solution. Structurally zero offsets in the capability allocator —
converter-target apparent power and ``\Delta V_{2,max}`` on any filter without
an explicit LCL midpoint — are therefore detected and stamped as constants
rather than lifted. The implementation likewise avoids an unnecessary `|U2|`
variable for `NoUnbalanceControl`. Studies using negative-sequence droop should
still include balanced cases in their solver-robustness audit and report any
initialization dependence, because a *near*-zero ``U_2`` remains poorly
conditioned even when the exactly-zero case has been removed.

## Compose with the physical inverter

```@example phase_control
using BMOPFTools: parse_bmopf

network = parse_bmopf("""
{"bus":{
  "grid":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"]},
  "poc":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"],
         "v_min":[180,180,180],"v_max":[270,270,270]}},
 "voltage_source":{"vs":{"bus":"grid","terminal_map":["a","b","c"],
   "v_magnitude":[245,215,230],"v_angle":[0.05,-2.15,2.0]}},
 "linecode":{"lc":{"R_series_1_1":0.05,"R_series_2_2":0.05,
   "R_series_3_3":0.05,"R_series_4_4":0.05}},
 "line":{"l1":{"bus_from":"grid","bus_to":"poc",
   "terminal_map_from":["a","b","c","n"],
   "terminal_map_to":["a","b","c","n"],"linecode":"lc","length":1}}}
"""; from_string=true)

inverter = AdvancedInverter(
    id="pv", bus="poc", phase_terminals=["a", "b", "c"], neutral="n",
    topology=:THREE_LEG, s_max=20e3, i_max=40.0,
    v_dc=700.0, c_dc=1.1e-3,
    r_filter=0.05, x_filter=0.15)

device = ControlledDevice(inverter, controller)
result = solve_controlled_inverter(
    network, device, request;
    solver_options=("max_iter" => 500, "tol" => 1e-8))

(status=result.termination_status,
 poc_power=result.plant.p_poc,
 converter_apparent_power=abs(result.converter_terminal.total_power),
 current_scale=result.control.current_scale,
 exact_smooth_current_error=result.exact_smooth_current_residual)
```

The converter-referred capability backoff remains publishable through PV
oversizing. The following table is regenerated by Documenter:

```@example phase_control
saturation_controller = SequenceController(AverageVoltageVoltVarWatt())
map((18e3, 22e3, 30e3)) do available
    r = solve_controlled_inverter(
        network, ControlledDevice(inverter, saturation_controller),
        InverterControlRequest(
            p_available=available, p_rated=available, q_scale=0.0);
        solver_options=("max_iter" => 500, "tol" => 1e-8))
    (p_available_W=available, status=r.termination_status,
     converter_apparent_power_VA=
         round(abs(r.converter_terminal.total_power); digits=2),
     current_scale=round(r.control.current_scale; digits=6))
end
```

`current_target=:converter` constrains the converter-side leg-current handles;
`current_target=:grid` constrains the grid-side current after the output filter.
The former most directly controls semiconductor current. The latter most
directly realizes the current seen by the network and is meaningful for LCL
filters, where shunt current separates the two locations.

The limiter reconstructs all three target currents and applies per-leg ratings,
never a positive-sequence or aggregate-current proxy. Its final plant-aware
backoff checks converter and grid currents separately, including LCL shunt
current, and checks converter-terminal apparent power. If `dv2_max` is declared,
it also backs off against the corresponding 2ω-power allowance
``2\omega C_{dc}V_{dc}\Delta V_{2,max}``. This prevents ordinary saturation
from becoming an infeasible command equality; the exact plant inequalities
remain authoritative backstops.

The plant-aware step uses a conservative triangle-bound scalar allocation. It
does not yet allocate switching-margin, PWM carrier-current, capacitor thermal,
or modulation headroom inside the controller; those hard limits can still make
a snapshot infeasible. Such runs must be retained and classified rather than
dropped from penetration statistics. Positive/negative-sequence service
priority is also future work; physical backoff can curtail both together after
the requested P/Q priority has been applied.

The plant remains responsible for filter KVL, internal voltage, modulation,
apparent power, conductor currents, losses, negative-sequence limits, 2ω ripple,
and DC-capacitor constraints. No parallel current-injection plant is introduced.

For simultaneous network snapshots with many controlled devices, use
[`ControlledInverterFleetSpec`](@ref) and
[`solve_controlled_inverter_fleet`](@ref). The [network-scale study
contract](../ibr/network_control_studies.md) explains native-IBR replacement,
objective semantics, tidy extraction, and verification obligations.

[`ControlledInverterResult`](@ref) contains:

- `plant`: the usual advanced-inverter result tuple;
- `control`: the solved smooth command and local measurements;
- `converter_terminal`: actual converter-terminal phase/sequence phasors and
  powers derived from the plant solution;
- `grid_phase_current`: actual post-filter grid-current phasors;
- `exact_control`: the exact law reevaluated at the solved phasor;
- `exact_smooth_current_residual`: maximum phase-current command difference;
- `bus` and `solve`: the network result and structured solve status.

Result extraction reads the stamped JuMP graph directly. `exact_control`
reevaluates the firmware law at the solved phasor and applies the same
plant-aware backoff.

!!! warning "What the exact/smooth residual does and does not bound"
    Both laws are evaluated at the **smooth model's own solved operating
    point**, and the exact law reuses that solve's converter voltage, converter
    current, and grid current. `exact_smooth_current_residual` therefore bounds
    the smoothing error of the controller algebra — PWL corners, smooth
    selectors, lifted magnitudes — given the plant solution. It does not bound
    the distance between the exact-law network equilibrium and the smooth-law
    network equilibrium, because a firmware controller iterating against the
    real network would settle at a different voltage.

    Closing that gap needs an independent oracle: a Picard/Gauss–Seidel
    iteration that evaluates [`evaluate_exact`](@ref) at the current voltage,
    imposes the resulting current on an otherwise uncontrolled power flow, and
    repeats to a fixed point, then reports ``\|\Delta V\|`` against the smooth
    equilibrium. That oracle is not implemented. Until it is, equilibrium-level
    agreement between the firmware law and the network surrogate is an
    assumption of every study on this page, not a tested property. Local
    Volt-var droop is known to have equilibrium existence, uniqueness, and
    convergence conditions of its own; see Farivar, Chen, and Low and Zhu and
    Liu in [IBR references](@ref ibr-references).

`converter_terminal` is the authoritative power record. It computes

```math
S_{conv}=\sum_\phi U_{conv,\phi}I_{conv,\phi}^*,\qquad
S_k=3U_{conv,k}I_{conv,k}^*.
```

Its zero-, positive-, and negative-sequence powers sum to total converter power.
The `control.sequence_power`, `control.total_power`, and `control.ripple_power`
fields instead describe the commanded current together with the measured POC
voltage. They are useful controller diagnostics, but are not converter powers
when an L or LCL filter introduces voltage drop, loss, or shunt current.

## Verification and validation

The canonical claim-to-test matrix, tolerances, OpenDSS scope, topology
reduction boundary, and publication gate live in
[Verification and benchmark cases](@ref ibr-verification). The controller tests
include dominant-conflict continuity, PWL-width refinement, all P/Q priorities,
rated-versus-available Volt-watt at partial irradiance, binding converter
apparent power, binding DC ripple, full stamped current phasors at converter and
grid targets, loss-versus-zero objective invariance, common-mode
(zero-sequence) invariance of the sequence-referred law, and sign safety of the
capability backoff when a per-leg limit is already exhausted. The P/Q-priority
regression sweeps the exact and smooth local laws over DC/AC ratios from 0.9
through 1.4, including a non-zero var request and the exact 1.0 transition. A
stamped saturated ratio of 1.1 is then solved for every priority in the
supported per-unit formulation. This crosses the capability boundary without
making the lightweight unit suite a broad Ipopt stress test. Raw-SI controller
stamping is outside the supported scope and is covered by a fail-fast API test.

### What OpenDSS does and does not validate

The OpenDSS case is intentionally balanced. Under balance, worst-phase,
average-voltage, and positive-sequence Volt-watt references coincide, while
`I2=0`. It therefore validates PWL convention, units, sign, available-power
scaling, rated-power scaling at partial irradiance, network/controller
fixed-point coupling, and the slope-to-saturation transition against an
independent implementation. It does **not** validate the
new negative-sequence law or three-leg feasibility under unbalance because
OpenDSS `InvControl` does not implement that law.

The BMOPFTools suite already validates its conventional single-phase, per-phase
four-leg, and averaged three-phase droop models against OpenDSS across deadband,
slope, and saturation regimes. PowerOptLab reuses that established oracle
pattern rather than duplicating the full upstream matrix.

Likewise, the four-leg/three-single-phase reduction is valid only for the
balanced fundamental-frequency AC case with identical phase plants and a slack
neutral. Three independent DC links do not reproduce a shared four-leg DC bus,
carrier-current correlation, capacitor heating, or neutral-leg loss. Under
unbalance, a four-leg inverter can exchange zero-sequence current and the simple
reduction is no longer the validation target.

### Publication gate

Before drawing control or hardware conclusions, each selected law must pass:

1. exact/smooth residual sweeps at every PWL knot and on both sides of each
   limiter transition;
2. balanced, pure-negative-sequence, mixed-sequence, and phase-regime-conflict
   cases;
3. global-angle, cyclic-phase, and common-mode metamorphic checks, plus the
   supported-unit boundary check;
4. conductor KCL, sequence-power, filter-power, switching, and capacitor
   residual checks from the physical plant;
5. multiple physically sensible Ipopt starts for claimed boundary points;
6. smoothing-width refinement demonstrating convergence of scientific
   observables; and
7. an averaged or switched time-domain comparison for representative
   current-limited and ripple-limited points.

OpenDSS is a fundamental-frequency network/control oracle, not an EMT or
DC-capacitor oracle. DC-link sizing claims still require the higher-fidelity
validation described in [Verification and benchmark cases](@ref
ibr-verification).

## Study boundary

`InverterControlRequest.p_available` currently represents non-negative PV active-
power availability. Battery operation is outside the scope of the present
control-law programme. Fleet-level replacement of native BMOPFTools IBRs belongs
in the next problem-builder layer; dataset fields should not be added to the
local law.

For hardware comparisons, retain at least converter/grid phase currents,
converter-terminal total and sequence powers, `current_scale`, `power_scale`,
`i_positive`, `i_negative`, `ripple`, `dv2`, `i_cap`, capacitor loss, switching
margin, curtailment, and the exact/smooth residual. Compare controls first with
fixed hardware, then perform outer sweeps over `i_max`, `c_dc`, and capacitor
ripple-current rating.
