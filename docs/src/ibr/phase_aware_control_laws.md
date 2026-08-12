# [Phase-aware local control laws](@id ibr-phase-aware-control-laws)

This note proposes steady-state local Volt-var and Volt-watt extensions for
unbalanced three-phase voltages. The controller may use the local RMS voltage
phasors and their relative angles, but no remote measurements or feeder-wide
optimisation. The immediate target is a three-leg, three-wire converter. A
four-leg converter can use the same framework with the zero-sequence channel
enabled.

The proposal deliberately separates two questions:

1. **What response would the voltage curves prefer?** Existing piecewise-linear
   primitives answer this question.
2. **What current can this converter actually produce?** Topology-aware
   algebraic limiters answer this question without pretending that three
   phase-power commands are independent. A numerical optimiser is useful as an
   offline oracle, but is not required in the deployed controller.

## Current behaviour in BMOPFTools

The present implementation is not an averaged three-leg controller:

- `PN_PER_PHASE` is the default Volt-var/Volt-watt reference for `FOUR_LEG` and
  `SINGLE_PHASE` IBRs;
- the `*_AVERAGED` references are explicit alternatives, and the legacy
  `voltage_aggregation="AVERAGE"` field can also request averaging; and
- Volt-var/Volt-watt is currently disabled for `THREE_LEG`, with a warning and
  fallback to static power boxes.

The built-in `THREE_LEG` injection model uses cyclic line-to-line branch
currents. [`AdvancedInverter`](@ref) instead represents physical phase-conductor
currents and imposes `Ia + Ib + Ic = 0`. The latter is the appropriate feasibility
oracle for developing a controller intended for a physical three-leg bridge.

## The unavoidable three-leg constraints

Let `U0`, `U1`, and `U2` and `I0`, `I1`, and `I2` be the RMS Fortescue voltage and
current phasors. For a three-leg, three-wire converter,

```math
I_0=0,\qquad I_a+I_b+I_c=0.
```

It therefore has four real fundamental-current degrees of freedom: the complex
positive- and negative-sequence currents. It can affect negative-sequence
voltage, but it cannot inject zero-sequence current or directly correct a
zero-sequence voltage. This is a physical boundary, not a control-design
limitation.

Using the sequence convention in `AdvancedInverter`, aggregate mean complex
power and the double-frequency oscillating-power phasor are

```math
S=3(U_1I_1^*+U_2I_2^*),
\qquad
\widetilde S=3(U_1I_2+U_2I_1).
```

These equations expose the principal tradeoff. Negative-sequence current is
the useful actuator for voltage-unbalance attenuation, but it interacts with
positive-sequence voltage to change phase-current peaks and DC-link 2ω power.
A balanced positive-sequence current is not ripple-free under unbalanced
voltage either: with `I2=0`, `S_tilde=3 U2 I1`.

Every candidate controller must ultimately satisfy all of the following, not
only an aggregate apparent-power circle:

- `I0=0` and each phase-conductor current limit;
- net DC active-power availability, including whether net import is permitted;
- converter-side apparent-power and loss limits;
- internal-voltage/filter KVL and the instantaneous switching hull;
- negative-sequence policy limits; and
- DC-bus 2ω voltage and capacitor RMS-current limits.

The controller should use converter-side `U_int` and `I_conv` for the last two
checks. POC powers alone omit the filter voltage drop and loss.

## Candidate laws

### 1. Worst-phase scalar droop

This is the lowest-risk baseline. Keep a balanced positive-sequence output, but
replace the mean-voltage input by direction-aware envelopes:

- Volt-watt uses the highest monitored phase voltage, so one overvoltage phase
  cannot be hidden by two normal phases;
- the injecting side of Volt-var uses the lowest phase voltage;
- the absorbing side of Volt-var uses the highest phase voltage; and
- if low- and high-voltage requests coexist, a declared priority or deadlock
  rule resolves the conflict rather than allowing their average to select the
  wrong direction.

A conservative smooth proxy can use `|U1| ± k|U2|` instead of hard minimum and
maximum operators. This law prevents the most harmful averaged decisions and
needs only firmware, but it does not actively compensate voltage unbalance.

### 2. Virtual-delta line-to-line droop

Measure

```math
U_{ab}=U_a-U_b,\quad U_{bc}=U_b-U_c,\quad U_{ca}=U_c-U_a,
```

apply the existing PWL curves to `|Uxy|/sqrt(3)`, and translate the desired
virtual branch currents into terminal currents with the delta incidence matrix.
The result satisfies `Ia + Ib + Ic = 0` by construction.

This is attractive when the product already senses line-to-line voltages: it
requires no neutral-voltage sensor and does not react to an unobservable and
uncontrollable zero-sequence voltage. It should be implemented in terminal
current space, or with a minimum-norm branch-current gauge. Treating three
cyclic branch powers as independently rated physical legs would introduce a
nonphysical circulating degree of freedom.

### 3. Closest-feasible per-phase droop

Evaluate the familiar curves independently to obtain preferred phase powers
`s_phi^0 = p_phi^0 + j q_phi^0`, then convert them to preferred currents,

```math
i_\phi^0=\left(s_\phi^0/u_\phi\right)^*.
```

Do **not** apply those powers directly. Project the preferred current vector
onto the converter feasible set:

```math
\mathop{\mathrm{minimize}}_{i_a,i_b,i_c}
\sum_\phi w_\phi|i_\phi-i_\phi^0|^2
\quad\text{subject to}\quad
i_a+i_b+i_c=0
```

together with current, net-power, switching, and DC-ripple constraints. With
only the KCL equality and equal weights, the projection has the simple closed
form `i = i0 - mean(i0)`. The constrained version is a small local allocator.

This is the most literal generalisation of per-phase Volt-var/Volt-watt. It can
use phase-asymmetric active power, reactive power, and power circulation to
counter unbalance. The realised phase powers will not in general lie exactly on
all three curves; the projection residual is therefore a required diagnostic,
not an error to hide. Because the fully constrained projection is an online
optimisation, this law is better retained initially as a research benchmark and
offline oracle than as the manufacturer-facing implementation.

### 4. Dual-sequence droop

Use the positive-sequence channel for ordinary aggregate voltage/power control
and reserve the negative-sequence channel for unbalance:

```math
I_1^*=F_1(U_1),\qquad
I_2^*=\operatorname{sat}\!\left(-k_2\widehat H_2^{-1}U_2\right),
```

where `H2` is the locally observed negative-sequence voltage sensitivity to
injected current. A fixed virtual admittance is the simplest implementation.
An online perturb-and-observe estimate can learn the useful angle using the same
local voltage/current sensors and measurement history, without communications.

The sensitivity angle matters. Pure negative-sequence reactive current is not
generally the best voltage actuator on a resistive LV connection, and an
incorrect high gain can worsen voltage or destabilise a weak-grid controller.
This law therefore needs current priority, gain scheduling, anti-windup, and a
dynamic impedance/stability study beyond the steady-state OPF model.

### 5. Ripple-cancelling negative-sequence current

For a three-wire converter, the unconstrained choice

```math
I_2^{ripple}=-\frac{U_2}{U_1}I_1
```

makes `S_tilde=0`. This is valuable as a reference policy and as one endpoint of
a tradeoff sweep. It is not automatically a voltage-unbalance controller: its
current angle is selected to cancel DC power pulsation and may oppose the angle
needed to attenuate `U2`.

A useful blended controller chooses `I2` by minimising

```math
w_v|I_2-I_2^{unbalance}|^2
+w_{dc}|U_1I_2+U_2I_1|^2
+w_i\max_\phi |I_\phi|^2,
```

subject to the exact converter limits. Varying `wv/wdc` produces a transparent
Pareto frontier between voltage unbalance, semiconductor utilisation, and
capacitor stress.

### 6. Sequence droop with algebraic limiters (recommended)

The manufacturer-facing controller can avoid a QP or MPC entirely. At each
phasor update:

1. Compute `U1` and `U2` with the fixed Fortescue transform. If the controller
   already receives phasors, this is only a few complex additions and constant
   multiplications.
2. Evaluate the existing Volt-var/Volt-watt PWL curves for the aggregate
   positive-sequence `P,Q` request. Retain a maximum-phase-voltage guard for
   Volt-watt so that a single high phase cannot be hidden.
3. Evaluate a new PWL **voltage-unbalance droop** on
   `eta=|U2|/max(|U1|, epsilon)`. Its output `kappa(eta)` is a
   negative-sequence admittance gain.
4. Give that admittance a fixed configured angle relative to `U2`,

   ```math
   I_2^v=-\kappa(\eta)e^{-j\phi_2}U_2,
   ```

   where `phi2` represents the expected local negative-sequence impedance
   angle. A product family can expose a small set of LV/resistive,
   mixed, and MV/inductive presets. Online impedance identification can remain
   an optional later feature. The admittance form goes continuously to zero at
   `U2=0`, needs no unit-phasor division, and is especially convenient for a
   smooth network model. A bounded current-magnitude droop can still be obtained
   with the downstream current limiter.
5. Optionally blend toward the ripple-cancelling target using one fixed setting,

   ```math
   I_2^{req}=(1-\lambda)I_2^v+\lambda I_2^{ripple},
   \qquad 0\le\lambda\le1.
   ```

6. Apply the ordinary product priority and saturation logic described below.

For fixed `I1`, the 2ω power limit is especially cheap to enforce. Since

```math
\widetilde S=3U_1(I_2-I_2^{ripple}),
```

`|S_tilde| <= S_tilde_max` is just a disk in the complex `I2` plane:

```math
|I_2-I_2^{ripple}|\le
\frac{\widetilde S_{max}}{3|U_1|}.
```

Clipping `I2_req` to this disk is one magnitude comparison and one scalar
rescaling. No optimisation solver is involved.

The phase currents are reconstructed directly,

```math
I_a=I_1+I_2,\qquad
I_b=\alpha^2I_1+\alpha I_2,\qquad
I_c=\alpha I_1+\alpha^2I_2.
```

If any phase exceeds its limit, a common scale
`gamma=min(1, Imax/max(|Ia|,|Ib|,|Ic|))` is a safe, very cheap fallback. A less
conservative implementation can use the same familiar priority logic as
positive/reactive current limiters: reserve either `I1` or `I2`, then find the
largest admissible scale on the other request with a fixed small number of
bisection steps. This is a one-dimensional limiter, not optimisation-based
control.

Two declared modes are preferable to hidden weights:

- **power priority** preserves the conventional positive-sequence Volt-var/
  Volt-watt request and gives unbalance control the remaining phase-current and
  ripple headroom; and
- **balance priority** reserves a declared negative-sequence current fraction,
  then curtails positive-sequence current if necessary.

In both modes, protection retains absolute priority: phase current, DC voltage,
capacitor current, and modulation margin; then net active-power availability;
then the selected power/balance service priority. The deployed code therefore
consists of PWL evaluation, a Fortescue transform, complex arithmetic, and
scalar saturation—the same broad implementation class as conventional droop
plus current limiting.

## Network-scale smooth formulation

Large-network simulation strengthens, rather than weakens, the case for the
sequence law. The deployed controller remains algebraic, so the OPF should stamp
those equations directly. It must not embed one optimisation problem, KKT
system, complementarity system, or integer mode selector per inverter.

### Smooth sequence measurement

The Fortescue transform is affine in rectangular voltage components and adds no
nonlinearity. Near balanced operation, avoid the singular derivative of `|U2|`
and the division by `|U1|` by defining

```math
\nu_2=\sqrt{|U_2|^2+\epsilon_v^2}-\epsilon_v,
\qquad
\nu_1=\sqrt{|U_1|^2+\epsilon_v^2},
\qquad
\eta=\frac{\nu_2}{\nu_1}.
```

This gives `eta=0` at exact balance, has finite derivatives, and preserves
rotational invariance. `epsilon_v` should be relative to the local voltage base,
not a fixed SI value shared across voltage levels.

The PWL curve `kappa(eta)` can then use BMOPFTools' existing smooth-ReLU/
softplus machinery. In rectangular form, multiplication by
`-kappa*exp(-j*phi2)` is only two smooth algebraic equalities for `I2`. It avoids
angle variables, `atan`, and a normalised `U2/|U2|` direction.

### Smooth worst-phase guards

The three phase-voltage magnitudes can reuse the existing registered magnitude
objects. Replace hard `max` and `min` in the model by pairwise smooth maximum and
minimum operators using the same voltage-relative smoothing policy as the PWL
curves. Because there are only three arguments, this adds constant work per
device. The firmware may retain exact comparisons; the smoothed OPF law should
be tested against it around every transition.

### Smooth protection behaviour

Hard current, apparent-power, switching, and ripple capability inequalities are
already smooth when written as squared norms. However, stamping an unconstrained
controller equality alongside those inequalities makes a stressed equilibrium
infeasible instead of reproducing the controller's saturation. The saturation
policy must therefore be part of the algebraic controller model.

For the simple common-scale fallback, first form unconstrained sequence commands
`I1_hat`, `I2_hat` and their phase currents. Define regularised magnitudes and

```math
M_I=\operatorname{smax}_\epsilon
  (|\widehat I_a|,|\widehat I_b|,|\widehat I_c|),
\qquad
\gamma_I=\operatorname{smin}_\epsilon
  \left(1,\frac{I_{max}}{M_I+\epsilon_i}\right).
```

Because oscillating power is linear in current for fixed voltage, an analogous
factor is

```math
\gamma_{dc}=\operatorname{smin}_\epsilon
\left(1,\frac{\widetilde S_{max}}
{|\widehat{\widetilde S}|_\epsilon}\right).
```

Set `gamma=smin(gamma_I,gamma_dc,...)` and
`(I1,I2)=gamma*(I1_hat,I2_hat)`. A conservative smooth minimum should remain no
larger than either argument. The exact squared-norm capability inequalities stay
in the model as backstops and diagnostics.

Power-priority and balance-priority limiters can be added later as smooth radial
clips of only the lower-priority sequence command. The common-scale law is a
better first large-network implementation because it has fixed structure,
bounded derivatives, no active-set logic, and an obvious firmware counterpart.

### Sparsity and numerical policy

Each controller touches only its own terminal voltage and current variables.
The Fortescue expressions, norms, two PWL evaluations, and limiter add `O(1)`
work per IBR; different devices still couple only through the existing bus KCL.
The Jacobian and Hessian therefore retain the network's block sparsity as device
count grows.

For scalable Ipopt studies:

- work in rectangular coordinates and per unit;
- cache and register sequence magnitudes instead of rebuilding repeated
  expression trees;
- use one documented smoothing scale for curve corners and separate,
  voltage-relative regularisation for complex norms;
- bound droop slopes and gains so that steep controls do not create a badly
  conditioned equilibrium problem;
- keep the model structure fixed when coefficients are parameterised;
- warm-start across time points and use continuation in control gain or
  smoothing when a feeder is difficult; and
- report the exact-versus-smooth controller residual as well as ordinary KCL
  and device-limit residuals.

The additional nonlinearities remain nonconvex, so Ipopt still finds a local
equilibrium/optimum. Large-scale validation should include multiple starts or
continuation for difficult cases and explicit scaling sweeps in both device
count and controller gain.

## DC capacitor implications

For a monolithic link with a DC source that contributes negligibly at 2ω, the
small-ripple relation already implemented by `AdvancedInverter` gives

```math
|D|=\frac{|\widetilde S|}{2\omega C_{dc}V_{dc}},\qquad
C_{dc,min}=\frac{|\widetilde S|}{2\omega V_{dc}\Delta V_{2,max}},
```

where `|D|` and `Delta V2,max` are peak sinusoidal amplitudes. The associated
capacitor-current RMS component is

```math
I_{2\omega,rms}=\frac{|\widetilde S|}{\sqrt2 V_{dc}}.
```

Consequently:

- a negative-sequence voltage controller can require a larger capacitor even
  though a three-leg converter carries no neutral current;
- a ripple-aware negative-sequence angle can reduce the required 2ω capacitance
  and RMS-current rating, at the cost of less voltage-unbalance attenuation or
  more phase current;
- capacitance and ripple-current rating are separate design constraints;
- the calculation must be combined with switching ripple, ESR versus frequency,
  DC-source impedance, hold-up energy, lifetime, and transient-control needs;
  and
- a finite-bandwidth DC source or DC/DC stage can share 2ω current, so the
  open-source formula is a conservative topology-specific limit, not a universal
  sizing rule.

For a four-leg converter, zero-sequence current flows in the fourth leg and does
not pass through a monolithic DC capacitor as fundamental neutral current. For a
split-link three-leg converter it does pass through the half-banks, creates
midpoint motion, and shares their thermal budget with 2ω and switching current.
Those are materially different capacitor-sizing problems.

## Bill-of-materials consequences

The recommended three-leg laws need no fourth leg, neutral conductor, split DC
link, or extra magnetic component. They need:

- positive/negative-sequence extraction in firmware, normally using the DSP and
  PLL-class processing already present;
- the existing phase-current and DC-voltage feedback; and
- either three phase-to-neutral voltage phasors or a line-to-line-only mode.

The line-to-line mode has the strongest no-new-sensor case. If phase-to-neutral
voltage is not already measured, adding a neutral-referenced sensing channel can
be a real BOM, isolation, and certification change. It still does not give a
three-leg bridge authority over zero-sequence voltage.

Larger capacitance is not intrinsically required by phase-aware control. It is
required only when the chosen negative-sequence policy leaves more `|S_tilde|`
than the existing link and ripple specification can tolerate. The algebraic
ripple-disk limiter can instead derate or redirect negative-sequence support
before hardware limits are crossed. A separate scalar thermal limiter can turn
`i_cap_max` into the corresponding allowable `S_tilde_max` after reserving the
measured or modelled switching-current component.

## Recommended development sequence

1. Add **worst-phase scalar droop** as the safe behavioural baseline.
2. Implement the **sequence droop with algebraic limiters**, initially with
   a fixed `phi2`, `lambda=0`, admittance-form `kappa(eta)`, and the smooth
   common-scale fallback.
3. Add the **ripple-disk clip** and balance-priority mode; both reuse the same
   sequence-current reference generator.
4. Retain the **closest-feasible per-phase projection** as an offline optimum
   and regression oracle. It quantifies how much performance the simple law
   leaves on the table without putting an optimiser in firmware.
5. Add **virtual-delta** as the sensor-minimal comparison for products that do
   not measure phase-to-neutral voltage.
6. Sweep the unbalance/ripple priority to publish capability frontiers for
   phase voltage extrema, `|U2|/|U1|`, phase-current peak, `|S_tilde|`, `dv2`,
   capacitor RMS current, converter loss, and switching margin.
7. Validate representative points in an averaged and switched EMT model before
   making dynamic stability or hardware-sizing claims.
8. Add network-size regressions that compare exact firmware evaluation with the
   smooth OPF law and record variable/constraint growth, iteration count,
   factorisation time, and memory as the number of controlled IBRs increases.

The first implementation should not simply enable the existing four-leg
per-phase equalities for `THREE_LEG`. That would confuse preferred curve values
with independently realisable phase powers and would bypass the topology, DC,
and modulation constraints that motivate this work.
