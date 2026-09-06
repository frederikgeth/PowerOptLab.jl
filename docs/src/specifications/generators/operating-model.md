# Generator voltage laws, controls, and capabilities

This page is the shared operating part of the [generalized-generator](generator.md)
and [source-generator](source.md) drafts. It specifies the same field meanings
and equations for either component, using its own definitions of port, conductor,
internal, and complete PCC power. Parts 1–5 are physical model definitions;
solver-specific realization is on the [implementation page](integration.md).

Use ``d`` for either device ID, ``r`` for its port count, and ``m`` for its external
conductor count. A field is present only when its corresponding assumption is
intended. Optional capabilities are inequalities, not implicit controllers.

## 1. Data model

The following are **flat fields on the component record**. They are not nested
Julia `voltage`, `control`, or `capability` objects. Symbols and dimension rules
are defined below. ✔ means required; other conditions are stated explicitly.

### Voltage-source law and voltage measurements

| Field | Type | Unit | Req. | Description |
|-------|------|:----:|:----:|-------------|
| `voltage_model` | string | – | ✔ | `NONE`, `FIXED_PHASOR`, `FIXED_MAGNITUDES`, `COMMON_MAGNITUDE`, or `PHASE_MAGNITUDES` |
| `v_magnitude` | number[] | V | b/c | Length ``r``; strictly positive internal EMF-template magnitudes |
| `v_angle` | number[] | rad | b/c | Length ``r``; internal template angles |
| `angle_offsets` | number[] | rad | d/e | Length ``r``; prescribed relative-angle pattern |
| `e_min`, `e_max` | number[] | V | | Length ``r``; nonnegative internal EMF magnitude limits; e requires strictly positive `e_min` |
| `v_min`, `v_max` | number[] | V | | Length ``r``; nonnegative **PCC port-voltage** magnitude limits, not phase-to-ground limits |
| `v_target` | number | V | | Strictly positive hard terminal-voltage target |
| `v_target_measurement` | string | – | with `v_target` | `POSITIVE_SEQUENCE` for three common-return or delta ports, or `PORT` for exactly one port |

Here b–e denote the laws in Part 4. `v_magnitude`/`v_angle` are forbidden in the
other modes; `angle_offsets` is forbidden in a–c. All required vectors have the
full port length. In c, adding the same constant to every template angle leaves
the feasible set unchanged; in b it changes the absolute source reference.
An explicit angle pattern is required for d/e even for one or three ports.
Examples are ``[0]``, ``[0,\pi]``, and ``[0,-2\pi/3,2\pi/3]``.

`e_min` and `e_max` retain per-port vector syntax in mode d, although every
magnitude is the same. Their intervals must intersect. In modes b/c, bounds
are consistency checks on the known template magnitudes. In mode a, optional
EMF bounds constrain the reconstructed internal voltage without adding a
voltage-source equality.

### Power setpoints and capability

| Field | Type | Unit | Req. | Description |
|-------|------|:----:|:----:|-------------|
| `power_setpoint_location` | string | – | with P/Q setpoint | `POC` or `INTERNAL` |
| `p_set`, `q_set` | number[] | W, var | | Length ``r``; per-port fixed active/reactive powers |
| `p_total_set`, `q_total_set` | number | W, var | | Fixed complete aggregate active/reactive powers |
| `power_limit_location` | string | – | with P/Q/S limit | `POC` or `INTERNAL`, independent of the setpoint location |
| `p_min`, `p_max` | number[] | W | | Length ``r``; per-port active-power bounds |
| `q_min`, `q_max` | number[] | var | | Length ``r``; per-port reactive-power bounds |
| `p_total_min`, `p_total_max` | number | W | | Complete aggregate active-power bounds |
| `q_total_min`, `q_total_max` | number | var | | Complete aggregate reactive-power bounds |
| `s_max` | number[] | VA | | Length ``r``; nonnegative per-port apparent-power ratings |
| `s_total_max` | number | VA | | Nonnegative magnitude limit on complete aggregate complex power |
| `cost_total` | number | currency/kWh | ✔ | Linear price on complete **PCC** active power; not affected by either location field |

For each setpoint select at most one representation: `p_set` or `p_total_set`,
and independently `q_set` or `q_total_set`. The draft follows the current runtime
capability by also allowing at most one member of each pair
`p_min`/`p_total_min`, `p_max`/`p_total_max`, `q_min`/`q_total_min`,
`q_max`/`q_total_max`, and `s_max`/`s_total_max`. Lower and upper sides can use
different representations, for example a total minimum and per-port maxima.
This is an explicit representational restriction, not a physical claim that
aggregate and phase ratings cannot coexist. Supporting both members of a pair
would need a corresponding API/schema extension.

P and Q may be positive or negative. Compare a lower and upper bound directly
only when they refer to the same quantity. A total P lower bound is not compared
entrywise against a vector of phase upper bounds. Any actual combination may
still be infeasible when assembled with the device and network equations.

### Current and sequence capability

| Field | Type | Unit | Req. | Description |
|-------|------|:----:|:----:|-------------|
| `i_max` | number[] | A | | Length ``m``; nonnegative external conductor ratings, including the neutral |
| `i_port_max` | number[] | A | | Length ``r``; nonnegative port-current ratings |
| `ig_max` | number | A | source only | Nonnegative star-to-earth current rating |
| `v_sequence_location` | string | – | with voltage sequence bound | `POC` or `INTERNAL` |
| `v_seq_min`, `v_seq_max` | number[] | V | | Length 3; nonnegative zero/positive/negative voltage-magnitude limits |
| `i_seq_min`, `i_seq_max` | number[] | A | | Length 3; nonnegative zero/positive/negative **port-current** magnitude limits |

Sequence fields require exactly three ordered ports with a common return, or
closed-delta windings ordered ab,bc,ca. Delta voltage/current sequences refer to
line-line winding voltages and winding currents; see the delta equations in
[Generalized generators](generator.md).
The order of the ports defines a,b,c; labels are not permuted automatically.
All three entries of a supplied sequence-bound vector are required, ordered
**zero, positive, negative**. A current lower bound is an optional service
obligation, not a default equipment capability. This series-only model has no
shunts, so internal phase and PCC port currents are the same; the current
sequence requires no separate measurement-location field.

A component-specific `v_min` concerns a port difference. It must not inherit the
phase-to-ground meaning of a bus's `v_min`. The input location is part of each
field's documented semantics, not inferred from similar names elsewhere.

## 2. Input symbols

| Field / definition | Symbol | Meaning |
|-------|:------:|-------|
| `v_magnitude`, `v_angle` | ``\textcolor{brown}{\overline{\mathbf E}_d}`` | Fixed complex internal template |
| `angle_offsets` | ``\textcolor{red}{\boldsymbol\phi_d}`` | Real relative-angle pattern |
| Template ratio / unit rotation | ``\textcolor{brown}{\beta_{d,k}}``, ``\textcolor{brown}{\gamma_{d,k}}`` | Defined below, not extra inputs |
| EMF and PCC voltage bounds | ``\textcolor{red}{\mathbf E_d^{\min/\max}}``, ``\textcolor{red}{\Delta\mathbf U_d^{\min/\max}}`` | Per-port magnitude limits |
| P/Q setpoints | ``\textcolor{red}{\widehat{\mathbf P}_d},\textcolor{red}{\widehat{\mathbf Q}_d}`` or ``\textcolor{red}{\widehat P_d},\textcolor{red}{\widehat Q_d}`` | Per-port or complete total targets |
| P/Q/S bounds | ``\textcolor{red}{P^{\min/\max}},\textcolor{red}{Q^{\min/\max}},\textcolor{red}{S^{\max}}`` | Each indexed quantity at its declared location |
| `v_target` | ``\textcolor{red}{V_d^\star}`` | Selected terminal-voltage magnitude target |
| Sequence bounds | ``\textcolor{red}{V_d^{\sigma,\min/\max}},\textcolor{red}{I_d^{\sigma,\min/\max}}`` | ``\sigma=0,1,2`` |

For b/c construct the template and its ratios; for d/e construct rotations:

```math
\overline E_{d,k}=v_{d,k}^{\mathrm{magnitude}}
\exp(jv_{d,k}^{\mathrm{angle}}),\qquad
\beta_{d,k}=\frac{\overline E_{d,k}}{\overline E_{d,1}},\qquad
\gamma_{d,k}=\exp\big(j(\phi_{d,k}-\phi_{d,1})\big).
```

The bar on ``\overline E`` denotes prescribed data, **not conjugation**; conjugation
uses ``*``. Strictly positive template magnitudes make the ratios well defined.

## 3. Variables and derived measurements

Each component defines complex ``\mathbf E_d``, ``\Delta\mathbf U_d``,
``\mathbf I_d``, ``\mathbf J_d``, ``\mathbf S_d^{port}``, ``\mathbf S_d^{emf}``,
``S_d^{poc}``, and ``S_d^{internal}``. No extra physical state is introduced here.

For either power location ``\ell\in\{\mathrm{POC},\mathrm{INTERNAL}\}``, define

```math
(\mathbf S_d^\ell,S_d^\ell)=
\begin{cases}
(\mathbf S_d^{\mathrm{port}},S_d^{\mathrm{poc}}),&\ell=\mathrm{POC},\\
(\mathbf S_d^{\mathrm{emf}},S_d^{\mathrm{internal}}),&\ell=\mathrm{INTERNAL},
\end{cases}.
```

```math
\begin{aligned}
\mathbf P_d^\ell&=\mathfrak R(\mathbf S_d^\ell),&\quad
\mathbf Q_d^\ell&=\mathfrak I(\mathbf S_d^\ell),\\
P_d^\ell&=\mathfrak R(S_d^\ell),&\quad Q_d^\ell&=\mathfrak I(S_d^\ell).
\end{aligned}
```

For a grounded source at `POC`, ``S_d^\ell`` is generally **not** the sum of
``\mathbf S_d^\ell``. The source page gives the exact correction.

## 4. Equality constraints

### a. No voltage-source equality — `NONE`

There is no source-shape equality on ``E_d``. The circuit reconstructs it from
bus voltages and currents. P/Q equality setpoints create a PQ operating law;
capability bounds alone leave dispatch and potentially phase allocation free.

### b. Full phasor — `FIXED_PHASOR`

```math
\textcolor{blue}{\mathbf E_d}=\textcolor{brown}{\overline{\mathbf E}_d}.
```

All ``2r`` real components are fixed internally. This supplies an absolute
orientation for the source, not automatically a fixed terminal voltage behind
nonzero impedance. Network P/Q is free unless separately constrained.

### c. Fixed magnitudes and relative angles — `FIXED_MAGNITUDES`

Use the canonical anchored equations

```math
E_{d,k}=\beta_{d,k}E_{d,1}\quad(k=2,\ldots,r),\qquad
E_{d,1}E_{d,1}^*=|\overline E_{d,1}|^2.
```

They are equivalent to ``E_d=\lambda\overline E_d`` with ``|\lambda|=1``:
only a common rotation remains. There are ``2(r-1)+1`` real source-law equations.
No separate global angle variable is required.

With nonzero impedance, this fixes **internal excitation**, not PCC magnitude.
Even at zero impedance it prescribes all terminal magnitudes and relative angles,
which is a stronger unbalanced assumption than a positive-sequence PV target.

### d. One common free magnitude — `COMMON_MAGNITUDE`

```math
E_{d,k}=\gamma_{d,k}E_{d,1}\quad(k=2,\ldots,r).
```

There are ``2(r-1)`` real source-law equations, leaving common magnitude and
rotation. Their conceptual form is ``E_{d,k}=\rho\exp(j(\theta+\phi_{d,k}))``
with ``\rho\ge0``. The canonical equations do not add ``\rho`` or ``\theta``
as optimization variables. If zero EMF is permitted by the bounds, its angle
has no physical meaning; the affine equations still have a defined zero limit.

### e. Separate phase/port magnitudes — `PHASE_MAGNITUDES`

Rotate each nonanchor EMF onto the anchor's desired ray:

```math
W_{d,k}=\gamma_{d,k}^*E_{d,k},\qquad
\mathfrak I(W_{d,k}E_{d,1}^*)=0\quad(k=2,\ldots,r).
```

The **orientation inequalities in Part 5 are also required**. Together with
strictly positive EMF magnitude minima, they give
``E_{d,k}=\rho_k\exp(j(\theta+\phi_{d,k}))``, ``\rho_k>0``.
There are ``r-1`` real angle equalities, leaving ``r`` magnitudes and one rotation.
The cross-product equality alone would also admit a phase reversed by ``\pi``.

### Fixed active/reactive power

At the separately declared setpoint location ``\ell_c``, each supplied field adds
its own equality:

```math
\mathbf P_d^{\ell_c}=\widehat{\mathbf P}_d\quad\text{or}\quad
P_d^{\ell_c}=\widehat P_d,\qquad
\mathbf Q_d^{\ell_c}=\widehat{\mathbf Q}_d\quad\text{or}\quad
Q_d^{\ell_c}=\widehat Q_d.
```

An omitted P or Q setpoint adds no equation. An aggregate target does not fix
phase allocation. Component law degrees of freedom do not alone prove that an
assembled network is square, nonsingular, or solvable; repeated references,
fixed powers, and ideal connections can add dependencies or conflicts.

### Sequence definitions and f. Sequence capability

For three common-return ports, let

```math
\alpha=\exp(j2\pi/3),\qquad
F=\frac13\begin{bmatrix}1&1&1\\1&\alpha&\alpha^2\\1&\alpha^2&\alpha\end{bmatrix},\qquad
F^{-1}=\begin{bmatrix}1&1&1\\1&\alpha^2&\alpha\\1&\alpha&\alpha^2\end{bmatrix}.
```

Define the sequences in the order zero, positive, negative:

```math
\begin{bmatrix}\Delta U_d^0\\\Delta U_d^1\\\Delta U_d^2\end{bmatrix}
=F\Delta\mathbf U_d,\qquad
\begin{bmatrix}E_d^0\\E_d^1\\E_d^2\end{bmatrix}=F\mathbf E_d,\qquad
\begin{bmatrix}I_d^0\\I_d^1\\I_d^2\end{bmatrix}=F\mathbf I_d.
```

These are amplitude-invariant Fortescue quantities, not power-invariant coordinates.
For a common-return generalized generator,

```math
I_{d,n}=-3I_d^0,\qquad
S_d^{\mathrm{poc}}=3\sum_{\sigma=0}^2\Delta U_d^\sigma(I_d^\sigma)^*.
```

For a source generator, add the source-page neutral/earth correction to the
last expression and use ``I_{d,n}=-3I_d^0-I_d^g``. Internal EMF power for either
component is ``3\sum_\sigma E_d^\sigma(I_d^\sigma)^*``.
Sequence bounds in Part 5 are an overlay on **any** mode a–e, not a sixth
mutually exclusive voltage law.

A diagonal sequence impedance may be converted to port coordinates by

```math
Z_{abc}=F^{-1}\operatorname{diag}(Z^0,Z^1,Z^2)F,\qquad
Z_{012}=FZ_{abc}F^{-1}.
```

An arbitrary asymmetric phase matrix need not become diagonal. A three-by-three
port transform is not a transform of the source's four-conductor primitive.

### Hard terminal PV target

A terminal target requires `COMMON_MAGNITUDE`, a P setpoint, and no Q setpoint.
Select the measured phasor

```math
V_d^{\mathrm{meas}}=
\begin{cases}
\Delta U_d^1,&\mathrm{POSITIVE\_SEQUENCE},\\
\Delta U_{d,1},&\mathrm{PORT}\ (r=1),
\end{cases}
\qquad
V_d^{\mathrm{meas}}(V_d^{\mathrm{meas}})^*=(V_d^\star)^2.
```

The familiar three-phase machine-inspired preset combines an aggregate P target,
this positive-sequence **terminal** magnitude target, and balanced internal
angle offsets. Common internal magnitude and angle can adjust. There is no
automatic Q-limit release, PV/PQ switching, saturation, droop, or frequency law.
A capability can make the hard target infeasible. Per-port P targets are allowed
by the prototype but can overconstrain this preset; their physical meaning and
feasibility must be assessed for the assembled study.

## 5. Inequality constraints

### Orientation domain for mode e

```math
\mathfrak R(W_{d,k}E_{d,1}^*)\ge0\quad(k=2,\ldots,r),\qquad
0<E_{d,k}^{\min}\le|E_{d,k}|\quad(k=1,\ldots,r).
```

The strict positivity is a condition on the **input minimum**; the solver enforces
the non-strict lower bound ``|E_k|\ge E_k^{min}``. No strict inequality on an
optimization variable is introduced.

### Cartesian variable bounds

For every explicit complex current variable ``z`` with a declared magnitude
rating ``h``, ``-h\le z^{\Re},z^{\Im}\le h`` is implied by its current circle.
It may be supplied as a variable box as specified on each component page. It is
not an independent engineering assumption. No arbitrary finite box is invented
when the physical rating is omitted.

### P/Q and apparent power at the capability location

At ``\ell_b=`` `power_limit_location`, per-port bounds are

```math
P_{d,k}^{\min}\le P_{d,k}^{\ell_b}\le P_{d,k}^{\max},\qquad
Q_{d,k}^{\min}\le Q_{d,k}^{\ell_b}\le Q_{d,k}^{\max},\qquad
(P_{d,k}^{\ell_b})^2+(Q_{d,k}^{\ell_b})^2\le(S_{d,k}^{\max})^2.
```

The corresponding total fields give

```math
P_d^{\min}\le P_d^{\ell_b}\le P_d^{\max},\qquad
Q_d^{\min}\le Q_d^{\ell_b}\le Q_d^{\max},\qquad
(P_d^{\ell_b})^2+(Q_d^{\ell_b})^2\le(S_d^{\max})^2.
```

Use each lower/upper side only when declared. ``|\sum_kS_k|`` is not
``\sum_k|S_k|``: `s_total_max` bounds the magnitude of the complete aggregate
complex power, not a sum of per-port ratings. A zero S rating requires P=Q=0.

### PCC and internal voltage magnitudes

In vector form, with componentwise inequalities,

```math
\Delta\mathbf U_d^{\min}\le|\Delta\mathbf U_d|\le\Delta\mathbf U_d^{\max},\qquad
\mathbf E_d^{\min}\le|\mathbf E_d|\le\mathbf E_d^{\max}.
```

The continuously differentiable forms are

```math
\Delta\mathbf U_d^{\min}\circ\Delta\mathbf U_d^{\min}
\le\Delta\mathbf U_d\circ\Delta\mathbf U_d^*
\le\Delta\mathbf U_d^{\max}\circ\Delta\mathbf U_d^{\max},\qquad
\mathbf E_d^{\min}\circ\mathbf E_d^{\min}
\le\mathbf E_d\circ\mathbf E_d^*
\le\mathbf E_d^{\max}\circ\mathbf E_d^{\max}.
```

### Port, conductor, earth, and sequence magnitudes

The component pages specify the separate port, conductor, and earth circles.
For each sequence ``\sigma=0,1,2``, select ``V_d^\sigma=\Delta U_d^\sigma``
at `POC` or ``V_d^\sigma=E_d^\sigma`` at `INTERNAL`, then impose

```math
(V_d^{\sigma,\min})^2\le V_d^\sigma(V_d^\sigma)^*\le(V_d^{\sigma,\max})^2,\qquad
(I_d^{\sigma,\min})^2\le I_d^\sigma(I_d^\sigma)^*\le(I_d^{\sigma,\max})^2.
```

Every sequence rating is independent of phase, neutral, and earth ratings; none
can replace the others. For example ``I^1=I^2=6`` A and ``I^0=0`` gives
``I_a=12`` A, so 8 A sequence ratings do not imply an 8 A phase rating.

### Exact boundaries and omitted constraints

For any magnitude-bounded phasor ``z`` with nonnegative limits ``l,h``, the
canonical zero-upper-limit realization is

```math
h=0:\quad z^{\Re}=0,\ z^{\Im}=0.
```

A lower limit of zero contributes no constraint. Equal positive limits give
``(z^{\Re})^2+(z^{\Im})^2=h^2`` once. A positive lower limit together with a zero
upper limit is invalid. These rules apply equally to phase, neutral, earth, EMF,
PCC-voltage, and sequence quantities. They avoid the vanishing gradient of a
squared-norm upper bound at zero. Fixed P/Q bounds similarly give one equality.

## 6. Objective contribution

The fields define feasible operation independently of the selected network
objective. Under the draft's linear total-price convention, the new components'
contribution to snapshot cost rate is

```math
\Phi_{\mathrm{new}}=\sum_{d\in\mathcal G_{\mathrm{new}}\cup\mathcal S_{\mathrm{new}}}
 c_d\frac{P_d^{\mathrm{poc}}}{1000}\quad[\mathrm{currency/h}].
```

For independent snapshots of duration ``\Delta t_\tau`` hours,

```math
\Phi_{\mathrm{new,horizon}}=\sum_\tau\Delta t_\tau
\sum_d c_d\frac{P_{d,\tau}^{\mathrm{poc}}}{1000}\quad[\mathrm{currency}].
```

Legacy elements keep their separately defined objective terms. Internal generation
cost, nonlinear fuel curves, distinct import/export prices, and inter-period
machine dynamics are not inferred from `cost_total`.
