# Generalized generators

A **generalized generator** exchanges currents through an oriented collection
of oriented two-terminal ports at one bus. Its internal voltage source is connected
to those ports through a linear impedance. Source-voltage laws, power controls,
and capability limits determine its operating region. Without an explicit power
control, it is not necessarily a PQ injection.

This is a [draft extension](index.md), distinct from the upstream `generator`
record. Parts 1–5 define its physical data and equations. The shared
[operating-model specification](operating-model.md) forms part of those definitions.

## 1. Data model

An entry of `generalized_generator`, keyed by string ID ``g``:

| Field | Type | Unit | Req. | Description |
|-------|------|:----:|:----:|-------------|
| `bus` | string | – | ✔ | Host bus ID ``i`` |
| `terminal_map` | string[] | – | ✔ | Distinct external conductor terminals, in order |
| `configuration` | string | – | ✔ | `WYE`, `SINGLE_PHASE`, `DELTA`, or `PORTS`, as defined below |
| `port_map` | integer[][] | – | for `PORTS` | Ordered pairs of 1-based conductor indices: outgoing, return |
| `R_series_p_q`, `X_series_p_q` | number | Ω | | Port-impedance entries; absent entries are zero |
| `voltage_model` and conditional voltage data | see shared table | V, rad | ✔ / conditional | Internal source law a–e |
| Power setpoints and limits, voltage/current limits | see shared table | SI | | Optional operating controls and capability |
| `cost_total` | number | currency/kWh | ✔ | Linear price of complete PCC active-power injection; explicit zero is permitted |

Let ``m`` be the length of `terminal_map` and ``r`` the number of ports:

| Configuration | Dimensions | Derived port pairs ``(f_k,t_k)`` |
|:--|:--|:--|
| `WYE` | ``m\ge2``, ``r=m-1`` | ``(1,m),\ldots,(m-1,m)``; the final conductor is the shared return |
| `SINGLE_PHASE` | ``m=2``, ``r=1`` | ``(1,2)``; either line-neutral or line-line, according to the actual terminals |
| `DELTA` | ``m=r=3`` | ``(1,2),(2,3),(3,1)``; winding order ab, bc, ca |
| `PORTS` | ``m\ge2``, ``r\ge1`` | Each row of `port_map` gives one ordered pair |

`port_map` is forbidden for the three derived configurations. For `PORTS`, both
indices lie in ``1,\ldots,m``, are unequal, and every declared conductor participates.
The oriented incidence matrix has column rank ``r`` for a forest. An ordered
closed triangle is also supported, with rank two and three winding currents.
Duplicate windings and other cycles are rejected. A `PORTS` triangle must have
the same cyclic orientation and ordering as `DELTA`.

A split-phase source uses `WYE` with two legs and a return, for example
`["L1","L2","N"]`. Its 180-degree displacement belongs to the **voltage law**;
it does not follow merely from a string or connection count.

The impedance is an ``r\times r`` **port constitutive matrix**, not an
``m\times m`` external-conductor matrix. It can have mutual coupling. To model
an explicit internal neutral lead, use the conductor primitive in
[Source generators](source.md), with open grounding if there is no earth path.
Do not Kron-reduce the network's neutral terminals or lines to supply this field.
Existing reduced internal-return equivalents need explicit provenance and upstream
review; their parameters do not describe separate neutral-lead current or loss.

## 2. Input symbols

| Data / construction | Symbol | Meaning |
|-------|:------:|-------|
| `terminal_map` | ``\textcolor{purple}{\mathbf{N}_g}`` | Ordered host-bus terminal labels |
| Derived port map | ``\textcolor{red}{\mathbf{C}_g}\in\mathbb R^{m\times r}`` | Oriented conductor–port incidence matrix |
| `R_series`, `X_series` | ``\textcolor{brown}{\mathbf{Z}_g}=\textcolor{red}{\mathbf{R}_g}+\textcolor{brown}{j}\textcolor{red}{\mathbf{X}_g}`` | Port series impedance |
| `cost_total` | ``\textcolor{red}{c_g}`` | Price of complete PCC active power |

The imaginary unit satisfies ``\textcolor{brown}{j}^2=-1``. Complex unknowns are
blue, real input parameters red, complex parameters brown, and terminal strings
purple; real unknowns are black. Colors may be dropped in longer derivations
without changing a symbol's type.

Construct the incidence matrix from the ordered pairs:

```math
[C_g]_{pk}=\begin{cases}
+1,&p=f_k,\\
-1,&p=t_k,\\
0,&\text{otherwise},
\end{cases}
\qquad \mathbf1_m^{\mathrm T}C_g=\mathbf0_r^{\mathrm T},\qquad
\operatorname{rank}(C_g)=\begin{cases}2,&\mathrm{DELTA},\\r,&\mathrm{forest}.\end{cases}
```

Physical impedance data must be passive at the represented frequency:

```math
H_g=\frac{Z_g+Z_g^{\mathrm H}}{2}\succeq0.
```

Passivity does not imply that ``Z_g`` is invertible. Reciprocity is not inferred
from missing data. A virtual control impedance is not automatically a physical
impedance with a dissipation ledger.

## 3. Variables

The bus supplies complex conductor voltages to common ground. Port currents are
the independent component unknowns:

```math
\textcolor{blue}{\mathbf{U}_g}
=\textcolor{blue}{\mathbf{U}_i}[\textcolor{purple}{\mathbf{N}_g}]
\in\mathbb C^m,\qquad
\textcolor{blue}{\mathbf{I}_g}\in\mathbb C^r.
```

``I_{g,k}`` flows outward through conductor ``f_k`` and returns through ``t_k``.
The conductor injection ``\textcolor{blue}{\mathbf{J}_g}\in\mathbb C^m``, port
voltage ``\textcolor{blue}{\Delta\mathbf{U}_g}\in\mathbb C^r``, and internal
EMF ``\textcolor{blue}{\mathbf{E}_g}\in\mathbb C^r`` are derived below. They need
not be additional optimization variables. Every phasor is RMS.

## 4. Equality constraints

### Port voltages, terminal injections, and network KCL

```math
\textcolor{blue}{\Delta\mathbf{U}_g}
=\textcolor{red}{\mathbf{C}_g}^{\mathrm T}\textcolor{blue}{\mathbf{U}_g},\qquad
\textcolor{blue}{\mathbf{J}_g}
=\textcolor{red}{\mathbf{C}_g}\textcolor{blue}{\mathbf{I}_g},\qquad
\mathbf1_m^{\mathrm T}\textcolor{blue}{\mathbf{J}_g}=0.
```

For the signed sum ``\kappa_{i,t}^{\mathrm{other}}`` of all other currents
**injected into** host terminal ``t``, the complete network balance is

```math
\kappa_{i,t}^{\mathrm{other}}
+\sum_{p:\,N_{g,p}=t}J_{g,p}=0.
```

Sum the contributions of all devices at a terminal. If the network residual
instead counts withdrawals positively, subtract ``J_g``. No additional return or
earth current is appended to this component's KCL contribution.
For a wye connection this specializes to

```math
J_{g,k}=I_{g,k}\ (k=1,\ldots,r),\qquad
J_{g,n}=-\sum_{k=1}^{r}I_{g,k}.
```

### Series drop and source law

```math
\textcolor{blue}{\mathbf{E}_g}
=\textcolor{blue}{\Delta\mathbf{U}_g}
+\textcolor{brown}{\mathbf{Z}_g}\textcolor{blue}{\mathbf{I}_g}.
```

Apply exactly the selected shared voltage law to ``E_g``, plus every declared
power control and terminal-PV target from [Operating model](operating-model.md).
This is an impedance equation, including when ``Z_g=0`` or is singular; an
admittance inverse is not part of the definition.

### Powers and losses

The conjugate ``*`` acts elementwise, ``\circ`` denotes the Hadamard product,
and ``\mathrm H`` denotes conjugate transpose:

```math
\mathbf S_g^{\mathrm{port}}=\Delta\mathbf U_g\circ\mathbf I_g^*,\qquad
\mathbf S_g^{\mathrm{emf}}=\mathbf E_g\circ\mathbf I_g^*,\qquad
\mathbf S_g^{\mathrm{terminal}}=\mathbf U_g\circ\mathbf J_g^*.
```

```math
S_g^{\mathrm{poc}}=\mathbf1_m^{\mathrm T}\mathbf S_g^{\mathrm{terminal}}
=\mathbf1_r^{\mathrm T}\mathbf S_g^{\mathrm{port}},\qquad
S_g^{\mathrm{internal}}=\mathbf1_r^{\mathrm T}\mathbf S_g^{\mathrm{emf}}.
```

For every named complex power, ``P=\mathfrak R(S)`` and ``Q=\mathfrak I(S)``.
The series complex absorption and physical dissipation are

```math
S_g^{\mathrm{series}}=\mathbf I_g^{\mathrm H}Z_g\mathbf I_g,\qquad
S_g^{\mathrm{internal}}-S_g^{\mathrm{poc}}=S_g^{\mathrm{series}},\qquad
P_g^{\mathrm{series}}=\mathbf I_g^{\mathrm H}H_g\mathbf I_g\ge0.
```

The imaginary part of ``S_g^{\mathrm{series}}`` is reactive absorption, which
can have either sign. These identities contain no core, friction, DC-link,
excitation-system, or mechanical loss model. Internal power is not shaft power.

## 5. Inequality constraints

### Cartesian variable bounds

If an explicit port-current rating ``I_{g,k}^{\max,\mathrm{port}}`` is supplied,
the corresponding redundant but useful variable boxes are

```math
-I_{g,k}^{\max,\mathrm{port}}\le I_{g,k}^{\Re}\le I_{g,k}^{\max,\mathrm{port}},\qquad
-I_{g,k}^{\max,\mathrm{port}}\le I_{g,k}^{\Im}\le I_{g,k}^{\max,\mathrm{port}}.
```

A conductor rating also implies boxes on ``J_{g,p}^{\Re},J_{g,p}^{\Im}``.
Those are bounds on derived expressions unless that current is itself an
independent variable. Boxes do not replace the current-magnitude circle.

### Engineering bounds

Apply all declared shared bounds, including

```math
\mathbf J_g\circ\mathbf J_g^*
\le \mathbf I_g^{\max,\mathrm{terminal}}\circ\mathbf I_g^{\max,\mathrm{terminal}},\qquad
\mathbf I_g\circ\mathbf I_g^*
\le \mathbf I_g^{\max,\mathrm{port}}\circ\mathbf I_g^{\max,\mathrm{port}}.
```

Their JSON fields are `i_max` and `i_port_max`, respectively. In particular, a
wye neutral rating constrains

```math
\left|\sum_{k=1}^r I_{g,k}\right|^2\le (I_{g,n}^{\max,\mathrm{terminal}})^2.
```

The shared specification gives the P/Q/S, PCC-voltage, EMF-voltage, and sequence
bounds, including the exact-zero cases. `ig_max` is forbidden for this component;
it has no earth port.

## 6. Implementation and ideal limits

The independent component variables can be only the ``2r`` real/imaginary current
components. Reconstruct ``J_g``, ``\Delta U_g``, and ``E_g`` as affine expressions.
When ``Z_g=0``, ``E_g\equiv\Delta U_g``: there are no internal voltage nodes,
voltage variables, or duplicate voltage equalities. The same elimination is exact
for a nonzero or singular matrix. The [implementation equations](integration.md)
give the full rectangular form and degree-two realization of the bounds.

The upstream wye dispatch model is the special case `WYE`, zero impedance,
`voltage_model="NONE"`, no power setpoint, and per-port capability bounds at
`"POC"`. Equal lower/upper powers give the same electrical fixed-PQ equations,
although upstream prose presently directs fixed injections to negative loads.
Conductor current ratings retain their upstream meaning. Cost equivalence needs
the conditions in [Migration](integration.md).

## 7. Data example

A three-port internal common-magnitude source with positive-sequence terminal PV
regulation and unequal conductor ratings:

```json
{
  "generalized_generator": {
    "machine": {
      "bus": "poc",
      "terminal_map": ["a", "b", "c", "n"],
      "configuration": "WYE",
      "R_series_1_1": 0.3, "X_series_1_1": 0.6,
      "R_series_2_2": 0.3, "X_series_2_2": 0.6,
      "R_series_3_3": 0.3, "X_series_3_3": 0.6,
      "voltage_model": "COMMON_MAGNITUDE",
      "angle_offsets": [0.0, -2.0943951023931953, 2.0943951023931953],
      "e_min": [190.0, 190.0, 190.0],
      "e_max": [270.0, 270.0, 270.0],
      "power_setpoint_location": "POC",
      "p_total_set": 1800.0,
      "v_target": 230.0,
      "v_target_measurement": "POSITIVE_SEQUENCE",
      "i_port_max": [35.0, 30.0, 32.0],
      "i_max": [35.0, 30.0, 32.0, 15.0],
      "cost_total": 0.15
    }
  }
}
```

This example specifies no Q controller or automatic current-limit transition.
Whether its hard voltage target is feasible depends on the assembled unbalanced
network. The same JSON conventions describe a line-line port with
`SINGLE_PHASE` and `terminal_map=["a","b"]`.

## 7. Closed delta and circulating current

For terminal order a,b,c and winding order ab,bc,ca,

```math
C=\begin{bmatrix}1&0&-1\\-1&1&0\\0&-1&1\end{bmatrix},\quad
\Delta U=\begin{bmatrix}U_a-U_b\\U_b-U_c\\U_c-U_a\end{bmatrix},\quad
J=\begin{bmatrix}I_{ab}-I_{ca}\\I_{bc}-I_{ab}\\I_{ca}-I_{bc}\end{bmatrix}.
```

```math
\mathbf1^{\mathrm T}\Delta U=0,\qquad
\mathbf1^{\mathrm T}J=0,\qquad
\mathbf1^{\mathrm T}E=\mathbf1^{\mathrm T}ZI.
```

These are consequences of the incidence and series equations, not extra KCL/KVL
rows. With ``I=I_\perp+I_{\mathrm{circ}}\mathbf1``, where
``\mathbf1^{\mathrm T}I_\perp=0``,

```math
I_{\mathrm{circ}}=\tfrac13\mathbf1^{\mathrm T}I,\qquad
J=CI_\perp,\qquad
S^{\mathrm{poc}}=\Delta U^{\mathrm T}I^*.
```

Circulation changes winding loading and series loss although it cancels from
line currents. No artificial condition ``I_{\mathrm{circ}}=0`` is imposed.
Controls, impedance, limits, or the objective determine it; an ideal source with
only terminal controls can have nonunique winding currents. Use a physical
circulating-current restriction when that is intended, and report nonuniqueness
otherwise. Delta has no neutral or earth port.

When ``\mathbf1^{\mathrm T}Z=0`` (including ``Z=0``), EMFs must sum to zero.
Fixed templates are validated against this requirement and only two complex
voltage equations are stamped. In a positive-magnitude angle law, the real
nullspace of ``[\Re(a);\Im(a)]``, ``a_k=\exp(j\phi_k)``, determines admissible
magnitude ratios. The implementation requires a unique positive ratio vector;
for 120-degree offsets all three magnitudes are equal. Thus mode e has two
voltage degrees of freedom in an ideal delta, not four. With nonzero loop
impedance, an unequal template can drive circulation and must not be rejected
merely because its EMFs do not sum to zero.

Sequence bounds use winding coordinates. With ``a=\exp(j2\pi/3)``,

```math
\Delta U^0=0,\quad I^0=I_{\mathrm{circ}},\quad
J^0=0,\quad J^1=(1-a)I^1,\quad J^2=(1-a^2)I^2.
```

A positive-sequence voltage target is a **line-line** RMS quantity. For a
balanced system, 400 V line-line corresponds to ``400/\sqrt3`` V phase-neutral;
zero-sequence terminal voltage is not observable through these winding voltages.
`i_port_max` and `i_seq_max` protect winding quantities; `i_max` protects the
external line conductors. A winding current rating is not a line current rating.
