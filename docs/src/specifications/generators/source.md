# Source generators

A **source generator** has an internal common star point, phase EMFs, a full
conductor series impedance, and a separately declared connection from the star to
common ground. Its external neutral and earth currents are independent return
paths. The phase EMFs may obey any of the shared source-voltage laws; a fixed
phasor is one choice, not the definition of the component.

This [draft extension](index.md) has its own `source_generator` collection. It
does not silently change the upstream ideal `voltage_source`. The
[shared operating model](operating-model.md) forms part of this component's
specification, using the source-specific power definitions below.

## 1. Data model

An entry of `source_generator`, keyed by string ID ``s``:

| Field | Type | Unit | Req. | Description |
|-------|------|:----:|:----:|-------------|
| `bus` | string | – | ✔ | Host bus ID ``i`` |
| `terminal_map` | string[] | – | ✔ | Outgoing phase/leg terminals, then the explicit neutral terminal |
| `configuration` | string | – | ✔ | `WYE`; one or more phase/leg ports with a shared internal star |
| `R_series_p_q`, `X_series_p_q` | number | Ω | | Full conductor-impedance entries, phases/legs then neutral; absent entries are zero |
| `grounding_model` | string | – | ✔ | `OPEN`, `IMPEDANCE`, or `IDEAL` |
| `r_ground`, `x_ground` | number | Ω | only with `IMPEDANCE` | Star-to-common-ground impedance; absent entries are zero |
| `voltage_model` and conditional voltage data | see shared table | V, rad | ✔ / conditional | Internal phase-to-star EMF law |
| Power setpoints and limits, voltage/current limits | see shared table | SI | | Optional controls and capability |
| `ig_max` | number | A | | Earth-current magnitude limit, distinct from the neutral limit |
| `cost_total` | number | currency/kWh | ✔ | Price of complete PCC conductor active-power injection |

With ``m=|\mathbf N_s|\ge2``, there are ``r=m-1`` phase/leg ports. The final
conductor is the neutral **by position**, regardless of its label. All terminal
names are distinct host-bus terminals; none denotes the common ground directly.
`port_map` is forbidden. Two legs can model split phase when the voltage law
explicitly sets a ``\pi`` relative displacement.

The series matrix is ``m\times m``. Its neutral row and column are retained,
including mutual coupling. It must satisfy the same Hermitian-part passivity
condition as the generalized generator. `r_ground` is nonnegative; `x_ground`
may have either sign. Both grounding parameters are forbidden for `OPEN` and
`IDEAL`, so ignored impedance data cannot conceal an intended connection.
`IMPEDANCE` with both entries zero is exactly the ideal limit, not an open path.

An absent neutral *impedance* is an ideal neutral lead. An absent neutral
*conductor* would be an open neutral and is **not** represented by this data model.
The explicit neutral terminal is required, even if the assembled solution has
zero neutral current. No bus neutral is grounded merely by using this component.

## 2. Input symbols

| Data | Symbol | Meaning |
|-------|:------:|-------|
| `terminal_map` | ``\textcolor{purple}{\mathbf N_s}`` | Ordered external terminals, final entry neutral |
| `R_series`, `X_series` | ``\textcolor{brown}{\mathbf Z_s^c}=\textcolor{red}{\mathbf R_s^c}+\textcolor{brown}{j}\textcolor{red}{\mathbf X_s^c}`` | Full conductor impedance |
| `r_ground`, `x_ground` | ``\textcolor{brown}{z_s^g}=\textcolor{red}{r_s^g}+\textcolor{brown}{j}\textcolor{red}{x_s^g}`` | Star-to-common-ground impedance |
| Common ground | ``U^g=0`` | Network-wide reference potential, not a conductor index |
| `ig_max` | ``\textcolor{red}{I_s^{g,\max}}`` | Earth-current rating |
| `cost_total` | ``\textcolor{red}{c_s}`` | Price of complete PCC active power |

Passivity means

```math
H_s^c=\frac{Z_s^c+(Z_s^c)^{\mathrm H}}2\succeq0,\qquad r_s^g\ge0.
```

## 3. Variables

The bus supplies the external conductor voltages; the component has ``m``
independent complex conductor-current unknowns:

```math
\textcolor{blue}{\mathbf U_s}
=\textcolor{blue}{\mathbf U_i}[\textcolor{purple}{\mathbf N_s}]
=\begin{bmatrix}\mathbf U_{s,\mathcal P}\\U_{s,n}\end{bmatrix},\qquad
\textcolor{blue}{\mathbf J_s}
=\begin{bmatrix}\mathbf I_s\\I_{s,n}\end{bmatrix}\in\mathbb C^m.
```

Here ``\mathbf I_s\in\mathbb C^r`` is the phase/leg current subvector. Every
entry of ``J_s`` is positive **outward from the source into its bus conductor**.
The star voltage ``\textcolor{blue}{U_{s,*}}``, phase-to-star EMFs
``\textcolor{blue}{\mathbf E_s}``, and earth current
``\textcolor{blue}{I_s^g}`` are derived quantities. ``I_s^g`` is positive from the
star toward common earth.

## 4. Equality constraints

### Internal circuit and terminal injection

```math
\textcolor{blue}{\mathbf U_s^{\mathrm{int}}}
=\begin{bmatrix}
\textcolor{blue}{\mathbf E_s}+\mathbf1_r\textcolor{blue}{U_{s,*}}\\
\textcolor{blue}{U_{s,*}}
\end{bmatrix},\qquad
\textcolor{blue}{\mathbf U_s^{\mathrm{int}}}-\textcolor{blue}{\mathbf U_s}
=\textcolor{brown}{Z_s^c}\textcolor{blue}{\mathbf J_s}.
```

Equivalently, without introducing internal voltage variables,

```math
U_{s,*}=U_{s,n}+(Z_s^c\mathbf J_s)_n,\qquad
E_{s,k}=U_{s,k}+(Z_s^c\mathbf J_s)_k-U_{s,*}
\quad(k=1,\ldots,r).
```

Stamp ``J_{s,p}`` once into its host terminal's injection-positive KCL:

```math
\kappa_{i,t}^{\mathrm{other}}
+\sum_{p:\,N_{s,p}=t}J_{s,p}=0.
```

The phase-to-PCC-neutral voltages used by port measurements are

```math
\Delta U_{s,k}=U_{s,k}-U_{s,n}.
```

They are generally different from the internal EMFs. Apply the selected shared
voltage law to ``E_s``, not to ``U_s`` or ``\Delta U_s``.

### Star KCL and grounding

```math
\mathbf1_r^{\mathrm T}\mathbf I_s+I_{s,n}+I_s^g=0,
\qquad I_s^g=-\mathbf1_m^{\mathrm T}\mathbf J_s.
```

The selected grounding law is

```math
\begin{array}{rll}
\mathrm{OPEN}:& I_s^g=0,&\text{no star-to-earth path},\\
\mathrm{IMPEDANCE}:&U_{s,*}=z_s^g I_s^g,&\text{finite or exact-zero impedance},\\
\mathrm{IDEAL}:&U_{s,*}=0,&\text{earth current follows star KCL}.
\end{array}
```

In particular, **do not impose** ``\mathbf1_m^{\mathrm T}J_s=0`` in the grounded
cases. The return through common earth closes conservation. For a diagonal
primitive the neutral drop is ``U_{s,*}-U_{s,n}=z_{s,n}I_{s,n}``; with mutual
coupling, the entire neutral row participates.

A network bus/shunt ground at the PCC neutral is a separate physical location
when the neutral lead has impedance. A coincident ideal source bond and ideal
bus bond with an ideal neutral lead describe the same bond twice. Their separate
currents become nonunique; the draft rejects that duplicate declaration. A
nonzero mutual neutral-row term must be retained when deciding whether the lead
is ideal.

### Powers and the neutral/earth correction

Define individual port, internal EMF, and conductor powers:

```math
\mathbf S_s^{\mathrm{port}}=\Delta\mathbf U_s\circ\mathbf I_s^*,\qquad
\mathbf S_s^{\mathrm{emf}}=\mathbf E_s\circ\mathbf I_s^*,\qquad
\mathbf S_s^{\mathrm{terminal}}=\mathbf U_s\circ\mathbf J_s^*.
```

The complete PCC and internal powers are

```math
S_s^{\mathrm{poc}}=\sum_{k=1}^{r}U_{s,k}I_{s,k}^*+U_{s,n}I_{s,n}^*,\qquad
S_s^{\mathrm{internal}}=\sum_{k=1}^{r}E_{s,k}I_{s,k}^*.
```

Substituting star KCL gives the measurement correction:

```math
S_s^{\mathrm{poc}}
=\sum_{k=1}^{r}(U_{s,k}-U_{s,n})I_{s,k}^*-U_{s,n}(I_s^g)^*
=\mathbf1_r^{\mathrm T}\mathbf S_s^{\mathrm{port}}-U_{s,n}(I_s^g)^*.
```

Thus a total PCC P/Q setpoint, bound, or price uses **complete conductor power**.
A per-port P/Q setpoint or bound uses the individual entries of ``S_s^{port}``.
Their sums are equivalent only when the complex correction vanishes.

### Complex absorption and physical loss

```math
S_s^{\mathrm{series}}=\mathbf J_s^{\mathrm H}Z_s^c\mathbf J_s,\qquad
S_s^{\mathrm{ground}}=U_{s,*}(I_s^g)^*,\qquad
S_s^{\mathrm{internal}}-S_s^{\mathrm{poc}}
=S_s^{\mathrm{series}}+S_s^{\mathrm{ground}}.
```

For impedance grounding,

```math
S_s^{\mathrm{ground}}=z_s^g|I_s^g|^2,\qquad
P_s^{\mathrm{series}}=\mathbf J_s^{\mathrm H}H_s^c\mathbf J_s\ge0,\qquad
P_s^{\mathrm{ground}}=r_s^g|I_s^g|^2\ge0.
```

Ground power is zero for `OPEN` (zero current) and `IDEAL` (zero voltage).
Common earth is at zero potential, so its exported port power is zero. A model
with nonzero earth potential ``U^g`` would instead have
``U_{s,*}-U^g=z_s^g I_s^g`` and exported earth power ``U^g(I_s^g)^*``.
That extended soil/earth network is not encoded here.

### Three-phase sequence identities

For exactly three ordered phase-to-neutral ports, the shared Fortescue transform
applies to ``\Delta U_s``, ``E_s``, and ``I_s``. In particular,

```math
I_s^0=\frac{I_{s,a}+I_{s,b}+I_{s,c}}3,\qquad
3I_s^0+I_{s,n}+I_s^g=0.
```

The neutral current is ``-3I_s^0`` only with no earth return. An earth rating is
not a substitute for a neutral-lead rating or a zero-sequence phase-current rating.

## 5. Inequality constraints

### Cartesian variable bounds

For a declared per-conductor rating ``I_{s,p}^{\max,\mathrm{terminal}}``,

```math
-I_{s,p}^{\max,\mathrm{terminal}}\le J_{s,p}^{\Re}\le I_{s,p}^{\max,\mathrm{terminal}},\qquad
-I_{s,p}^{\max,\mathrm{terminal}}\le J_{s,p}^{\Im}\le I_{s,p}^{\max,\mathrm{terminal}}.
```

Port-current ratings can tighten the corresponding phase-variable boxes.
Earth-current boxes, if used, are bounds on the derived ``I_s^g`` expression.
They are implied by its circle; the implementation need not add them separately.

### Engineering bounds

Apply every declared shared bound using the source-specific powers above.
Conductor, port, and earth ratings are distinct:

```math
\mathbf J_s\circ\mathbf J_s^*
\le\mathbf I_s^{\max,\mathrm{terminal}}\circ\mathbf I_s^{\max,\mathrm{terminal}},\qquad
\mathbf I_s\circ\mathbf I_s^*
\le\mathbf I_s^{\max,\mathrm{port}}\circ\mathbf I_s^{\max,\mathrm{port}},\qquad
|I_s^g|^2\le(I_s^{g,\max})^2.
```

For `ig_max=0`, use ``(I_s^g)^{\Re}=(I_s^g)^{\Im}=0``. With `OPEN` those
constraints already follow from the grounding equality and are not repeated.
A declared limit on a fixed quantity is checked for consistency, not silently
removed. A capability-limited fixed-phasor source may make the network infeasible;
it is not an unlimited slack.

## 6. Implementation and relation to the ideal source

Only ``2m`` real conductor-current variables are needed in addition to existing
bus voltages. Reconstruct star potential, EMFs, and earth current using Part 4.
There are no internal voltage nodes, including at zero impedance. The shared
[implementation page](integration.md) writes all rectangular equations explicitly.

With zero conductor impedance, ideal grounding, and a fixed internal phasor,
``U_{s,n}=0`` and ``U_{s,k}=E_{s,k}``. This recovers the **voltage behavior** of
the existing ideal source. It does not generally recover its stated neutral-current
rule: the dedicated source allows a separate earth return. Do not claim complete
terminal-current equivalence or copy ``I_n=-\sum I_k`` into this limit.

Reference-angle sufficiency is an assembled-network property. A fixed internal
phasor constrains its own absolute orientation; `NONE` and rotating laws do not.
Multiple sources, island references, and the upstream exactly-one-source policy
need explicit problem-level treatment when this draft is adopted.

## 7. Data example

A three-phase source with finite neutral lead and star grounding:

```json
{
  "source_generator": {
    "supply": {
      "bus": "poc",
      "terminal_map": ["a", "b", "c", "n"],
      "configuration": "WYE",
      "R_series_1_1": 0.2, "X_series_1_1": 0.5,
      "R_series_2_2": 0.25, "X_series_2_2": 0.45,
      "R_series_3_3": 0.18, "X_series_3_3": 0.55,
      "R_series_4_4": 0.1, "X_series_4_4": 0.05,
      "grounding_model": "IMPEDANCE",
      "r_ground": 0.8, "x_ground": 0.1,
      "voltage_model": "FIXED_PHASOR",
      "v_magnitude": [235.0, 225.0, 232.0],
      "v_angle": [0.02, -2.13, 2.05],
      "i_max": [40.0, 40.0, 40.0, 20.0],
      "ig_max": 10.0,
      "cost_total": 0.2
    }
  }
}
```

To open its star-earth connection, set `grounding_model` to `OPEN` and omit
`r_ground`/`x_ground`. For an ideal bond use `IDEAL` and omit those parameters.
Neither change creates or removes a PCC bus-ground declaration.
