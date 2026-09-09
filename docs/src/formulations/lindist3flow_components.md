# LinDist3Flow component model

This page is the normative mathematical specification of PowerOptLab's
LinDist3Flow formulation. Its organization follows the BMOPFTools component
specifications: data model, input symbols, variables, equality constraints,
inequality constraints, and implementation. The model accepts only
Kron-reduced networks. An explicit-neutral BMOPF case may be passed to the
public builder because it first calls [`kron_reduce_bmopf`](@ref) on a copy.

## 1. Data model and units

The public input is normalized BMOPF JSON in SI units, irrespective of solver
coordinates. Supported top-level component families are shown below.

| Component | Supported BMOPF data | SI input units |
|-----------|----------------------|:--------------:|
| `bus` | retained phase terminals; `v_min`, `v_max` | V |
| `voltage_source` | fixed `v_magnitude`, `v_angle`; P/Q/S/I bounds; linear `cost` | V, rad, W, var, VA, A, currency/kWh |
| `line` / `linecode` | full series R/X matrix; S/I ratings | Ω, VA, A |
| `load` | `SINGLE_PHASE`, `WYE`, `DELTA`; constant P, constant Z, or ZIP with zero I fractions | W, var, V |
| `generator` | channel P/Q boxes, S/I ratings, linear `cost` | W, var, VA, A, currency/kWh |
| `shunt` | fixed full G/B matrix | S |
| `transformer/single_phase` | fixed ideal ratio | V/V |
| `transformer/single_phase_autotransformer` | fixed ideal ANSI A/B ratio | – |
| `transformer/open_delta_regulator` | fixed ideal two-unit bank, ABBC/BCAC/CABA | –, VA, A |

Line shunts, series losses, variable taps, and every component not listed here
are rejected by [`check_l3f_applicability`](@ref); they are not silently
omitted. The full exclusion list is on the [formulation overview](lindist3flow.md).

`L3FOptions(per_unit=true, s_base=S_b)` is the default. BMOPFTools' public
classic scaling preparation supplies the working copy and bases

```math
V_{b,i},\qquad S_b,\qquad
Z_{b,i}=V_{b,i}^2/S_b,\qquad
I_{b,i}=S_b/V_{b,i},\qquad
Y_{b,i}=S_b/V_{b,i}^2.
```

In per-unit coordinates, voltage magnitudes divide by ``V_{b,i}``, powers by
``S_b``, impedances by ``Z_{b,i}``, currents by ``I_{b,i}``, and admittances by
``Y_{b,i}``. Costs are transformed so the physical objective is invariant.
`per_unit=false` uses unit bases. In both cases the caller's dictionary,
reference state, result voltages, result powers, and nonlinear replay are SI.

## 2. Common input symbols

Let ``i,j`` index buses, ``\phi,\psi`` retained terminals, ``\ell:i\to j`` an
edge oriented away from its unique source, and ``k`` a device power channel.

| Symbol | Meaning |
|:------:|---------|
| ``\bar{\boldsymbol v}_i`` | fixed nonzero reference phasors |
| ``\boldsymbol Z_\ell`` | full phase series-impedance matrix |
| ``\boldsymbol D_d`` | real channel-to-terminal incidence matrix |
| ``\boldsymbol T_t`` | fixed parent-to-child voltage map |
| ``\boldsymbol H(D,\bar v)`` | channel-power to terminal-power map |

For each terminal pair, the fixed-angle affine cross-voltage closure is

```math
\widehat{v_\phi v_\psi^*}
=c_{\phi\psi}+a_{\phi\psi}w_\phi+b_{\phi\psi}w_\psi,
\quad
a_{\phi\psi}=\frac{\bar v_\psi^*}{2\bar v_\phi^*},\quad
b_{\phi\psi}=\frac{\bar v_\phi}{2\bar v_\psi},
```

```math
c_{\phi\psi}=\bar v_\phi\bar v_\psi^*
-a_{\phi\psi}|\bar v_\phi|^2-b_{\phi\psi}|\bar v_\psi|^2.
```

It is exact at the reference point. This is a fixed-coefficient affine closure,
not a Taylor-series component model. For a real incidence row ``d_k``, define

```math
\widehat{|d_kv|^2}=\sum_{\phi,\psi}d_{k\phi}d_{k\psi}
\Re\!\left(\widehat{v_\phi v_\psi^*}\right)
=c_{d,k}+a_{d,k}^{\mathsf T}w.
```

The fixed-reference power-allocation map is

```math
H(D,\bar v)=\operatorname{diag}(\bar v)D^{\mathsf T}
\operatorname{diag}(D\bar v)^{-1}.
```

Thus terminal complex power is ``s^{term}=Hs^{ch}``. Both real and imaginary
parts are stamped as real affine expressions. `SINGLE_PHASE`/`WYE` use identity
rows after Kron reduction; `DELTA` uses the relevant phase-pair incidence rows.

## 3. Variables

The formulation has no voltage-angle, current, tap, binary, or integer decision
variables.

| Variable | Domain | Meaning |
|:--------:|:------:|---------|
| ``w_{i\phi}\ge0`` | real | squared phase-to-ground voltage magnitude |
| ``p_{\ell\phi},q_{\ell\phi}`` | real | lossless receiving-end branch terminal power |
| ``p^g_{dk},q^g_{dk}`` | real | generator channel injection |
| ``p^s_{dk},q^s_{dk}`` | real | voltage-source terminal injection |

Load and shunt powers are affine expressions, not independent variables.

## 4. Equality constraints

### Buses and voltage sources

A fixed source at bus ``i`` imposes

```math
w_{i\phi}=|\bar v_{i\phi}|^2.
```

Source active and reactive powers remain free injections subject to their
declared boxes and conic ratings. Exactly one fixed source per radial island is
required.

### Series lines

For ``\ell:i\to j``, define

```math
\Gamma_{\phi\psi}=\bar v_{i\phi}/\bar v_{i\psi},\qquad
M_\ell=2\Re(\overline Z_\ell\odot\Gamma),\qquad
N_\ell=-2\Im(\overline Z_\ell\odot\Gamma).
```

The lossless multiphase voltage drop is

```math
\boldsymbol w_j=\boldsymbol w_i-M_\ell\boldsymbol p_\ell
-N_\ell\boldsymbol q_\ell.
```

The same ``p_\ell+jq_\ell`` enters the child and leaves the parent. No loss
variable or physics slack is present.

### Loads

For load channel ``k``, constant power is

```math
p^d_k=P^{nom}_k,\qquad q^d_k=Q^{nom}_k.
```

Constant impedance and pure ZP ZIP loads use

```math
p^d_k=P^{nom}_k\left(\alpha_{P,k}
+\alpha_{Z,k}\frac{\widehat{|d_kv|^2}}{(V^{nom}_k)^2}\right),
```

```math
q^d_k=Q^{nom}_k\left(\beta_{P,k}
+\beta_{Z,k}\frac{\widehat{|d_kv|^2}}{(V^{nom}_k)^2}\right).
```

Constant impedance is the special case ``\alpha_Z=\beta_Z=1``. A missing ZIP
coefficient family defaults to constant power, consistent with BMOPFTools.
Nonzero ``\alpha_I`` or ``\beta_I`` is rejected because it would require a
non-affine magnitude law or another approximation. Terminal absorption is
``H(D_d,\bar v_i)(p^d+jq^d)``.

### Generators

Each generator channel injects ``p^g_k+jq^g_k``. Its terminal injection is

```math
s^{g,term}=H(D_g,\bar v_i)(p^g+jq^g).
```

This covers retained grounded-wye/single-phase and delta connections with the
same connection map as loads.

### Fixed shunts

For a fixed terminal admittance matrix ``Y``, the absorbed terminal power is

```math
s^{sh}_\phi=\sum_\psi
\widehat{v_\phi v_\psi^*}\,Y_{\phi\psi}^*.
```

It is affine in ``w`` through the common cross-voltage closure.

### Fixed ideal transformers and regulators

Let ``v_j=T_tv_i`` be the fixed parent-to-child map. Each child squared voltage
is

```math
w_{j\phi}=\widehat{|T_{t,\phi:}v_i|^2}.
```

If ``s_t`` denotes child-side terminal power, the lossless parent-side terminal
power is

```math
s_t^{parent}=H(T_t,\bar v_i)s_t.
```

For `single_phase`, ``T_t`` is the inverse fixed nominal/tap ratio. For
`single_phase_autotransformer`, the executable BMOPFTools convention at the
pinned dependency revision is ANSI A: declared `tap_ratio` directly, ANSI B:
its reciprocal, followed by inversion to obtain ``T_t``. PowerOptLab follows
that executable convention. The prose/schema convention in BMOPFTools should
be reconciled upstream if it states the opposite.

For `open_delta_regulator`, ``v_{from}=Av_{to}`` and ``T=A^{-1}``. With effective
ratios ``r_1,r_2`` and ABBC connection,

```math
A=\begin{bmatrix}r_1&1-r_1&0\\0&1&0\\0&1-r_2&r_2\end{bmatrix}.
```

BCAC and CABA are cyclic permutations. The WYE and closed-delta matrices from
Bazrafshan, Gatsis, and Zhu are available through `regulator_gain_matrix`, but
are not stamped as invented BMOPF component subtypes.

### Nodal power balance

For every retained bus terminal, complex injections equal complex absorptions:

```math
s_i^{source}+s_i^{generator}+\sum_{\ell:\,\ell\to i}s_\ell
-\sum_{\ell:\,i\to\ell}s_\ell^{parent}
-s_i^{load}-s_i^{shunt}=0.
```

Lines have ``s_\ell^{parent}=s_\ell``; transformers use the fixed map above.
The implementation stamps the real and imaginary parts separately.

## 5. Inequality constraints

### Voltage and box bounds

Retained bus limits impose

```math
(V^{min}_{i\phi})^2\le w_{i\phi}\le(V^{max}_{i\phi})^2.
```

Generator and source channel boxes impose their declared ``p_min/p_max`` and
``q_min/q_max`` directly.

### Native second-order-cone bounds

Every apparent-power rating is represented exactly as

```math
\left\|\begin{bmatrix}p\\q\end{bmatrix}\right\|_2\le S^{max}.
```

With the fixed-reference current policy, a channel current rating is

```math
\left\|\begin{bmatrix}p\\q\end{bmatrix}\right\|_2
\le I^{max}|\bar u|,
```

where ``\bar u=D\bar v`` for a connected device and ``\bar u=\bar v_\phi``
for a line/source terminal. This is an explicit reference-based approximation,
not an outer linearization of the cone.

For an open-delta unit whose coil spans terminals ``a,b`` but whose terminal
power variable is at terminal ``a``, its winding nameplate is converted by the
fixed reference ratio:

```math
\left\|\begin{bmatrix}p_a\\q_a\end{bmatrix}\right\|_2
\le S^{max}_{coil}\frac{|\bar v_a|}{|\bar v_a-\bar v_b|}.
```

Both declared winding-current sides are likewise converted with the appropriate
fixed terminal voltage. No exponential, power, or other exotic cone appears;
no SOC is replaced by a polyhedral outer approximation.

## 6. Objective and implementation

`objective=:cost` minimizes the sum of per-channel linear energy-cost
coefficients times generator and source active power. `:source_import` minimizes
total source active injection, and `:feasibility` uses a zero objective. The
reported objective is in the physical SI interpretation in both coordinate
modes.

The build sequence is:

1. Parse/copy BMOPF JSON and optionally Kron-reduce explicit neutrals.
2. Run applicability checks and orient each radial island from its source.
3. Construct or validate the fixed SI reference phasors.
4. If requested, use BMOPFTools' public classic scaling preparation to create a
   per-unit working copy and bases.
5. Build the independent PowerOptLab affine/SOC JuMP model.
6. Extract all published voltages and powers back to SI; optionally replay the
   fixed dispatch with BMOPFTools' nonlinear power flow.

`L3FBuild.network` and `L3FBuild.reference` are always SI.
`L3FBuild.working_network`, `L3FBuild.working_reference`, and `L3FBuild.bases`
expose the optimization-coordinate snapshot for inspection. The result's
`formulation` dictionary records `per_unit`, `s_base`, `working_units`, and
`result_units`.
