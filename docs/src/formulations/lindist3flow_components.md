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
| `transformer/wye_delta`, `transformer/delta_wye` | fixed ideal three-phase bank, three retained conductors per side | V, VA, A |
| `transformer/open_delta_regulator` | fixed ideal two-unit bank, ABBC/BCAC/CABA | –, VA, A |

Ratings follow BMOPF's declared shapes: `s_rating` is a scalar nameplate, while
`i_max_from`/`i_max_to` are per-conductor arrays — two entries for an open-delta
bank, one for a single-conductor device, where a bare scalar is also accepted.

Line shunts, series losses, variable taps, and every component not listed here
are rejected by [`check_l3f_applicability`](@ref); they are not silently
omitted. The full exclusion list is on the [formulation overview](lindist3flow.md).

`L3FOptions(unsupported=:lower)` widens the *input* vocabulary without changing
any equation on this page: switches, capacitors, line shunts, and transformer
leakage and no-load admittance are rewritten into the components above, exactly.
`unsupported=:approximate` additionally substitutes load laws and taps, which
does change the problem. Both are described under
[widening the admissible input](lindist3flow.md#Widening-the-admissible-input).

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

!!! note "Voltage-source ratings are restated locally"
    BMOPFTools' preparation at the pinned revision scales a source's
    `v_magnitude`, active/reactive bounds and `cost`, but leaves `s_max` and
    `i_max` in SI. PowerOptLab recomputes both from the caller's SI network on
    every build, so a source nameplate binds identically in either coordinate
    mode. Without that, a rating would be silently ``S_b`` times too loose in the
    default per-unit coordinates. A regression test pins both modes together.

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

This is the first-order Taylor expansion of
``\sqrt{w_\phi w_\psi}\,e^{j(\bar\theta_\phi-\bar\theta_\psi)}`` in the
squared-magnitude coordinates, evaluated about ``\bar w``. It is exact at the
reference point, and exact for every ``w`` when ``\phi=\psi`` — which is why a
grounded-wye device incurs no closure error at all. For a real incidence row
``d_k``, define

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
rows after Kron reduction, so ``H=I`` and the split is exact; `DELTA` uses the
phase-pair incidence rows, and freezing ``H`` at ``\bar v`` is then a second
approximation alongside the cross-voltage closure.

Two invariants hold for **any** reference and are worth relying on:

```math
\sum_\phi (Hs)_\phi=\sum_k s_k,
\qquad
\operatorname{rank}H=\operatorname{rank}D .
```

Total complex power is therefore never distorted by the choice of reference —
only its distribution across terminals is. The rank identity means a closed
three-channel delta (``\operatorname{rank}D=2``) has a one-complex-dimensional
family of channel powers ``s=c\,\bar u`` with identical terminal injection, so
its published per-channel dispatch is determined only up to that circulating
component.

## 3. Variables

The formulation has no voltage-angle, current, tap, binary, or integer decision
variables.

| Variable | Domain | Meaning |
|:--------:|:------:|---------|
| ``w_{i\phi}\ge0`` | real | squared phase-to-ground voltage magnitude |
| ``p_{\ell\phi},q_{\ell\phi}`` | real | sending-end branch terminal power; identical to the receiving end under the lossless balance |
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

For `single_phase`, ``T_t`` is the inverse of ``N=(V^{nom}_{from}/V^{nom}_{to})
\cdot\texttt{tap}``, matching BMOPFTools' ``N = N_0\cdot\texttt{tap}``.

For `single_phase_autotransformer`, BMOPFTools' `_autotransformer_neff` gives
``n_{eff}=a`` for ANSI type A and ``n_{eff}=1/a`` for type B under
``V_{from}=n_{eff}V_{to}``; ``T_t=n_{eff}^{-1}``. This is the self-consistent
reading of `tap_ratio` as a regulated/source ratio: a type-B unit with
``a=1.05`` raises the regulated side by 5 %. The prose table in BMOPF's
regulator specification states the reciprocal pairing and is the side that needs
correcting upstream.

For `wye_delta` and `delta_wye` the coil relation is
``D\boldsymbol v_\Delta=g\boldsymbol v_Y`` with
``g=\sqrt3\,V^{nom}_\Delta/V^{nom}_Y``. Since ``\operatorname{rank}D=2`` the
map exists only with the delta winding upstream, where ``T_t=D/g``; the reverse
orientation is closed by the zero-zero-sequence gauge ``T_t=g\,D^{+}`` under
`unsupported=:approximate`. See
[delta-wye and wye-delta banks](lindist3flow.md#Delta-wye-and-wye-delta-banks).

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

A device rating is enforced once per conductor. Because a supported transformer
is ideal and lossless, its from- and to-side terminal powers coincide, so
`i_max_from` and `i_max_to` both constrain the same branch variable through
their own side's reference voltage. For an open-delta bank the two sides differ
and are constrained separately; the shared phase, which carries both units'
current, is deliberately unrated, matching BMOPF's two-element declaration.

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
coefficients times generator and source active power, following BMOPF's
``\sum_g\sum_k (c^g_k/1000)\,p^g_k`` in currency/h. `:source_import` minimizes
total source active injection in W, and `:feasibility` uses a zero objective.
The reported objective is in the physical SI interpretation in both coordinate
modes: BMOPFTools scales `cost` by ``S_b`` when it builds the per-unit working
copy, which exactly cancels the ``1/S_b`` on the power variable.

Because series losses are omitted, both priced objectives understate the true
cost of serving a load — on IEEE 37 by about 1.95 % of feeder load. They are
sound for comparing dispatches under one formulation, not as absolute costs.

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
