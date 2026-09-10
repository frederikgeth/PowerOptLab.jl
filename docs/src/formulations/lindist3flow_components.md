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
| `bus` | retained phase terminals; `v_min`, `v_max`, pair-ordered `vpp_min`, `vpp_max`; grounded-reduction `vpn_*` aliases | V |
| `voltage_source` | fixed `v_magnitude`, `v_angle`; P/Q/S/I bounds; linear `cost` | V, rad, W, var, VA, A, currency/kWh |
| `line` / `linecode` | full series R/X matrix; S/I ratings | Ω, VA, A |
| `load` | `SINGLE_PHASE`, `WYE`, `DELTA`; constant P, constant Z, or ZIP with zero I fractions | W, var, V |
| `generator` / restricted `ibr` | channel P/Q boxes, S/I ratings, signed fixed PF, aggregate P/Q/S bounds, fixed-P/free-Q voltage target, linear `cost` | W, var, VA, A, V, currency/kWh |
| `shunt` | fixed full G/B matrix | S |
| `transformer/single_phase` | fixed ideal ratio | V/V |
| `transformer/center_tap` | fixed coupled 1→2 split-phase ratio; HV side upstream | V, VA, A |
| `transformer/single_phase_autotransformer` | fixed ANSI A/B ratio; connection-aware leakage and from-side exciting branch under `:lower` | Ω, S, VA, A |
| `transformer/wye_delta`, `transformer/delta_wye` | fixed three-phase bank, three retained conductors per side; wye-referred leakage and external winding-2 shunt under `:lower` | V, Ω, S, VA, A |
| `transformer/open_delta_regulator` | fixed ideal two-unit bank, ABBC/BCAC/CABA | –, VA, A |
| L3F-local `grounded_wye_wye`, `delta_delta`, `closed_delta_regulator` | fixed ideal three-phase maps; `delta_delta` is experimental | V, A |

Ratings follow BMOPF's declared shapes: `s_rating` is a scalar nameplate, while
`i_max_from`/`i_max_to` are per-conductor arrays — two entries for an open-delta
bank, one for a single-conductor device, where a bare scalar is also accepted.
After neutral reduction a center tap has one `i_max_from` entry and two
`i_max_to` entries; its scalar `s_rating` belongs to the aggregate HV coil.
For Yd/Dy banks, `s_rating` is total bank VA and is enforced only on the three
wye coils, at `s_rating/3`; the delta terminals do not receive a duplicate
nameplate cone.

Line shunts, series losses, variable taps, and every component not listed here
are rejected by [`check_l3f_applicability`](@ref) unless an explicitly
documented canonical lowering is selected; they are not silently omitted. The
full exclusion list is on the [formulation overview](lindist3flow.md).

`L3FOptions(unsupported=:lower)` widens the *input* vocabulary without changing
the one-shot S-W model class: switches, capacitors, line shunts, supported
connection-aware transformer equivalents, and restricted fixed-control IBRs
are rewritten into the components above. These are canonical L3F lowerings,
not exact AC equivalences.
`unsupported=:approximate` additionally substitutes load laws and taps, and
admits the explicitly reported delta-delta common-mode projection, which
does change the problem. Both are described under
[widening the admissible input](lindist3flow.md#Widening-the-admissible-input).
`unsupported=:permissive` further removes well-formed operational constraints
that lack an S-W representation. Each removal is recorded with its original
value, and the result marks the solved network semantics as projected; this
does not certify feasibility for the supplied nonlinear network.

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
grounded-wye device incurs no cross-phase closure error. For a real incidence row
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
rows after Kron reduction, so ``H=I`` and the channel split is exact; `DELTA` uses the
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
| ``p_{\ell\phi},q_{\ell\phi}`` | real | lossless series or transformer child-side power, oriented from the source-side parent toward the child; endpoint terminal totals used by ratings are affine expressions derived from it |
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

Under `unsupported=:approximate`, each current fraction is replaced by its
nominal-voltage tangent, half Z and half P. An exponential term with any finite
``\gamma`` becomes ``\alpha_Z=\gamma/2`` and
``\alpha_P=1-\gamma/2`` (independently for P and Q). Coefficients are retained
verbatim without normalization or clamping, so values outside the usual
exponent range can produce negative Z or P coefficients. The replacement
matches value and slope at nominal voltage; it is not a global load-law claim.

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

For a Kron-reduced `center_tap`, `v_nom_to` is the per-leg voltage and the
retained maps have one HV terminal and two LV hot legs. With the centre tap
merged into ideal ground,

```math
T_{ct}=\begin{bmatrix}1/N\\-1/N\end{bmatrix},\qquad
N=(V^{nom}_{from}/V^{nom}_{to})\cdot\texttt{tap}.
```

The sign is the series-aiding winding polarity: the two terminal reference
phasors are 180 degrees apart. The existing power map gives
``H(T_{ct},\bar v)=[1\;1]``, hence

```math
s^{parent}_{ct}=s_1+s_2.
```

This is a coupled three-winding device, not two independent transformers. Its
fixed map is currently defined only with the declared HV/from winding upstream;
power may still reverse through the resulting source-oriented branch.

For `single_phase_autotransformer`, BMOPFTools' `_autotransformer_neff` gives
``n_{eff}=a`` for ANSI type A and ``n_{eff}=1/a`` for type B under
``V_{from}=n_{eff}V_{to}``; ``T_t=n_{eff}^{-1}``. This is the self-consistent
reading of `tap_ratio` as a regulated/source ratio: a type-B unit with
``a=1.05`` raises the regulated side by 5 %. The prose table in BMOPF's
regulator specification states the reciprocal pairing and is the side that needs
correcting upstream.

For `wye_delta` and `delta_wye` the coil relation is
``D\boldsymbol v_\Delta=g\boldsymbol v_Y`` with
nominal ``g_0=\sqrt3\,V^{nom}_\Delta/V^{nom}_Y``. The native from-winding tap
gives ``g=g_0/tap`` for Yd and ``g=g_0tap`` for Dy. Since ``\operatorname{rank}D=2`` the
map exists only with the delta winding upstream, where ``T_t=D/g``; the reverse
orientation is closed by the zero-zero-sequence gauge ``T_t=g\,D^{+}`` under
`unsupported=:approximate`. This pseudo-inverse also projects the upstream wye
voltage through ``DD^+=I-\mathbf1\mathbf1^T/3``; it is therefore an upstream
zero-sequence projection as well as a downstream gauge choice. See
[delta-wye and wye-delta banks](lindist3flow.md#Delta-wye-and-wye-delta-banks).

The local `grounded_wye_wye` map is diagonal. The local `delta_delta` map
projects away upstream common mode and chooses zero downstream common mode, so
it is available only under `:approximate`; its neutral and zero-sequence
current-path consistency is not certified. `closed_delta_regulator` uses the
three-unit matrix returned by `regulator_gain_matrix("CLOSED_DELTA", ...)`.
These local subtypes are ideal-only. `delta_delta` and `closed_delta_regulator`
reject `s_rating` because no coil/unit nameplate convention is declared;
side-terminal `i_max` remains available. Nonlinear replay is unavailable for
these local-only subtype names.

For `open_delta_regulator`, ``v_{from}=Av_{to}`` and ``T=A^{-1}``. With effective
ratios ``r_1,r_2`` and ABBC connection,

```math
A=\begin{bmatrix}r_1&1-r_1&0\\0&1&0\\0&1-r_2&r_2\end{bmatrix}.
```

BCAC and CABA are cyclic permutations. The WYE and closed-delta matrices from
Bazrafshan, Gatsis, and Zhu are also reachable directly through
`regulator_gain_matrix`, and back the `closed_delta_regulator`,
`grounded_wye_wye` and `delta_delta` banks. Those three subtype names are
**L3F-local**: BMOPF does not define them, so a network using one is not
schema-valid and cannot be replayed through the nonlinear power flow. See
[L3F-local bank subtypes](lindist3flow.md#L3F-local-bank-subtypes).

### Nodal power balance

For every retained bus terminal, complex injections equal complex absorptions:

```math
s_i^{source}+s_i^{generator}+\sum_{\ell:\,\ell\to i}s_\ell
-\sum_{\ell:\,i\to\ell}s_\ell^{parent}
-s_i^{load}-s_i^{shunt}=0.
```

Lines have ``s_\ell^{parent}=s_\ell``; transformers use the fixed map above.
At a child balance the same lossless flow enters with the opposite device-side
sign. Endpoint shunts remain separate terms in nodal balance.
The implementation stamps the real and imaginary parts separately.

## 5. Inequality constraints

### Voltage and box bounds

A rating is enforced at the physical endpoint where BMOPF declares it. Lines
and lowered closed switches receive constraints at both input endpoints. Let
``e`` be an endpoint and define its sign from the source-oriented topology,
not from the input's `from`/`to` labels:

```math
\sigma_e=\begin{cases}+1,&e=\operatorname{parent}(\ell),\\
-1,&e=\operatorname{child}(\ell),\end{cases}
\qquad
s^{end}_{\ell e\phi}=\sigma_e s^{series}_{\ell\phi}+s^{sh}_{\ell e\phi}.
```

The endpoint is still keyed by its original `from` or `to` label, even when
that label is opposite to the oriented traversal. Thus line pi-shunt power is
inside the corresponding rating cone even though canonical lowering represents
the shunt as a standalone element.

Retained bus limits impose

```math
(V^{min}_{i\phi})^2\le w_{i\phi}\le(V^{max}_{i\phi})^2.
```

Generator and source channel boxes impose their declared ``p_min/p_max`` and
``q_min/q_max`` directly.

### Native second-order-cone bounds

For compactness define the two native cone templates

```math
\mathcal S(p,q;S):\quad [S,p,q]\in\mathcal Q_3
\quad\Longleftrightarrow\quad p^2+q^2\le S^2,
```

```math
\mathcal I(p,q;w,I):\quad [w,I^2/2,p,q]\in\mathcal Q_4^r
\quad\Longleftrightarrow\quad p^2+q^2\le wI^2.
```

The builder stamps the algebraically identical, numerically conditioned form

```math
[I w/\bar V,\ I\bar V/2,\ p,\ q]\in\mathcal Q_4^r,
```

where ``\bar V>0`` is the magnitude of the fixed voltage reference for that
terminal or winding channel. The reciprocal scaling of the first two rotated-
cone axes cancels exactly, so ``\bar V`` does not alter the feasible set. It
keeps raw-SI cone coordinates on the apparent-power scale and is unnecessary
but harmless in per unit.

The model uses the following complete set of conic bounds.

**Generators.** For every channel ``k`` with declared ratings,

```math
\mathcal S(p^g_k,q^g_k;S^{max}_k),\qquad
\mathcal I(p^g_k,q^g_k;\widehat{|d_kv_i|^2},I^{max}_k).
```

The current cone therefore uses the live affine fixed-angle winding-voltage
closure, not ``|d_k\bar v|^2``.

**Voltage sources.** For each source terminal ``\phi``,

```math
\mathcal S(p^s_\phi,q^s_\phi;S^{max}_\phi),\qquad
\mathcal I(p^s_\phi,q^s_\phi;w_{i\phi},I^{max}_\phi).
```

The source equality fixes ``w_{i\phi}=|\bar v_{i\phi}|^2``, but the same native
rotated-SOC template is retained.

**Lines and lowered closed switches.** For every conductor ``\phi`` and both
physical endpoints ``e``, using the endpoint total defined above,

```math
\mathcal S(\Re s^{end}_{\ell e\phi},\Im s^{end}_{\ell e\phi};S^{max}_{\ell\phi}),
```

```math
\mathcal I(\Re s^{end}_{\ell e\phi},\Im s^{end}_{\ell e\phi};
           w_{e\phi},I^{max}_{\ell\phi}).
```

A rating declared directly on a line takes precedence; otherwise the linecode
rating is inherited. A lowered switch copies both `s_max` and `i_max` onto its
zero-impedance line. A lowered pi shunt contributes ``s^{sh}_{\ell e\phi}`` at
its own endpoint.

Both endpoint cones are stamped even for a shunt-free line. In an exact AC
series branch,

```math
|S_{parent}|^2/w_{parent}=|S_{child}|^2/w_{child}=|I^{series}|^2.
```

The lossless model instead reuses one series ``S`` while allowing the endpoint
voltages to differ, so ``|S|/\sqrt{w_{parent}}`` and
``|S|/\sqrt{w_{child}}`` need not agree. Requiring both current cones takes the
tighter endpoint-derived value within the LinDist3Flow surrogate. This does not
certify AC ampacity and differs
from BMOPFTools' nonlinear branch builder, which omits its separate to-side
current cone on a shunt-free line. Both ``S^{max}`` cones are also retained at
all lines and lowered switches for one uniform endpoint contract.

**Ordinary single-phase transformers.** Let
``s^{end}_{te\phi}`` be power entering the transformer at original endpoint
``e\in\{from,to\}``: the oriented parent expression is ``Hs_t``, the child
expression is ``-s_t``, and any lowered exciting shunt is added at the original
to-side endpoint. The scalar nameplate is applied at both endpoints,

```math
\mathcal S(\Re s^{end}_{te\phi},\Im s^{end}_{te\phi};S^{rating}_t),
\qquad e\in\{from,to\},
```

while a current rating is stamped only on a side that declares it:

```math
\mathcal I(\Re s^{end}_{te\phi},\Im s^{end}_{te\phi};
           w_{e\phi},I^{max}_{te\phi}).
```

**Single-phase autotransformers.** The scalar `s_rating` bounds bare
from-side through power only, matching BMOPFTools' autotransformer advantage
contract. The exciting branch is at the original from terminal. Its power is
added for `i_max_from`, while `i_max_to` bounds bare output series current.
Thus the nameplate excludes exciting power even though input ampacity includes
the exciting current.

Leakage lowering preserves the original physical endpoint bus and terminal map
for these cones; a synthetic internal-bus voltage is never substituted for the
declared endpoint voltage. The exciting shunt remains on the original external
winding-2 bus and is included in the ordinary transformer's endpoint-current
expression.

**Center-tap transformers.** Let ``s_1,s_2`` be the two child-leg powers and
``s_\Sigma=s_1+s_2`` the aggregate HV-coil power. BMOPFTools' scalar nameplate
is imposed once on that coil,

```math
\mathcal S(\Re s_\Sigma,\Im s_\Sigma;S^{rating}_{ct}).
```

The retained primary and leg-current bounds are live-voltage rotated SOCs,

```math
\mathcal I(\Re s_\Sigma,\Im s_\Sigma;w_h,I^{max}_{from}),
\qquad
\mathcal I(\Re s_k,\Im s_k;w_k,I^{max}_{to,k}),\quad k=1,2.
```

The center-tap `i_max_to` convention bounds the bare leg-winding currents, so
its exciting-shunt current is not added to either leg cone. The live voltages
are still those at the original external LV endpoints after leakage lowering.

Under `unsupported=:lower`, the three-winding star leakage becomes one
from-side series arm carrying ``s_\Sigma`` and two identical to-side arms
carrying ``s_1,s_2``. Consequently the explicit line/ideal-transformer network
stamps exactly the lossless LinDistFlow equations

```math
w_k={w_h\over N^2}
-2\left({r_hP_\Sigma+x_hQ_\Sigma\over N^2}
       +r_\ell P_k+x_\ell Q_k\right),\qquad k=1,2.
```

The shared primary term is the retained coupling between the legs. The no-load
admittance is materialized once across winding 2 (LV leg 1), following the
BMOPFTools/OpenDSS convention; it is not duplicated on leg 2.

**Yd/Dy banks.** BMOPF's scalar ``S^{rating}_{bank}`` is a total-bank
nameplate. It is divided equally across the three wye coils and is not copied
onto the delta terminals:

```math
\mathcal S(p_{Y,k},q_{Y,k};S^{rating}_{bank}/3),
\qquad k=1,2,3.
```

Declared side-specific current limits retain BMOPF's terminal/bushing-current
semantics on both wye and delta sides:

```math
\mathcal I(\Re s^{end}_{te\phi},\Im s^{end}_{te\phi};
           w_{e\phi},I^{max}_{te\phi}),
\qquad e\in\{from,to\}.
```

In particular, a delta-side current cone uses live phase-to-ground terminal
``w_{e\phi}``; it is not an internal delta-coil current bound.

**Open-delta regulators.** For unit ``k`` spanning terminals ``a,b``, choose
the unit's non-shared terminal ``a`` as prescribed by the ABBC/BCAC/CABA
connection. At each original side ``e``, the fixed-reference conversion from
coil nameplate to terminal apparent power is

```math
\mathcal S(p_{tea},q_{tea};
S^{rating}_{t}|\bar v_{ea}|/|\bar v_{ea}-\bar v_{eb}|).
```

Its current rating instead uses the live terminal voltage:

```math
\mathcal I(p_{tea},q_{tea};w_{ea},I^{max}_{te,k}).
```

The common terminal, which carries both units' current, receives no separate
cone because BMOPF declares two unit ratings rather than three conductor
ratings. No exponential, power, or other exotic cone appears, and no SOC is
replaced by a polyhedral outer approximation.

## 6. Objective, extraction, and implementation

`objective=:cost` minimizes the sum of per-channel linear energy-cost
coefficients times generator and source active power, following BMOPF's
``\sum_g\sum_k (c^g_k/1000)\,p^g_k`` in currency/h. `:source_import` minimizes
total source active injection in W, and `:feasibility` uses a zero objective.
The reported objective is in the physical SI interpretation in both coordinate
modes: BMOPFTools scales `cost` by ``S_b`` when it builds the per-unit working
copy, which exactly cancels the ``1/S_b`` on the power variable.

Both priced objectives omit the effect of series losses; on IEEE 37 the omitted
source import is about 1.95 % of feeder load. Their values and rankings apply
only to the stated surrogate. With general cost coefficients they are neither
absolute AC costs nor guaranteed bounds or rankings for the nonlinear AC
problem.

### Result extraction

The stable result is deliberately a small semantic projection of the JuMP
solution, not a reconstructed AC state. Numerical decision-variable values are
published only when the solve is optimal; otherwise they are `NaN`. In SI mode
they are copied directly. In per-unit mode the extractor applies only

```math
w^{SI}_{i\phi}=V_{b,i}^2w^{pu}_{i\phi},\qquad
(p^{SI},q^{SI})=S_b(p^{pu},q^{pu}).
```

It then derives

```math
v^{mag}_{i\phi}=\sqrt{\max(0,w^{SI}_{i\phi})}.
```

The `max` protects the square root against a tiny negative solver-tolerance
artifact; the separately published `w` is the unclamped scaled value. The
extractor also attaches fixed or discrete metadata rather than new electrical
solutions:

- `reference_angle` is ``\arg\bar v`` and is not optimized;
- `parent`, `child`, terminal maps, and `reversed_from_input` describe
  the source-oriented topology;
- `effective_ratio_from_to` is recomputed from transformer input data.

For a line, published `p` and `q` are the lossless series-flow variables. For a
transformer, they are the child-side variables ``s_t`` used by the fixed power
map; a center tap therefore returns two entries in `terminal_map_child` order.
Neither is converted back to input `from → to` orientation, and the affine
endpoint totals used inside rating cones are not published separately. Loads,
shunts, currents, complex voltage phasors, losses, and residual corrections are
not synthesized during extraction.

In per-unit mode `objective=:source_import` is multiplied by ``S_b``. The cost
objective is already physical because its prepared coefficient scaling cancels
the power-variable scaling, while the feasibility objective remains zero.

Canonical preprocessing is also not undone: `L3FBuild.network` and
`L3FResult.network` are the SI Kron-reduced, lowered, or projected snapshots,
and synthetic `_l3f_...` buses and branches remain in the published topology.
An optional nonlinear replay fixes generator dispatch on a deep copy and fills
the `validation` dictionary; it does not overwrite or polish any linear result.

### Build sequence

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
