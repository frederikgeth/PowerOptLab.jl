# Generator implementation and upstream migration

The component and shared operating pages define the physical feasible set.
This page documents an equivalent rectangular realization, the mapping to the
current PowerOptLab prototype, and the changes needed before upstream adoption.
The proposed JSON collections are **not accepted by an implemented PowerOptLab
JSON importer**, and these pages do not supply a released upstream JSON Schema.

## 1. Rectangular electrical equations

For a complex phasor ``z=z^{\Re}+jz^{\Im}``, each complex equality represents
two real equalities. Real incidence matrices do not mix real and imaginary parts.

### Generalized generator

With ``Z=R+jX`` and the component's incidence matrix ``C``,

```math
\Delta U^{\Re}=C^{\mathrm T}U^{\Re},\qquad
\Delta U^{\Im}=C^{\mathrm T}U^{\Im},\qquad
J^{\Re}=CI^{\Re},\qquad J^{\Im}=CI^{\Im},
```

```math
D^{\Re}=RI^{\Re}-XI^{\Im},\qquad
D^{\Im}=XI^{\Re}+RI^{\Im},\qquad
E^{\Re}=\Delta U^{\Re}+D^{\Re},\qquad
E^{\Im}=\Delta U^{\Im}+D^{\Im}.
```

Each independent current component enters the network balance through

```math
\kappa_{i,t}^{\mathrm{other},\Re}+\sum_{p:N_p=t}J_p^{\Re}=0,\qquad
\kappa_{i,t}^{\mathrm{other},\Im}+\sum_{p:N_p=t}J_p^{\Im}=0.
```

These are affine equalities. Quantities on the left of the preceding definitions
may remain expressions; none requires a duplicated voltage variable.

### Source generator

Use the full conductor current ``J``, with phase subvector ``I`` and final neutral
entry ``J_n=I_n``:

```math
D^{\Re}=R^cJ^{\Re}-X^cJ^{\Im},\qquad
D^{\Im}=X^cJ^{\Re}+R^cJ^{\Im},
```

```math
U_*^{\Re}=U_n^{\Re}+D_n^{\Re},\qquad
U_*^{\Im}=U_n^{\Im}+D_n^{\Im},\qquad
E_k^{\Re}=U_k^{\Re}+D_k^{\Re}-U_*^{\Re},\qquad
E_k^{\Im}=U_k^{\Im}+D_k^{\Im}-U_*^{\Im},
```

```math
\Delta U_k^{\Re}=U_k^{\Re}-U_n^{\Re},\qquad
\Delta U_k^{\Im}=U_k^{\Im}-U_n^{\Im},\qquad
(I^g)^{\Re}=-\sum_{p=1}^mJ_p^{\Re},\qquad
(I^g)^{\Im}=-\sum_{p=1}^mJ_p^{\Im}.
```

For ``z^g=r^g+jx^g``, impedance grounding becomes

```math
U_*^{\Re}=r^g(I^g)^{\Re}-x^g(I^g)^{\Im},\qquad
U_*^{\Im}=x^g(I^g)^{\Re}+r^g(I^g)^{\Im}.
```

Ideal grounding is ``U_*^{\Re}=U_*^{\Im}=0``. Open grounding is
``(I^g)^{\Re}=(I^g)^{\Im}=0``. These equations are well defined at exact zero
impedance. The conductor KCL equations use ``J`` exactly as above, including the
independent neutral component.

### Powers, loss, and the grounded-source correction

For any voltage/current pair ``v,i`` at the same measurement location,

```math
P(v,i)=v^{\Re}i^{\Re}+v^{\Im}i^{\Im},\qquad
Q(v,i)=v^{\Im}i^{\Re}-v^{\Re}i^{\Im}.
```

Apply this definition to ``(\Delta U_k,I_k)``, ``(E_k,I_k)``, and ``(U_p,J_p)``
for port, EMF, and conductor powers. Complete PCC P/Q sums all conductor terms;
internal P/Q sums the EMF terms. For a source,

```math
P^{\mathrm{poc}}=\sum_kP^{\mathrm{port}}_k
-U_n^{\Re}(I^g)^{\Re}-U_n^{\Im}(I^g)^{\Im},
```

```math
Q^{\mathrm{poc}}=\sum_kQ^{\mathrm{port}}_k
-U_n^{\Im}(I^g)^{\Re}+U_n^{\Re}(I^g)^{\Im}.
```

Let ``L=I`` for the port primitive or ``L=J`` for the conductor primitive, and
let ``D=ZL`` be the corresponding drop. Then

```math
P^{\mathrm{series}}=\sum_p(D_p^{\Re}L_p^{\Re}+D_p^{\Im}L_p^{\Im}),\qquad
Q^{\mathrm{series}}=\sum_p(D_p^{\Im}L_p^{\Re}-D_p^{\Re}L_p^{\Im}),
```

```math
P^{\mathrm{ground}}=U_*^{\Re}(I^g)^{\Re}+U_*^{\Im}(I^g)^{\Im},\qquad
Q^{\mathrm{ground}}=U_*^{\Im}(I^g)^{\Re}-U_*^{\Re}(I^g)^{\Im}.
```

Ground terms are zero for the generalized generator. For either component,

```math
P^{\mathrm{internal}}-P^{\mathrm{poc}}-P^{\mathrm{series}}-P^{\mathrm{ground}}=0,\qquad
Q^{\mathrm{internal}}-Q^{\mathrm{poc}}-Q^{\mathrm{series}}-Q^{\mathrm{ground}}=0.
```

These are useful independent physical residuals, not extra equality rows when
the circuit already implies them. In particular, a loss variable must not be
allowed to increase arbitrarily to satisfy another constraint or objective.

## 2. Rectangular operating equations

### Fixed phasor and rotating shapes

Mode b uses

```math
E_k^{\Re}=\overline v_k\cos\overline\theta_k,\qquad
E_k^{\Im}=\overline v_k\sin\overline\theta_k.
```

For a fixed complex shape ratio ``b_k=b_k^{\Re}+jb_k^{\Im}`` (``\beta_k`` in c,
``\gamma_k`` in d), the anchored equations are

```math
E_k^{\Re}=b_k^{\Re}E_1^{\Re}-b_k^{\Im}E_1^{\Im},\qquad
E_k^{\Im}=b_k^{\Im}E_1^{\Re}+b_k^{\Re}E_1^{\Im}\quad(k>1).
```

Mode c additionally uses ``(E_1^{\Re})^2+(E_1^{\Im})^2=\overline v_1^2``.
Mode d uses the common-magnitude bounds instead. Their sine/cosine coefficients
are evaluated from **data**, not nonlinear functions of optimization variables.

For mode e, write ``\gamma_k=c_k+js_k`` and rotate

```math
W_k^{\Re}=c_kE_k^{\Re}+s_kE_k^{\Im},\qquad
W_k^{\Im}=c_kE_k^{\Im}-s_kE_k^{\Re},
```

```math
W_k^{\Im}E_1^{\Re}-W_k^{\Re}E_1^{\Im}=0,\qquad
W_k^{\Re}E_1^{\Re}+W_k^{\Im}E_1^{\Im}\ge0\quad(k>1).
```

Positive per-port magnitude minima complete its angle domain. Both expressions
may be divided by one fixed positive voltage scale squared. No tangent, atan,
all-pairs angle equations, or independent absolute-angle pins are required.

### Sequence measurements and terminal PV

For any three-port vector ``z`` and ``F=F^{\Re}+jF^{\Im}``,

```math
z_{012}^{\Re}=F^{\Re}z^{\Re}-F^{\Im}z^{\Im},\qquad
z_{012}^{\Im}=F^{\Im}z^{\Re}+F^{\Re}z^{\Im}.
```

Use ``z=\Delta U,E,I`` as required. Each entry is affine. The PV equality is

```math
(V^{\mathrm{meas},\Re}/V^\star)^2
+(V^{\mathrm{meas},\Im}/V^\star)^2=1.
```

P/Q control equalities and capability boxes use the quadratic power expressions
above at the declared locations; this does not introduce fourth-degree terms.

## 3. Bound realization and exact elimination

For a magnitude-bounded affine complex expression ``z`` and a positive upper
bound ``h``, use

```math
(z^{\Re}/h)^2+(z^{\Im}/h)^2\le1.
```

For positive lower bound ``l``, use

```math
(z^{\Re}/l)^2+(z^{\Im}/l)^2\ge1.
```

Zero upper bounds become the two affine equations ``z^{\Re}=z^{\Im}=0``;
zero lower bounds add nothing. Equal positive bounds give one normalized circle
equality. Explicit current variables also receive their implied rectangular
boxes. Identical affine norms share the intersection of their bounds, including
sign-reversed single-phase returns. Identical affine equalities are imposed once.

If an apparent-power bound is applied after substituting bilinear P/Q, directly
squaring those expressions would produce quartic terms. For a positive rating
``\overline S``, the prototype instead introduces normalized real auxiliary
variables ``\widehat p,\widehat q``:

```math
\widehat p=P(v,i)/\overline S,\qquad
\widehat q=Q(v,i)/\overline S,\qquad
-1\le\widehat p,\widehat q\le1,\qquad
\widehat p^2+\widehat q^2\le1.
```

Each equality and inequality is degree at most two. For a total rating, replace
``P(v,i),Q(v,i)`` by the corresponding complete totals. For ``\overline S=0``,
use P=Q=0 directly. If a same-location control already fixes a bounded power,
substitute and check that known value instead of adding a duplicate control
constraint. These auxiliary variables are an **implementation choice**; the
foundational apparent-power equation remains the one in the operating page.

### Shape-dependent sequence constraints

For a known rotating shape ``E=bE_1``, let ``f=Fb``. Then

```math
E^\sigma=f_\sigma E_1,\qquad |E^\sigma|=|f_\sigma|\,|E_1|.
```

If ``f_\sigma=0``, its zero upper bound is already satisfied and a positive lower
bound conflicts with the shape. Otherwise its sequence bounds become bounds
on the one anchor magnitude:

```math
V^{\sigma,\min}/|f_\sigma|\le|E_1|
\le V^{\sigma,\max}/|f_\sigma|.
```

Intersect those intervals with all ``E_k^{\min}/|b_k|`` and
``E_k^{\max}/|b_k|``. With fixed template magnitude, check the intervals directly.
This avoids appending dependent zero-sequence equations to an already balanced
source law. It does not remove a physical limit or change the feasible set.

For mode e with exact internal sequence zeros indexed by ``\mathcal K``, write
``E=e^{j\theta}\operatorname{diag}(e^{j\phi})\rho`` with real positive
``\rho``. The reduction operates on

```math
M\rho=0,\qquad
M=\begin{bmatrix}
\mathfrak R(F[\mathcal K,:]\operatorname{diag}(e^{j\phi}))\\
\mathfrak I(F[\mathcal K,:]\operatorname{diag}(e^{j\phi}))
\end{bmatrix}.
```

If its nullspace is one-dimensional with a strictly positive shape ``a``, normalize
``a_1=1`` and use ``b_k=a_k e^{j(\phi_k-\phi_1)}`` and the preceding affine
shape equations. The current prototype rejects exact-zero combinations without
such a unique positive shape. This is a **prototype support restriction**;
the canonical angle/sequence equations themselves also define other cases.
No general Jacobian-rank or global-solvability certificate is claimed.

### Working units

JSON and published results stay in SI. With local engine bases satisfying
``S_b=V_bI_b`` and ``Z_b=V_b/I_b``, solver coordinates are

```math
\widetilde U=U/V_b,\quad\widetilde E=E/V_b,\quad
\widetilde I=I/I_b,\quad\widetilde Z=Z/Z_b,\quad
\widetilde P=P/S_b,\quad\widetilde Q=Q/S_b.
```

Convert ratings and targets by the corresponding base before building their rows.
Bases, starts, tolerances, redundant variable boxes, and optimizer choice belong
to the implementation. They must not introduce undeclared physical limits or
claim dynamic stability, uniqueness, or global optimality.

## 4. Mapping the draft data to PowerOptLab

This is an explicit translation contract for a **future** importer. It does not
claim that `parse_bmopf` or a current JSON Schema recognizes the new collections.

| Draft field(s) | Current Julia API | Translation |
|:--|:--|:--|
| Collection key, `bus` | `id`, `bus` | Copy strings |
| Generalized `terminal_map`, `configuration`, `port_map` | `GeneralizedGenerator.connections` | Construct ordered outgoing/return terminal-string pairs |
| Source `terminal_map` | `SourceGenerator.phase_terminals`, `neutral` | All but final entry; final entry |
| `R_series_p_q`, `X_series_p_q` | `impedance` | Materialize the full port or conductor matrix with omitted entries zero |
| `grounding_model`, `r_ground`, `x_ground` | `grounding` | `OPEN` → `nothing`; `IDEAL` → `0.0`; `IMPEDANCE` → `r_ground + im*x_ground` |
| `voltage_model` | `GeneratorVoltageLaw.mode` | Corresponding lowercase symbol |
| `v_magnitude`, `v_angle` | `GeneratorVoltageLaw.phasor` | `v_magnitude .* cis.(v_angle)` |
| `angle_offsets` | `GeneratorVoltageLaw.angles` | Copy full vector |
| `e_min`, `e_max` | `GeneratorVoltageLaw.magnitude_min/max` | Copy vectors, except common mode uses max of minima / min of maxima |
| `power_setpoint_location` | `GeneratorControl.power_location` | Lowercase symbol; absent when no P/Q setpoint can map to unused default `:poc` |
| `p_set` / `p_total_set`, `q_set` / `q_total_set` | `GeneratorControl.p/q` | Vector / scalar; omission → `nothing` |
| `v_target`, `v_target_measurement` | `voltage_target`, `voltage_metric` | Copy target; `POSITIVE_SEQUENCE` → `:positive_sequence`; `PORT` → `:phase` |
| `power_limit_location` | `GeneratorCapability.power_location` | Lowercase symbol; absent when no P/Q/S limit can map to unused default `:poc` |
| `p_min/max` / `p_total_min/max` | `GeneratorCapability.p_min/max` | Per-side vector / scalar |
| `q_min/max` / `q_total_min/max` | `GeneratorCapability.q_min/max` | Per-side vector / scalar |
| `s_max` / `s_total_max` | `GeneratorCapability.s_max` | Vector / scalar |
| `i_max` | `GeneratorCapability.terminal_i_max` | External conductor vector, permuted if needed |
| `i_port_max` | `GeneratorCapability.i_max` | Port-current vector |
| `v_min`, `v_max` | `GeneratorCapability.voltage_min/max` | PCC port-voltage vectors |
| `ig_max` | `GeneratorCapability.earth_i_max` | Source-only scalar |
| `v_sequence_location` | `GeneratorSequenceLimits.location` | Lowercase symbol; current-only bounds may use the unused default `:poc` |
| `v_seq_min/max`, `i_seq_min/max` | `GeneratorSequenceLimits.voltage_min/max`, `current_min/max` | Copy vectors in 0,1,2 order; attach to `GeneratorCapability.sequence` |
| `cost_total` | `cost` | Scalar PCC price |

The component key is unique within its JSON collection. The prototype's combined
device list requires unique IDs across both new component types. If the same
string occurs in both collections, an importer must qualify internal IDs by
collection and retain the reversible correspondence; it must not merge or drop
one of the records. Direct string copying applies when the IDs are already unique.

JSON `number` includes integer-valued numbers. A Julia importer must explicitly
convert physical numeric scalars and arrays to `Float64` before constructing the
prototype's scalar-or-vector bound fields; array indices remain integers. It must
not rely on a JSON library preserving a decimal spelling as a floating-point type.

For general `PORTS` data, the prototype reports distinct outgoing terminals first,
then previously unseen return terminals. If that differs from draft `terminal_map`,
an importer must permute **conductor** ratings and exported conductor results.
It must not permute the port impedance or port quantities, whose order remains
the order of the port pairs. Source conductor order already agrees directly.

The API permits scalar broadcasts for some current/voltage bounds and convenient
scalar/diagonal impedance input. The canonical draft JSON uses full-length arrays
and explicit matrix entries; shorthand constructors are not additional schema
encodings. `voltage_scale` is an implementation-only numerical option.

Published result correspondence is:

| Mathematical measurement | `GeneratorResult` field |
|:--|:--|
| ``U``, ``\Delta U``, ``E`` | `terminal_voltage`, `port_voltage`, `emf` |
| ``J``, ``I`` | `terminal_current`, `port_current` |
| ``I^g``, ``U_*`` | `earth_current`, `star_voltage` (zero for generalized generator) |
| ``P^{poc},Q^{poc}`` | `p`, `q` |
| ``P^{internal},Q^{internal}`` | `p_internal`, `q_internal` |
| ``P^{series},P^{ground}`` | `series_loss`, `ground_loss` |
| ``\Delta U^{012},E^{012},I^{012}`` | `voltage_sequence`, `emf_sequence`, `current_sequence`, or `nothing` |
| Active-power conservation residual | `power_balance_error` |

Unpublished solver outcomes expose NaN measurements under the package's status
contract. A result serialization schema is not proposed by these tables.

## 5. Baseline compatibility and issues found in the upstream review

The revision-pinned [source register](index.md) identifies the reviewed files.
The following differences must be addressed explicitly when deriving upstream PRs.

| Existing contract or statement | Draft treatment and required action |
|:--|:--|
| `generator` means dispatchable WYE, with fixed injection described as a negative load | Keep legacy records unchanged. New source laws and fixed controls are explicit opt-in assumptions; do not redefine old data by inference. |
| Generator power arrays are per phase; `i_max` is per conductor and its neutral entry can be optional | Draft power arrays are per ordered port, with separately named total fields. Draft `i_max` has all ``m`` entries if present. A legacy wye phase-only current vector maps exactly to `i_port_max`, since its phase conductor currents equal its port currents; a full phase-plus-neutral vector maps to `i_max`. Do not invent a missing neutral rating. |
| `cost` is a per-phase vector; the objective explicitly rejects a scalar | Use a new `cost_total` field. Equal phase prices map exactly for a closed-return generator. Unequal prices cannot be collapsed without changing the objective. |
| Source PCC power is often accounted for by summing phase-to-neutral port powers | With earth return, use the complete conductor sum. For equal phase price ``c``, the old phase-priced cost minus total-priced cost is ``c\,\mathfrak R(U_n(I^g)^*)/1000``; it need not vanish. |
| The source fixes bus terminal voltages and implicitly grounds its neutral | A finite-impedance source fixes internal EMF only. Retain ordinary PCC voltage bounds and make the internal-star ground explicit. |
| Source neutral return is stated as ``-\sum I_k`` | Retain that statement only for the closed-return model. The new source requires ``I_n=-\sum I_k-I^g``. Ideal-voltage equivalence does not establish neutral-current equivalence. |
| Grounding page says generators never ground, only buses declare optional grounding, and sources always connect to ground | Add the source-generator exception and its internal grounding location, including `OPEN`; keep ordinary generalized generators without earth. |
| Bus text describes KCL exceptions at grounded/source-fixed terminals while also describing free balancing currents | State one explicit KCL convention: retain current balance with the appropriate physical free balancing current. Never omit balance at a rated custom source because its bus happens to contain source metadata. |
| Notation prose mixes into-element and into-bus current descriptions; the displayed bus residual subtracts generation | Define orientation at each interface. These drafts make component injection positive and state how to change sign in a withdrawal-positive residual; any upstream editorial correction must preserve existing equations' intended signs. |
| Source page permits exactly one source | Reference sufficiency and multiple source placement are problem-level scope decisions, not consequences of introducing the component data type. |
| Explicit neutral/no Kron reduction is a foundational principle | Keep network neutrals and conductor matrices explicit. Use the source primitive for an internal neutral lead. Do not present a reduced-neutral port impedance as an explicit neutral-conductor model. |
| Objective refers to source/IBR cost arrays, but the reviewed source field table does not list `cost` | Reconcile source fields and objective in an upstream issue; do not infer an undocumented legacy default while translating. |

The local builder already handles native/custom ownership, retires unused native
currents, preserves native bus limits when source physics is replaced, and adds
custom complete-PCC injection/cost ledgers. These are integration requirements,
not new physical equations. They must remain covered by regression tests as the
engine and upstream specifications evolve.

## 6. Suggested upstream contribution units

1. **Generalized port-source primitive and shared operating data.** Propose the
   new component identity or explicit subtype, field tables, topology rules,
   voltage laws a–e, capability overlay f, measurements, and examples. Explain
   that the canonical circuit is shared while source/control assumptions differ.
2. **Dedicated neutral/earth source.** Add its conductor primitive, explicit
   grounding, power correction, and ratings; update grounding, bus KCL/limits,
   source-reference scope, and objective semantics together.
3. **Schema and migration for each normative increment.** The upstream contribution
   guide locates the canonical schema in
   [distribution-system-opt/dsopt-schema](https://github.com/distribution-system-opt/dsopt-schema).
   Coordinate the prose/examples PR with its schema PR and version policy.
   Do not tag the draft JSON as conforming to an old schema.
4. **Explanatory implementation and evidence.** Move rectangular derivations,
   exact-limit realization, measurement examples, and benchmark evidence into
   implementation/tutorial sections while keeping solver options outside the
   foundational data contract.

The upstream guide asks normative proposals to begin with an issue or discussion.
That is a future upstream contribution step, not a prerequisite for drafting
these files locally. No issue, PR, schema release, or publication has been created.

Before proposing a schema, encode conditional required fields, dimensions,
configuration/port incidence rules, enum sets, mutually exclusive per-port/total
representations, finite numeric values, positive magnitude domains, passivity,
and bus-terminal references. Several of these require semantic validation beyond
basic JSON Schema. New terms should be added to upstream notation/nomenclature
and glossary: port EMF, source star, PCC complete power, earth return, power
measurement location, and voltage-source shape.

## 7. Evidence boundary

The [scientific model](../../components/generalized_generator.md) and
[literature register](../../components/generalized_generator_evidence.md) discuss
prior work and physical limitations. The [fourteen-case network tutorial](../../tutorials/generalized_generator_tradeoffs.md)
compares model choices under unbalance. Unit tests independently check circuit
power/grounding identities, every source law and capability field, SI/per-unit
conversion, exact zeros, singular impedances, connection types, and native/custom
ownership. These support the prototype's declared scope, not universal device
fidelity or comparative speed claims.

Open source-neutral leads, soil networks, dynamical machine
states, physical shunts/filters, automatic PV/PQ switching, and controller limit
priority remain outside the implemented model. A mathematical feasible point is
not proof that a particular machine or converter can realize the selected control
freedoms or withstand the declared unbalance.


### Delta constraint lowering

The incidence equations above also apply to the rank-two delta matrix. Keep
three complex winding currents to retain circulation. Never invert ``C`` or
replace winding currents with line currents. For an exact ideal loop
``\mathbf1^{\mathrm T}Z=0``, validate template closure, impose only the first
two fixed-phasor equations (or only the second-to-first complex ratio), and
reconstruct the third winding voltage. For mode e, combine ideal-loop closure
and any exact internal sequence zeros in the real magnitude-space nullspace
before forming ratios. This removes dependent equations at construction time.
No impedance inverse, internal voltage variable, or extra cycle row is needed.


An exact sequence-current restriction can also follow from a delta voltage law
and its series matrix. Before adding affine delta equalities, the implementation
checks dependence within this component's small row set using normalized
elimination (relative coefficient tolerance ``10^{-12}``). It retains the original
scaled physical rows in the optimization model. This is local algebraic cleanup,
not a proof of full-network Jacobian rank or numerical robustness for every case.
