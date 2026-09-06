# Generalized generators: scientific model

> **Kind:** Component model · **Maturity:** research prototype · **Direction:** forward · **Temporal:** single-snapshot, fundamental frequency

PowerOptLab implements optional generator components on BMOPFTools' staged
engine. The [implementation/API guide](generalized_generator_api.md) distinguishes
the implemented modes and connections from later extensions in this scientific
design. The [network tutorial](../tutorials/generalized_generator_tradeoffs.md)
runs the model choices on the same loaded, unbalanced four-wire feeder.
The research question is: **when do assumptions about source stiffness, voltage
regulation, sequence response, and neutral connection change a network decision?**
Relevant decisions include voltage operability, unbalance mitigation, operating
envelopes, and inference of source impedance from terminal observations.

Start with a multiconductor voltage source behind a linear series impedance,
including an exact zero-impedance case. Include a dedicated source-generator
variant with an explicit neutral grounding point and earth-current port. Compose connection topology, source
constraints, operating controls, and capability limits separately. Retain the
advanced IBR as a distinct component with converter-specific physics.

The [worked examples](../tutorials/generalized_generator_models.md) explain the
distinctions numerically. The [evidence register](generalized_generator_evidence.md)
records the reviewed literature, its scope, and the boundary between established
models and this proposal. Review date: 6 September 2026.

## Scientific positioning

Voltage-behind-impedance models and unbalanced PV formulations are established.
A directly relevant precedent specifies total active power and positive-sequence
terminal voltage magnitude, with a positive-sequence internal EMF and passive
negative/zero-sequence responses. See Fernandes et al., Section 3.1.2 and Fig. 3
([2019](https://doi.org/10.1049/iet-gtd.2018.6176)).

Recent work still develops phasor models for differing negative-sequence controls:
Yang et al. compare augmented rectangular power-flow models against PSCAD
([2024](https://doi.org/10.1049/gtd2.13108)). Thus the useful contribution is a
well-specified, reproducible family of models and evidence about which assumptions
matter, rather than a claim that generalized generator equations are new.

Candidate research outputs are:

- quantify the decision error from replacing finite source impedance with an
  ideal source, or terminal positive-sequence regulation with fixed phase voltages;
- find cases where independent sequence bounds overestimate realizable phase or
  neutral capability, and compare the resulting dispatch against a declared controller;
- distinguish identifiable source/grounding hypotheses using controlled injections;
- compare equivalent formulations for exact ideal connections and relative-angle
  constraints using residuals, Jacobian rank, solve reliability, and sensitivities.

These are proposed experiments. Neither novelty nor superiority is established
by this focused review.

## Four independent modeling choices

| Layer | Specifies | Examples |
|:--|:--|:--|
| Connection | Physical terminals, returns, winding orientation, grounding | phase-neutral, phase-phase, wye, delta, center tap |
| Electrical primitive | Relation between internal and terminal voltages/currents | ideal connection, full series impedance, later explicit shunt/two-port |
| Source and operating law | Admissible EMF shape and what determines operation | fixed phasor, rotating template, fixed P/Q, terminal voltage regulation |
| Capability | Permissible operating region at a named location | phase/neutral current, P/Q, EMF, sequence limits |

Angle reference selection belongs to the assembled problem. A physical fixed
phasor can supply that reference, but an arbitrary coordinate choice must not
silently become a stiff three-phase source or a neutral-ground connection.

## Electrical model and sign conventions

Use RMS fundamental phasors. Positive current and positive complex power mean
injection from the device into the network. Let ``v`` collect external conductor
potentials and let the real incidence matrix ``C`` map them to oriented port
voltages. Each column of ``C`` has a +1 at the outgoing terminal and a -1 at its
return; a physical ground can be represented explicitly with known zero potential.
For port currents ``i`` and port voltages ``u``,

```math
u=C^\mathsf{T}v,\qquad j=Ci,\qquad
S_{\mathrm{poc}}=\sum_k u_k\overline{i_k}.
```

Stamp ``j`` once into network KCL. This construction conserves complex power:
``\sum_t v_t\overline{j_t}=S_{\mathrm{poc}}``. It covers arbitrary terminal
pairs without assuming three phases or an earth return. Internal winding nodes
and explicit grounding branches extend the same circuit construction. The matrix
includes **all** physical ports, including earth where connected. Its projection
onto just the network's phase/neutral conductors need not have zero current sum.

For a series-only primitive in compatible independent port coordinates,

```math
e-u=Zi,\qquad Z=R+\mathrm{j}X,
```

where ``e`` is internal source voltage. Its rectangular stamp is linear:

```math
e_r-u_r=Ri_r-Xi_i,\qquad e_i-u_i=Xi_r+Ri_i.
```

There is no need to invert ``Z``. Full matrices allow mutual coupling and unequal
sequence impedances. For a physical passive impedance require
``H=(Z+Z^\mathsf{H})/2\succeq0`` within a documented numerical tolerance;
check reciprocity separately when it is part of the declared physical model.
Do not equate this single-frequency passivity check with dynamic stability.

The source and terminal power ledgers must remain distinct:

```math
S_{\mathrm{int}}=\sum_k e_k\overline{i_k},\qquad
S_{\mathrm{int}}-S_{\mathrm{poc}}=i^\mathsf{H}Zi,\qquad
P_{\mathrm{loss}}=i^\mathsf{H}Hi\geq0.
```

Name the location of every power setpoint, bound, and cost. Terminal active
power is not shaft power: copper, core, mechanical, and excitation losses require
their own assumptions. A virtual impedance is a control relation and must not
automatically enter the physical dissipation ledger.

### Neutral impedance with one return path

For a wye source with phase currents ``i`` and a series neutral, define
``K=[I;-\mathbf{1}^\mathsf{T}]``. If ``Z_c`` is the primitive conductor impedance,
the phase-to-neutral drop uses ``Z=K^\mathsf{T}Z_cK``. In the diagonal case,

```math
Z=\operatorname{diag}(z_a,z_b,z_c)+z_n\mathbf{1}\mathbf{1}^\mathsf{T}.
```

This reduction assumes a single corresponding neutral return, with no intermediate
shunt/earth paths. Otherwise retain the actual internal circuit. A physical
ground bond, a finite grounding impedance, and a floating star point are different
models. Never silently fix the PCC neutral to ground.

### Dedicated source generator with neutral grounding

`SourceGenerator` shares the electrical/source-law primitives and has an explicit internal neutral/star node ``v_*`` and a grounding
connection. Its outward external conductor currents are
``j=(I_a,I_b,I_c,I_n)``; its earth current ``I_g`` is positive from the source
neutral toward earth. Conservation is

```math
I_a+I_b+I_c+I_n+I_g=0.
```

Consequently **do not constrain** ``I_a+I_b+I_c+I_n=0`` in this model. A nonzero
sum is supplied or returned through the grounding port, not a violation of KCL.
For three-phase sequence quantities the correct relationship is
``3I_0+I_n+I_g=0``. Do not reuse a helper that sets ``I_n=-3I_0`` here.

Let ``e=(E_a,E_b,E_c)`` be internal EMFs measured to the source star point, and
let ``v=(V_a,V_b,V_c,V_n)`` be network terminal potentials. The four-conductor
series primitive is

```math
v_{\mathrm{int}}=\begin{bmatrix}v_*+E_a\\v_*+E_b\\v_*+E_c\\v_*\end{bmatrix},
\qquad v_{\mathrm{int}}-v=Z_c j.
```

For diagonal ``Z_c``, the neutral relation is ``v_*-V_n=z_n I_n``. Full mutual
coupling remains allowed without substituting a phase-current sum for ``I_n``.
The grounding connection is a separately declared choice:

| Grounding connection | Equations | Meaning |
|:--|:--|:--|
| Open | ``I_g=0`` | No earth return at this component |
| Finite impedance | ``v_*-V_g=z_g I_g`` | Ground electrode/bond equivalent at the declared frequency |
| Ideal bond | ``v_*=V_g``; ``I_g`` determined by KCL | Exactly grounded star point, no tiny impedance |

Use ``V_g=0`` for the ideal remote-earth reference. If ground potential rise or
shared soil paths matter, expose an actual earth-network terminal instead. A
PCC neutral-ground bond and an internal source-star bond can produce different
currents when ``z_n\ne0``; the data model must identify the grounding location.
Omitting a network neutral lead is also different from an ideal neutral lead.

Maintain an earth-current handle even for an ideal bond, preferably as a
reconstructed KCL expression when no separate variable is necessary. Stamp
network conductor injections, local internal KCL, and the earth return exactly
once. Do not combine an implicit engine ground injection and a second custom
ground branch for the same physical bond. With zero conductor impedance,
internal potentials alias the corresponding network potentials; earth current
still permits a nonzero four-conductor sum.

For remote earth at zero, the power identity for this series-plus-grounding
topology is

```math
S_{\mathrm{emf}}=\sum_{k\in\{a,b,c\}} E_k\overline{I_k},\qquad
S_{\mathrm{emf}}-\sum_{t\in\{a,b,c,n\}} V_t\overline{I_t}
=j^\mathsf{H}Z_cj+v_*\overline{I_g}.
```

For a finite grounding impedance the final term is ``z_g|I_g|^2``. For a nonzero
earth terminal potential also account for exported earth-port power
``V_g\overline{I_g}``. Report phase, neutral, and earth currents independently,
with separate lead/bond ratings and losses. A source preset can default to a fixed
internal phasor, but grounding can also compose with the other source laws.

For this source, define PCC power using the complete conductor sum above. A
phase-to-neutral power sum alone omits a term when both ``V_n\ne0`` and
``I_g\ne0``:

```math
S_{\mathrm{poc}}=\sum_{k\in\{a,b,c\}}(V_k-V_n)\overline{I_k}
-V_n\overline{I_g}.
```

Name this measurement convention in power setpoints and results. Include any
explicit earth-terminal power separately when its potential is not zero.

### Exact ideal connections

When the **entire effective series primitive is exactly zero**, use ``e\equiv u``
as expression aliases. Create no additional internal voltage variables or network
buses and no duplicate equality ``e=u``. Internal measurements remain available
through those aliases. The implementation goes further: it reconstructs
``e=u+Zi`` as affine expressions even for nonzero impedance. This retains exact
internal measurements and quadratic power expressions without introducing local
voltage variables or mutating the engine's network bus dictionary.

A singular nonzero matrix is still valid in the IVR stamp. A zero diagonal entry
alone does not justify merging nodes if mutual terms remain. Eliminate only
structurally proven ideal connections; arbitrary nullspace directions need not
correspond to physical node mergers. Never replace zero by a tiny impedance.
Treat a user-requested near-zero reduction as a separate approximation with an
error budget. A nonzero neutral impedance also prevents the all-zero reduction.

For a PQ injection with terminal setpoints and no internal limits or losses of
interest, the internal voltage is merely reconstructible as ``e=u+Zi``. Adding
an unconstrained internal voltage cannot change the terminal feasible region.
This should be an explicit simplification, not evidence that impedance is
irrelevant to voltage-controlled sources.

## Mapping the requested source modes

Let ``n`` be the number of independent, nonzero source-port phasors and let
``\phi_k`` be specified relative phase offsets. The degrees of freedom below
count source voltage only, before circuit coupling, active limits, power control,
or any global angle reference.

| Requested mode | Proposed descriptive name | Source-voltage law | Real voltage freedoms |
|:--|:--|:--|:--|
| a | No voltage law | ``e`` unrestricted by a source law; omit it if unused | ``2n`` if retained |
| b | Fixed phasor | ``e_k=\bar e_k`` | 0 |
| c | Rotating fixed template | ``e_k=\bar m_k\exp(\mathrm{j}(\delta+\phi_k))`` | 1 |
| d | Rotating common magnitude | ``e_k=m\exp(\mathrm{j}(\delta+\phi_k))`` | 2 |
| e | Rotating phase magnitudes | ``e_k=m_k\exp(\mathrm{j}(\delta+\phi_k))`` | ``n+1`` |
| f | Sequence capability overlay | Bounds on named voltage/current sequence components | Depends on underlying source law and active bounds |

For d, allow specified winding ratios ``\rho_k>0`` as an extension:
``e_k=\rho_k m\exp(\mathrm{j}(\delta+\phi_k))``. Equal ratios give exactly the
requested common magnitude. For three-phase positive sequence use
``\phi=(0,-2\pi/3,2\pi/3)``. Unequal magnitudes with these offsets generally
have nonzero negative and zero sequences; fixed 120-degree separation does not
by itself imply a balanced phasor.

Mode a becomes **PQ** only with specified P and Q. Bounds on P/Q instead produce
a dispatchable injection for OPF. Specify whether powers are aggregate or per
winding; aggregate P/Q alone may leave phase allocation free.

Mode b is a stiff reference source **at the location of its fixed phasor**.
With nonzero impedance, the internal phasor is fixed while PCC voltages can vary.
Its P/Q are solved subject to declared capability. Fixing its active power as well
can overconstrain a power-flow problem. Angle reference and mismatch allocation
are distinct roles; distributed slack is a separate problem-level policy.

Mode c with fixed P is a useful idealized PV-like preset, but it is not the unique
unbalanced PV generalization. At zero impedance it fixes every terminal magnitude
and all relative angles, which assumes more regulation than fixing just positive
sequence magnitude. Behind impedance it fixes internal excitation magnitude,
which does **not** maintain terminal voltage.

### Add an explicit terminal PV preset

For a simple round-rotor synchronous-machine-inspired model, use mode d with a
balanced internal EMF, appropriate sequence impedance, and

```math
P_{\mathrm{poc}}=P^\star,\qquad |U_1|^2=(V_1^\star)^2.
```

Internal magnitude and angle can then adjust. Negative/zero-sequence terminal
responses follow the circuit and grounding. This is the useful PV precedent in
Fernandes et al.; it is not a full machine model. Other valid presets can regulate
an average magnitude, a selected phase, or line-to-line voltage. State the sensor
and actuator assumptions instead of naming every version simply `PV`.

Per-phase magnitude regulation would require suitable independent actuators and
enough remaining freedoms. It should not be attributed to an ordinary common
field winding without a physical model supporting it.

Reactive, current, or excitation limits require a declared policy. In fixed-control
power flow, a hard voltage target may become infeasible, or a documented limiter
can release it and enter a saturated mode. MATPOWER's optional reactive-limit
enforcement provides a familiar PV-to-PQ example
([manual](https://matpower.app/manual/matpower/ACPowerFlow.html)). In OPF, a
voltage or excitation setting can instead be an explicit decision variable.
Do not let an unspecified optimizer choose the controller response. Droop,
limit priority, hysteresis, and island frequency balance are later control/problem
extensions, not implied by a free angle or magnitude.

## Sequence capability and its limits

For ordered three-phase, phase-to-the-same-neutral RMS quantities, choose

```math
T=\frac13\begin{bmatrix}1&1&1\\1&a&a^2\\1&a^2&a\end{bmatrix},
\quad a=\exp(\mathrm{j}2\pi/3),\quad
x_{012}=Tx_{abc},\quad Z_{abc}=T^{-1}Z_{012}T.
```

Use distinct ``Z_0,Z_1,Z_2`` when a diagonal sequence model is appropriate.
For a general asymmetric phase primitive, ``TZT^{-1}`` need not be diagonal.
The transform is a coordinate change, not a reason to decouple the physics.

Mode f should compose with a-e. Specify whether limits apply to ``E_s`` internally,
``U_s`` at the PCC, or source/PCC current when shunts make those currents different.
Typical constraints include

```math
\underline V_s^2\leq |U_s|^2\leq\overline V_s^2,\qquad
|I_s|^2\leq\overline I_s^2,\qquad s\in\{0,1,2\}.
```

Make current lower bounds optional and zero by default. A minimum sequence
injection requirement is a control/service obligation, not a universal capability
requirement. When any upper magnitude bound is exactly zero, impose the two
linear equalities on its real/imaginary parts, not ``x_r^2+x_i^2\leq0`` with a
zero gradient at the solution. Do not also stamp a redundant zero lower bound.

Retain phase current ratings and neutral current ratings alongside sequence
limits. The phase currents combine all sequences. For the chosen normalization,
``I_n=-3I_0`` for a wye port with one neutral return and no earth path. For the
dedicated grounded source use ``I_n=-3I_0-I_g`` instead. A three-wire connection
with no earth return requires zero external zero-sequence current; it does not prohibit internal delta
circulating current. Line-current and winding-current limits are different.

If a voltage unbalance factor is used, specify ``|U_2|/|U_1|`` and a strictly
positive lower bound on ``|U_1|``. Prefer the quadratic relation
``|U_2|^2\leq\epsilon^2|U_1|^2`` to division, with the denominator domain retained.
Independent magnitude bounds alone do not enforce current priority, converter
DC/ripple feasibility, thermal dynamics, or fault ride-through behavior.
Recent asymmetric-disturbance controller research makes that distinction relevant
([Liu et al., 2026](https://doi.org/10.1109/TPEL.2025.3632684)).

## How much output impedance is enough?

The minimum useful implementation is a constant full series matrix at a declared
frequency, with scalar/diagonal and sequence-data constructors. This captures
source stiffness, neutral drop, coupling, losses, and a useful approximate
machine sequence response. EPRI's source and generator documentation illustrates
both Thevenin modeling and distinct negative-sequence behavior
([source](https://opendss.epri.com/Modeling.html),
[generator dynamics](https://opendss.epri.com/GeneratorDynamicsModel.html)).

Choose the model order for the study:

| Intended study | Useful next detail | What a constant series matrix omits |
|:--|:--|:--|
| Fundamental-frequency dispatch and stiffness sensitivity | Full series matrix and correct grounding | Voltage-dependent losses or hidden internal limits |
| Filter currents, internal ratings, reactive shunt exchange | Explicit shunt or LCL/two-port circuit | Source and PCC currents can differ |
| Machine excitation capability or saliency | Machine-specific algebraic laws, field/stator limits, saturation | ``d/q`` anisotropy and operating-point dependence |
| Electromechanical transients and faults | Consistent transient/subtransient parameters plus rotor/exciter/governor states | Time evolution and limiter transitions |
| Converter impedance stability or harmonics | Frequency-dependent coupled impedance/admittance or dynamic model | Control bandwidth, frequency coupling, resonances |

Steady-state, transient, and subtransient reactances are not interchangeable.
Store parameter regime, frequency, base values, and provenance. Rotating-frame
machine models explicitly include field/damper dynamics and saturation
([MathWorks equations](https://www.mathworks.com/help/sps/ref/synchronousmachinesalientpole.html)).
Sequence/dq impedance analysis for converters also has frequency-coupling issues
([Rygg et al.](https://arxiv.org/abs/1605.00526)). These are reasons to preserve
an electrical-primitive interface, not to implement every model order now.

A linear circuit can sometimes be reduced to a terminal Thevenin equivalent for
a fixed operating condition. That reduction does not preserve unrecorded internal
current limits or controller measurements automatically. Keep LCL and DC physics
in the advanced IBR until a demonstrated shared use warrants extraction.

## Topology and numerical formulation

Use terminal labels and oriented winding pairs as the core data; phase names and
sequence transforms are optional metadata. Single-phase line-neutral and
line-line ports are two-terminal instances. Split phase uses two leg-to-center
ports with a 180-degree template and the actual center-tap return. Generic
multiphase templates use their declared offsets, not a hard-coded 120 degrees.

Winding voltages must satisfy circuit cycles. For a closed delta,
``u_{ab}+u_{bc}+u_{ca}=0``. Arbitrary fixed unequal magnitudes at exact 120-degree
offsets cannot generally be imposed on these three winding voltages. The source
DOF table assumes independent ports; rank/cycle validation must precede applying
it to delta or other redundant winding descriptions. Reject unsupported topology
in the first implementation rather than silently reducing it to wye.

For relative-angle source laws, avoid tangent equations and all-pairs redundant
angle constraints. A rectangular implementation can use a shared complex
``w=c+\mathrm{j}s`` with ``c^2+s^2=1`` and
``e_k=m_k\exp(\mathrm{j}\phi_k)w``. Positive magnitudes preserve phase orientation.
For common-magnitude mode d, use one free complex phasor ``z`` and
``e_k=\rho_k\exp(\mathrm{j}\phi_k)z`` with bounds on ``|z|``; this avoids an
unnecessary amplitude/angle factorization. For c, ``e_k=\bar e_k w``.
Do not impose these parameterizations and their equivalent pairwise laws together.

Fix one angle per independent rotational symmetry of the assembled energized
model, usually one per coupled island, unless a physical reference already fixes
it. Electrically decoupled phase subsystems can have additional symmetries. A
rectangular reference ``\Im(z)=0,\ \Re(z)>0`` is implemented using a declared
positive operating lower bound. Do not pin the angle of every generator: relative
rotor/source angles carry physical power transfer. Floating potential/common-mode
freedoms are distinct from angle freedoms and require their own circuit treatment.

The attached Geth-Pacaud-Heidari paper motivates careful reference choices,
phase-consistent initialization, scaling, and removal of redundant constraints.
Its Table 2 fixes absolute angles in `3vafix`; requested mode e leaves a common
rotation free until a separate reference is supplied. Its `vafixseq` retains more
relative-angle freedom than e. Its 3-by-3 benchmarks do not validate explicit
neutral or split-phase extensions. See the evidence register for exact scope and
normalization notes.

Use SI inputs and the engine's voltage/current bases internally. Compute
``Z_b=V_b/I_b`` from those bases; do not guess whether the supplied power base is
per-phase or three-phase. Initialize positive-sequence sources with the expected
rotation and respect transformer vector groups. Report physical-unit residuals
as well as solver termination. Local feasibility is not a stability or global
optimality certificate.

## Integration with BMOPFTools and the advanced IBR

The inspected dependency is BMOPFTools commit
`5b51d2f361dab91bd7c16711019584407da79ed8`, pinned in both project files.
The installed source confirms:

- `OpfBuildSpec` / `OpfDeviceBuilder` support per-component ownership for both
  `:generator` and `:voltage_source`; the public callback receives `(ctx, ids)`.
- Native generators provide power/current capability and KCL stamping. Native
  voltage sources fix voltage phasors and ground their source-bus neutral.
  A custom source must implement its own declared grounding semantics.
- Current variables are declared before device ownership is resolved. Native
  results, costs, initializations, and diagnostics may still consult native data.
  Ownership of a physics callback alone does not complete result integration.
- `opf_model`, `opf_bases`, `opf_object`, semantic current/voltage keys, and
  `add_terminal_injection!` provide the public integration surface.

The implemented `GeneralizedGenerator <: AbstractDevice` and
`SourceGenerator <: AbstractDevice` types live in
`src/components/`, using the existing `validate_device`, `stamp_device!`,
`link_device!`, and `extract_device` lifecycle. The vocabulary below describes
the architectural separation; exact constructor names and fields are in the
[API guide](generalized_generator_api.md):

```text
GeneralizedGenerator
  id, terminal_connection, electrical_primitive
  source_law, operating_law, capability
  parameter_provenance

SourceGenerator
  id, terminal_connection, conductor_primitive
  source_neutral, grounding_connection, optional_earth_terminal
  source_law, operating_law, capability, parameter_provenance

SourceLaw = NoVoltageLaw | FixedPhasor | RotatingTemplate |
            CommonMagnitude | PhaseMagnitudes
OperatingLaw = Dispatchable | FixedPQ | FixedPVoltageTarget
Capability = phase + neutral + power + internal-voltage + sequence limits
```

Keep the network solve and objective in the existing problem layer. A new device
can be added through the staged hook; replacement of an existing native record
uses the typed ownership seam. Avoid a second generator-specific network solver.

Reuse native current handles only when their orientation, arity, and semantics
match. Otherwise retire unused placeholders through public APIs and provide
authoritative custom outputs/costs, following the controlled-inverter ownership
pattern. Prevent double KCL stamping, duplicate generation cost, and native
placeholder results masquerading as physical measurements. Audit native source
warm starts and solution validators for assumptions about fixed PCC voltage.

The advanced IBR already has primitive conductor impedance and neutral-drop
logic, Fortescue helpers, internal voltages, and a balanced common-magnitude EMF
mode. Extract shared port maps, transforms, linear drop stamps, magnitude-limit
lowering, unit conversion, and measurement/result conventions only with regression
evidence. Preserve IBR switching, DC-link, topology-specific feasibility, ripple,
loss, and control machinery. Sharing a voltage-behind-impedance primitive does not
make the generalized generator a converter capability model.

## Implementation increments and acceptance evidence

The first implementation covers all five voltage laws, optional sequence
capability, hard terminal-PV control, open/finite/ideal source-star grounding,
independent single-phase/split-phase/wye port connections, and stateless
multi-period composition. The saturation policy in increment 3, closed delta
cycles, explicit earth-network terminals, higher-order primitives and IBR
refactoring remain future work. The table preserves the development/evidence
plan rather than implying that every item below has been delivered.

| Increment | Deliverable | Required evidence before claiming support |
|:--|:--|:--|
| 1 | Wye/single-phase ports; dedicated source with neutral/earth ports; exact ideal/full series primitive; a-e; P/Q and conductor/ground limits | Native PQ/reference equivalence, internal/PCC/earth power balance, nonzero four-conductor current sum with earth return, analytic voltage drops, correct zero-impedance variable count, mixed native/custom ownership and costs |
| 2 | Sequence transforms/limits; terminal positive-sequence PV preset | Pure-sequence normalization, phase overload counterexample, grounding/neutral cases, fixed-excitation versus terminal-regulation comparison |
| 3 | Fixed-control saturation behavior and source uncertainty studies | Published limit policy, PV-to-limited transition, objective-invariance checks where the response is intended to be unique, physical residuals and local rank |
| 4 | Split phase and further connection presets | Center-tap KCL, line-line versus leg voltage checks, delta cycle constraints and winding/line current distinction |
| 5 | Shared primitive extraction and optional higher-order devices | Unchanged advanced-IBR results over its existing regression suite; a named research case requiring each extra physical detail |

Design the port representation for arbitrary configurations in increment 1,
even if support is initially narrower. Keep f as a composable capability object
from the start; implementation can follow the basic source laws.

Validation should include positive/negative/zero-sequence excitation, floating
versus grounded neutral, open/finite/ideal source grounding, internal versus PCC
ground bond, independent neutral/earth current reporting, finite versus exact-zero impedance, singular nonzero
matrices, island references, conflicting ideal sources, active bounds, and SI/pu
equivalence. For sequence-zero limits, check the nondegenerate linear stamp.
For equivalent formulations compare local Jacobian rank after reference removal,
not merely solver success. A zero-impedance collapse changes structure; derivatives
with respect to impedance should use a declared fixed-structure formulation.

Compare against BMOPFTools' native overlapping models and independent circuit
calculations first. Use OpenDSS for matching finite-impedance cases with documented
model/connection/limit choices; its nominal ideal source uses a tiny impedance,
so it is not an exact ideal-source oracle. Use a machine or EMT model only when
testing claims about that device's controlled equilibrium or dynamics.

For each published experiment preserve the source law, regulated quantity,
actuator allocation, grounding, impedance regime, capability location, sequence
normalization, solver versions, bases, tolerances, starts, ownership manifest,
physical residuals, and study scripts. Separate parameter assumptions, analytic
identities, numerical case evidence, and untested hypotheses.
