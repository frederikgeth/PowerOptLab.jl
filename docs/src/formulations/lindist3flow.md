# LinDist3Flow BMOPF prototype

PowerOptLab provides an independent LinDist3Flow model builder for normalized
BMOPF JSON. It does not depend on PowerModelsDistribution and does not turn the
BMOPFTools nonlinear OPF engine into a multi-formulation engine. BMOPFTools is
reused at the boundaries: parsing, public SI/per-unit coordinate preparation,
and nonlinear replay; PowerOptLab's
[`kron_reduce_bmopf`](@ref) supplies the neutral-reduced network representation.
The complete component equations are specified in the
[LinDist3Flow component model](lindist3flow_components.md), using the same
data/symbols/variables/equalities/inequalities/implementation organization as
the BMOPFTools mathematical-model documentation.

```julia
options = L3FOptions(validate_nonlinear=true, per_unit=true, s_base=1e6)
result = solve_l3f_opf(network, Clarabel.Optimizer; options)
```

BMOPF JSON is always supplied in SI and results are always returned in SI.
`per_unit=true` (the default) changes only the optimization coordinates;
`per_unit=false` builds the identical formulation in raw SI. PowerOptLab uses
BMOPFTools' public classic-base preparation to scale a private working copy,
then constructs its own JuMP model. It does not invoke the nonlinear OPF
component builders. SI/per-unit equivalence is regression-tested across lines,
ZP loads, shunts, generators, sources, and fixed regulators.

An explicit-neutral input is reduced on a deep copy by default. Pass
`L3FOptions(kron_reduce=false)` to require callers to supply an already reduced
case. `check_l3f_applicability` performs the same boundary checks without
constructing a JuMP model. Unsupported data is represented by typed findings
with stable codes; it is never silently dropped.

## Implemented slice

Version `0.1-prototype` builds a continuous affine LP or SOCP for radial AC islands with
exactly one fixed-voltage source per island. It supports neutral-reduced series
lines without line shunts, grounded-wye, single-phase and delta constant-power,
constant-impedance, and pure ZP ZIP loads, constant-power generators, fixed bus
shunts, ideal fixed-ratio single-phase
transformers and ANSI A/B autotransformer regulators, retained phase-to-ground
voltage-magnitude bounds, linear per-channel energy costs, and reverse power
flow. Apparent-power bounds are native second-order cones. Ampacity bounds use
the declared fixed reference voltage and are also native second-order cones.
Fixed ideal three-phase `open_delta_regulator` banks are supported on three
retained phase terminals, including AB/CB and its two cyclic permutations.

Radiality is checked **per conductor**, not per bus: the graph whose nodes are
`(bus, terminal)` pairs must be a forest. A three-unit single-phase regulator
bank on one bus pair — the wye-bank topology of the IEEE 13, 34 and 123 feeders
— is therefore admissible, because its units occupy disjoint conductors. Two
devices sharing a conductor between the same buses remain non-radial.

### Where the line equation comes from

A series branch obeys ``v_j = v_i - Z i`` exactly. Forming the outer product and
keeping its diagonal,

```math
\operatorname{diag}(v_jv_j^H)=\operatorname{diag}(v_iv_i^H)
-2\operatorname{diag}\!\big(\Re(Ziv_i^H)\big)
+\operatorname{diag}(Zii^HZ^H).
```

LinDist3Flow drops the last term — the loss term, second order in current — and
substitutes ``i_\psi^*=s_\psi/v_{i\psi}`` with the *reference* phasors in the
voltage ratio. Writing ``\Gamma_{\phi\psi}=\bar v_\phi/\bar v_\psi``, the
result is affine in ``w=|v|^2`` and the branch powers:

```math
w_j = w_i - M p_{ij} - N q_{ij},\qquad
M = 2\Re(\overline Z \odot \Gamma),\quad
N = -2\Im(\overline Z \odot \Gamma).
```

``p,q`` are sending-end branch powers. Because the dropped term is exactly what
distinguishes the two ends, the same variable serves as the receiving-end power
and the nodal balance is lossless. On IEEE 37 the omitted losses are about
**1.95 % of feeder load** (see [accuracy](#Accuracy-and-when-it-degrades)).

### The linearization point

The default reference is the **flat, no-load profile**: the source phasors are
propagated outward through the radial topology, applying each transformer or
regulator ratio but **no line drop**. Every bus therefore sits at its nominal
magnitude with the source's angles, which is the classical LinDist3Flow
linearization point.

A better reference makes the closure tighter. Pass an explicit
`L3FReferenceState` or a phasor dictionary — for instance one taken from a
converged nonlinear power flow — as the `reference` argument.
`reference_policy` decides how the two interact: `:auto` prefers a supplied
reference, `:explicit` demands one, and `:source_propagated` always uses the
propagated profile. Source values, missing or zero phasors, and topology maps
are validated before model construction. Angles are coefficient data, not
decision variables.

The reference actually used is published as
`result.formulation["reference_provenance"]` and a content hash under
`reference_hash`, so a result can be traced to the point it was linearized at.

`cross_voltage_coefficients`, `winding_voltage_coefficients`,
`connection_power_map`, and `line_drop_coefficients` are pure coefficient
oracles. The connection map
implements ``H=\operatorname{diag}(\bar v)D^T
\operatorname{diag}(D\bar v)^{-1}`` for phase-to-phase and delta
devices. For each voltage-dependent load channel, the physical winding voltage
is represented by the same fixed-angle affine closure,
``\widehat{|Dv|^2}=c_D+a_D^T w``. Pure ZP laws are therefore exact affine
functions of that closure:

```math
p=p_{nom}\left(\alpha^P+\alpha^Z\widehat{|Dv|^2}/v_{nom}^2\right),\qquad
q=q_{nom}\left(\beta^P+\beta^Z\widehat{|Dv|^2}/v_{nom}^2\right).
```

Scalar or per-channel BMOPF coefficients and `v_nom` values are accepted.
Fitted ZIP coefficients are used verbatim and are not normalized. A coefficient
family with no fields defaults to constant power, following BMOPFTools semantics.

### The three approximations, stated plainly

The formulation is affine in ``w`` by construction, but it is a *linearization*
and it is worth being exact about where the error enters.

1. **Fixed-angle cross-voltage closure.** ``\widehat{v_\phi v_\psi^*}`` is the
   first-order Taylor expansion of ``\sqrt{w_\phi w_\psi}\,e^{j\Delta\bar\theta}``
   in ``(w_\phi,w_\psi)`` about the reference. It is exact at the reference
   point and exact for all ``w`` when ``\phi=\psi`` — so a **grounded-wye**
   device, whose incidence matrix is the identity, incurs no error here at all.
   A **delta** device does.
2. **Frozen channel-to-terminal split.** ``H`` is built from reference phasors,
   so the allocation of a device's channel powers across its terminals does not
   move with voltage. The *total* complex power is conserved exactly for any
   reference — ``\sum_\phi(Hs)_\phi=\sum_k s_k`` identically — but the
   per-phase split drifts. Again exact for wye, approximate for delta and for
   the open-delta regulator bank.
3. **Omitted series losses.** The nodal balance carries the same ``p+jq`` into
   the child and out of the parent.

Claims that the formulation contains "no approximation" apply only to the
algebra downstream of these three, not to the physics.

!!! note "Closed-delta channel powers are not uniquely determined"
    A closed delta's map ``H`` has rank 2: channel powers of the form
    ``s=c\,(D\bar v)`` are a circulating component the terminals cannot see.
    The published `pg`/`qg` for a three-channel delta device are therefore fixed
    only up to that component, even though its terminal injection is unique.

Fixed regulator and transformer settings preserve affine physics. Adjustable
tap intervals are rejected: there are no integer taps, McCormick envelopes, or
continuous voltage-ratio relaxations. Nonzero transformer leakage and no-load
admittance must be represented as separate supported network elements.

## Fixed regulator banks

`regulator_gain_matrix` implements the WYE, closed-delta, and open-delta gain
matrices from Table I of Bazrafshan, Gatsis, and Zhu, in that paper's own
convention ``v_n = A_{nm} v_m`` with ``n`` the source side — that is,
``v_{from}=A v_{to}``. ANSI type A uses the declared `tap_ratio` directly and
type B uses its reciprocal.

!!! note "The BMOPF regulator spec table has the ANSI mapping backwards"
    BMOPFTools' executable `_autotransformer_neff` gives ``n_{eff}=1/a`` for
    type B and ``n_{eff}=a`` for type A under ``V_{from}=n_{eff}V_{to}``, which
    is the self-consistent reading of `tap_ratio` as a regulated/source ratio.
    The prose table in the BMOPF regulator specification states the reciprocal
    pairing and should be corrected upstream. PowerOptLab follows the
    implementation.

The matrices are tested by the physical relation each bank must produce — that a
wye unit scales its own phase-to-ground voltage, that each closed-delta row is
an affine combination, and that an open-delta unit scales exactly its own
line-to-line voltage while leaving the common phase as the gauge — rather than
by restating the implementation's entries.

The optimizer currently stamps only `open_delta_regulator`, because that is the
only three-phase regulator-bank component in the BMOPF schema. Its voltage law
uses the fixed matrix ``A^{-1}``, and its terminal-power law uses the same
fixed-reference complex power transformation as other multi-terminal devices.
Two-winding apparent-power and reference-current ratings become native SOC
constraints; note that the shared (common) phase of an open-delta bank carries
both units' current and is deliberately unrated, matching the two-element
`i_max_from`/`i_max_to` shape BMOPF declares for the subtype. Closed-delta and three-unit WYE banks remain coefficient oracles
until BMOPF provides an unambiguous component representation; PowerOptLab does
not invent private JSON subtypes.

## Literature and OpenDSS evidence

Two feeder regressions run against OpenDSS as a nonlinear oracle.

### Modified IEEE 13

Transcribed from the IEEE PES Distribution Test Feeder Working Group archive and
cross-checked against the maintained OpenDSS deck; the impedance matrices are
Ω/1000 ft and reproduce the published Ω/mile values exactly on multiplication by
5.28. This feeder is the one whose regulator is a **bank of three single-phase
units on a single bus pair** — radial per conductor, three parallel edges on the
bus graph — so it is the direct test of the conductor-level radiality rule. It
also covers what IEEE 37 does not: buses retaining only one or two phases, a
single delta leg spanning two terminals, constant-impedance loads in both wye
and delta connection, an in-line transformer, and fixed shunt capacitors derived
from the two capacitor banks.

The modifications are each forced by the supported slice and mirrored in the
OpenDSS deck, so both describe identical physics: the 115 kV source and
substation transformer are dropped in favour of a nominal 4.16 kV source at bus
650; the regulator bank and XFM-1 are ideal at fixed taps 1.0625 / 1.0500 /
1.06875; and loads 692 and 611, constant-current in the original, become
constant power. Everything else is the published feeder.

Because the OpenDSS regulator bank is built as ideal single-phase transformers
at the same fixed taps rather than being bypassed, the regulated bus matches to
better than 1e-4 pu and the bank itself is independently validated — not just
the feeder downstream of it.

### IEEE 37

The larger regression independently transcribes the IEEE 37 topology and the
four impedance matrices used by [Bazrafshan, Gatsis, and Zhu
(PSCC 2018)](https://arxiv.org/abs/1901.04566). It applies the repository's
documented half-adjacent delta-to-grounded-wye conversion, removes line shunts,
uses the published fixed effective open-delta ratios `[0.9062, 0.9062]`
(stored as their reciprocals in a BMOPF Type-B component), and builds the
result as a 37-bus, 35-line BMOPF dictionary. The fixture is cross-checked
against the maintained [DSS-Extensions IEEE 37
deck](https://github.com/dss-extensions/electricdss-tst/tree/master/Version8/Distrib/IEEETestCases/37Bus).
The pinned source commits are recorded beside the test data transcription.

The tests separate three claims:

- Table I regulator matrices and fixed-ratio voltage/power transformations are
  exact algebraic oracles.
- The transformed feeder total, ``0.9828+j0.4804`` pu on a 2.5 MVA base, is a
  published data oracle. The public `MultiphaseVRs` workbook plus its conversion
  function yields phase totals ``[0.3636,0.2732,0.3460] +
  j[0.1774,0.1342,0.1688]`` pu, not the different per-phase split printed in
  the paper; the regression records the reproducible repository result rather
  than concealing that discrepancy.
- A nonlinear OpenDSS replay uses the same ideal regulator-secondary phasors,
  constant-power grounded-wye loads, series impedances, and zero line shunts.
  It checks all 108 downstream phase-voltage magnitudes against LinDist3Flow.
  OpenDSS's discrete `RegControl` and finite transformer impedance are excluded
  deliberately so this is a test of the continuous fixed-regulator contract.

Table III's import values (including the reported open-delta nonlinear value
1.0351 pu) are not equality assertions here. That OPF includes series losses,
variable tap selection, and a different transformer representation, whereas
this formulation is lossless and accepts fixed taps. Treating those objectives
as interchangeable would be a misleading test.

## Accuracy and when it degrades

The two error sources with a measurable size are the omitted losses and the
fixed-angle closure. On the IEEE 37 regression:

| Quantity | IEEE 37 | Modified IEEE 13 |
|---|---|---|
| Omitted series losses | 48.0 kW, **1.95 %** of load | 103.6 kW, **2.98 %** of load |
| Max phase-voltage error vs OpenDSS | **≈1.2 %** of nominal (108 terminals) | **0.51 %** of nominal (35 terminals) |
| RMS phase-voltage error vs OpenDSS | **≈0.6 %** | **0.37 %** |
| Signed mean error | positive | **+0.30 %** |

The signed mean is the informative one: every deviation on IEEE 13 is an
*overestimate*. Omitting series losses can only make the linear model optimistic
about voltage, so a sign flip in that column would mean a different error source
had appeared. The regression asserts the sign, not just the magnitude.

These sit inside the published LinDist3Flow envelopes. Sankur, Dobbe, Stewart,
Callaway and Arnold ([arXiv:1606.04492](https://arxiv.org/abs/1606.04492))
report under 1 % maximum relative voltage error on modified IEEE 13 and 37
feeders and an OPF objective within 0.2 % of the nonlinear optimum; Table I of
[arXiv:2210.08550](https://arxiv.org/abs/2210.08550) reports 0.007–0.01 pu
maximum deviation from a Z-bus power flow on IEEE 13, 0.008–0.02 pu on IEEE 123,
and 0.008–0.06 pu on the IEEE 8500-node feeder. The regression brackets its
OpenDSS comparison on **both** sides, so a change that makes the comparison
vacuous fails as loudly as one that makes it inaccurate.

Accuracy degrades where the underlying assumptions weaken:

- **High R/X and heavy loading.** The dropped loss term grows with ``|i|^2``, so
  long, heavily loaded feeders lose accuracy fastest — and always in the same
  direction: source import is understated.
- **Large voltage deviation from the reference.** The closure is a first-order
  expansion about ``\bar v``. [`l3f_reference_from_powerflow`](@ref) makes it
  exact at a real operating point — but see the caveat below before assuming
  that improves the answer.
- **Strong unbalance.** The fixed-angle assumption freezes the inter-phase angle
  differences at their reference values; delta devices and the open-delta bank
  additionally freeze their channel-to-terminal split.

Because losses are omitted, `objective=:cost` and `:source_import` systematically
understate the true cost of serving a load. Treat them as comparative rather
than absolute.

### A better reference is not automatically a better answer

[`l3f_reference_from_powerflow`](@ref) removes the closure error by linearizing
at a converged operating point. It does nothing about the omitted losses, and
the two errors can partially cancel. Measured against the same nonlinear power
flow:

| Case | Flat reference | Power-flow reference |
|---|---|---|
| IEEE 37 (wye, constant power) | 0.081 % max | 0.077 % max |
| IEEE 37 at 4x load | 2.19 % max | 2.08 % max |
| 5-bus delta constant-impedance feeder | 0.046 % max | 0.071 % max |
| the same feeder at 3x load | 0.359 % max | **0.558 % max** |

The last two rows are not a defect. Omitting losses makes LinDist3Flow
*overestimate* voltages; linearizing a ZP load at nominal makes it overestimate
the load's power draw, which pushes voltages back down. On a delta ZP feeder
those two biases partly cancel, and a reference that removes the second one
leaves the first standing alone. The signed error is positive under both
references — the flat profile is simply closer by accident, not by construction.

So: use a power-flow reference when you need the closure itself to be faithful —
delta devices, voltage-dependent loads, strong unbalance, or a linearization you
intend to differentiate through. Do not assume it lowers voltage error. IEEE 37
gains about 5 % because its all-wye constant-power loads barely exercise the
closure at all; the only reference dependence left there is ``\Gamma`` in the
line-drop coefficients.

## Widening the admissible input

`L3FOptions(unsupported=...)` controls what happens to data outside the
supported vocabulary. The important point is that "unsupported" covers three
different situations, and collapsing them into one lenient switch would hide
which one you are in.

| | What it is | Example | Accuracy cost |
|---|---|---|---|
| **Lowering** | outside the component vocabulary, inside the mathematical class | transformer leakage | **none** |
| **Missing feature** | expressible in the existing closure, simply not written yet | `vpp` / sequence limits | none, once implemented |
| **Projection** | genuinely destroys information | adjustable tap | real, and not quantifiable from inside |

The policy is a ladder:

| `unsupported` | Behaviour |
|---|---|
| `:reject` (default) | Nothing is rewritten. Unsupported data is an error, as the formulation's contract promises. |
| `:lower` | Exact re-representations only, reported as `L.L3F.*` at severity `:info`. **The solved model is the same physics.** |
| `:approximate` | Also the lossy projections, reported as `A.L3F.*` at severity `:warning`. **The solved model is a different problem.** |

Every rewrite appears in the applicability report, and
`result.formulation["unsupported_policy"]` and `["lowered"]` record what was in
force, so a result can always be traced to what was actually solved.

### Exact lowerings

Each of these replaces a component with supported components that reproduce the
same two-port or shunt behaviour. There is no approximation, and the regression
suite proves it by solving each case twice — once from the richer component,
once from a hand-written equivalent — and requiring agreement to solver
tolerance.

- **Closed switch → zero-impedance line.** A closed ideal switch is a branch
  with no series drop. An open switch is removed; if that leaves a subnetwork
  unenergized, the island check reports it rather than this pass guessing.
- **Fixed capacitor → fixed shunt.** ``B = q_{rated}/v_{nom}^2`` per coil, and
  for a delta bank ``D^{\mathsf T}\operatorname{diag}(b)D`` — the terminal
  admittance matrix of the same three coils. BMOPF capacitors carry no switching
  state, so nothing is assumed.
- **Line shunt → terminal shunts.** BMOPF already declares the from- and to-side
  halves separately, so moving each onto its own bus restates the same π model.
  Linecode entries are per unit length and scale with `length`.
- **Transformer leakage and no-load admittance → series line and shunt.** A
  winding leakage is a series impedance in the coil's own coordinates and the
  no-load admittance is a shunt across the winding-2 coil. Introducing one
  internal bus per non-zero winding and stamping ordinary line and shunt
  elements reproduces the two-port exactly. The internal buses and branches
  appear in the result under `_l3f_` names, which is what makes the rewrite
  auditable rather than hidden.

That last one is the substantive one: it turns "ideal transformers only" into
"any transformer of a supported connection", at no cost in fidelity.

### Projections

- **Constant-current, ZIP with a current fraction, and exponential loads → ZP.**
  A term ``(V/V_{nom})^\gamma`` is matched in value and first derivative at
  ``V=V_{nom}`` by ``\alpha_P+\alpha_Z(V/V_{nom})^2`` with
  ``\alpha_Z=\gamma/2``, ``\alpha_P=1-\gamma/2``. Exact for ``\gamma=0`` and
  ``\gamma=2``; for ``\gamma=1`` it is the familiar half-to-Z, half-to-P split
  of a constant-current term. The error is second order in the voltage deviation
  from nominal.
- **Adjustable tap → fixed tap.** The declared operating tap when it lies inside
  the interval, otherwise the midpoint. This removes a decision variable: the
  answer is feasible *for that setting*, not optimal over the range.
- **Unassessed bus limits are dropped.** `vpp_*`, `vpos_*`, `vneg_max`,
  `vzero_max` and `vm_unbalance_max` are removed with a warning. The solved
  problem is a **relaxation** and its solution may violate them. These are the
  "missing feature" row of the table above rather than a true impossibility —
  each is affine in ``w`` through the same cross-voltage closure the model
  already forms, so implementing them properly is the right eventual fix.

!!! warning "The replay does not measure the projection error"
    Under `:approximate` the snapshot the nonlinear replay runs on is the
    *projected* network, so both sides describe the same substituted physics.
    `validation["replayed_network"]` is `"projected"` in that case and
    `"as_supplied"` otherwise.

### What is never projected

Projection must not invent physics. These stay errors at every policy level,
because there is no defensible substitution:

- **`wye_delta`, `delta_wye`, `center_tap`, `n_winding` transformers.** Each
  needs a real voltage map; treating one as a per-conductor ratio would destroy
  the phase shift and the zero-sequence blocking. These are a missing feature —
  the maps are fixed matrices of exactly the form ``T`` already takes — not an
  impossibility.
- **Meshed islands, multiple or missing sources.** Choosing a spanning tree or
  a slack would change which problem is being solved.
- **IBR component models.** Replacing a control law with a free P/Q box would
  make the OPF optimistic in exactly the dimension the IBR exists to constrain.
- **DC subsystems.**

## Diagnostic codes

Every rejection is a typed [`L3FFinding`](@ref) with a stable code. Errors make
the report inapplicable and cause [`build_l3f_opf`](@ref) to raise
[`L3FInapplicableError`](@ref); warnings and `:info` findings never block a
build. The prefix carries the severity: `E.` error, `W.` warning, `L.` an exact
lowering (`:info`), `A.` a lossy projection (`:warning`). A test enforces that
each code below has a reachable case.

| Code | Severity | Meaning |
|---|:--:|---|
| `E.L3F.EMPTY_NETWORK` | E | The network declares no buses. |
| `E.L3F.BUS_INVALID` | E | A bus entry is not an object. |
| `E.L3F.BUS_UNKNOWN` | E | A device references an undeclared bus. |
| `E.L3F.TERMINAL_MAP_INVALID` | E | A terminal map is empty, duplicated, or names a conductor the bus does not declare. |
| `E.L3F.DEVICE_ARITY` | E | Channel counts and terminal maps disagree, or a required per-channel vector is missing. |
| `E.L3F.DEVICE_DATA_INVALID` | E | Device data is non-finite or violates its own ordering (`p_min > p_max`, and similar). |
| `E.L3F.VOLTAGE_BOUND_INVALID` | E | A bus `v_min`/`v_max` is non-finite, negative, mis-sized, or inverted. |
| `E.L3F.LIMIT_INVALID` | E | A rating is non-positive, non-finite, or the wrong shape for its device. |
| `E.L3F.LIMIT_UNSUPPORTED` | E | A bus declares a sequence or phase-to-phase limit the formulation does not assess. |
| `E.L3F.CONNECTION_UNSUPPORTED` | E | A device configuration outside `WYE`/`SINGLE_PHASE`/`DELTA` (sources: wye only). |
| `E.L3F.LOAD_MODEL_UNSUPPORTED` | E | A load model other than constant power, constant impedance, or ZIP. |
| `E.L3F.ZIP_CURRENT_UNSUPPORTED` | E | A ZIP load with a nonzero current fraction; not affine in squared voltage. |
| `E.L3F.LINE_MATRIX_INVALID` | E | A line has no impedance source, both inline and linecode data, or a malformed matrix. |
| `E.L3F.LINE_SHUNT_UNSUPPORTED` | E | A line or linecode declares shunt admittance. |
| `E.L3F.SHUNT_INVALID` | E | A shunt admittance is non-finite or declares an entry outside its terminal-map arity. |
| `E.L3F.TRANSFORMER_UNSUPPORTED` | E | A transformer subtype outside the supported three. |
| `E.L3F.TRANSFORMER_RATIO_INVALID` | E | A ratio is missing, non-positive, non-finite, or the wrong arity. |
| `E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED` | E | Nonzero leakage or no-load admittance on a device modelled as ideal. |
| `E.L3F.ADJUSTABLE_TAP_UNSUPPORTED` | E | A tap interval with `min < max`; taps are fixed data here. |
| `E.L3F.TOPOLOGY_NOT_RADIAL` | E | Two devices share a conductor between the same buses, or a conductor component contains a cycle. |
| `E.L3F.SOURCE_MISSING` | E | An energized island has no voltage source. |
| `E.L3F.MULTIPLE_SOURCES` | E | An energized island has more than one. |
| `E.L3F.REFERENCE_MISSING` | E | A reference is required but absent, incomplete, or a source does not cover its root bus. |
| `E.L3F.REFERENCE_ZERO_WINDING` | E | A reference phasor or winding voltage is zero, so a coefficient is undefined. |
| `E.L3F.EXPLICIT_NEUTRAL_UNSUPPORTED` | E | An explicit neutral with `kron_reduce=false`. |
| `E.L3F.GROUNDED_TERMINAL_RETAINED` | E | A perfectly grounded terminal survived reduction; `kron_reduce_bmopf` eliminates only conductors declared as neutrals. |
| `E.L3F.KRON_REDUCTION_FAILED` | E | The reduction itself raised. |
| `E.L3F.NEUTRAL_REDUCTION_UNDECLARED` | E | `require_neutral_provenance` is set but no provenance or explicit reference was supplied. |
| `E.L3F.TIME_SERIES_UNSUPPORTED` | E | A component references a time series; resolve to a snapshot first. |
| `E.L3F.CONTROL_PROFILE_UNSUPPORTED` | E | Top-level control-profile or time-series tables must be resolved first. |
| `E.L3F.SWITCH_UNSUPPORTED` | E | Switches are outside the slice. |
| `E.L3F.COMPONENT_UNSUPPORTED` | E | Capacitors, IBR component models, and other unsupported families. |
| `E.L3F.DC_SUBSYSTEM_UNSUPPORTED` | E | Any DC subsystem table. |
| `E.L3F.CAPACITOR_INVALID` | E | A capacitor could not be lowered: non-positive `v_nom`, non-finite `q_rated`, or a connection the incidence map does not cover. |
| `W.L3F.COST_MISSING` | W | `objective=:cost` but a dispatchable unit declares no `cost`; it is priced at zero and the optimum may be non-unique. |
| `L.L3F.SWITCH_LOWERED` | I | Closed switch represented exactly as a zero-impedance line. |
| `L.L3F.SWITCH_OPEN_REMOVED` | I | Open switch removed; any subnetwork it alone energized is now a separate island. |
| `L.L3F.CAPACITOR_LOWERED` | I | Fixed capacitor represented exactly as a shunt. |
| `L.L3F.LINE_SHUNT_LOWERED` | I | A declared π half moved onto its own bus as a shunt. |
| `L.L3F.TRANSFORMER_LEAKAGE_LOWERED` | I | Winding leakage represented exactly as a series line through an internal bus. |
| `L.L3F.TRANSFORMER_NO_LOAD_LOWERED` | I | No-load admittance represented exactly as a shunt across the to-side coil. |
| `A.L3F.LOAD_LAW_PROJECTED` | W | A constant-current, ZIP-with-current, or exponential law projected onto its ZP tangent at `v_nom`. |
| `A.L3F.ADJUSTABLE_TAP_PROJECTED` | W | A tap interval collapsed to one fixed setting; the optimizer no longer selects the tap. |
| `A.L3F.BUS_LIMIT_DROPPED` | W | A bus limit the formulation does not assess was removed; the solved problem is a relaxation. |

## Result contract

[`solve_l3f_opf`](@ref) returns an [`L3FResult`](@ref). Every numeric field is
SI regardless of `per_unit`, and is `NaN` when the solve was not optimal.

| Path | Contents |
|---|---|
| `buses[bus][terminal]` | `"w"` (V²), `"vm"` (V), `"reference_angle"` (rad, coefficient data — not a solved angle) |
| `lines[id]`, `transformers[id]` | `"p"`, `"q"` (W, var), `"parent"`, `"child"`, `"terminal_map_parent"`, `"terminal_map_child"`, `"reversed_from_input"` |
| `transformers[id]` also | `"subtype"`, `"effective_ratio_from_to"` |
| `generators[id]`, `sources[id]` | `"pg"`, `"qg"` (W, var), `"terminal_map"` |
| `objective` | currency/h for `:cost`, W for `:source_import`, `0.0` for `:feasibility` |
| `formulation` | `"name"`, `"version"`, `"problem_class"`, `"per_unit"`, `"s_base"`, `"working_units"`, `"result_units"`, `"series_losses"`, `"voltage_angles"`, `"unsupported_policy"`, `"lowered"`, `"reference_provenance"`, `"reference_hash"` |
| `validation` | see below |

!!! warning "Branch flows are oriented parent → child"
    `lines[id]["p"]` is positive when power flows from `"parent"` toward
    `"child"`, which is the direction away from the island's source — **not**
    necessarily BMOPF's `bus_from → bus_to`. When `"reversed_from_input"` is
    `true` the sign is opposite to the input orientation, and the entries are
    indexed by `"terminal_map_parent"` rather than `terminal_map_from`.

Under `unsupported=:lower` or `:approximate` the result also contains the
internal buses and branches the lowering introduced, under `_l3f_`-prefixed
names. They are real model elements, not bookkeeping, and the corresponding
`L.L3F.*` finding names each one.

`validation["status"]` is `"replayed"`, `"failed"`, or `"not_run"`. It reports
whether the nonlinear replay *ran*, never whether the linear answer was accurate
— the omitted losses guarantee some difference. Pass `voltage_tolerance` to
[`solve_l3f_opf`](@ref) or [`validate_l3f_solution`](@ref) to have the margin
judged; the result then also carries `"within_tolerance"`. `"physical_limits"`
is always `"unassessed"`.

## Deliberate exclusions

The first implementation rejects meshed islands, multiple or missing sources,
line shunts, nonideal or adjustable transformers/regulators, switches,
controllable capacitors, IBR component models, DC subsystems, constant-current
loads, ZIP loads with a nonzero current fraction, exponential loads, time-series
controls, and sequence-voltage limits. Series losses are omitted.

Within the affine coordinates the model is exact: there is no artificial physics
slack, no integer variable, no non-SOC cone, and no polyhedral outer
approximation of a cone. The coordinates themselves rest on the three
approximations above.

When nonlinear validation is enabled, `solve_l3f_opf` fixes the optimized
generator dispatch in the reduced snapshot and calls BMOPFTools power flow. The
reported voltage error is a replay comparison, not a certificate that excluded
physical limits were satisfied.

One upstream gap is compensated here rather than worked around: BMOPFTools'
per-unit preparation at the pinned revision scales a voltage source's bounds and
cost but leaves `s_max` and `i_max` in SI. PowerOptLab restates both from the
caller's SI network on every build, so a source nameplate binds identically
under `per_unit=true` and `per_unit=false`.

## Staged use

```julia
report = check_l3f_applicability(network)
is_l3f_applicable(report) || foreach(println, report.findings)

build = build_l3f_opf(network, Clarabel.Optimizer;
    options=L3FOptions(validate_nonlinear=false))
@assert l3f_model_class(build) in (:LP, :QP, :SOCP)
optimize!(build.model)
```

The `L3FBuild` object exposes semantic variable and constraint dictionaries so
future refinements can add limits or objectives without reaching through JuMP
names. Use `solve_l3f_opf` for the stable result contract.

## API

```@docs
L3FOptions
L3FFinding
L3FApplicabilityReport
L3FInapplicableError
L3FReferenceState
CrossVoltageCoefficients
AffineScalarCoefficients
ConnectionPowerMap
LineDropCoefficients
L3FBuild
L3FResult
is_l3f_applicable
cross_voltage_coefficients
evaluate_cross_voltage
winding_voltage_coefficients
evaluate_affine
connection_power_map
line_drop_coefficients
regulator_gain_matrix
check_l3f_applicability
build_l3f_opf
l3f_reference_from_powerflow
l3f_model_class
validate_l3f_solution
solve_l3f_opf
```
