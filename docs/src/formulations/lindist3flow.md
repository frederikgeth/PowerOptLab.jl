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
lines without shunts under the default `unsupported=:reject` policy; `:lower`
adds the documented endpoint-shunt representation. It also supports
grounded-wye, single-phase and delta constant-power,
constant-impedance, and pure ZP ZIP loads, constant-power generators, fixed bus
shunts, ideal fixed-ratio single-phase
transformers, ideal Yd/Dy banks, and ANSI A/B autotransformer regulators, retained phase-to-ground
voltage-magnitude bounds, linear per-channel energy costs, and reverse power
flow. Apparent-power bounds are native second-order cones. Ampacity bounds use
the live squared terminal or affine winding voltage and are native rotated
second-order cones.
Fixed ideal three-phase `open_delta_regulator` banks are supported on three
retained phase terminals, including AB/CB and its two cyclic permutations.

Radiality is checked **per conductor**, not per bus: the graph whose nodes are
`(bus, terminal)` pairs must be a forest. A three-unit single-phase regulator
bank on one bus pair — the wye-bank topology of the IEEE 13, 34 and 123 feeders
— is therefore admissible, because its units occupy disjoint conductors. Two
devices sharing a conductor between the same buses remain non-radial.

### Where the line equation comes from

The underlying AC series branch obeys ``v_j = v_i - Z i``. LinDist3Flow forms
the outer product and keeps its diagonal,

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

``p,q`` are lossless series-flow variables oriented from the source-side parent
to the child. Because the dropped term is exactly what distinguishes the two
ends of the series path, the same variable serves both nodal balances with
opposite signs. Endpoint terminal power can nevertheless differ when a local
pi shunt is present; the rating cones include that shunt explicitly. On IEEE 37 the omitted losses are about
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

The BMOPF bus field `va_nom` is nominal angle-difference centering metadata,
not an absolute operating-point phasor. It therefore does not select or alter
the L3F reference. Supply explicit phasors, or use source propagation (including
the declared transformer vector groups), when choosing the fixed-angle closure.

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
``\widehat{|Dv|^2}=c_D+a_D^T w``. Pure ZP laws are therefore affine functions
of the chosen closure (exact as functions of that approximation, not an exact
AC load-flow law):

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

The canonical formulation is therefore never an exact AC physics model. A
phase-to-ground wye channel is the special case where the cross-phase closure
and channel split introduce no additional error; the line loss and fixed-angle
assumptions still apply to the surrounding network.

!!! note "Closed-delta channel powers are not uniquely determined"
    A closed delta's map ``H`` has rank 2: channel powers of the form
    ``s=c\,(D\bar v)`` are a circulating component the terminals cannot see.
    The published `pg`/`qg` for a three-channel delta device are therefore fixed
    only up to that component, even though its terminal injection is unique.

Fixed regulator and transformer settings preserve affine physics. Adjustable
tap intervals are rejected: there are no integer taps, McCormick envelopes, or
continuous voltage-ratio relaxations. Nonzero transformer leakage and no-load
admittance must be represented as separate supported network elements.

## Delta-wye and wye-delta banks

`wye_delta` and `delta_wye` are supported as ideal three-phase banks on three
retained conductors per side, the wye star point being the eliminated
conductor. Their coil relation is

```math
D\,\boldsymbol v_{\Delta}=g\,\boldsymbol v_{Y},\qquad
g=\sqrt3\,\frac{V^{nom}_{\Delta}}{V^{nom}_{Y}},
```

where ``D`` is the delta incidence and ``\boldsymbol v_Y`` is measured against
the grounded star point. This is BMOPFTools' executable convention:
`wye_delta` uses ``n_{eff}=\sqrt3/N`` and `delta_wye` uses
``n_{eff}=N\sqrt3``, with ``N=V^{nom}_{from}/V^{nom}_{to}`` and both nominals
quoted phase-to-neutral. The two spellings collapse to the same statement, and
the ``\sqrt3`` is exactly the difference between a coil that spans a
line-to-line voltage and a nominal quoted phase-to-neutral.

### Orientation is part of the model

``D`` has **rank 2**, and that single fact determines how the component
behaves. The relation fixes the wye voltages given the delta ones, but not the
reverse: a delta winding neither imposes nor carries a zero-sequence terminal
voltage.

- **Delta winding upstream** — ``T=D/g`` as an ideal connection map. The downstream wye voltages
  are fully determined, and their zero-sequence component is zero, which is the
  correct behaviour of an ideal bank whose zero-sequence impedance is zero.
  This is the ordinary substation and service-transformer arrangement, and it
  works at the default `unsupported=:reject`.
- **Wye winding upstream** — the delta terminal voltages are determined only up
  to a common offset. `E.L3F.DELTA_ORIENTATION_UNSUPPORTED` by default. Under
  `unsupported=:approximate` the pseudo-inverse ``T=g\,D^{+}`` selects the
  minimum-norm solution, which is the one with **zero zero-sequence voltage at
  the delta bus**, and `A.L3F.DELTA_ZERO_SEQUENCE_GAUGE` records the
  assumption. It is an assumption, not physics: in a phase-to-ground
  formulation the delta bus's ground reference actually comes from capacitive
  coupling this formulation has already discarded.

Because the map is not diagonal, the two sides of the bank carry different
terminal power. `i_max_from` and `i_max_to` are therefore applied to their
declared physical endpoints through live-voltage rotated SOCs rather than to
whichever side the traversal happened to reach. Yd/Dy `s_rating` retains its
BMOPF total-bank interpretation and is split across the three wye coils.

### Angles are coefficients, and the vector group is a real choice

The formulation has **no angle variables** — only ``w=|v|^2`` and device powers.
Reference phasors ``\bar v`` are fixed data, computed once by propagating the
source through the topology, and a Yd/Dy map contributes its ±30° to them. So
the bank's phase shift is encoded in the fixed coefficients and published as
`reference_angle`; it is reference data, not a solved quantity.

``\bar v`` enters only through **same-bus** ratios and conjugate products —
``\Gamma_{\phi\psi}=\bar v_\phi/\bar v_\psi``, the closure's
``\bar v_\phi\bar v_\psi^*``, and
``H=\operatorname{diag}(\bar v)D^{\mathsf T}\operatorname{diag}(D\bar v)^{-1}``,
where the leading rotation cancels the one from the inverse. Every one is
invariant under ``\bar v\to e^{j\theta}\bar v`` at that bus, so **rotating
each bus's reference independently changes nothing** (verified to 3e-16). A
constant phase shift therefore cannot affect a solution here, and an
angle-difference constraint between buses is inexpressible.

That invariance is what makes the ±30° harmless — but it is *not* a licence to
treat vector groups as interchangeable:

!!! warning "A balanced comparison cannot see the vector group"
    BMOPFTools pairs delta coil ``k`` with wye phase ``k``
    (``v_{\Delta,k}-v_{\Delta,k+1}=n_{eff}v_{Y,k}``). Writing the delta winding
    in OpenDSS's **default** node order (`buses=[d.1.2.3 …]`) pairs them
    differently. Under a **balanced** reference the two differ by a uniform 60°
    rotation and are indistinguishable in every solved quantity. Under an
    **unbalanced** reference they are not: the pairing is a cyclic relabelling
    of which delta pair drives which wye phase, so the magnitudes themselves
    move — about 11 % in the case the regression pins.

    To reproduce BMOPFTools' group in OpenDSS, write the delta winding as
    `buses=[d.2.3.1 …]`. The regression cross-checks against that deck with an
    unbalanced source, where magnitudes agree to solver tolerance, and separately
    pins the default order as materially different so the distinction cannot be
    "fixed" away later.

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
Two-winding apparent-power ratings become native SOCs and current ratings use
live-voltage rotated SOCs; note that the shared (common) phase of an open-delta bank carries
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
| **Canonical lowering** | outside the component vocabulary, inside the L3F component class | switch or single-phase transformer impedance | preserves the L3F approximation |
| **Missing feature** | expressible in the existing closure, simply not written yet | `vpp` / sequence limits | none, once implemented |
| **Projection** | genuinely destroys information | adjustable tap | real, and not quantifiable from inside |

The policy is a ladder:

| `unsupported` | Behaviour |
|---|---|
| `:reject` (default) | Nothing is rewritten. Unsupported data is an error, as the formulation's contract promises. |
| `:lower` | Canonical L3F-preserving lowerings only, reported as `L.L3F.*` at severity `:info`. The surrounding model remains a fixed-angle, lossless approximation. |
| `:approximate` | Also experimental projections, reported as `A.L3F.*` at severity `:warning`. The projected model is a different problem. |

Every rewrite appears in the applicability report, and
`result.formulation["unsupported_policy"]` and `["lowered"]` record what was in
force, so a result can always be traced to what was actually solved.

### Canonical lowerings

Each of these replaces a component with supported components while preserving
the canonical L3F data and equations. This does not remove the fixed-angle
closure or omitted-series-loss approximation.

- **Closed switch → zero-impedance line.** A closed ideal switch is a branch
  with no series drop. An open switch is removed; if that leaves a subnetwork
  unenergized, the island check reports it rather than this pass guessing.
- **Fixed capacitor → fixed shunt.** ``B = q_{rated}/v_{nom}^2`` per coil, and
  for a delta bank ``D^{\mathsf T}\operatorname{diag}(b)D`` — the terminal
  admittance matrix of the same three coils. BMOPF capacitors carry no switching
  state, so nothing is assumed.
- **Line shunt → terminal shunts.** BMOPF already declares the from- and to-side
  halves separately, so moving each onto its own bus preserves the canonical π
  data.
  Linecode entries are per unit length and scale with `length`.
- **Single-phase transformer leakage and no-load admittance → series line and
  shunt.** A phase-to-ground winding leakage is lowered in its own coordinates,
  and no-load admittance is placed across the to-side coil. Connection-aware
  Yd/Dy, autotransformer, and other coupled winding semantics remain
  unsupported; they are not replaced by diagonal phase lines.

The lowering is auditable through `_l3f_` internal elements, but it does not
turn the canonical approximation into an exact transformer or AC network model.

### Projections

- **Experimental constant-current, ZIP with a current fraction, and exponential
  loads → ZP.**
  A term ``(V/V_{nom})^\gamma`` is matched in value and first derivative at
  ``V=V_{nom}`` by ``\alpha_P+\alpha_Z(V/V_{nom})^2`` with
  ``\alpha_Z=\gamma/2``, ``\alpha_P=1-\gamma/2``. This is exact only for the
  endpoint exponent classes within the projection; for ``\gamma=1`` it is the
  familiar half-to-Z, half-to-P split of a constant-current term. For other
  exponents it is a tangent projection, reported with an experimental warning.
- **Adjustable tap → fixed tap.** The declared operating tap when it lies inside
  the interval, otherwise the midpoint. This removes a decision variable: the
  answer is feasible *for that setting*, not optimal over the range.
- **Unassessed voltage and angle limits are never dropped.** Bus `vpn_*`,
  `vn_max`, sequence/unbalance limits, phase-to-phase limits, and line
  `va_diff_*` remain applicability errors under every policy until they are
  stamped explicitly.

!!! warning "The replay does not measure the projection error"
    Under `:approximate` the snapshot the nonlinear replay runs on is the
    *projected* network, so both sides describe the same substituted physics.
    `validation["replayed_network"]` is `"projected"` in that case and
    `"as_supplied"` otherwise.

### What is never projected

Projection must not invent physics. These stay errors at every policy level,
because there is no defensible substitution:

- **`center_tap` and `n_winding` transformers.** Each needs a real voltage map;
  treating one as a per-conductor ratio would destroy the phase shift and the
  zero-sequence blocking. These remain a missing feature — the maps are fixed
  matrices of exactly the form ``T`` already takes — rather than an
  impossibility. (`wye_delta` and `delta_wye` were in this list and are now
  [implemented](#Delta-wye-and-wye-delta-banks).)
- **Meshed islands, multiple or missing sources.** Choosing a spanning tree or
  a slack would change which problem is being solved.
- **IBR component models.** Replacing a control law with a free P/Q box would
  make the OPF optimistic in exactly the dimension the IBR exists to constrain.
- **DC subsystems.**

## Diagnostic codes

Every rejection is a typed [`L3FFinding`](@ref) with a stable code. Errors make
the report inapplicable and cause [`build_l3f_opf`](@ref) to raise
[`L3FInapplicableError`](@ref); warnings and `:info` findings never block a
build. The prefix carries the severity: `E.` error, `W.` warning, `L.` a canonical
lowering (`:info`), `A.` an experimental projection (`:warning`). A test enforces that
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
| `E.L3F.LIMIT_UNSUPPORTED` | E | A bus or line declares a phase-to-neutral, neutral, phase-to-phase, sequence/unbalance, or angle limit the formulation does not assess. |
| `E.L3F.CONNECTION_UNSUPPORTED` | E | A device configuration outside `WYE`/`SINGLE_PHASE`/`DELTA` (sources: wye only). |
| `E.L3F.LOAD_MODEL_UNSUPPORTED` | E | A load model other than constant power, constant impedance, or ZIP. |
| `E.L3F.ZIP_CURRENT_UNSUPPORTED` | E | A ZIP load with a nonzero current fraction; not affine in squared voltage. |
| `E.L3F.LINE_MATRIX_INVALID` | E | A line has no impedance source or has a malformed selected series matrix. |
| `E.L3F.LINE_IMPEDANCE_SOURCE` | E | A line declares both inline series impedance and a linecode; BMOPF requires exactly one source. |
| `E.L3F.LINE_SHUNT_UNSUPPORTED` | E | A line or linecode declares shunt admittance. |
| `E.L3F.LINE_SHUNT_INVALID` | E | A selected line-shunt matrix declares an entry outside its terminal-map arity. |
| `W.L3F.LINE_SHUNT_ASYMMETRIC` | W | A lowered line-shunt matrix is asymmetric; it is preserved, but is not the reciprocal passive-line form normally expected. |
| `E.L3F.SHUNT_INVALID` | E | A shunt admittance is non-finite or declares an entry outside its terminal-map arity. |
| `E.L3F.TRANSFORMER_UNSUPPORTED` | E | A transformer subtype outside the supported set. |
| `E.L3F.TRANSFORMER_RATIO_INVALID` | E | A ratio is missing, non-positive, non-finite, or the wrong arity. |
| `E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED` | E | Nonzero leakage or no-load admittance on a device modelled as ideal. |
| `E.L3F.DELTA_ORIENTATION_UNSUPPORTED` | E | A Yd/Dy bank has its wye winding facing the source, leaving the delta terminals undetermined in their common component. |
| `A.L3F.DELTA_ZERO_SEQUENCE_GAUGE` | W | That orientation accepted under `:approximate` by assuming zero zero-sequence voltage at the delta bus. |
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
| `L.L3F.SWITCH_LOWERED` | I | Closed switch represented in the canonical L3F model as a zero-impedance line. |
| `L.L3F.SWITCH_OPEN_REMOVED` | I | Open switch removed; any subnetwork it alone energized is now a separate island. |
| `L.L3F.CAPACITOR_LOWERED` | I | Fixed capacitor represented in the canonical L3F model as a shunt. |
| `L.L3F.LINE_SHUNT_LOWERED` | I | A declared π half moved onto its own bus as a shunt. |
| `L.L3F.TRANSFORMER_LEAKAGE_LOWERED` | I | Single-phase winding leakage represented in the canonical L3F model as a series line through an internal bus. |
| `L.L3F.TRANSFORMER_NO_LOAD_LOWERED` | I | Single-phase no-load admittance represented in the canonical L3F model as a shunt across the to-side coil. |
| `A.L3F.LOAD_LAW_PROJECTED` | W | A constant-current, ZIP-with-current, or exponential law projected onto its ZP tangent at `v_nom`. |
| `A.L3F.ADJUSTABLE_TAP_PROJECTED` | W | A tap interval collapsed to one fixed setting; the optimizer no longer selects the tap. |

## Result contract

[`solve_l3f_opf`](@ref) returns an [`L3FResult`](@ref). Electrical result fields
are SI regardless of `per_unit`, and decision-variable values are `NaN` when
the solve was not optimal. Dimensionless ratios, angles in radians, and metadata
retain the units stated below.

| Path | Contents |
|---|---|
| `buses[bus][terminal]` | `"w"` (V²), `"vm"` (V), `"reference_angle"` (rad, coefficient data — not a solved angle) |
| `lines[id]`, `transformers[id]` | model branch `"p"`, `"q"` (W, var), `"parent"`, `"child"`, `"terminal_map_parent"`, `"terminal_map_child"`, `"reversed_from_input"` |
| `transformers[id]` also | `"subtype"`, `"effective_ratio_from_to"` |
| `generators[id]`, `sources[id]` | `"pg"`, `"qg"` (W, var), `"terminal_map"` |
| `objective` | currency/h for `:cost`, W for `:source_import`, `0.0` for `:feasibility` |
| `formulation` | `"name"`, `"version"`, `"problem_class"`, `"per_unit"`, `"s_base"`, `"working_units"`, `"result_units"`, `"series_losses"`, `"voltage_angles"`, `"unsupported_policy"`, `"lowered"`, `"reference_provenance"`, `"reference_hash"` |
| `validation` | see below |

!!! warning "Branch flows are oriented parent → child"
    For lines, `p` and `q` are the lossless series-flow variables and are
    positive from `"parent"` toward `"child"`, away from the island's source —
    **not** necessarily BMOPF's `bus_from → bus_to`. Transformer `p` and `q`
    are the corresponding child-side variables used by its fixed power map.
    When `"reversed_from_input"` is true the orientation is opposite to the
    input declaration. Endpoint totals used by rating cones, including local
    shunts, are affine model expressions and are not separately returned.

Under `unsupported=:lower` or `:approximate` the result also contains the
internal buses and branches the lowering introduced, under `_l3f_`-prefixed
names. They are real model elements, not bookkeeping, and the corresponding
`L.L3F.*` finding names each one.

Extraction rescales `w` by ``V_b^2`` and every published power by ``S_b`` in
per-unit mode, then reports `vm = sqrt(max(0,w))`. The published `w` itself is
not clamped. `reference_angle` comes from the fixed reference, transformer
ratios are recomputed from input data, and no phasors, endpoint flows, currents,
losses, or residual corrections are reconstructed. Canonical Kron reduction,
lowering, and projection are not undone. See the
[normative extraction contract](lindist3flow_components.md#Result-extraction).

`validation["status"]` is `"replayed"`, `"failed"`, or `"not_run"`. It reports
whether the nonlinear replay *ran*, never whether the linear answer was accurate
— the omitted losses guarantee some difference. Pass `voltage_tolerance` to
[`solve_l3f_opf`](@ref) or [`validate_l3f_solution`](@ref) to have the margin
judged; the result then also carries `"within_tolerance"`. `"physical_limits"`
is always `"unassessed"`.

## Deliberate exclusions

The first implementation rejects meshed islands, multiple or missing sources,
nonideal or adjustable transformers/regulators outside the justified
single-phase lowering, switches and line shunts unless `:lower` is selected,
controllable capacitors, IBR component models, DC subsystems, constant-current
loads, ZIP loads with a nonzero current fraction, exponential loads, time-series
controls, and sequence-voltage limits. Series losses are omitted.

Within the affine coordinates there is no additional relaxation: there is no
artificial physics slack, no integer variable, no non-SOC cone, and no
polyhedral outer approximation of a cone. The coordinates themselves still
rest on the three canonical approximations above.

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
