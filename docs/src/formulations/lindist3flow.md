# LinDist3Flow BMOPF prototype

PowerOptLab provides an independent LinDist3Flow model builder for normalized
BMOPF JSON. It does not depend on PowerModelsDistribution and does not turn the
BMOPFTools nonlinear OPF engine into a multi-formulation engine. BMOPFTools is
reused at the boundaries: parsing and nonlinear replay; PowerOptLab's
[`kron_reduce_bmopf`](@ref) supplies the neutral-reduced network representation.

```julia
options = L3FOptions(validate_nonlinear=true)
result = solve_l3f_opf(network, Clarabel.Optimizer; options)
```

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

For an oriented line, the model uses terminal powers ``p,q`` and squared voltage
magnitudes ``w``:

```math
w_j = w_i - M p_{ij} - N q_{ij},\qquad
M = 2\Re(\overline Z \odot \Gamma),\quad
N = -2\Im(\overline Z \odot \Gamma),\quad
\Gamma_{\phi\psi}=\bar v_\phi/\bar v_\psi.
```

The source phasors are propagated through the radial topology to form the fixed
coefficient reference. An explicit `L3FReferenceState` or phasor dictionary may
instead be supplied. Source values, missing or zero phasors, and topology maps
are validated before model construction. Angles are coefficient data, not
decision variables.

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

Fixed regulator and transformer settings preserve affine physics. Adjustable
tap intervals are rejected: there are no integer taps, McCormick envelopes, or
continuous voltage-ratio relaxations. Nonzero transformer leakage and no-load
admittance must be represented as separate supported network elements.

## Fixed regulator banks

`regulator_gain_matrix` implements the WYE, closed-delta, and open-delta gain
matrices from Table I of Bazrafshan, Gatsis, and Zhu, with the convention
``v_{from}=A v_{to}``. Under the BMOPF convention, ANSI type A uses the
declared `tap_ratio` directly and type B uses its reciprocal. This pure
coefficient oracle covers all
three published bank configurations and is tested against their tabulated
matrices.

The optimizer currently stamps only `open_delta_regulator`, because that is the
only three-phase regulator-bank component in the BMOPF schema. Its voltage law
uses the fixed matrix ``A^{-1}``, and its terminal-power law uses the same
fixed-reference complex power transformation as other multi-terminal devices.
Two-winding apparent-power and reference-current ratings become native SOC
constraints. Closed-delta and three-unit WYE banks remain coefficient oracles
until BMOPF provides an unambiguous component representation; PowerOptLab does
not invent private JSON subtypes.

## Literature and OpenDSS evidence

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

## Deliberate exclusions

The first implementation rejects meshed islands, multiple or missing sources,
line shunts, nonideal or adjustable transformers/regulators, switches,
controllable capacitors, IBR component models, DC subsystems, constant-current
loads, ZIP loads with a nonzero current fraction, exponential loads, time-series
controls, and sequence-voltage limits. Series losses are omitted. Pure ZP is
affine without approximation; the formulation contains no Taylor series,
artificial physics slack, integer variable, non-SOC cone, or polyhedral
approximation of a cone.

When nonlinear validation is enabled, `solve_l3f_opf` fixes the optimized
generator dispatch in the reduced snapshot and calls BMOPFTools power flow. The
reported voltage error is a replay comparison, not a certificate that excluded
physical limits were satisfied.

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
cross_voltage_coefficients
winding_voltage_coefficients
connection_power_map
line_drop_coefficients
regulator_gain_matrix
check_l3f_applicability
build_l3f_opf
l3f_model_class
validate_l3f_solution
solve_l3f_opf
```
