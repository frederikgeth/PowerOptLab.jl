# LinDist3Flow BMOPF prototype

PowerOptLab provides an independent LinDist3Flow model builder for normalized
BMOPF JSON. It does not depend on PowerModelsDistribution and does not turn the
BMOPFTools nonlinear OPF engine into a multi-formulation engine. BMOPFTools is
reused at the boundaries: parsing and nonlinear replay; PowerOptLab's
[`kron_reduce_bmopf`](@ref) supplies the neutral-reduced network representation.

```julia
options = L3FOptions(validate_nonlinear=true)
result = solve_l3f_opf(network, Ipopt.Optimizer; options)
```

An explicit-neutral input is reduced on a deep copy by default. Pass
`L3FOptions(kron_reduce=false)` to require callers to supply an already reduced
case. `check_l3f_applicability` performs the same boundary checks without
constructing a JuMP model. Unsupported data is represented by typed findings
with stable codes; it is never silently dropped.

## Implemented slice

Version `0.1-prototype` builds a continuous affine LP for radial AC islands with
exactly one fixed-voltage source per island. It supports neutral-reduced series
lines without shunts, grounded-wye or single-phase constant-power loads,
grounded-wye or single-phase generators with per-terminal P/Q boxes, retained
phase-to-ground voltage-magnitude bounds, linear per-phase energy costs, and
reverse power flow.

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
`connection_power_map`, `line_drop_coefficients`, and
`regular_polygon_coefficients` are pure coefficient oracles. The connection map
implements ``H=\operatorname{diag}(\bar v)D^T
\operatorname{diag}(D\bar v)^{-1}``, allowing later delta and transformer work
to share independently tested algebra without admitting those devices yet.

## Deliberate exclusions

The first implementation rejects meshed islands, multiple or missing sources,
line shunts and ratings, transformers, switches, shunts/capacitors, delta
devices, IBRs, DC subsystems, non-constant-power loads, time-series controls,
and sequence or phase-to-phase limits. Series losses are omitted. Options for
fixed loss correction and droop intentionally fail until their equations and
validity checks are implemented.

When nonlinear validation is enabled, `solve_l3f_opf` fixes the optimized
generator dispatch in the reduced snapshot and calls BMOPFTools power flow. The
reported voltage error is a replay comparison, not a certificate that excluded
physical limits were satisfied.

## Staged use

```julia
report = check_l3f_applicability(network)
is_l3f_applicable(report) || foreach(println, report.findings)

build = build_l3f_opf(network, Ipopt.Optimizer;
    options=L3FOptions(validate_nonlinear=false))
@assert l3f_model_class(build) == :LP
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
RegularPolygonCoefficients
L3FBuild
L3FResult
cross_voltage_coefficients
winding_voltage_coefficients
connection_power_map
line_drop_coefficients
regular_polygon_coefficients
check_l3f_applicability
build_l3f_opf
l3f_model_class
validate_l3f_solution
solve_l3f_opf
```
