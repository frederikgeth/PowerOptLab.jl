# Generalized generators: implementation and API

> **Kind:** Component model · **Maturity:** research prototype · **Direction:** forward · **Temporal:** single snapshot or independent snapshots

The implementation adds optional generator components to the BMOPFTools engine.
It supports all five source-voltage laws in the [scientific model](generalized_generator.md),
sequence capability overlays, a hard terminal-PV target, and a dedicated source
with separate neutral and earth currents. The [numerical tutorial](../tutorials/generalized_generator_tradeoffs.md)
compares these choices on one loaded, unbalanced four-wire feeder.
The [data-model specification drafts](../specifications/generators/index.md) give
upstream-style field tables, complete equations, and an explicit proposed JSON-to-API
mapping for later Task Force contributions.

## A first device

```julia
using PowerOptLab, BMOPFTools, JuMP, Ipopt

generator = GeneralizedGenerator(
    id="machine", bus="poc",
    connections=[("a","n"), ("b","n"), ("c","n")],
    impedance=generator_sequence_impedance(0.3+0.2im, 0.2+0.6im, 0.1+0.25im),
    voltage=GeneratorVoltageLaw(:common_magnitude;
        magnitude_min=190.0, magnitude_max=270.0),
    control=GeneratorControl(p=3000.0, voltage_target=230.0),
    capability=GeneratorCapability(i_max=30.0,
        terminal_i_max=[30.0,30.0,30.0,15.0]))

# net contains this bus and an appropriate system angle reference.
result = solve_generator_opf(net, [generator]; s_base=1e4)
if solve_status(result).publishable
    g = result.devices["machine"]
    println(g.emf, g.port_voltage, g.p, g.q)
end
```

The impedance here is an effective three-port impedance; its zero-sequence
value must already contain the intended return-path contribution. Do not supply
a four-conductor matrix to this constructor. For an explicit source neutral lead
and grounding connection, use `SourceGenerator` instead.

## Units, ordering, and power measurements

Inputs and results use RMS fundamental-frequency volts, amperes, ohms, watts, and
vars. `cost` is currency/kWh, so the single-snapshot cost objective is currency/h.
`voltage_scale` supplies a positive numerical scale and a fallback voltage start;
it neither fixes a voltage nor creates a bound. Frequency dependence is the
responsibility of the supplied impedance data; these components contain no
frequency-state or harmonic model.

The generalized device uses the order of `connections` for port quantities. Its
terminal order is first all distinct outgoing terminals, then any new return
terminals. A wye example therefore reports `a,b,c,n`. `SourceGenerator` always
orders phases first and its explicit neutral last. Sequence quantities use the
declared three-port order as a,b,c; labels are not automatically permuted.

Power-control/capability scalars refer to aggregate power. Vectors refer to the
individual oriented ports. The two sides of a power bound can have different
shapes: scalar aggregate P minimum and vector phase P maxima constrain different
quantities. For a grounded source, scalar PCC power is
``\sum_t V_t\overline{I_t}`` over phases **and neutral**. Vector PCC powers are
``(V_k-V_n)\overline{I_k}``; their sum differs from the complete conductor power
when both neutral voltage and earth current are nonzero. The scientific model
derives the correction. `power_location=:internal` refers to EMF power, including
the physical series/ground losses upstream of the PCC, not shaft power.

## Selecting a voltage law

| Mode | Required data | Meaning | Typical accompanying control |
|:--|:--|:--|:--|
| `:none` | None | No internal voltage-source constraint | Fixed P/Q, or an OPF dispatch region |
| `:fixed_phasor` | Complex `phasor` vector | Internal magnitude and absolute angle fixed | P/Q solved by the network |
| `:fixed_magnitudes` | Complex `phasor` vector | Fixed template up to one common rotation | Fixed aggregate P |
| `:common_magnitude` | Optional scalar magnitude bounds, relative angles | One common free magnitude and rotation | Aggregate P/Q, or P and terminal voltage |
| `:phase_magnitudes` | Positive magnitude minima, relative angles | Independent magnitudes and common rotation | E.g. aggregate P and three phase Q targets |

The phasor template can be unbalanced. Fixed 120-degree phase separation alone
does not enforce balanced magnitudes. Default relative offsets are 0 for one
port and 0,-120,+120 degrees for three ports. Other port counts require explicit
`angles`; use `[0.0, pi]` for two split-phase legs. The magnitude variables are
physical freedoms, not automatically regulator settings. An objective can
select among remaining freedoms; that does not make the selected point a
controller response.

For an isolated system supply one appropriate angle reference per independent
rotational symmetry. A full fixed-phasor source supplies an absolute reference;
the other modes do not. The builder does not automatically identify islands or
pin every generator. Explicit bounds with conflicting ideal references produce
an error or an infeasible solve. No epsilon impedance is inserted to resolve
inconsistent ideal sources.

## Grounded source

```julia
source = SourceGenerator(
    id="source", bus="poc", phase_terminals=["a","b","c"], neutral="n",
    impedance=ComplexF64[0.2+0.5im, 0.2+0.5im, 0.2+0.5im, 0.1+0.05im],
    grounding=0.8+0.1im,
    voltage=GeneratorVoltageLaw(:fixed_phasor;
        phasor=230 .* cis.([0.0,-2pi/3,2pi/3])),
    capability=GeneratorCapability(
        i_max=40.0, terminal_i_max=[40.0,40.0,40.0,20.0], earth_i_max=10.0))
```

The four impedance entries are phase leads then the neutral lead, not three
sequence impedances. A full four-by-four matrix includes mutual coupling.
`grounding=nothing` opens the earth path; zero creates an ideal internal-star
bond; finite impedance connects the star to remote earth at 0 V. A separate
network PCC ground can coexist with a finite neutral lead. The same ideal bond
must not also be declared as an engine perfect ground when the neutral lead is
ideal: that creates redundant return paths and nonunique current sharing, so
validation rejects it. The star potential and earth current are expressions,
including at ideal limits.

Results expose `terminal_current[end]`, `earth_current`, and `star_voltage`
separately. The identity is ``3I_0+I_n+I_g=0``, not ``I_n=-3I_0``. Do not use
the engine's retired native source-current placeholders to reconstruct these
quantities.

## Capability and numerical formulation

Every limit is optional. Current limits apply to the named physical current and
are never inferred from power ratings. Sequence limits use order zero, positive,
negative and an explicit voltage location:

```julia
limits = GeneratorCapability(
    p_min=0.0, p_max=15e3, q_min=-8e3, q_max=8e3, s_max=16e3,
    i_max=[30.0,30.0,30.0], terminal_i_max=[30.0,30.0,30.0,12.0],
    voltage_min=200.0, voltage_max=255.0,
    sequence=GeneratorSequenceLimits(location=:poc,
        voltage_max=[10.0,255.0,5.0], current_max=[4.0,30.0,3.0]))
```

These are capability inequalities. They do not prescribe negative-sequence
current phase, current-limit priority, droop, DC ripple, or fault response.
Only a three-port common-return connection can use standard sequence quantities.
Single-phase and split-phase devices use their port/conductor quantities.

Following the [BMOPFTools formulation principles](https://frederikgeth.github.io/BMOPFTools.jl/dev/opf/):

- The series primitive uses impedance form without an inverse. Internal EMF is
  an affine expression ``e=u+Zi``; the grounded source similarly reconstructs
  its star from the neutral lead. **No internal voltage variables are created**
  for either zero or nonzero impedance. This is exact elimination, not an
  approximate Thevenin reduction. Singular nonzero matrices remain supported.
- Fixed source laws use affine equalities. Fixed-template and common-magnitude
  laws use a spanning set of complex ratios to the first port; they add no
  separate global angle or amplitude variable. Mode e uses normalized quadratic
  cross products and orientation inequalities, with declared positive voltage
  minima excluding the undefined-angle domain. No tangent or all-pairs angle
  constraints are used.
- Internal sequence voltage bounds on a fixed template are checked directly.
  On a common-magnitude shape they tighten its one magnitude interval. Exact
  internal sequence zeros in mode e are reduced in real magnitude space when
  they determine one positive shape, then use the affine template equations.
  Unsupported or inconsistent shapes are rejected. This avoids appending
  dependent zero-sequence equations to an already balanced source law.
- Derived powers stay quadratic expressions. Apparent-power bounds deliberately
  lift normalized P and Q to variables in [-1,1], then impose a quadratic circle.
  No quartic power-magnitude row is introduced.
- Positive magnitude limits use normalized circles. Explicit current variables
  also receive the corresponding rectangular bounds. Exactly zero magnitude
  upper bounds become linear real/imaginary equalities or variable fixes;
  redundant zero lower bounds add nothing. Equal positive magnitude limits use
  one equality rather than two opposing inequalities. Identical affine norms
  share the intersection of their limits, including a single-phase device's
  outgoing and return current. Power targets substitute into same-location
  capability checks so a matching capability does not repeat a control row.
  Equal P/Q capability bounds use one equality per quantity. Identical affine
  equalities are stamped once, including open grounding plus a zero earth-current
  rating.
- Numerical bases come from `opf_coordinate_bases(ctx,bus)`, including local
  power bases. Initial currents follow the declared P/Q targets and the engine's
  voltage start; otherwise they use a deterministic small injection. These
  choices do not select or certify a particular power-flow branch.

The formulation is polynomial of degree at most two, but nonconvex. Magnitude
lower bounds, power equalities, and relative-angle equations remain nonconvex.
The implementation validates known domains and topology restrictions; it does
not claim a general Jacobian-rank certificate or a solve-speed improvement over
every alternative formulation. Physical power-balance residuals and normalized
solve status are reported separately.

## Engine composition and replacement

`stamp_device!(ctx, device)` works inside the existing `model_hook!`. It returns
a named tuple of model-unit expressions: `ur/ui` (PCC port voltage), `er/ei`
(internal EMF), `ir/ii` (port current), `jr/ji` (all external currents), `gr/gi`
(earth current), power/loss expressions, local `bases`, and the registered
`constraints`. It is deliberately possible to add research constraints to these
handles. Do not interpret their working coordinates as SI.

`build_generator_model` returns the engine `context` and handles before KCL and
optimization. Call `enforce_kcl!` exactly once before solving. Its hook receives
`(ctx, handles)` and can supply an explicit angle reference, an objective, or
additional problem equations. `solve_generator_opf` is the convenience composition
of those same steps; it contains no second network solver.

```julia
result = solve_generator_opf(net, [generator];
    replacements=Dict("machine" => (:generator,"old_generator")),
    objective=:cost, s_base=1e4,
    solver_options=Pair{String,Any}["tol"=>1e-9, "max_iter"=>2000])
```

A replacement transfers only the named native physics to PowerOptLab through
`OpfDeviceBuilder`. Other native devices remain active. Native placeholder currents
are retired; custom results and costs are authoritative. The returned native
result removes replaced records and adds a `custom_injection` PCC power ledger.
The input network is not mutated. Native costs/limits on a replaced device are
**not copied into the new device**: specify its intended capability and cost.

The pinned engine infers source-fixed bus terminals from native source metadata
even when the source builder is replaced. During the public operational-limit
stage the wrapper temporarily excludes replaced source metadata from the prepared
working network, restoring it before device construction. This preserves the
native bus limits without copying the engine's bus formulations. Source metadata
still seeds network scaling/initialization. This compatibility behavior is tested
against the pinned dependency; it is not a promise about arbitrary future pins.

`objective=:cost` includes native generation and custom PCC cost. `:loss`
minimizes **custom-device** series plus ground losses in kW, not total network
loss. `:feasibility` sets a zero objective. A stateless generator can also appear
in `solve_multiperiod_opf`: its costs are duration-weighted and its extracted
snapshots and custom injection ledger are retained. There is no inter-period
machine dynamic state or ramp law.

## Validation coverage and boundaries

The unit tests use phase-magnitude and phase-angle unbalance, unequal phase
injections, and neutral displacement. They exercise every voltage mode for both
device types, scalar/vector P/Q limits at both power locations, P/Q targets,
apparent-power lifts, PCC and internal voltage limits, positive/negative/zero
sequence capabilities, conductor/earth current limits, and deliberate infeasible
limits. Independent circuit calculations check full-matrix drops, neutral/earth
sharing, and complete power ledgers. Other tests cover SI/per-unit bases, native
replacement/cost equivalence, zero/singular impedances, single-phase/split-phase,
multi-period composition, domain rejection, and exact-zero constraint structure.

Supported connections are independent terminal pairs (a forest) and ordered
closed delta windings `[("a","b"),("b","c"),("c","a")]`. Delta ratings and
sequence quantities use winding coordinates; `terminal_i_max` limits line
currents. See the [delta equations](../specifications/generators/generator.md).
Open neutral leads on `SourceGenerator`, explicit soil/earth-network
terminals, physical shunts/LCL filters, machine saliency/saturation, and automatic
PV/PQ limiter switching are not implemented. The hard PV target reports
infeasibility if the capability cannot support it. The advanced IBR remains a
separate component; its DC, switching, ripple and hardware model is not replaced.

## API

```@docs
GeneratorVoltageLaw
GeneratorControl
GeneratorSequenceLimits
GeneratorCapability
GeneralizedGenerator
SourceGenerator
generator_sequence_impedance
GeneratorResult
build_generator_model
solve_generator_opf
GeneratorOPFResult
```
