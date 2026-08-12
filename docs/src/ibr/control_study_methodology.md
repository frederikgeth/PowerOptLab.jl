# [Inverter-control study methodology](@id ibr-control-study-methodology)

This chapter turns the candidate [phase-aware local control laws](@ref
ibr-phase-aware-control-laws) into an implementation and experimental programme.
The goal is to compare control laws on LV networks with many PV inverters, then
quantify how each law changes AC conductor-current requirements,
DC-link capacitance, capacitor ripple-current rating, losses, curtailment, and
voltage quality.

The controller is deliberately simple enough to plausibly implement in product
firmware. The network representation is a smooth algebraic surrogate of the same
law, designed for sparse nonlinear solvers such as Ipopt.

## Questions and estimands

The study should answer two different questions. They must not be collapsed into
one comparison.

1. **Fixed-hardware performance:** with identical `Imax`, `Vdc`, `Cdc`, filter,
   PWM, and capacitor-current ratings, how do the controllers differ in voltage
   quality, energy export, losses, curtailment, and limit activation?
2. **Required-hardware performance:** what minimum current and DC-link ratings
   let each controller deliver a declared service envelope over the same set of
   network scenarios?

The first is a fair control comparison. The second is a control/hardware
co-design comparison. A law that performs better only by consuming more phase
current or DC ripple should not be reported as an unconditional improvement.

Primary response variables are:

- maximum and minimum phase-to-neutral voltage;
- negative-sequence voltage ratio `|U2|/|U1|` and, when observable, zero sequence;
- PV energy exported and curtailed;
- converter-terminal total and zero/positive/negative-sequence complex power;
- converter and filter losses;
- per-phase fundamental and total RMS current, plus `|I1|` and `|I2|`;
- `|S_tilde|`, `dv2`, capacitor 2ω RMS current, carrier RMS current, weighted
  thermal current, and capacitor ESR loss;
- internal-voltage and switching margin; and
- solver status, residuals, iterations, factorisation time, and memory.

## Software layers

The implementation should have four layers with one-way dependencies.

### 1. Physical plant

[`AdvancedInverter`](@ref) remains the source of truth for converter-side
voltages and currents, filters, topology, switching feasibility, losses, and
DC-link stress. Control laws must constrain its existing handles; they must not
reimplement a simplified current-injection plant beside it.

The present programme treats non-negative PV availability as the DC-side input.
Battery schedules, state of charge, and bidirectional active-power semantics are
out of scope; they must not complicate the first control-law comparison.

### 2. Controller primitives

Add `src/components/inverter_controls.jl` containing no network-study logic:

```julia
abstract type AbstractInverterControlLaw end
abstract type AbstractPositiveSequencePolicy end
abstract type AbstractUnbalancePolicy end
abstract type AbstractLimiterPolicy end
abstract type AbstractCurrentTarget end

struct PiecewiseLinearLaw
    x::Vector{Float64}
    y::Vector{Float64}
end

struct SequenceController{P,U,L,T} <: AbstractInverterControlLaw
    positive::P
    unbalance::U
    limiter::L
    current_target::T
    power_voltage_floor::Float64
end
```

The concrete first policies should be:

- `AverageVoltageVoltVarWatt` — the legacy mean-magnitude comparator;
- `PositiveSequenceVoltVarWatt` — ``|U_1|`` Volt-var with a maximum-phase
  Volt-watt guard;
- `WorstPhaseVoltVarWatt` — high/low phase envelopes and declared conflict
  resolution for conventional PWL Volt-var/Volt-watt;
- `NegativeSequenceAdmittanceDroop` — PWL `kappa(eta)`, fixed `phi2`, and an
  optional fixed ripple blend;
- `NoUnbalanceControl` — balanced-current baseline; and
- `CommonScaleLimiter` — the first smooth, fixed-structure protection law.

`ConverterCurrentTarget` and `GridCurrentTarget` select which physical current
the algebraic command constrains. In both cases the plant retains converter- and
grid-side per-leg ratings. The final plant-aware backoff includes the solved LCL
shunt-current offset and checks both locations, converter-terminal apparent
power, and a configured `dv2_max` before imposing the target equality.

Watt/var/proportional P/Q priority is implemented. Positive-sequence-power versus
negative-sequence-balance priority, virtual-delta, and closest-feasible
projection policies can be added without changing the plant or study interfaces.

Every law needs two evaluators with the same field semantics:

```julia
evaluate_exact(control, measurement, request, ratings)
stamp_smooth_control!(ctx, control, request, inverter, inverter_handles)
```

The exact evaluator is a pure numeric firmware oracle with exact PWL corners,
comparisons, and saturation. The smooth evaluator adds JuMP expressions and
equalities. Their output records should use the same names so exact-versus-smooth
residuals are automatic.

### 3. Plant/controller composition

Do not add every experimental control field to `AdvancedInverter`. Preserve its
backward-compatible plant API and compose it with a typed decorator:

```julia
struct ControlledDevice{D<:AbstractDevice,C<:AbstractInverterControlLaw} <:
       AbstractDevice
    device::D
    controller::C
end
```

The device lifecycle becomes:

```julia
inverter = inverter_spec(controlled)
plant_handles = stamp_device!(ctx, inverter)
inv_handles = inverter_handles(controlled, plant_handles)
control_handles = stamp_smooth_control!(
    ctx, controlled.controller, request, inverter, inv_handles)
```

`inverter_handles` and `inverter_spec` are small PowerOptLab traits for the
controlled wrapper. `extract_device` returns plant results, a control diagnostic
record, and converter-terminal phasors and powers.

`solve_controlled_inverter` is a controlled power flow, not an OPF over local
control actions: controller equalities determine the command. Its default loss
objective only selects remaining plant allocation freedom; `selection_objective=:zero`
is provided for objective-invariance tests. Reported capacitor loss can therefore
be a minimum-loss allocation within the plant abstraction, not a unique dynamic
prediction.

The first implementation supports `ControlledDevice{AdvancedInverter}` only.
This composition boundary keeps local-law data out of the physical plant.

Control parameters are immutable law data. Time-varying inputs such as PV
availability, enable flags, or firmware modes belong
to a separate `ControlRequest`/schedule passed by the problem builder. This keeps
one controller configuration reusable across thousands of snapshots and avoids
placing arbitrary functions inside reproducibility-critical device structs.

### 4. Problem and study orchestration

Add `src/problems/inverter_control_fleet.jl` for formulation-level work:

- stamp many controlled devices into one BMOPFTools context;
- bind per-snapshot availability/request data;
- choose the network objective without changing the local law;
- enforce KCL once all devices are present;
- extract device-local and network-wide results; and
- expose a single-snapshot solve suitable for warm starts and outer parallelism.

Keep dataset-specific ingestion, scenario generation, plots, and batch output in
`scripts/studies/inverter_controls/`. The library API should consume prepared
BMOPF networks and typed device/request collections rather than know the schema
of one research dataset.

Independent PV snapshots should be solved independently, warm-started along a
time or penetration trajectory, and parallelised outside Ipopt. Do not put a
year of unrelated operating points into one NLP unless a genuine inter-temporal
state is later introduced.

## BMOPFTools integration

The current public APIs are sufficient for most of the implementation:

- `build_opf_model` / `enforce_kcl!` / `extract_result` provide the staged
  lifecycle;
- `OpfBuildSpec` and per-component `OpfDeviceBuilder`s can replace selected
  native `ibr` entries while retaining all other network devices;
- public bus-voltage keys expose the rectangular terminal voltages;
- `add_terminal_injection!` supports custom KCL contributions;
- `OpfModelKey`, constraint registration, extension state, and result extractors
  support semantic diagnostics; and
- coefficient providers support fixed-structure per-snapshot availability and
  control parameters.

For datasets that already contain native BMOPF IBRs, prefer per-component
`OpfDeviceBuilder` replacement over silently stamping a second inverter at the
same bus. The study adapter maps each selected native IBR id to its
`ControlledDevice`; unselected IBRs remain native. Builder ownership then appears
in the BMOPFTools construction manifest.

Advanced-inverter parameters not present in the BMOPF schema belong to the
PowerOptLab study configuration during this research phase. Do not expand the
BMOPF schema until the control and hardware interfaces have stabilised.

## BMOPFTools dependencies

### Public smooth PWL expression API (complete)

The pinned BMOPFTools revision exposes stable exact/smooth numeric evaluation and
JuMP curve construction backed by its context-cached softplus operators.

BMOPFTools now exposes the supported context API used here:

```julia
opf_piecewise_linear_expression(
    ctx, input, breakpoints, values;
    epsilon)

piecewise_linear_value(input, breakpoints, values; epsilon=nothing)
```

The context method reuses the context's selected `softplus` mode and operator
cache, validates strictly increasing knots, and accepts inputs in working units.
The numeric method provides the exact firmware oracle when `epsilon=nothing` and
the corresponding smooth numerical oracle when a positive width is supplied.

This removes downstream duplication and makes the intended extension use of
BMOPFTools' smooth curves explicit. PowerOptLab's CI and documentation pin the
tested upstream commit.

### Public cached bus sequence expressions (useful, not blocking)

Add semantic keys and a helper that returns cached rectangular `U0`, `U1`, and
`U2` expressions for a three-phase terminal ordering. BMOPFTools already forms
similar private expressions for bus sequence limits, but downstream controllers
cannot request or reuse them directly.

PowerOptLab can compute the affine Fortescue transform from public voltage keys,
so this PR need not block the prototype. It becomes worthwhile when many
controllers, limits, and diagnostics at the same bus would otherwise build
duplicate expressions or use incompatible terminal ordering.

### PRs to defer

Do not initially upstream:

- the new control-profile schema;
- advanced-inverter plant fields;
- controller priority modes;
- hardware-sizing study APIs; or
- generic DC-coupled custom-IBR replacement changes.

Those interfaces are still research questions. A DC-coupling seam may become
necessary if a later implementation replaces native IBRs that use BMOPFTools'
`dc_bus` or `dc_link_coupled` features, but it is not required for the PV studies.

## Experimental methodology

### Stage A — algebra and controller equivalence

Test pure numeric helpers before constructing a network:

- Fortescue forward/inverse round trips and phase-order failures;
- exact PWL values at, between, and outside knots;
- smooth PWL error and derivatives around every knot;
- zero-`U2` regularisation and rotational invariance;
- reconstruction of phase currents from `I1`/`I2` with `I0=0`;
- ripple identity `S_tilde=3(U1 I2 + U2 I1)`;
- ripple-cancelling target and common-scale limiter; and
- exact versus smooth outputs across dense synthetic phasor grids.

### Stage B — one-inverter circuit tests

Use the existing balanced, pure-sequence, and strongly unbalanced
`AdvancedInverter` fixtures. For each controller:

- recover the legacy balanced solution when `U2=0`;
- verify that the negative-sequence droop has the intended sign and angle;
- trigger phase-current, ripple, capacitor-current, and switching limits one at
  a time;
- run both converter-current and grid-current targets through an LCL filter and
  verify both physical per-leg ratings;
- test SI/per-unit invariance;
- test smoothing and sample-grid refinement; and
- evaluate the exact firmware law at the solved point.

### Stage C — feeder equilibrium and scaling

Start with deterministic feeders and increase controlled-device count while
holding the electrical case comparable. Record:

- variables, equalities, inequalities, Jacobian/Hessian nonzeros;
- model-build, derivative, factorisation, and total solve time;
- Ipopt iterations and restoration events;
- peak memory;
- convergence under flat, previous-snapshot, and continuation starts; and
- exact/smooth controller residuals.

Use device counts that expose the scaling trend rather than selecting only one
large successful case. Preserve failed runs and their diagnostics.

### Stage D — matched penetration study

For each feeder and random seed, reuse the same customer placement, phase
assignment, load/PV trace, and weather sample across every
control law. Vary at least:

- PV penetration and inverter oversizing;
- phase allocation and load unbalance;
- feeder impedance and `R/X` range;
- negative-sequence gain and configured `phi2`;
- controller priority and smoothing;
- `Imax`, `Vdc`, `Cdc`, capacitor ripple-current rating, and DC-source 2ω
  impedance; and
- credible load/PV forecast errors.

Use paired differences between laws for each identical scenario. Report
distributions and tails, not only averages: hardware is usually selected by a
worst credible condition while energy and loss outcomes accumulate over time.

## Rating and capacitor-sizing method

### Current rating

For each operating point, retain converter-side and grid-side fundamental phase
currents plus the carrier audit. The semiconductor requirement is based on the
maximum physical converter-leg total RMS current; grid conductors and filter
elements have their own per-phase ratings. Neither may be replaced by `|I1|`,
aggregate apparent power, or negative sequence alone. Keep short-duration
overload and long-term thermal ratings distinct.

First compare laws at a common current rating. Then find the minimum rating that
satisfies a declared service criterion using an outer parameter sweep or
bisection around the network solve. This is easier to audit than making the
nameplate rating a decision variable inside the nonconvex OPF.

### DC-link capacitance

For a monolithic link with negligible source contribution at 2ω, the diagnostic
post-processing requirement is

```math
C_{2\omega,req}=
\max_t\frac{|\widetilde S_t|}
{2\omega V_{dc}\Delta V_{2,allow}}.
```

Use this formula for an unconstrained first pass. When `Cdc` changes dispatch,
switching feasibility, or limiter activation, repeat the network solve in an
outer `Cdc` sweep/bisection; post-processing alone is then insufficient.

Capacitance is not the capacitor's RMS-current rating. Separately calculate

```math
I_{2\omega,rms}=\frac{|\widetilde S|}{\sqrt2V_{dc}}
```

and combine it with switching current using frequency-dependent ESR/thermal
weights. Report the capacitance requirement, ripple-current requirement, loss,
and temperature/lifetime assumption separately. Also retain switching-ripple,
hold-up-energy, transient, tolerance, and ageing constraints; the selected bank
must satisfy the maximum of all capacitance requirements and every thermal
requirement.

For split links, retain each half-bank and midpoint mechanism. Do not apply the
monolithic three-leg formula to neutral-current heating.

### Capability frontiers

The most informative result is not one “minimum capacitor” number. For each
control law and penetration level, construct a frontier over:

- phase current rating;
- capacitance and capacitor-current rating;
- voltage-unbalance attenuation;
- worst-phase voltage compliance;
- curtailed energy; and
- converter/capacitor loss.

Use monotone outer sweeps where the physics is monotone and retain the complete
grid where it is not. This exposes whether a law buys voltage quality with
silicon, capacitance, curtailment, or loss.

## Reproducibility and validation gates

Every published scenario should retain:

- input dataset/version and scenario seed;
- controller type and all curve knots, smoothing scales, gains, angles, and
  priority settings;
- topology, filter, DC-source, PWM, and hardware ratings;
- BMOPFTools build manifest and research provenance;
- exact/smooth controller residuals;
- KCL, device, switching, current, and capacitor residuals; and
- solver options, starts, status, iterations, and timing.

The reviewed reference environment for the initial controller suite was Julia
1.12.6, JuMP 1.31.1, Ipopt 1.15.0, and BMOPFTools commit
`4c0ec8b9c947eea5cd94966f32d2c97f65115b87`. Examples and publication runs
should set `tol` explicitly (the initial cases use `1e-8`) rather than rely on a
solver-version default, and should record all feasibility tolerances separately
from controller smoothing widths.

Before drawing hardware conclusions, validate representative balanced,
unbalanced, current-limited, and ripple-limited points against an averaged model
and a switched EMT model. The network-wide NLP supplies breadth; the higher-
fidelity subset supplies confidence in the mechanisms used for sizing.

## Implementation sequence

Steps 1--5 are implemented by the public BMOPFTools smooth-PWL API and the
PowerOptLab phase-aware control component. The exact and smooth evaluators,
three voltage-reference policies, converter/grid current targets, dominant
conflict policy, converter-terminal power extraction, three-leg plant
composition, SI/per-unit regression, and extracted exact-versus-smooth current
residual are covered by unit tests. The remaining sequence starts with fleet
construction; batteries are outside this programme.

1. **Complete:** submit the BMOPFTools public smooth-PWL API PR.
2. **Complete:** add PowerOptLab exact sequence/curve/limiter numeric helpers and tests.
3. **Complete:** add `SequenceController`, `ControlledDevice{AdvancedInverter}`, and semantic
   control handles.
4. **Complete:** stamp the smoothed worst-phase and negative-sequence admittance laws into one
   advanced inverter; verify exact/smooth agreement.
5. **Complete:** add continuous normalised conflict handling, rated/available
   Volt-watt bases, P/Q priority, plant-side converter/grid current and
   apparent-power backoff, optional `dv2_max` backoff, and binding-limit tests.
6. Add multi-device single-snapshot construction using per-component
   `OpfDeviceBuilder` replacement for dataset IBRs.
7. Add tidy device/network result extraction and scaling benchmarks.
8. Run fixed-hardware PV penetration studies.
9. Add outer current/capacitor sizing sweeps and capability frontiers.
