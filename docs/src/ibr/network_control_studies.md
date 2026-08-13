# Network-scale control studies

`solve_controlled_inverter_fleet` is the single-snapshot bridge between the
phase-aware controller/advanced-inverter model and prepared BMOPFTools network
datasets. It is intended for comparative studies of local control laws,
per-leg current requirements, and DC-link stress across many independent LV
network snapshots.

## Scientific scope

The solver answers a deliberately narrow question: given one network state,
one physical inverter design per selected IBR, one locally measured
voltage-phasor control law, and one available-power request, what simultaneous
network and inverter operating point satisfies all algebraic equations and
physical limits?

It is a controlled power flow, not model-predictive control and not a local
optimization. Every controller is a fixed smooth computation graph. Ipopt is
used to solve the coupled network/device equations at scale; it does not choose
the control policy. Independent time, penetration, and uncertainty snapshots
belong in an outer study loop so failures and limit activation remain visible.

The initial fleet implementation supports selected dataset IBR records with
`THREE_LEG` or `FOUR_LEG` AC port descriptions, each exposing three ordered
phase conductors. Their PowerOptLab replacements are three-leg
`AdvancedInverter`s. Native DC-network-coupled IBRs are rejected because
replacing their AC formulation would otherwise remove native DC balance
physics.

Selecting a native `FOUR_LEG` record therefore means deliberately studying a
three-leg replacement at the same four-wire connection; it is not a claim that
the hardware topologies are equivalent. A native `THREE_LEG` record is accepted
only when the replacement has `neutral=nothing` and has neither an LCL midpoint
shunt nor a POC shunt. These conditions preserve the native line-to-line port:
the replacement has zero summed phase current and introduces no phase-to-ground
or neutral current path.

## Replacement, not parallel stamping

The dataset continues to identify an inverter by `network["ibr"][id]`. A
`ControlledInverterFleetSpec` maps that same identifier to a
`ControlledDevice` and an `InverterControlRequest`:

```julia
spec = ControlledInverterFleetSpec(
    Dict("pv_17" => ControlledDevice(advanced_inverter, controller)),
    Dict("pv_17" => InverterControlRequest(
        p_available=6.4e3, p_rated=8e3, q_scale=8e3)))

result = solve_controlled_inverter_fleet(network, spec)
```

PowerOptLab assigns each selected `(:ibr, id)` to a per-component
`OpfDeviceBuilder`. BMOPFTools therefore omits the native IBR equations and KCL
injection for that identifier while retaining native physics for every
unselected IBR. The advanced plant adds the only live terminal injection.

BMOPFTools declares its semantic native IBR current variables before ownership
is resolved. PowerOptLab fixes those unused placeholders to zero. Published
results then follow one unambiguous ownership rule:

- `result.devices[id]` contains each replaced, controlled advanced inverter;
- `result.network["ibr"]` contains only unselected native IBRs; and
- `result.build_manifest.component_owners[(:ibr, id)]` records `:PowerOptLab`
  or `:BMOPFTools`.

This contract is stronger than inferring replacement from a coincident bus or
from a near-zero duplicate injection. A selected native record with tight
native power bounds cannot constrain the replacement plant. The number of
placeholders fixed is derived from BMOPFTools' resolved dataset neutral labels,
not assumed to equal the replacement's three phases. This also closes every
unused variable when a non-conventional neutral spelling makes BMOPFTools
declare an additional semantic IBR-current pair.

One upstream qualification remains: if BMOPFTools does not recognise that
neutral spelling, its native `:network_cost` construction also treats the
terminal as a fourth phase. A costed record must consequently provide four cost
entries or BMOPFTools rejects it while indexing the objective. This affects
native and replaced IBRs alike; PowerOptLab closes the extra unused current
pair, but does not reinterpret native cost data. Prefer BMOPFTools-recognised
neutral metadata/spelling for production datasets.

## Input contract

For each selected identifier, construction validates before model mutation:

1. the device and request dictionaries have identical nonempty key sets;
2. the dictionary key equals the wrapped advanced-inverter identifier;
3. the identifier exists in the network's native `ibr` collection;
4. native and replacement bus identifiers agree;
5. phase order agrees exactly, and a `FOUR_LEG` record places the matching
   neutral after the three phases; and
6. a `THREE_LEG` record is paired only with a neutral-free, shunt-free
   replacement; and
7. the native IBR has no DC-network coupling.

Phase order is a physical convention, not presentation metadata: it determines
the symmetrical-component transform and the reconstructed phase currents.
Silently sorting phase labels would change the control law.

The controller formulation is currently supported only with `per_unit=true`.
Configuration and extracted results nevertheless remain in SI units. The
advanced inverter must use `pwm_strategy=:NONE`; calibrated PWM reserve terms
may be supplied explicitly, but the separate PWM outer iteration is not part of
the fleet solve.

## Objective semantics

Local controller equalities determine their current commands. The objective
only selects freedom left elsewhere in the algebraic model:

- `selection_objective=:loss` (default) minimizes selected inverter
  semiconductor and capacitor losses;
- `:zero` is an objective-invariance diagnostic; and
- `:network_cost` applies BMOPFTools' native generation-cost expression to
  dispatchable unselected network devices.

There is intentionally no silent weighted sum of economic cost and electrical
loss: the quantities have no universal commensurate scale. For comparative
control-law studies, pin exogenous loads, generation, and unselected IBR
setpoints in each prepared snapshot, then use `:loss`. If dispatch of native
devices is part of the experiment, state and use the economic objective
explicitly. Results from different objective semantics are not interchangeable.

## Result tables and aggregation

`ControlledInverterFleetResult` retains both the device-local results and the
ordinary BMOPFTools network result, plus a defensive copy of the study
specification used for provenance. `solve_status(result)` is the authoritative
publication gate. Do not drop infeasible or non-publishable snapshots; retain
their scenario identifiers and classify the failure mechanism.

The bus snapshot is canonicalized once per fleet result and shared by
`result.network["bus"]` and every `result.devices[id].bus`. Treat extracted
results as immutable: mutating that shared dictionary changes every view of the
same snapshot.

Two dependency-free extraction helpers produce table-ready named tuples:

```julia
device_rows = controlled_inverter_rows(result)
phase_rows = controlled_inverter_phase_rows(result)
```

`device_rows` has one row per inverter and scalar SI columns for POC and
converter power, losses, voltage extrema and sequences, VUF, requested power,
limiter scales, maximum converter/grid current, sequence current, 2ω ripple,
DC-link ripple voltage, capacitor current, and the exact-versus-smooth current
residual. `phase_rows` has one row per inverter and phase, retaining rectangular
and magnitude voltage, converter-current, and grid-current values. Both helpers
sort identifiers lexically, independent of dictionary insertion order.

The conventional `voltage_unbalance_factor` is ``|U_2|/|U_1|``. It is distinct
from `regularized_voltage_unbalance`, the controller's finite-denominator
quantity. Keeping both prevents the controller's numerical regularization from
being mistaken for a reported power-quality metric.

`solve_diagnostics(result)` supplies lightweight fleet aggregates, including
total controlled POC P/Q, total selected-inverter losses, the largest
converter/grid phase current, and the largest exact-versus-smooth current
residual. Study scripts should add scenario keys—time, feeder, penetration,
random seed, and hardware/control configuration—before writing rows.

## Verification obligations

The unit tests exercise discriminative rather than merely feasible cases:

- selected native records have near-zero native ratings while their replacement
  plants deliver kilowatts, detecting accidental parallel native stamping;
- selected and unselected IBRs coexist, with ownership and result channels
  checked explicitly;
- average-voltage and worst-phase Volt-watt laws share one unbalanced network
  state and produce materially different requests;
- output ordering and aggregate identities are deterministic; and
- each solution satisfies three-leg zero-sum current, per-leg current limits,
  converter-terminal complex-power reconstruction, filter active-power balance,
  and DC-link power balance.

These are software and model-consistency checks, not validation of a hardware
controller. The study methodology still requires a stratified subset of
balanced, unbalanced, current-limited, and ripple-limited cases against an
independent averaged model and a switched EMT model. OpenDSS is suitable as a
network-voltage oracle for cases whose inverter injection can be frozen, but it
cannot independently validate the advanced converter, DC capacitor, or the
simultaneous smooth controller equations. Use it to isolate network stamping;
use averaged/EMT models for converter and control-law validation.

## Study workflow

For a publishable comparison, freeze a factorial design before inspecting
results: feeder/scenario, PV penetration and phase allocation, irradiance,
source-voltage unbalance, load realization, hardware rating, current target,
positive-sequence policy, negative-sequence policy, and limiter settings. Use
the same snapshot and hardware across controller variants.

Report at least:

- voltage compliance and phase extrema;
- conventional VUF and negative-sequence current;
- energy curtailed relative to the declared available-power trajectory;
- empirical quantiles and worst cases of converter- and grid-side per-leg
  current;
- frequency and identity of binding physical limits;
- 2ω ripple power, DC voltage ripple, capacitor RMS/thermal current, and
  capacitor stored energy ``C_{dc}V_{dc}^2/2``; and
- solver status, exact-versus-smooth residual, and sensitivity to smoothing
  widths on a representative subset.

Current and capacitor sizing conclusions require outer parameter sweeps. A
feasible point at one installed rating does not identify the minimum rating,
and a smooth-NLP failure is not by itself a physical infeasibility certificate.

## API

```@docs
ControlledInverterFleetSpec
ControlledInverterFleetResult
solve_controlled_inverter_fleet
controlled_inverter_rows
controlled_inverter_phase_rows
```
