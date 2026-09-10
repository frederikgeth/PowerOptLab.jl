# Electric vehicles, equipment and sessions

> **Scope:** deterministic fixed-assignment AC charging, smooth NLP or NLP with
> complementarity. **Evidence:** analytical cases and regression tests; no
> fleet-wide calibration or field-performance claim.

A charging study needs to distinguish **what the driver needs**, **what the
vehicle accepts**, **what the equipment can supply**, and **what the controller
can command**. A departure-energy constraint addresses only the first question.
PowerOptLab represents the first three as a scheduling feasible set. It does not
yet simulate a pilot-response controller.

## Start with the question

| Question | Model and additional evidence |
|:--|:--|
| Can a parked vehicle obtain its requested energy? | PE energy account, session, equipment ratings and AC network |
| Which phase or conductor is limiting? | Explicit connection, per-conductor ratings and phase policy |
| Does charging taper change deadline feasibility? | Calibrated acceptance envelope and time-grid refinement |
| Is V2G profitable over repeated operation? | Equal or economically valued terminal state, losses, tariffs; aging for lifetime claims |
| Does the proposed dispatch track real equipment? | Requires a response/replay model and measured pilot/charging data |
| Can a fleet promise a network service? | Requires disaggregation, uncertainty and delivery-duration validation |

The [literature notes](literature.md) distinguish research findings from our
implementation choices. More detailed battery physics is useful only when the
data support it; see [PE versus IVQ](../tutorials/battery_storage_models.md).

## Object responsibilities

| Object | Owns | Does not imply |
|:--|:--|:--|
| `EV` | Capacity/bounds, AC-referred onboard limits, efficiency, phase capability, optional acceptance curve | Location or a reusable initial state for multiple visits |
| `EVSE` | One fixed AC outlet, conductor mapping, equipment P/Q/S/I ratings, operational availability | An onboard converter, a DC cabinet, or independent Q support while empty |
| `ChargingSession` | Fixed EV–EVSE association, occupancy, initial energy, departure target, V2G/Q permission | A driving model or an optimized assignment |
| `EVDevice` | Compact legacy combination of energy and network-port parameters | Equipment identity, occupancy conflict detection, trip energy |

AC conversion occurs onboard the vehicle; DC charging requires offboard
conversion and potentially shared power electronics. The present `EVSE` is
explicitly an **AC outlet envelope**. It does not double-count an onboard
converter by stamping another inverter at the outlet. Background and source:
[NREL EVSE technical overview](https://docs.nrel.gov/docs/fy21osti/78085.pdf).

## A complete fixed session

```julia
using PowerOptLab
car = EV(id="car", energy_max=40e3,
    onboard_charge_max=7e3, eff_charge=0.9)
outlet = EVSE(id="outlet", bus="poc", s_max=22e3, i_max=32.0)
visit = ChargingSession(id="visit", ev=car, evse=outlet,
    available=[false,true,true,false],
    energy_init=10e3, departure_energy=16.3e3)
result = solve_multiperiod_opf(nets, [visit];
    time_grid=TimeGrid([1.0,0.5,0.5,1.0]))
d = result.dispatch["visit"]
d.energy_wh, d.current_magnitude_a, d.diagnostics
```

Here `nets` contains four parsed AC networks with bus `poc`. A self-contained
network and executable version are supplied in `examples/ev_charging.jl` and
[the session tutorial](sessions.md).

`available` in a session means **occupied**, including periods when its EVSE is
out of service. The EVSE's separate mask controls electrical availability.
Sessions must occupy one contiguous interval; their departure is the end of the
last occupied interval. `departure_period` can be supplied explicitly but must
agree with that boundary. `energy_init` is a horizon-start stored-energy value;
energy is held before arrival. No driving energy is inferred.

Different vehicles can use the same outlet successively. Overlapping occupancy,
conflicting specifications with the same EVSE id, and repeated sessions for one
EV id are rejected. The last restriction prevents accidental state resets; it
is not a general mobility model. Distinct outlet ids at the same bus are allowed.

## Permissions and capability intersections

Charging cannot exceed either the onboard or outlet power rating. The AC
current and apparent-power constraints apply in addition. V2G requires a
positive vehicle discharge limit, positive EVSE discharge limit, and
`allow_v2g=true`. Q is zero unless `allow_reactive_power=true`; then the vehicle
and EVSE Q intervals are intersected. Measured charger reactive consumption is
not the same as controllable reactive capability.

For multiple phases, the default is equal phase P and equal phase Q. This is
**equal power**, not balanced current under unequal voltages. Independent phase
redistribution requires `phase_policy=:independent` on both the EV and EVSE.
Set `phase_count` on the EV to match the fixed outlet connection. There is no
phase-switching decision.

## Numerical and compatibility contract

- Powers and ratings are W/var/VA; energy is Wh; current is A; duration is hours.
- `energy_wh` is the explicit PE state name. `soc` retains the same Wh values
  for compatibility. Neither is IVQ's dimensionless charge-based SoC.
- `StorageDevice` and `EVDevice` now default to `operation=:relaxed`, a finite
  normalized charge/discharge product bound. This changes the old feasible set.
- `operation=:complementarity` supplies exact mathematical mode exclusion to an
  MPCC backend. `:independent` explicitly restores the outer relaxation.
- All multi-phase PE devices require `i_max`. Legacy single-phase devices may
  omit `i_max` and `s_max`, retaining a declared power-only abstraction.
- A disconnected port has every phase-current component fixed to zero.
- `solve_status(result).publishable` is a solver-status contract, not a physical
  accuracy or global optimality certificate. Inspect physical diagnostics.

See [equations and solver choices](model.md), [mode misconceptions](modes.md),
[approximation and flexibility](approximations.md), and [the roadmap](roadmap.md).
The [workplace feeder tutorial](workplace.md) demonstrates equipment-dependent
cost and deadline feasibility. Existing studies should read the
[compatibility and migration review](migration.md).
