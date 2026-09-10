# Compatibility review and migration

This release changes the PE port's feasible set. Reproducing an old numerical
study requires explicit modeling choices, not just the same input data and
solver version. `EVDevice` remains available; new studies can use `EV`, `EVSE`
and `ChargingSession` without moving conversion losses into the EVSE.

| Previous behavior | Current behavior | Migration decision |
|:--|:--|:--|
| Independent charge/discharge variables | Default `operation=:relaxed`, normalized product tolerance `1e-8` | Use `:complementarity` for exact MPCC semantics. Use `:independent` only when intentionally studying the former outer relaxation. |
| Free redistribution across phases | Equal phase P and equal phase Q | Set `phase_policy=:independent` explicitly if the hardware/model permits it. A session requires both vehicle and equipment permission. |
| Aggregate power bounds alone on multiphase ports | `i_max` is required | Supply a defensible per-conductor rating in A; do not infer it from aggregate kW without a voltage and phase-sharing assumption. |
| An inactive EV constrained only in aggregate P/Q | Every phase-current component is fixed to zero | No option restores nonzero cancelling currents on an unplugged device. |
| EV terminal inventory specified only by a floor | Optional `energy_final` equality | Match initial/final inventory when isolating recurring arbitrage economics. |
| `soc` in stored-energy units | `energy_wh` plus the existing `soc` alias | Both are Wh, not dimensionless electrochemical SoC. |

A 10 kW/10 kW port with normalized product tolerance `1e-8` can still have
simultaneous power of order 1 W before numerical error is added. Calling this
mode “relaxed” is intentional. Choose its tolerance from a physical budget and
inspect diagnostics; it is not exact mode exclusion.

## Introducing equipment identities

A session intersects vehicle and outlet power limits, applies outlet I/S limits,
and gates discharge/reactive support through explicit session permissions.
`EVSE.available` means operational; `ChargingSession.available` means occupied.
An outage cannot release an occupied outlet for another driver.

Use one EVSE identity for the same physical outlet. Different specifications
under the same id and overlapping occupancy are rejected. One EV id may have
only one session in a solve: repeated visits need a persistent mobility state,
which is not yet implemented. This prevents apparently feasible schedules that
silently reset stored energy between visits.

These are fixed AC installations. `phase_terminals` and `phase_count` are inputs,
not phase-switching variables. An empty outlet has no independent power
injection or reactive service. DC cabinets need a distinct conversion model.

## Extension authors and result consumers

Use keyword constructors for device specifications. Built-in port and linking
handles are implementation objects; the PE linking handle now contains
`energy` and `durations_h` so extraction can report physical energy residuals.
Custom device lifecycle methods retain their dispatch signatures. Consumers
should use `result.dispatch[id]` rather than unpacking internal handles.

`build_multi_context(...; configure! = callback)` configures a new model before
optimizer attachment, allowing complementarity bridges. With an existing model,
the caller retains ownership of optimizer attachment. `optimizer=nothing` still
supports building an unattached model.

Published dispatch follows the existing strict solve-status contract. Rejected
raw candidates remain unpublished, even if a subset of residuals looks small.
Conversely, a publishable local solve still requires physical residual checks.
Neither status nor these regressions certify global optimality.
