# Tutorial: what changes when an outlet has an identity?

Consider a workplace with two successive visitors and one all-day employee.
Each car has a 7 kW onboard charging limit. The visitors share one 7 kVA outlet,
but its conductor rating is 16 A. The employee uses a separate 32 A outlet.
Does a 7 kW nameplate make these charging opportunities equivalent?

## A deliberately interpretable feeder

`examples/ev_workplace.jl` supplies a complete three-phase AC network and the
sessions. Source magnitudes are 245/215/230 V, each phase has a 0.2 ohm resistive
service conductor, and office loads are 0.5/1.0/1.5 kW. The visitor outlet is on
the weaker 215 V phase; the employee is on the 245 V phase.

The neutral is ideal and equipotential. Thus phase-to-ground and phase-to-neutral
voltage coincide here; this example cannot establish neutral displacement,
neutral heating or general four-wire feeder behavior. The parameters are
illustrative, not a measured workplace or a standards-compliance example.

There are four one-hour intervals with source prices 0.30/0.10/0.10/0.30 in
currency/kWh. The morning visitor occupies intervals 1–2, the afternoon visitor
3–4, and the employee 1–4. Each visitor needs 5.5 kWh of **stored energy**; the
employee needs 10 kWh. All charging efficiencies are 0.9. Initial energy is
10 kWh, and all capacities are 40 kWh.

```@example workplace
using PowerOptLab
include(joinpath(dirname(pathof(PowerOptLab)), "..", "examples", "ev_workplace.jl"))
base = EVWorkplaceExample.study()
@assert solve_status(base.result).publishable
EVWorkplaceExample.summary(base)
```

The outlet is occupied continuously, but its identity links two different
vehicles and deadlines. There is no shared vehicle energy state between them.
The model checks occupancy conflicts and intersects each vehicle's capability
with the same equipment limits.

## Nameplate kW does not determine charging cost

The morning visitor cannot obtain all its required energy in cheap interval 2:
current is limited even though onboard and outlet power ratings are 7 kW. Some
charging must occur in expensive interval 1. The afternoon visitor faces the
same issue in intervals 3–4. No extra objective penalty is needed to expose this.

```@example workplace
upgraded = EVWorkplaceExample.study(weak_current=32.0)
@assert solve_status(upgraded.result).publishable
@assert base.result.dispatch["morning"].p_charge[1] > 1000
@assert upgraded.result.dispatch["morning"].p_charge[1] < 1
@assert upgraded.result.objective < base.result.objective
(base.result.objective, upgraded.result.objective)
```

Only the outlet current rating changed. All requests, prices, efficiencies and
power nameplates stayed fixed. The feasible sets are nested, so the global
optimum cannot worsen after this upgrade; these are local NLP solutions, and
we additionally check the observed cost difference. Reported cost includes
**office demand and AC line losses**, not just EV charging. The difference is an
operating-cost comparison, not an investment appraisal or a universal upgrade
benefit. Adding upgrade costs would answer a different question.

At Q=0, current magnitude is ``p/|V|``. Equal power at unequal voltages therefore
requires unequal currents. Actual currents and voltage drops come from the AC
model, rather than converting 7 kW into amperes using an assumed 230 V.

## A necessary bound can prove a deadline impossible

Raise each visitor's request to 8 kWh. A 7 kW-only view suggests that two hours
are ample: ``0.9\times7\times2=12.6`` kWh. But this fixture bounds outlet voltage
by 250 V. Even using that optimistic voltage throughout, a 16 A outlet satisfies

```math
\Delta E \le 2\times0.9\times\min(7000,250\times16)
          =7200\text{ Wh}<8000\text{ Wh}.
```

This is an analytical necessary bound on the original model, independent of
whether a local NLP solver can diagnose infeasibility. Actual weak-phase voltage
and line losses only make this optimistic bound less attainable in this case.

```@example workplace
@assert EVWorkplaceExample.visitor_energy_bound() == 7200
impossible = EVWorkplaceExample.study(visitor_gain=8e3)
@assert !solve_status(impossible.result).publishable
possible = EVWorkplaceExample.study(weak_current=32.0,visitor_gain=8e3)
@assert solve_status(possible.result).publishable
(impossible.result.termination_status, possible.result.termination_status)
```

A failed local solve alone would not be a proof. Here the bound supplies the
proof for 16 A, while the 32 A case supplies a numerically feasible trajectory
for the changed installation. Published results are also checked for departure
shortfall, energy conservation, current limits and exact disconnection in the
regression tests.

## What this study does and does not establish

This example motivates equipment identity, actual conductor ratings, fixed phase
connections, consecutive occupancy and individual deadlines. It assumes known
arrivals, continuous charging commands and constant efficiencies. A real study
needs measured charger response, availability uncertainty, installation data,
and appropriate feeder validation; see the [literature notes](literature.md).
Shared DC cabinets, discrete pilots and driving between visits require the
additional models described in the [roadmap](roadmap.md).
