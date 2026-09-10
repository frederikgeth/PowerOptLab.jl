# Tutorial: equipment ratings, energy requests and occupied outlets

## Experiment 1: a 22 kVA outlet and a 7 kW vehicle

**Misconception:** the equipment nameplate determines the charging power.

Use an ideal 230 V source so the device accounting has an analytical answer.
This removes feeder losses deliberately; feeder hosting studies must restore
the actual network. The executable example defines a 40 kWh vehicle, a 7 kW
onboard charger, 90% AC-to-stored-energy efficiency, a 22 kVA outlet with a 32 A
conductor limit, and one hour of occupancy split into two half-hour intervals.

```@example ev_sessions
using PowerOptLab
include(joinpath(dirname(pathof(PowerOptLab)), "..", "examples", "ev_charging.jl"))
r = EVChargingExamples.ratings()
d = r.dispatch["visit"]
@assert solve_status(r).publishable
@assert isapprox(d.energy_wh[end], 16_300; atol=0.1)
(d.p_charge, d.energy_wh, r.objective)
```

The candidate charging limits are 7 kW onboard, 22 kW from the outlet's active
rating, 22 kVA apparent power, and ``230\times32=7360`` W at unity PF. The
onboard limit binds. The battery gains

```math
0.9\times7000\times(0.5+0.5)=6300\text{ Wh}.
```

At 0.20 currency/kWh the grid energy costs 1.40, before any separately modeled
charges. The example stores the source price as 0.20 currency/kWh. BMOPFTools
converts its W-valued power expression to kW when forming the cost rate.
Do not divide the price by 1000 a second time.

## Experiment 2: a request is not a battery state

**Misconception:** requesting 6.3 kWh from an outlet guarantees a 6.3 kWh battery
increase.

In this example 7 kWh crosses the AC port and 6.3 kWh enters the stored-energy
account. `departure_energy=16.3e3` is the **absolute stored-energy floor**, given
10 kWh initially. It is neither an AC energy request nor a requested increment.
For AC-only session data, battery capacity and initial stored energy may be
unknown. Do not invent them and then call the result measured SoC.

The state indices are boundaries: occupancy is in intervals 2 and 3, so the
vehicle leaves at the end of interval 3, represented by `energy_wh[4]`. The
last horizon interval is disconnected and holds that state. Equipment outages
also hold this PE state but do not release the outlet's occupancy.

Try changing `eff_charge` to 0.8 while retaining the same deadline. The required
AC energy becomes 7.875 kWh, which cannot pass a 7 kW onboard charger in one
hour. The infeasibility is physical for this model; it should not be repaired by
relaxing the departure condition silently.

## Experiment 3: two vehicles, one outlet

**Misconception:** two independently feasible sessions are jointly feasible.

```@example ev_sessions
shared = EVChargingExamples.successive_visits()
@assert solve_status(shared).publishable
[(id, x.energy_wh) for (id,x) in shared.dispatch]
```

Both sessions reference the same outlet id. Their occupancy does not overlap,
so each receives its own energy account while the outlet's capability is reused.
Changing both occupancy masks to `[true,false]` is rejected before optimization.
Two distinct outlet ids on the same bus would instead require electrical
network or group limits to represent a shared upstream bottleneck.

A disconnected socket and a disabled socket occupied by a vehicle have different
assignment semantics. The EVSE outage mask disables electrical exchange without
making the parking space available.

## Experiment 4: zero net power is not isolation

**Misconception:** ``P=Q=0`` forces a three-phase device to be disconnected.

For nonzero phase voltages one can choose ``S_a=+2300`` VA,
``S_b=-2300`` VA and ``S_c=0`` and form nonzero currents
``I_\phi=\overline{S_\phi/V_{\phi n}}``. The aggregate complex power vanishes,
but two conductors still carry current. Such a fictitious device can redistribute
load across phases despite being reported as unplugged.

The regression in `test/evse_tests.jl` attempts to impose a 10 A real conductor
current on an unplugged three-phase EV. The corrected model rejects that
candidate; unconstrained disconnected solves return zero conductor currents in
both SI and per-unit coordinates. The connected test uses unequal source
voltages to check equal **power** sharing without assuming equal currents.

This also explains why an aggregate kVA circle needs conductor limits. When
reporting phase-balancing benefits, state whether the equipment supports
independent phase currents, equal powers, or another control law.
