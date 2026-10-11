# Compact EVDevice model

> **Kind:** Component model · **Maturity:** promotion candidate · **Direction:** forward · **Temporal:** inter-temporal

[`EVDevice`](@ref) combines stored energy, an AC network port, availability and
one departure requirement in a single specification. It remains available for
compact studies and compatibility. For explicit vehicle and outlet identities,
use [`EV`](@ref), [`EVSE`](@ref) and [`ChargingSession`](@ref) as described in the
[electric vehicle overview](index.md) and [session guide](sessions.md).

## A compact charging schedule

An [`EVDevice`](@ref) is a storage device that is only controllable while plugged
in and must reach a target energy by departure:

```julia
ev = EVDevice(
    id               = "ev1",
    bus              = "bus1",
    p_charge_max     = 20e3,               # W  (V1G: no discharge)
    energy_max       = 40e3,               # Wh
    energy_init      = 10e3,               # Wh
    available        = [true, true, true, false],  # unplugged in period 4
    departure_energy = 30e3,               # Wh required …
    departure_period = 3,                  # … by the end of period 3
)
```

While a period is unavailable the charger is idle and the state of charge is
held. Set `p_discharge_max > 0` for bidirectional (V2G) operation, letting the
vehicle discharge into expensive periods subject to its departure requirement.

See the [EV API](api.md) for the full field list of [`EVDevice`](@ref).

## Shared physics and migration

The device uses the [PE storage energy balance](../components/devices.md).
Its availability and departure constraints are described in the
[scientific model](model.md). Read the [operating-mode guide](modes.md) before
using V2G, and the [compatibility guide](migration.md) when reproducing studies
that used earlier phase-sharing or charge/discharge defaults.
