# Electric vehicle API

For equations, limitations and numerical semantics see [the model guide](model.md).

The shared PE battery constructor is documented in the
[battery storage API](../components/battery_storage_api.md).
The [multi-period scheduling API](../api/multiperiod.md) documents the solver
that accepts `EVDevice` and `ChargingSession` objects.

```@docs
EVDevice
EV
EVSE
ChargingSession
```
