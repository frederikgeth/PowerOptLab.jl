# Storage & EV devices

> **Kind:** Component model · **Maturity:** promotion candidate · **Direction:** forward · **Temporal:** inter-temporal (state of charge)

Devices are the reusable building blocks that
[`solve_multiperiod_opf`](@ref) stamps into each network snapshot. A device is
described by an SI-valued struct; the package handles the per-unit scaling, the
current-injection stamping into the engine's KCL, and the inter-temporal
state-of-charge linking.

## Model

Per snapshot, a device contributes a per-phase current injection `(cr, ci)` added
to the engine's KCL accumulators, so it draws or supplies real physical current.
Its aggregate AC injection is split into nonnegative charge and discharge power:

```math
P^{\text{inj}} = p^{\text{d}} - p^{\text{c}}, \qquad p^{\text{c}}, p^{\text{d}} \ge 0
```

with `pᶜ` the power drawn to charge and `pᵈ` the power delivered when
discharging (discharge is positive injection). The state of charge integrates
these with one-way efficiencies across a period of length ``\Delta t``:

```math
E_{t+1} = E_t + \left(\eta^{\text{c}}\, p^{\text{c}}_t - \frac{p^{\text{d}}_t}{\eta^{\text{d}}}\right)\Delta t,
\qquad E^{\min} \le E_{t+1} \le E^{\max}.
```

Round-trip loss (``\eta < 1``) makes simultaneous charging and discharging
suboptimal, so the split stays physical without an explicit complementarity
constraint.

This is an economic, not algebraic, exclusion: at unit efficiency or under
unusual negative-price objectives the split can be degenerate. Studies that need
an explicit operating-mode guarantee should add a complementarity/disjunctive
mode formulation rather than interpreting both nonnegative variables as a
certified physical mode.

The solve validates finite nonnegative power limits, ordered energy bounds,
efficiencies in `(0, 1]`, terminal existence in every snapshot, and an EV
availability entry for every interval before constructing the optimization
model.

## Battery / storage

```julia
using PowerOptLab
bat = StorageDevice(
    id            = "bat",
    bus           = "bus1",
    p_charge_max  = 40e3,     # W
    p_discharge_max = 40e3,   # W
    energy_max    = 100e3,    # Wh
    energy_init   = 40e3,     # Wh
    eff_charge    = 0.95,
    eff_discharge = 0.95,
    cyclic        = true,     # terminal SOC == initial SOC
)
```

Use `energy_final` to pin the terminal energy to a specific value instead of the
cyclic default, and `q_min`/`q_max` to allow reactive support (default is unity
power factor). Multi-phase inverters are supported via `phase_terminals`; the
`neutral` terminal is the return (`nothing` if referenced directly to ground).

## Electric vehicles

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

See the API reference for the full field list of [`StorageDevice`](@ref) and
[`EVDevice`](@ref).
