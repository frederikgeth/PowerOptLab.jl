# Power-energy battery storage

> **Kind:** Component model · **Maturity:** promotion candidate · **Direction:** forward · **Temporal:** inter-temporal (state of charge)

[`StorageDevice`](@ref) describes battery storage in power-energy (PE) space.
Start with the [battery storage overview](battery_storage.md) or the
[PE versus IVQ comparison](../tutorials/battery_storage_models.md) when choosing
model fidelity.

PE devices are the reusable building blocks that
[`solve_multiperiod_opf`](@ref) stamps into each network snapshot. A device is
described by an SI-valued struct; the package handles the per-unit scaling, the
current-injection stamping into the engine's KCL, and the inter-temporal
state-of-charge linking.

All package devices subtype [`AbstractDevice`](@ref) and use the same lifecycle:
[`validate_device`](@ref), [`stamp_device!`](@ref), [`link_device!`](@ref), and
[`extract_device`](@ref). New experimental devices can implement that interface
without adding another type switch to the multi-period builder.

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

Efficiencies do **not** guarantee exclusive charging/discharging. Bidirectional
PE devices now default to `operation=:relaxed`, imposing
`(p_charge/p_charge_max)*(p_discharge/p_discharge_max) <= 1e-8`.
This is a finite smooth NLP relaxation with reported mode residuals, not exact
physical exclusion. `operation=:complementarity` uses a native MPCC pair with
an external solver such as CCOpt. `:independent` explicitly restores the outer
relaxation for research comparisons. Zero-rated directions are fixed to zero.
See [the counterexample and solver tutorial](../ev/modes.md).

Optional `s_max`, `i_max` and `neutral_current_max` add aggregate apparent-power,
phase-current and neutral-current limits. Multi-phase PE ports require `i_max`.
`phase_policy=:equal_power` imposes equal P and Q on each phase; `:independent`
explicitly permits redistribution within the conductor limits. Disconnected
ports fix **every** conductor-current component, not just aggregate power.

The solve validates finite nonnegative power limits, ordered energy bounds,
efficiencies in `(0, 1]`, terminal existence in every snapshot, and an EV
availability entry for every interval before constructing the optimization
model.

## Usage

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

EV scheduling reuses this energy balance and port model, with availability and
departure requirements. See the [electric vehicle overview](../ev/index.md)
and the [compact EVDevice guide](../ev/legacy_device.md) for EV-specific usage.

## Dedicated vehicle, equipment and session objects

The [session guide](../ev/sessions.md) explains [`EV`](@ref), [`EVSE`](@ref) and
[`ChargingSession`](@ref), including fixed AC outlet assignments, capability and
permission intersections, and occupancy conflicts.

## Results and API

PE results expose `energy_wh` alongside the compatible Wh-valued `soc`, phase
currents/powers and independent physical residuals in `diagnostics`. Inspect
those and [`solve_status`](@ref) before interpreting a schedule. For the full
storage field list see the [battery storage API](battery_storage_api.md).
