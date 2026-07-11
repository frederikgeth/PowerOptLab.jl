# PowerOptLab.jl

[![CI](https://github.com/frederikgeth/PowerOptLab.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/frederikgeth/PowerOptLab.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://frederikgeth.github.io/PowerOptLab.jl/dev/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](LICENSE.md)

A staging ground for **component models** and **problem specifications** built on
top of the [BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl) reference
optimal-power-flow engine.

BMOPFTools ships one small, correct four-wire rectangular current–voltage OPF and
deliberately refuses to become a "model zoo". PowerOptLab is where the zoo lives:
experimental devices and alternative problem formulations that reuse the engine's
device physics, per-unit handling, and result extraction **through its public
extension seams** (`model_hook!` / `solution_hook!` and the staged
`build_opf_model` / `enforce_kcl!` / `generation_cost` / `extract_result` API) —
without forking the engine. Anything here that matures into accepted practice can
later be folded back into the BMOPF spec.

## What's inside

| Capability | Entry point | Reuses |
|---|---|---|
| **Storage / battery** devices with state of charge | [`StorageDevice`](src/devices.jl) | `model_hook!` current injection + KCL |
| **EV charging** (V1G / V2G) with availability & departure energy | [`EVDevice`](src/devices.jl) | storage device + per-period masking |
| **Multi-period OPF** co-optimising many snapshots | [`solve_multiperiod_opf`](src/multiperiod.jl) | staged API (one shared model) |
| **State estimation** (weighted least squares) | [`solve_state_estimation`](src/state_estimation.jl) | bounds-free physics + custom objective |
| **Dynamic operating envelopes** (DER export limits) | [`solve_operating_envelope`](src/operating_envelope.jl) | operational bounds + fairness objective |
| **Advanced inverter** (internal-node IBR prototype) | [`AdvancedInverter`](src/advanced_inverter.jl) | internal node + filter / EMF / losses / ripple |

## Examples

### Battery arbitrage across two prices

```julia
using PowerOptLab
# `nets` are BMOPFTools net dicts, one per period, with a time-varying slack
# import price (voltage_source `cost`).
bat = StorageDevice(id="bat", bus="bus1",
                    p_charge_max=40e3, p_discharge_max=40e3,   # W
                    energy_max=100e3, energy_init=40e3,        # Wh
                    eff_charge=0.95, eff_discharge=0.95, cyclic=true)
res = solve_multiperiod_opf(nets, [bat]; dt_h=1.0)
res.dispatch["bat"].p_net    # SI discharge (+) / charge (−) per period
res.dispatch["bat"].soc      # SI state of charge (Wh) at each step boundary
```

### EV charging with a departure deadline

```julia
ev = EVDevice(id="ev1", bus="bus1", p_charge_max=20e3,        # V1G (charge only)
              energy_max=40e3, energy_init=10e3,
              available=[true, true, true, false],            # unplugged in period 4
              departure_energy=30e3, departure_period=3)      # ≥30 kWh by end of period 3
res = solve_multiperiod_opf(nets, [ev]; dt_h=1.0)
```

Set `p_discharge_max > 0` for bidirectional (V2G) charging.

### State estimation from noisy measurements

```julia
using PowerOptLab
meas = [
    Measurement(kind=:vmag, bus="bus1", value=979.5, sigma=2.0),   # volts
    Measurement(kind=:pinj, bus="bus1", value=-20_000.0, sigma=400.0),  # watts
    # …
]
# `net` is a physics-only net (buses, lines, source; no operational limits).
se = solve_state_estimation(net, meas)
se.bus["bus1"]["1"]["vm"]   # estimated voltage magnitude (V)
se.residuals                # per-measurement (measured, estimated, residual, normalized)
```

### Dynamic operating envelope (DER export limits)

```julia
using PowerOptLab
cps  = [ConnectionPoint(id="der1", bus="bus1", export_max=10e3),   # W
        ConnectionPoint(id="der2", bus="bus2", export_max=10e3)]
# `nets` are per-interval snapshots (differing baseline loads). Each interval's
# envelope respects that interval's voltage/thermal limits.
env = solve_operating_envelope(nets, cps; fairness=:equal)  # or :sum, :proportional
env.envelope["der1"]   # allocated export limit per interval (W)
env.total_export       # total allocated across connection points, per interval
```

`:equal` allocates the same limit to every point (equitable); `:sum` maximises
the total (efficient, but may starve electrically weaker points); `:proportional`
maximises `sum(log(pₑ))` — a middle ground where no point is starved but stronger
points still get more.

### Advanced inverter (prototype internal-node IBR)

A more detailed IBR than the engine's built-in current-injection model: an
explicit internal AC node behind the converter, with an output filter, internal
EMF / DC-link modulation bounds, grid-forming operation, converter losses, and
double-frequency ripple limits (following the BMOPFTools
[IBR model extensions](https://github.com/frederikgeth/BMOPFTools.jl/blob/main/docs/ibr_model_extensions.md)
design doc). Every feature is opt-in.

```julia
using PowerOptLab
inv = AdvancedInverter(id="inv", bus="poc", s_max=5e3,   # VA converter rating
                       r_filter=0.2, x_filter=0.5,        # Ω output filter
                       p_loss_fixed=20.0, a_loss=0.3, c_loss=0.02)  # 3-term loss
r = solve_advanced_inverter(net, inv; objective=:min_loss, p_set=3e3)
r.p_poc    # power delivered to the grid (W)
r.p_dc     # DC-link power = p_conv + p_loss  (non-branching)
r.v_int_mag  # internal EMF magnitude per phase (V)
```

## Development setup

BMOPFTools is not yet registered, so develop it from a local checkout:

```julia
using Pkg
Pkg.activate(".")
Pkg.develop(path="../BMOPFTools.jl")   # adjust to your checkout
Pkg.instantiate()
Pkg.test()
```

`Manifest.toml` is intentionally not committed (library convention); the
environment is reproducible from `Project.toml` plus the dev-path above.

## License

Dual-licensed:

- **Source code** — BSD-3-Clause. See [LICENSE.md](LICENSE.md).
- **Data files** (test cases, benchmark networks, and any dataset shipped with
  the package) — Creative Commons Attribution 4.0 International (CC BY 4.0). See
  [LICENSE-DATA.md](LICENSE-DATA.md).

This mirrors the BMOPFTools convention: code under a permissive software license,
data under CC BY so cases can be shared and cited with attribution.
