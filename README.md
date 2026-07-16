> [!WARNING]  
> This project is currently ongoing rapid development and may have breaking changes made directly to main. Use at your own risk until further notice.
> Note that the bar for quality here is research contributions, not production grade code.

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

Contributions are organised by *what layer of the engine they extend* — see
[Concepts](docs/src/concepts.md).

**Component models** (`src/components/`) — new network elements via `model_hook!` / `solution_hook!`:

| Capability | Entry point | Reuses |
|---|---|---|
| **Storage / battery** devices with state of charge | [`StorageDevice`](src/components/devices.jl) | `model_hook!` current injection + KCL |
| **EV charging** (V1G / V2G) with availability & departure energy | [`EVDevice`](src/components/devices.jl) | storage device + per-period masking |
| **Advanced inverter** (internal-node IBR prototype) | [`AdvancedInverter`](src/components/advanced_inverter.jl) | internal node + filter / EMF / losses / ripple |

**Problem specifications** (`src/problems/`) — new formulations via the staged API:

| Capability | Entry point | Reuses |
|---|---|---|
| **Multi-period OPF** co-optimising many snapshots | [`solve_multiperiod_opf`](src/problems/multiperiod.jl) | staged API (one shared model) |
| **Legacy WLS state estimation** | [`solve_state_estimation`](src/problems/state_estimation.jl) | bounds-free physics + custom objective |
| **Constrained NLLS state estimation** | [`solve_sparse_state_estimator`](src/problems/constrained_state_estimation.jl) | compiled four-wire residual/constraint model |
| **Parameter estimation** (calibrate line lengths / taps) | [`solve_parameter_estimation`](src/problems/parameter_estimation.jl) | shared parameters across snapshots + WLS objective |
| **Dynamic operating envelopes** (active import/export capacity) | [`solve_operating_envelope`](src/problems/operating_envelope.jl) | operational bounds + scenarios + fairness policy |

**Bespoke algorithms** (`src/algorithms/`) — custom solution methods, currently
including [HELM power flow](docs/src/algorithms/helm.md).

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

### Calibrate line lengths / taps from smart-meter data

```julia
using PowerOptLab
# `nets` are per-snapshot physics nets (known source + transformers; uncertain
# lines omitted; no loads). `meas[t]` are that snapshot's noisy (P,Q,|V|) readings
# as `Measurement`s. Taps use the engine's native free-tap; lines are re-stamped
# with a shared free length. Follows Vanin, Geth et al., arXiv:2506.04949.
r = solve_parameter_estimation(nets, meas;
        lines=[CalibLine(id="l1", bus_from="src", bus_to="b1",
                         r_per_length=0.4, x_per_length=0.25)],
        taps =[CalibTap(id="t1", tap_min=0.9, tap_max=1.1)],
        objective=:wls)       # or :wlav (robust, bad-data-rejecting)
r.line_length["l1"]   # estimated line length
r.tap["t1"]           # estimated tap multiplier τ
r.residual_rms        # RMS voltage misfit (V); near the meter noise floor if the fit is good
```

Diverse loads across multiple time steps are what make the shared parameters
identifiable — a single snapshot cannot separate a long line from a heavy load.
`per_unit=true` and `per_unit=false` give identical estimates.

### Dynamic operating envelopes

```julia
using PowerOptLab
cps = [ConnectionPoint(id="der1", bus="bus1", export_max=10e3),
       # Bind realistic PV/batteries to an existing BMOPFTools IBR so its
       # mandatory Q-V law and converter limits remain active.
       ConnectionPoint(id="der2", bus="bus2", ibr_id="pv2", export_max=10e3)]

env = solve_operating_envelope(nets, cps;
    direction=:export,
    fairness=FairnessPolicy(kind=:max_min, normalization=:capacity),
    security=:bound_point)  # use :corners for every zero/full box corner

env.envelope["der1"]  # positive directional capacity per interval (W)
env.total_capacity
env.diagnostics
```

The solver supports export or import, forecast/model scenario groups, weighted
and normalized fairness policies, and explicit bound-point versus all-corner
security semantics. Loads retain snapshot P/Q; connection-bound IBRs retain
their prescribed Volt-VAr/Volt-Watt laws. Add a BMOPFTools STATCOM to a copy of
the network to quantify its impact on the active-power envelope.

### Advanced inverter (prototype internal-node IBR)

A more detailed IBR than the engine's built-in current-injection model: an
explicit internal AC node behind the converter, with an output filter, internal
EMF / DC-link modulation bounds, grid-forming operation, converter losses, and
double-frequency ripple limits (following the BMOPFTools
[IBR model extensions](https://github.com/frederikgeth/BMOPFTools.jl/blob/main/docs/ibr_model_extensions.md)
design doc). For three phases it carries the exact **switching-polytope**
feasible-region models of the 3-leg, 4-leg, and split-DC-link topologies (exact
time-sampled DC-utilisation limits, 2ω bus-ripple derating, neutral-current
limits). Every feature is opt-in; runs in SI or per-unit.

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

BMOPFTools is not yet registered. For development, use a local checkout; CI and
documentation builds pin commit `c8df8f353637a46a70f60fa95c2d2184e7475d39`
so their dependency source is reproducible:

```julia
using Pkg
Pkg.activate(".")
Pkg.develop(path="../BMOPFTools.jl")   # adjust to your checkout
Pkg.instantiate()
Pkg.test()
```

`Manifest.toml` is intentionally not committed (library convention); the
package compat pins BMOPFTools 0.1.0 exactly, while CI pins the tested source
commit shown above.

## License

Dual-licensed:

- **Source code** — BSD-3-Clause. See [LICENSE.md](LICENSE.md).
- **Data files** (test cases, benchmark networks, and any dataset shipped with
  the package) — Creative Commons Attribution 4.0 International (CC BY 4.0). See
  [LICENSE-DATA.md](LICENSE-DATA.md).

This mirrors the BMOPFTools convention: code under a permissive software license,
data under CC BY so cases can be shared and cited with attribution.
