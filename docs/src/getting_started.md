# Your first PowerOptLab study

This guide installs the package, builds a complete small network, schedules a
battery across two electricity prices, and checks the returned result. It uses
synthetic data and the standard Ipopt dependency; no external dataset or optional
solver is needed. The Julia examples on this page execute during documentation
builds.

## Install from a checkout

Use Julia 1.10 or later. Clone the repository, enter it, and install the pinned
dependencies:

```sh
git clone https://github.com/frederikgeth/PowerOptLab.jl.git
cd PowerOptLab.jl
julia --project=. scripts/instantiate_pinned.jl
julia --project=.
```

If you already have a checkout, run the last two commands from its root.
The setup script installs the BMOPFTools revision declared in `Project.toml`;
BMOPFTools is unregistered. It replaces any local development-path override in
this project's manifest. The first installation and compilation take longer
than subsequent runs.

Paste the following Julia blocks into that REPL in order. You can also save them
together in a script and run `julia --project=. your_script.jl`.

## Build two network snapshots

The network is a 230 V single-phase circuit with an explicit phase and neutral
line. Only the source neutral is grounded; the load-bus neutral is a network
variable. A 2 kW load is present in both one-hour intervals. The only changed
input is the source's price, in currency/kWh.

```@example first_run
using PowerOptLab
using BMOPFTools: parse_bmopf

function feeder(price)
    parse_bmopf("""
    {
      "bus": {
        "grid": {"terminal_names": ["1", "n"],
                 "perfectly_grounded_terminals": ["n"]},
        "house": {"terminal_names": ["1", "n"]}
      },
      "voltage_source": {
        "supply": {"bus": "grid", "terminal_map": ["1"],
                   "v_magnitude": [230.0], "v_angle": [0.0],
                   "cost": [$price]}
      },
      "linecode": {
        "cable": {"R_series_1_1": 0.05, "R_series_2_2": 0.05}
      },
      "line": {
        "service": {"bus_from": "grid", "bus_to": "house",
                    "terminal_map_from": ["1", "n"],
                    "terminal_map_to": ["1", "n"],
                    "linecode": "cable", "length": 1.0}
      },
      "load": {
        "demand": {"bus": "house", "terminal_map": ["1", "n"],
                   "configuration": "SINGLE_PHASE",
                   "p_nom": [2000.0], "q_nom": [0.0]}
      }
    }
    """; from_string=true)
end

nets = [feeder(0.30), feeder(0.10)]
grid = TimeGrid([1.0, 1.0])
nothing # hide
```

Terminal labels are strings and the terminal-map order matches the linecode
matrix order. Network powers are in W, voltages in V and impedances in Ω.
This small single-phase example introduces the workflow; it is not a study of
three-phase unbalance or an equipment-calibrated feeder.

## Add storage and solve

The battery can charge or discharge at 1 kW, starts at 2 kWh, and must finish at
the same stored energy. Its efficiencies refer to the AC-to-stored-energy
boundary. The solve includes line losses and couples the intervals through the
energy balance.

```@example first_run
battery = StorageDevice(
    id="battery", bus="house",
    p_charge_max=1000.0, p_discharge_max=1000.0,
    energy_max=4000.0, energy_init=2000.0,
    eff_charge=0.95, eff_discharge=0.95,
    cyclic=true, i_max=10.0, s_max=1000.0,
)

result = solve_multiperiod_opf(
    nets, [battery]; time_grid=grid,
    solver_options=(tol=1e-8,),
)

status = solve_status(result)
(termination=status.termination_status, publishable=status.publishable)
```

Ipopt solves a nonconvex nonlinear model. A published local solution is useful
for this example, but its status does not establish global optimality.
The default `operation=:relaxed` bounds simultaneous charging and discharging
by a finite normalized product tolerance; inspect the physical residuals as
well as the status.

## Check the status before reading dispatch

```@example first_run
if !status.publishable
    error("No published dispatch: $(status.termination_status)")
end

dispatch = result.dispatch["battery"]
(
    net_injection_w=round.(dispatch.p_net; digits=2),
    stored_energy_wh=round.(dispatch.energy_wh; digits=2),
    horizon_cost=round(result.objective; digits=4),
)
```

Positive `p_net` means injection into the network (discharging); negative means
charging. Expect discharge in the expensive interval and recharge in the cheap
one. There are two power values and three energy values: the initial state and
the end of each interval. `energy_wh` is stored energy in Wh; the compatible
`soc` field has the same units, rather than being a fraction.

On a failed solve, the wrapper retains status information and masks unpublished
numerical results with `NaN`. Changing the solver settings without checking that
status can make a failed candidate look like a schedule.

## Check physical residuals and a network voltage

The following checks use illustrative budgets of 0.01 Wh for energy and 0.1 W
for power. Choose budgets from the accuracy required by your own study.

```@example first_run
diagnostics = dispatch.diagnostics
@assert maximum(abs.(diagnostics.energy_balance_wh)) < 0.01
@assert abs(diagnostics.terminal_energy_error_wh) < 0.01
@assert maximum(abs.(diagnostics.ac_power_balance_w)) < 0.1
@assert maximum(diagnostics.power_limit_violation_w) < 0.1
@assert maximum(abs.(diagnostics.simultaneous_power_w)) < 0.1
@assert dispatch.p_net[1] > 0 && dispatch.p_net[2] < 0

# Each snapshot reports conductor potentials in SI units.
bus = result.snapshots[1]["bus"]["house"]
phase_voltage = complex(bus["1"]["vr"], bus["1"]["vi"])
neutral_voltage = complex(bus["n"]["vr"], bus["n"]["vi"])
load_voltage = abs(phase_voltage - neutral_voltage)
@assert 220.0 < load_voltage < 240.0

(
    maximum_energy_residual_wh=maximum(abs.(diagnostics.energy_balance_wh)),
    maximum_power_bound_violation_w=maximum(diagnostics.power_limit_violation_w),
    maximum_mode_residual_w=maximum(abs.(diagnostics.simultaneous_power_w)),
    load_voltage_v=round(load_voltage; digits=3),
)
```

These checks validate this returned schedule under its declared model. The
mode residual is the absolute smaller charge/discharge power. Solver tolerances
can produce a small negative value for a nominally nonnegative power; the
power-bound check measures that violation separately, without clipping dispatch.
The
voltage check is an observation at the solved point; this example does not impose
a voltage operating band or prove feasibility under other loads or prices.
Read [Concepts](concepts.md) for the shared units and result contracts.

## Choose your next study

| Your question | Next guide |
|:--|:--|
| How do storage assumptions affect dispatch? | [Battery storage](components/battery_storage.md) and [multi-period OPF](problems/multiperiod.md) |
| Can an EV meet its departure target through a given outlet? | [Electric vehicles](ev/index.md) |
| How do source impedance, grounding and voltage regulation affect the result? | [Generators](components/generators.md) |
| Do converter topology, limits or local controls change the operating point? | [Inverter-based resources](ibr/index.md) |
| What states or model parameters are supported by evidence? | [State and model inference](estimation/studies.md) |
| What can a declared set of operating-envelope checks establish? | [DOE in 15 minutes](tutorials/doe_getting_started.md) |

For power-flow branches, static checks or neutral reduction, use the
[power-flow and operability overview](algorithms/index.md). For numerical
representations and solver tolerances, use
[formulation studies](formulations/index.md). The [home page](index.md) maps all
six study topics, and the [API index](api.md) locates their callable interfaces.
