# A mixed-controller Australian feeder

<!-- study-executable -->

One network solve can contain mean-voltage, worst-phase, positive-sequence and
per-phase four-leg VV+VW controllers. Each site owns its plant, controller,
available-power request and operating mode. Unselected native devices retain
their network equations. Run this example from a repository checkout; it uses
the prepared Australian feeder and the declared fixture, which passes the
rated-capability gate (`scripts/australian_pv/rated_capability.jl`).

```julia
using PowerOptLab, Test
include(abspath("scripts", "australian_pv", "coordinated_engineering.jl"))

condition = first(engineering_conditions())
base = parse_bmopf(joinpath(ROOT, "networks", "LV3_55bus.bmopf.json"))
network = voltage_aware_network(base, condition)
mixed = engineering_fleet(network, condition)

@test Set(values(mixed.labels)) == Set(MIXED_CONTROLLERS)
@test length(mixed.spec.devices) == 12
@test Set(r.p_available for r in values(mixed.spec.requests)) == Set((1250., 3000., 5000.))
```

The study factory constructs the public composition
`ControlledInverterFleetSpec(Dict(id => ControlledDevice(plant, controller)),
Dict(id => request))`. The three-leg laws use `SequenceController` with
`AverageVoltageVoltVarWatt`, `WorstPhaseVoltVarWatt`, or
`PositiveSequenceVoltVarWatt`; the four-leg sites use `PerPhaseVoltVarWatt`.
The mean/worst/sequence distinction is voltage sensing. These three-leg sites
retain balanced converter current and request no negative-sequence service.
Four-leg sites retain equal phase active powers and may redistribute Q.

```julia
local_result = solve_voltage_aware_inverter_fleet(network, mixed.spec;
    operating_modes=mixed.modes, s_base=VA_BASE, solver_options=ENGINEERING_OPTIONS)
@test solve_status(local_result).publishable
@test all(r.passed for r in values(physical_feasibility(local_result)))
@test mixed_thermal_replay(network, local_result).accepted

# Local controllers do not enforce the 253 V planning limit.
local_point = engineering_point(local_result, condition; voltage_max=nothing)
local_export = sum(row.p_poc_W for row in voltage_aware_inverter_rows(local_result))

# Same reactive service: each site's exact local volt-var request and volt-watt
# ceiling, evaluated at the OPF's own voltages. Without the 253 V row this must
# reproduce the local fixed point.
same = solve_engineering_equal_service(network, mixed; local_point, voltage_max=nothing)
@test same.service_converged
agreement = engineering_agreement(local_point,
    engineering_point(same, condition; coordinated=true), mixed)
@test agreement["accepted"]

# Volt-var released: Q free, 253 V planning limit enforced.
released = solve_engineering_coordinated(network, mixed;
    local_point, minimum_p_W=local_export)
@test released.model_candidate_accepted
@test released.temporary_neutral_constraints_remaining == 0
```

The agreement check is the main test that simulation and optimisation describe
the same physics: two separately built formulations, one equilibrium. The
released OPF answers a different question (what Q freedom buys) and needs
network-operator agreement under AS/NZS 4777.2. With the declared fixture, all
sites in this nominal condition are limited by PV availability, so the released
OPF gains only converter and filter losses.

| Path | Supported scope | Explicit boundary |
|---|---|---|
| Implicit controlled PF | Mixed fixed-graph controller laws, fundamental three/four-leg models | Nonlinear DC energy and per-phase staged policy use the voltage-aware path |
| Voltage-aware local PF | Mixed mean/worst/sequence and per-phase four-leg policies, energized/disconnected sites, reactive-priority stages | Quasi-static fixed-point equilibrium; failures issue no replacement command |
| Constrained OPF plant | Three/four-leg and split-link fundamental models; user-defined objectives/service constraints | No built-in global optimality or firmware/dynamic guarantee |
| Embedded coupled PWM thermal allocation | Centered monolithic three/four-leg L filters, proportional passive R/L, open carrier-frequency source | LCL, split-link PWM and finite carrier-source dynamics are rejected |
| Nonlinear DC energy | Monolithic three/four-leg links, constant input power and fitted mean loss | Finite small-signal source and split-bank dynamics require separate models |
| Current/peak and combined voltage envelopes | Declared physical limits under the supported periodic carrier model | Expert RMS allowances do not imply a peak waveform |
| Switched EMT replay | Parameter-identical fixed-POC averaged/switched snapshots | Does not validate feeder/controller dynamics or all-time transients |

All-duty thermal bounds can reject a point that a tighter waveform-specific
model could accept. Retain that distinction when comparing controllers. The
current campaign, historical negative controls and restoration commands are
listed in the Australian study evidence index.
