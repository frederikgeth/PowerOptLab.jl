# A mixed-controller Australian feeder

<!-- study-executable -->

One network solve can contain mean-voltage, worst-phase, positive-sequence and
per-phase four-leg VV+VW controllers. Each site owns its plant, controller,
available-power request and operating mode. Unselected native devices retain
their network equations. Run this example from a repository checkout; it uses
the prepared Australian feeder and its fixed hardware rather than resizing
devices to make the comparison pass.

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

local_point = engineering_point(local_result, condition)
local_export = sum(row.p_poc_W for row in voltage_aware_inverter_rows(local_result))
coordinated = solve_engineering_coordinated(network, mixed;
    local_point, minimum_p_W=local_export)
@test coordinated.model_candidate_accepted
@test coordinated.temporary_neutral_constraints_remaining == 0
@test last(coordinated.attempts).stage === :full_capability
@test last(coordinated.attempts).weight == 0.
```

The coordinated solve releases the VV/VW commands and maximizes total active
export with free Q. Hardware, availability, physical current capabilities and
network limits stay fixed. This compares two service choices; it cannot show
that coordination improves export while delivering the same reactive service.
The final solve removes every temporary neutral row and uses no regularizer.
Only raw solver success plus physical/network/replay acceptance can publish it.
Accepted-stage primal/dual warm starts and exception-safe row cleanup live in
package solver infrastructure; the neutral-release ordering and dispatch policy
remain choices of this case study.

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
