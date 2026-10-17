# Phase-aware control API

The [control guide](../components/inverter_controls.md) shows how to specify
and solve the local control laws. Their equations are derived in the
[phase-aware control chapter](../ibr/phase_aware_control_laws.md), and the
[validation guide](../ibr/verification.md) describes the reference calculations
and numerical tests. [Network studies](../ibr/network_control_studies.md) cover
multiple inverters and comparisons of equipment ratings.

```@docs
AbstractInverterControlLaw
AbstractPositiveSequencePolicy
AbstractUnbalancePolicy
AbstractLimiterPolicy
AbstractCurrentTarget
ConverterCurrentTarget
GridCurrentTarget
PiecewiseLinearLaw
WorstPhaseVoltVarWatt
AverageVoltageVoltVarWatt
PositiveSequenceVoltVarWatt
NoUnbalanceControl
NegativeSequenceAdmittanceDroop
negative_sequence_virtual_impedance
CommonScaleLimiter
SequencePriorityLimiter
SequenceController
InverterControlMeasurement
InverterControlRequest
InverterControlRatings
InverterControlResult
ConverterTerminalResult
ControlledDevice
ControlledInverterResult
evaluate_exact
dc_power_offset
sensing_terminal
evaluate_smooth
stamp_smooth_control!
solve_controlled_inverter
inverter_spec
inverter_handles
```
