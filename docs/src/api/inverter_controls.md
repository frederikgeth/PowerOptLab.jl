# Phase-aware control API

Use the [control modeling guide](../components/inverter_controls.md) for
construction and executable examples, the
[control-law derivation](../ibr/phase_aware_control_laws.md) for equations, and
[verification](../ibr/verification.md) for the scope of the evidence.
[Fleet and hardware studies](../ibr/network_control_studies.md) document the
network-level interfaces.

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
CommonScaleLimiter
SequenceController
InverterControlMeasurement
InverterControlRequest
InverterControlRatings
InverterControlResult
ConverterTerminalResult
ControlledDevice
ControlledInverterResult
evaluate_exact
evaluate_smooth
stamp_smooth_control!
solve_controlled_inverter
inverter_spec
inverter_handles
```
