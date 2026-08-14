# Closed-loop evidence diagnostics

The inverter-control evaluators are deliberately stateless: they map a local
voltage measurement and an operating request to a current command. For a
network study, that map can be embedded in a reduced feeder model, a
quasi-static re-solve, or a measured/discrete-time experiment. This page
provides two lightweight diagnostics for that outer loop:

1. `fixed_point_oracle` iterates a supplied map and reports convergence,
   residuals, and short repeated-state cycles.
2. `screen_fixed_point_gain` computes a finite-difference Jacobian and reports
   its spectral radius and induced infinity norm.

These are evidence-generation tools, not formal dynamic stability proofs. A
fixed point can be locally attractive while a different initial condition
fails, and a local gain screen says nothing about response time, delays,
sampling, ramps, protection state, or saturation outside the linearisation
point. Those effects must be represented in the supplied map when they are
material to the question.

## Inverter-specific campaign adapters

`inverter_control_fixed_point_oracle` wraps the exact three-phase controller
law around a supplied affine feeder sensitivity. It makes the reference
operating voltage explicit and iterates

```math
v_{k+1} = v_{ref} + Z\,[i_{exact}(v_k)-i_{exact}(v_{ref})].
```

This is the smallest useful equilibrium oracle for a campaign: it can expose
cycling, multiple initial-condition basins, and the distance from a smooth
equilibrium before a full feeder re-solve is commissioned. `Z` is still a
study input, not an inferred network property; use a full network callback
when voltage sensitivity varies materially over the operating range.

When a complete `Dict` network and `ControlledDevice{AdvancedInverter}` are
available, `solve_inverter_control_network_fixed_point` performs that higher-
fidelity callback internally. It solves the smooth controlled case first, then
rebuilds the physical plant at every exact controller iterate with the selected
converter- or grid-side current fixed. A non-publishable plant solve is
retained as an equilibrium failure; it is not converted into a false
convergence result.
The returned and iterated voltages use the controller's declared sensing
reference: phase-to-neutral when `neutral` is present, otherwise
phase-to-ground.

For a selected fleet, use `solve_controlled_inverter_fleet_network_fixed_point`.
All exact current targets are fixed in one re-solve, so voltage coupling between
controlled devices and unselected native IBRs is preserved. The deterministic
`controlled_inverter_network_fixed_point_rows` output is intended for campaign
tables and includes smooth-versus-exact voltage/current and P/Q differences,
iteration, cycle, and publishability fields.

Use `solve_controlled_inverter_fleet_multistart` when one initial condition is
not enough evidence. Starts are named and sorted deterministically, and the
smooth fleet solve is reused. `controlled_inverter_fleet_multistart_rows`
reports convergence, cycles, residuals, and final voltage spread relative to
the first converged start, plus the maximum pairwise spread across all
converged starts. Non-converged iterates are not labelled as equilibria and
receive `NaN` spread fields. A nonzero spread is evidence of initial-condition
dependence; it is not by itself a proof of multiple physical equilibria.

`inverter_control_network_voltage_sensitivity` estimates the feeder matrix
directly by perturbing the selected physical current target and rebuilding the
plant. It records plus/minus `SolveStatus` values for every column and leaves
failed columns as `NaN`; this makes current-limit boundaries visible in the
loop-gain evidence instead of silently extrapolating through them.
The result uses a 6×6 phase-rectangular matrix for four-leg devices and a 4×4
positive/negative-sequence matrix for three-leg converter-current targets,
because independent phase-current perturbations violate the three-leg
zero-neutral-current constraint.
The reduced screen intentionally projects sensed phase voltage onto `U₁/U₂`;
if a study's controller is materially sensitive to common-mode voltage, retain
that channel in a separate phase-ground study rather than treating the 4×4
screen as a complete six-coordinate model.

For a selected fleet, `controlled_inverter_fleet_network_voltage_sensitivity`
performs the same audit with all selected current targets fixed simultaneously.
Its `6N × 6N` four-leg or `4N × 4N` three-leg converter-current matrix retains
cross-device coupling, rather than assembling
independent one-device sensitivities. Device blocks and perturbation columns
are ordered by sorted device id, and failed columns remain `NaN` with their
plus/minus statuses retained. Feed the resulting matrix to
`controlled_inverter_fleet_loop_gain` together with the fleet's local voltage
measurements to obtain a block-diagonal controller Jacobian composed with the
physical network response.
`controlled_inverter_fleet_network_sensitivity_rows` converts the matrix and
its per-column solve provenance to deterministic long-form records for CSV,
Arrow, or `DataFrame` campaign outputs.

The `InverterControlScalingAudit` returned by
`inverter_control_scaling_audit` records the SI starts/scales for voltage,
sequence voltage, current, apparent power, and the priority-capacity
auxiliary, together with the rating-relative smoothing widths. Record it with
each campaign manifest. In particular, the capacity auxiliary is scaled to
the residual capability implied by `priority_headroom_fraction`, not to the
full apparent-power rating.

## Generic fixed-point oracle

The callable passed to `fixed_point_oracle` accepts and returns a real vector.
The iteration is

```math
x_{k+1} = (1-\alpha)x_k + \alpha F(x_k), \qquad 0 < \alpha \le 1.
```

Convergence is declared when the update residual is below an absolute-plus-
relative tolerance. The optional trajectory includes the initial state and
all iterates, making it suitable for small reproducibility plots or tabular
evidence. A repeated state within `cycle_window` is reported as a cycle rather
than being silently labelled a failed solve.

```julia
map(x) = feeder_resolve_or_measurement_map(x)
oracle = fixed_point_oracle(
    map, initial_voltage;
    relaxation=0.8,
    max_iterations=100,
    atol=1e-8,
    rtol=1e-8,
)
oracle.converged, oracle.residual_norm, oracle.iterations
```

The map may be a reduced sensitivity model or a complete small feeder
re-solve. The oracle does not assume that the vector is voltage, nor does it
assign physical time to an iteration.

## Local loop-gain screen

For a map `F`, `screen_fixed_point_gain(F, x₀)` evaluates

```math
J_F(x_0) = \frac{\partial F}{\partial x}(x_0).
```

The reported `spectral_radius` is the largest eigenvalue magnitude. A value
below `threshold` is labelled `local_contractive`; this is explicitly the
unrelaxed Picard screen and the `margin` is `threshold - spectral_radius`.
`maximum_real_eigenvalue` and `continuous_time_margin` apply the ideal
first-order response model, while `alpha_max` gives the largest update-to-
response-time ratio predicted by the relaxed discrete iteration. These are
local numerical screens, not certificates for a nonlinear or hybrid controller.

For a four-leg controller adapter, the state is six real rectangular
coordinates in declared a-b-c order:

```text
(Re Va, Im Va, Re Vb, Im Vb, Re Vc, Im Vc).
```

`inverter_control_current_jacobian` differentiates the exact evaluator. A
caller supplies a 6×6 real `voltage_sensitivity` matrix, `Z`, that maps a
phase-current perturbation to a phase-voltage perturbation. The adapter then
screens the local loop

```math
J_{loop} = Z\,\frac{\partial i}{\partial v}.
```

For a three-leg converter-current target, use the reduced four-coordinate
positive/negative-sequence basis `(Re U₁, Im U₁, Re U₂, Im U₂)` and the
corresponding 4×4 network sensitivity. Independent phase-current perturbations
are not physically realizable because the converter enforces zero neutral
current. Grid-current targets may have a different realizable subspace when
phase-to-ground shunts or an LCL are present; the sensitivity result records
the coordinate system selected for the target.

The units of `Z` must match the phasor convention used by the supplied
measurement and evaluator (normally V RMS per A RMS in rectangular
coordinates). Its sign and reference direction are part of the caller's
network convention; physical network sensitivities returned by
`inverter_control_network_voltage_sensitivity` and its fleet counterpart use
the network's positive-current injection convention and phase-to-neutral
voltages whenever the inverter declares a neutral.

```julia
controller = SequenceController(PositiveSequenceVoltVarWatt(...))
measurement = InverterControlMeasurement(phase_voltage)
request = InverterControlRequest(
    p_available=p_available, p_rated=p_rated, q_scale=q_scale)
ratings = InverterControlRatings(s_max=s_max, i_max=i_max)

J_i = inverter_control_current_jacobian(
    controller, measurement, request, ratings)
screen = inverter_control_loop_gain(
    controller, measurement, request, ratings, Z)
screen.spectral_radius, screen.maximum_real_eigenvalue,
screen.continuous_time_margin, screen.alpha_max
```

The adapter uses the exact piecewise controller law, so a finite-difference
point on a hard curve corner can be one-sided or numerically ambiguous. For
evidence, record the finite-difference step and perturbation point, and use
`evaluate_smooth` or a smooth feeder map when a differentiable surrogate is
the intended object of study. The smooth surrogate is a numerical modelling
choice; it does not change the firmware interpretation of the exact law.

## Recommended evidence record

For each operating point, retain:

- the controller, topology, current-target, and limiter configuration;
- the voltage phasors and request/ratings used for the evaluation;
- the supplied feeder sensitivity or outer-loop map provenance;
- `iterations`, `residual_norm`, and any detected `cycle_period`;
- finite-difference step, spectral radius, induced norm, threshold, and
`local_contractive` unrelaxed Picard result, `maximum_real_eigenvalue`,
`continuous_time_margin`, and `alpha_max` response-rate screen;
- whether the point is inside a curve corner, current limit, or apparent-power
  limit.

Across a campaign, compare these diagnostics with the physical quantities in
the network result: per-leg current, converter/grid terminal current, sequence
power, zero/negative-sequence voltage, and DC-link ripple. This keeps the
diagnostic layer reusable for general inverter studies while providing an
auditable path to standards-oriented evidence when the study assumptions are
fixed separately.

## API

```@docs
FixedPointIterationResult
FixedPointGainScreen
InverterControlScalingAudit
InverterControlFixedPointResult
InverterControlNetworkFixedPointResult
fixed_point_oracle
finite_difference_jacobian
screen_fixed_point_gain
inverter_control_scaling_audit
inverter_control_current_jacobian
inverter_control_fixed_point_oracle
solve_inverter_control_network_fixed_point
ControlledInverterFleetNetworkFixedPointResult
solve_controlled_inverter_fleet_network_fixed_point
controlled_inverter_network_fixed_point_rows
ControlledInverterFleetMultiStartResult
solve_controlled_inverter_fleet_multistart
controlled_inverter_fleet_multistart_rows
InverterControlNetworkSensitivityResult
inverter_control_network_voltage_sensitivity
ControlledInverterFleetNetworkSensitivityResult
controlled_inverter_fleet_network_voltage_sensitivity
controlled_inverter_fleet_network_sensitivity_rows
controlled_inverter_fleet_loop_gain
inverter_control_loop_gain
```
