# Concepts

PowerOptLab studies four-wire distribution networks in which measurements,
parameters or operating conditions are uncertain. The network equations use RMS
voltage and current phasors, with an explicit neutral conductor. They are based
on the [BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl)
current–voltage formulation.

The intended research sequence is to estimate the network state and parameters,
identify additional measurements or operating-point changes that reduce the
uncertainty, and test operating decisions against the remaining possible models.
The present solvers address parts of that sequence. Joint topology and parameter
identification and automatic experiment design remain research objectives.
The [research programme](research_program.md) describes those objectives.

The nonlinear operating-envelope solvers test specified scenarios and utilisation
points. The static operability methods also provide sufficient fixed-point
conditions for supported network models. These establish different results:
finite scenario tests do not establish feasibility over a continuous uncertainty
set, and a fixed-point condition applies only to its specified model and region.

## Models, formulations and solution methods

A component model describes a device, a problem formulation specifies the
network calculation, and a solution method defines the numerical procedure.
These distinctions help identify what changes between two studies. The source
folders follow the same three categories; the documentation groups them by
[study topic](index.md#Study-topics).

| Category | Contribution | Programming interface | Source |
|---|---|---|---|
| **Component model** | a new network *element* (device physics) | `model_hook!` / `solution_hook!` | `src/components/` |
| **Problem specification** | a new *objective + variable/constraint structure* over the whole network | staged `build_opf_model` / `optimize_opf!` / `generation_cost` / `extract_result` | `src/problems/` |
| **Solution method** | a numerical procedure for solving the equations | custom solve loops around the staged API | `src/algorithms/` |

- A **component model** adds the electrical equations and ratings of a device.
  Its terminal currents enter network KCL; storage models also introduce a state
  equation between time intervals.
- A **problem formulation** defines an objective and constraints over the network
  model. Examples include OPF, state estimation and multi-period scheduling.
- A **solution method** determines how those equations are solved. HELM, for
  example, constructs a power series from the network admittance matrix.

## Study assumptions and development status

The metadata at the start of a model guide describes its development status and
intended use:

- **Maturity** — *proof of concept* (demonstrates a narrow idea), *prototype*
  (test-backed but experimental), *experimental* (an evolving model/interface
  with explicitly documented validation limits), *research prototype* (intended for declared
  studies rather than operational use), or *promotion candidate* (stable enough
  to propose for the BMOPF spec).
- **Direction** — *forward* (given parameters, find an operating point:
  dispatch/OPF), *inverse* (given measurements, find the state or parameters),
  or *hierarchical* (different decision-makers or response levels).
- **Temporal structure** — *single-snapshot* vs *inter-temporal* (state-of-charge
  coupling, shared parameters, or per-interval envelopes spanning many snapshots).

- **Uncertainty/security treatment** — deterministic, explicit scenarios,
  tested utilization points, or globally certified. These labels describe
  different levels of assurance. State which one was used when reporting results.

## Units and solution status

Inputs and results are expressed in physical units, including W, var, Wh, V and A.
The OPF-based
models use the BMOPFTools bases obtained with `opf_bases(ctx)` for per-unit
scaling. Solvers that support both `per_unit=true` and `per_unit=false` retain
the same physical units at their interfaces. Numerical scaling, convergence and
the local solution reached can differ. Coupled phase-aware inverter-control
modelling currently requires per-unit coordinates; its parameters and results
remain SI-valued. Check each solver's supported formulation.

The solvers use the following common definitions:

- [`AbstractDevice`](@ref) describes how a device is checked, added to the
  network equations, linked across intervals and reported in the results;
- [`TimeGrid`](@ref), which gives every snapshot an explicit positive duration;
- [`AbstractMeasurement`](@ref) describes a measurement value and its uncertainty;
- [`build_multi_context`](@ref), which builds ordered snapshots into one JuMP
  model, plus [`solve_status`](@ref) and [`solve_diagnostics`](@ref) for
  convergence information and numerical checks.

Results from the JuMP solvers retain a [`SolveStatus`](@ref), including the
termination and primal feasibility statuses. Check
`solve_status(result).publishable` before reading the reported solution. This
flag applies the solver-convergence criteria used by the package; physical
residuals and operating limits still need to be checked for the study.

When composing an external OPF, `solve_status(model)` applies the same contract
to the solved JuMP model. Pass that status to `extract_device(device, handles,
status)` so an unsuccessful solve cannot be presented as a valid dispatch.

See [Contributing](contributing.md) for adding models and methods, documenting
their validation and proposing established work for inclusion in BMOPFTools.

## Capability status

The following tables summarise the development status of the component models
and network calculations. Each linked guide describes its supported equations
and validation.

### Component models

| Capability | Entry point | Maturity |
|---|---|---|
| [Generalized generators](components/generators.md) | [`GeneralizedGenerator`](@ref), [`SourceGenerator`](@ref) | research prototype |
| [Battery storage](components/battery_storage.md) with state of charge | [`StorageDevice`](@ref) | promotion candidate |
| [EV charging](ev/index.md) (V1G / V2G) with availability & departure energy | [`EVDevice`](@ref), [`ChargingSession`](@ref) | promotion candidate |
| [Advanced inverter](ibr/index.md) (circuit-aware IBR) | [`AdvancedInverter`](@ref) | experimental |
| [Phase-aware inverter controls](components/inverter_controls.md) (local three-leg laws) | [`SequenceController`](@ref) | research prototype |
| [IVQ battery](components/ivq_battery.md) (current–voltage–charge model) | [`IVQBattery`](@ref) | prototype |

### Problem formulations

| Capability | Entry point | Direction | Maturity |
|---|---|---|---|
| [Multi-period OPF](problems/multiperiod.md) co-optimising many snapshots | [`solve_multiperiod_opf`](@ref) | forward | promotion candidate |
| [Legacy WLS state estimation](problems/state_estimation.md) | [`solve_state_estimation`](@ref) | inverse | prototype |
| [Constrained NLLS state estimation](problems/constrained_state_estimation.md) | [`solve_sparse_state_estimator`](@ref) | inverse | prototype |
| [Post-OPF operability](problems/operability.md) (native static snapshot diagnostics and sufficient certificates) | [`check_opf_operability`](@ref) | forward | research prototype |
| [Parameter estimation](problems/parameter_estimation.md) (line lengths / taps) | [`solve_parameter_estimation`](@ref) | inverse | prototype |
| [Inverse Carson reconstruction](problems/inverse_carson.md) (compatible overhead constructions) | [`solve_inverse_carson`](@ref) | inverse | prototype |
| [Dynamic operating envelopes](tutorials/doe_getting_started.md) (active import/export capacity) | [`solve_operating_envelope`](@ref) | forward | research prototype |
| [Bilevel PV/tap POC](problems/bilevel.md) (DiffOpt lower-level response) | [`solve_bilevel_pv_tap`](@ref) | hierarchical | proof of concept |
