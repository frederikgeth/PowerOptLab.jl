# PowerOptLab.jl

PowerOptLab provides models and numerical methods for studying **four-wire
distribution networks with uncertain states and parameters**.

## Start here

The [getting-started example](getting_started.md) builds a small feeder,
schedules a battery and checks the solution. It uses synthetic data and Ipopt.
The example is run automatically when the documentation is built.

For an existing checkout, install the pinned dependencies from its root:

```sh
julia --project=. scripts/instantiate_pinned.jl
```

Then choose a study topic below. For a complete worked experiment, follow the
[battery-model comparison](tutorials/battery_storage_models.md) or
[four-wire state-estimation study](tutorials/constrained_nlls_state_estimation.md).
Both include figures, expected results and executable checks.

## Study topics

Choose a chapter according to the quantities you want to estimate, optimise or
check.

| Study topic | Question | Start here |
|:--|:--|:--|
| State and parameter estimation | Estimate voltages, line lengths, taps or overhead-line configurations from measurements and network data. | [Estimation methods](estimation/studies.md) |
| Device models and scheduling | Select generator, battery and EV models, and calculate charging and discharging schedules. | [Device models](components/index.md) |
| Inverter models and controls | Study voltage control, converter losses, conductor currents and DC-link requirements. | [Inverter studies](ibr/index.md) |
| Dynamic operating envelopes | Calculate import or export limits and test them over specified operating conditions. | [DOE in 15 minutes](tutorials/doe_getting_started.md) |
| Power flow and voltage operability | Compare power-flow solutions, examine static voltage margins and assess neutral reduction. | [Power-flow methods](algorithms/index.md) |
| Numerical formulations and solvers | Compare mathematical representations, approximation errors and solver tolerances. | [Formulation studies](formulations/index.md) |

The [API reference](api.md) lists functions, arguments and returned quantities.
The [symbol index](symbol_index.md) provides an alphabetical lookup.
[Concepts](concepts.md#Capability-status) summarises the development status of
the models and methods.

## Research context

The research programme addresses three questions:

1. Which network states and parameter values are consistent with measurements
   and network records?
2. Which additional measurements or changes in operating point would distinguish
   the plausible models while respecting voltage and equipment limits?
3. Which operating decisions remain feasible under those models?

The package builds on the
[BMOPFTools](https://github.com/frederikgeth/BMOPFTools.jl) reference
current–voltage OPF formulation. It uses the same network equations, with an
explicit neutral conductor and per-unit scaling, and adds models through the
public model-building functions.

The current package provides separate estimation, optimisation and checking
methods. Joint network-model identification and automatic experiment design
remain research objectives. The nonlinear DOE solvers also lack a general
feasibility certificate for continuous uncertainty and utilisation sets.
[Concepts](concepts.md) explains the modelling assumptions and conventions; the
[research programme](research_program.md) sets out the longer-term questions.

## Extending the package

See [Contributing](contributing.md) for adding a model or solution method.
Models with sufficiently established interfaces and validation may be proposed
for inclusion in BMOPFTools.
