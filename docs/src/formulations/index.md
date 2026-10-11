# Function formulations

Volt-var and volt-watt characteristics can be represented by smooth
approximations, exact piecewise-linear constraints, complementarity conditions
or convex-hull relaxations. The choice affects the optimisation problem, its
numerical conditioning and the error relative to the specified characteristic.

This chapter compares those representations for bounded continuous
piecewise-linear (PWL) curves in physical units. It also covers approximation
error, solver tolerances and their use in the inverter models.

## Worked comparisons

| Your task | Worked guide |
|:--|:--|
| Model volt-var or volt-watt without tying the data to a solver | [Volt-var and volt-watt formulations](compilation.md) |
| Exploit voltage bounds or distinguish tracking from a power limit | [Voltage bounds and constraints](bounds_and_relations.md) |
| Compare encodings, scales, starts and solver settings reproducibly | [Configurable experiments](experiments.md) |
| Select smoothing in a physical feeder/inverter model | [Physical inverter example](controllers.md) |
| Understand approximation error, derivatives and relaxation witnesses | [Mathematics and analytic reference](theory.md) |
| Translate an engineering accuracy requirement into a numerical budget | [Physical error and conditioning](error_budgets.md) |
| Interpret OpenDSS tolerances and choose a whole-model solver | [Tolerances and solver choice](tolerances_and_solvers.md) |
| Implement a different smooth approximation | [Smoothing-family extension](extensions.md) |

Start with the control-characteristic example, the effect of voltage bounds,
and the numerical experiments. The feeder example shows how formulation choices
affect an inverter study. The analytical reference examines approximation error
and stationary solutions. The examples are run when the documentation is built.
The [API reference](api.md) gives the function arguments, and the
[references](references.md) identify the relevant literature and software.

## Choose the modelling interface

| Interface | Modeler's responsibility | Returned object |
|:--|:--|:--|
| `VoltVarWattIntent` and `VoltVarWattEncoding` | Separate physical curves/sensing from per-curve encoding choices | Reusable intent and encoding data |
| `formulate_control_curve!` / `formulate_pwl!` | Request a graph output and couple it to the rest of the model | `PWLFormulationHandle` |
| `formulate_control_relation!` / `formulate_pwl_relation!` | Declare `:equal`, `:upper` or `:lower` for an existing output, with valid input bounds | `PWLRelationHandle` with the selected lowering |
| `lower_positive_policy` | Select smooth curves within the supported physical controller equations | Controller policy for the advanced IBR model |
| `FormulationCase` | Build a fresh model and declare observations/custom metrics | Reproducible experiment definition |

A curve is reusable; its lowering depends on the occurrence. A volt-watt upper
limit over a concave portion may reduce exactly to linear inequalities even when
tracking the same curve needs a nonlinear or discrete representation. Scalar
convexity does not establish convexity of the network model.

## Interpreting the formulations

- **Specified characteristic:** the prescribed PWL function, before approximation.
- **Smoothed characteristic:** equality to a smooth approximation; its error budget is a
  function-value bound, not a network equilibrium or DOE safety guarantee.
- **Exact constraints:** a mathematically equivalent formulation on the
  declared domain; finite solver tolerances still require candidate checks.
- **Outer relaxation:** a superset of the intended feasible set. A graph hull may
  admit points off the curve; some projected inequality relations are exact.
- **Solution checks:** residuals at a returned point, separate from solver success,
  stationarity, global optimality and physical safety.

The library specializes scalar occurrences using declared bounds. It does not
infer general network bounds, automatically convert whole AC controllers to
MIP/MPCC/conic models, or supply global nonlinear error certificates. Discontinuous
switching, hysteresis and dynamic controller state require explicit models.
