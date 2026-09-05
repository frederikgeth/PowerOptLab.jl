# Function formulations

PowerOptLab separates **what a control law means** from **how its equations are
encoded**. Define a bounded continuous piecewise-linear (PWL) curve in physical
units, specify whether an occurrence tracks that curve or imposes a limit, and
choose a numerical representation and a solver for the resulting whole model.

The toolkit supports smooth approximations, exact segment graphs,
complementarity graphs, and convex-hull relaxations. It also provides physical
error bounds, purpose-specific magnitude approximations, configurable experiment
runners, and integration with the advanced inverter models.

## Choose a starting point

| Your task | Worked guide |
|:--|:--|
| Model volt-var or volt-watt without tying the data to a solver | [Control intent and compilation](compilation.md) |
| Exploit voltage bounds or distinguish tracking from a power limit | [Bounds, equalities and limits](bounds_and_relations.md) |
| Compare encodings, scales, starts and solver settings reproducibly | [Configurable experiments](experiments.md) |
| Select smoothing in a physical feeder/inverter model | [Physical inverter example](controllers.md) |
| Understand approximation error, derivatives and relaxation witnesses | [Mathematics and analytic reference](theory.md) |
| Translate an engineering accuracy requirement into a numerical budget | [Physical error and conditioning](error_budgets.md) |
| Interpret OpenDSS tolerances and choose a whole-model solver | [Tolerances and solver choice](tolerances_and_solvers.md) |
| Implement a different smooth approximation | [Smoothing-family extension](extensions.md) |

For a first encounter, follow compilation → bounds → experiments. Power-system
modelers can then use the physical inverter example; OR researchers can use the
analytic reference and solver guide to examine feasible-set and stationarity
questions. The examples execute during documentation builds. The
[API reference](api.md) lists constructors and options; the
[references](references.md) connect the mathematics to primary literature and
external software.

## Choose the modeling interface

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

## Read results with the right meaning

- **Canonical curve:** the prescribed PWL function, before approximation.
- **Surrogate graph:** equality to a smooth approximation; its error budget is a
  function-value bound, not a network equilibrium or DOE safety guarantee.
- **Exact graph or relation:** a mathematically equivalent encoding on the
  declared domain; finite solver tolerances still require candidate checks.
- **Outer relaxation:** a superset of the intended feasible set. A graph hull may
  admit points off the curve; some projected inequality relations are exact.
- **Candidate audit:** residuals at a returned point, separate from solver success,
  stationarity, global optimality and physical safety.

The library specializes scalar occurrences using declared bounds. It does not
infer general network bounds, automatically convert whole AC controllers to
MIP/MPCC/conic models, or supply global nonlinear error certificates. Discontinuous
switching, hysteresis and dynamic controller state require explicit models.
