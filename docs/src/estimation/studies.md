# State and parameter estimation

The appropriate estimation method depends on which quantities are unknown and
which measurements are available. PowerOptLab provides state estimation,
calibration of selected line lengths and taps, and overhead-line reconstruction.

| Estimation problem | Required data and assumptions | Guide |
|:--|:--|:--|
| Estimate bus voltages from noisy measurements. | Measurements and their uncertainties; known topology and parameters | [State estimation](index.md) |
| Estimate line lengths or transformer taps using several operating points. | Measurements at each operating point; fixed topology, conductor connections and per-length R/X | [Parameter calibration](../problems/parameter_estimation.md) |
| Identify overhead-line configurations consistent with sequence impedances. | Sequence impedances and candidate conductor geometries and materials | [Inverse Carson reconstruction](../problems/inverse_carson.md) |

## Selecting a method

For state estimation, begin with the [formulation comparison](comparison.md),
then the [modelling tutorial](../tutorials/constrained_nlls_state_estimation.md).
Use [observability](observability.md) and
[likelihood, loss and priors](maximum_likelihood.md) to interpret the result.
The [state estimation API](../api/state_estimation.md) lists measurement types,
solvers and returned diagnostics.

The [calibration guide](../problems/parameter_estimation.md) contains a complete
example in which line lengths and taps are common to several snapshots. The
[calibration API](../api/parameter_estimation.md) describes the input parameters
and estimated quantities.

The [line reconstruction tutorial](../tutorials/inverse_carson_reconstruction.md)
tests candidate conductor arrangements and constructs the corresponding
four-wire line models. The [line reconstruction API](../api/inverse_carson.md)
describes the inputs, fitted models and parameter profiles.

## Interpreting the estimates

A state estimate is conditional on the assumed network parameters. A
calibration result is a local fit to the specified measurements. Line
reconstruction can return several configurations that fit the same sequence
impedances. When using these results in an
[operating-envelope study](../tutorials/doe_scenario_design.md), carry the
remaining uncertainty into the network scenarios. Joint estimation across these
problems is a research objective described in the
[research programme](../research_program.md).
