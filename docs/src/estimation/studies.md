# State and model inference

Choose an inference workflow from what is uncertain and the evidence available.
State estimation, parameter calibration and line reconstruction have different
inputs and identify different quantities.

| Your question | Evidence and fixed assumptions | Start here |
|:--|:--|:--|
| What network state is supported by noisy measurements? | Telemetry with declared uncertainties; known topology and parameters | [State estimation](index.md) |
| Which shared line lengths or transformer taps fit multiple operating points? | Snapshot measurements; fixed topology, nominal per-length R/X and terminal maps | [Parameter calibration](../problems/parameter_estimation.md) |
| Which supplied overhead constructions agree with sequence line data? | Sequence impedance observations and declared geometry/material candidates | [Inverse Carson reconstruction](../problems/inverse_carson.md) |

## Follow a study path

For state estimation, begin with the [formulation comparison](comparison.md),
then the [modeling tutorial](../tutorials/constrained_nlls_state_estimation.md).
Use [observability](observability.md) and
[likelihood, loss and priors](maximum_likelihood.md) to interpret the result.
The [state estimation API](../api/state_estimation.md) lists the constructors,
compiled-model interfaces, solvers and diagnostics.

The [calibration guide](../problems/parameter_estimation.md) contains a complete
shared-parameter example. The [calibration API](../api/parameter_estimation.md)
defines its selected lengths/taps and result fields.

The [line reconstruction tutorial](../tutorials/inverse_carson_reconstruction.md)
shows candidate screening and model materialization; the
[line reconstruction API](../api/inverse_carson.md) lists the observation,
candidate and profile interfaces.

## Keep the inference boundary visible

These are separate workflows. A fitted state, a local parameter optimum and a
compatible construction are different forms of evidence. Their models and
uncertainty assumptions must remain explicit when preparing a downstream
[operating-envelope study](../tutorials/doe_scenario_design.md).
The [research program](../research_program.md) describes the proposed integrated
model-forensics workflow; it is not an existing joint inference API.
