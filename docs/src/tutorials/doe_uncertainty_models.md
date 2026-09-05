# Constructing and comparing DOE uncertainty models

> **Question:** which uncertainty representation matches the physical evidence,
> and how can alternatives be compared without changing the DOE claim?
>
> **Prerequisite:** [DOE scenario design and held-out
> evaluation](doe_scenario_design.md)
>
> **Scope:** uncertainty construction and sensitivity analysis; no uncertainty
> model is validated automatically

The uncertainty source, mathematical model, sampled values, and materialized
network are different objects:

```text
measurements or engineering records
        → uncertainty model
        → reproducible samples or candidates
        → materialized network scenarios
        → allocation and held-out evaluation
```

Collapsing those steps hides assumptions. For example, a Gaussian covariance
does not explain physical bounds; a list of impedance candidates does not imply
probabilities; and a bootstrap does not establish stationarity.

## 1. Choose the representation from the evidence

| Available evidence | Useful representation | Principal question to justify |
|---|---|---|
| Forecast or DSSE residual history | Row or moving-block residual bootstrap | Are residuals transferable, ordered, and approximately stationary over the study period? |
| Estimated mean and covariance | Gaussian or box-conditioned Gaussian | Is the distributional form defensible, and how are nonphysical values excluded? |
| Engineering tolerances | Box, norm ball, or designed stress grid | Do the bounds represent credible simultaneous errors or independent worst cases? |
| Alternative conductor/topology/phase records | Categorical model candidates | Which candidates remain plausible, and are weights evidential or merely priorities? |
| Rare operational concerns | Named stress scenarios | What question does each stress case test? |

Do not call all five rows “Monte Carlo uncertainty.” Their semantics and
permitted claims differ.

## 2. Draw and materialize explicitly

```julia
using PowerOptLab

include("scripts/cases/doe_uncertainty_tutorial.jl")

samples = sample_doe_box_truncated_gaussian_uncertainty(
    [5000.0, 5000.0],
    [500.0^2 0.4 * 500.0^2;
     0.4 * 500.0^2 500.0^2];
    parameter_names=("p1_W", "p2_W"),
    lower=[0.0, 0.0],
    upper=[Inf, Inf],
    count=100,
    seed=20260905,
    draw_budget=1000,
    metadata=Dict("model" => "synthetic teaching assumption"))

function apply_load_sample!(network, sample)
    network["load"]["d1"]["p_nom"][1] = sample.parameters["p1_W"]
    network["load"]["d2"]["p_nom"][1] = sample.parameters["p2_W"]
    return nothing
end

load_scenarios = materialize_doe_scenarios(
    doe_uncertainty_feeder(5000.0, 5000.0),
    samples,
    apply_load_sample!;
    dataset_id="tutorial-load-model-v1",
    materializer_id="apply-load-sample-v1",
    role=:calibration,
    source="synthetic box-conditioned Gaussian")
```

The draw is reproducible and its support rule is explicit. This does not
validate Gaussianity, covariance calibration, the bounds, or the network
transformation. Those remain scientific responsibilities.

Use [`sample_doe_empirical_residual_bootstrap`](@ref) when actual residual rows
are available. Resample all correlated quantities together. Moving blocks can
retain short-range serial structure, but their block length, seasonality, and
transferability still need evidence.

## 3. Compare models with a matched design

The comparison target is not “which model gives the smallest envelope?” A more
conservative answer can arise from a wider uncertainty set, a different norm,
more dimensions, different controls, or a genuinely more influential physical
quantity.

Hold these fixed:

- feeder, base operating points, participant set, and units;
- allocation objective and normalization;
- control information structure;
- utilization search and replay procedure;
- calibration/test separation;
- evaluation scenarios and pairing identities; and
- the physical interpretation of uncertainty magnitude or coverage.

Then evaluate one issued envelope or one consistently re-issued method under
each labelled construction:

```julia
comparison = compare_doe_uncertainty_models(
    ["model_a" => scenarios_a,
     "model_b" => scenarios_b],
    cps,
    issued;
    scales=(0.5, 0.75, 1.0),
    roles=:test,
    utilizations=:corners,
    control_policy=IssuePlusLocalLaws(),
    pairing=:auto)

comparison.rows
comparison.pairwise_rows
comparison.diagnostics["pairings"]
```

Pairing reveals which same-case outcomes change classification. It does not
prove that the uncertainty model caused the change or that either model is
correct.

## 4. Interpret the impedance-versus-load result carefully

[Liu, Braslavsky, and Mahdavi
(2023)](https://doi.org/10.35833/MPCE.2023.000653) report a much larger capacity
reduction for one impedance-uncertainty construction than for one load-
uncertainty construction on their Australian network. This is a strong
motivation for parameter-uncertainty research, but the norms, dimensions,
radii, and physical quantities are not commensurate.

A scientifically useful extension should derive both perturbation models from
measurement or estimation evidence and compare them on common held-out network
realizations. Until then, say “larger in the reported case,” not “impedance
uncertainty universally dominates load uncertainty.”

## Research exercises

1. Compare clipped Gaussian samples with correctly box-conditioned samples.
   Identify which moments and boundary masses change.
2. Repeat a residual bootstrap with block lengths 1, 6, and 48. State what
   temporal assumption each choice makes.
3. Treat two line models first as equally weighted candidates and then as
   unweighted stress cases. Explain why the network scenarios are unchanged but
   the interpretation is different.

## What this does not prove

- Reproducible sampling is not uncertainty-model validation.
- Physical bounds do not create calibrated probabilities.
- A norm radius in one parameter space is not directly comparable with the
  same percentage in another.
- Paired finite outcomes do not establish a causal model effect.
- A robust feasible region for exogenous uncertainty is not automatically a
  robust participant-utilization box.

Continue with [Statistical validation of DOE
claims](doe_statistical_validation.md). Full covariance, bootstrap, materializer,
and model-sensitivity examples remain in the [advanced uncertainty
laboratory](doe_uncertainty_coverage.md).
