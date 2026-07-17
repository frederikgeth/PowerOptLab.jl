# HELM versus nonlinear power flow

> **Audience:** power-system researchers · **Scope:** four-wire distribution power flow, solvability diagnostics, and cross-method validation.

The Holomorphic Embedding Load-flow Method (HELM) is not merely another nonlinear iteration. It constructs the branch connected to an energized no-load network and evaluates it by analytic continuation. This gives it useful, independent numerical evidence beside the Ipopt-based nonlinear `solve_pf`, but a finite-order divergent series is not a proof of voltage collapse.

The trade-off is explicit: the present embedding supports constant-power and constant-impedance load components only. It cannot replace the nonlinear engine for arbitrary ZIP/exponential loads, IBRs, generators, controls, or OPF limits.

## 1. Two solvers, different evidence

| Question | HELM `solve_pf_helm` | Nonlinear `solve_pf` |
|---|---|---|
| Operating point | no-load germ, power series, Padé continuation | iterative nonlinear feasibility/optimization solve |
| Initialization | none | numerical initialization and local convergence matter |
| Branch | continuously connected to no-load state | feasible local point reached by the solve |
| Failed solve | exposes Padé/coefficient diagnostics | exposes nonlinear solver status |
| Collapse study | `singularity_estimate` is a heuristic | loading continuation is separate |
| v1 load model | constant-P and constant-Z | full engine-supported operational model |

Agreement on common physics is strong cross-method evidence because the numerical mechanisms differ. Disagreement should trigger a model/branch audit, not a vote for a preferred solver.

## 2. HELM's non-obvious construction

HELM introduces a complex loading parameter ``s``. Sources stay fixed; supported loads scale from the energized no-load germ at ``s=0`` to the study case at ``s=1``.

```math
s=0:\ \text{energized no-load germ},\qquad s=1:\ \text{requested operating point}.
```

The germ is not a flat-voltage approximation. In an unbalanced four-wire feeder it includes the actual source boundary and no-load neutral equilibrium. Each voltage and ideal-coupling current is a holomorphic series in ``s``. One LU factorization of the augmented nodal system serves every coefficient order; Wynn-epsilon Padé continuation evaluates the series at ``s=1``.

```julia
using PowerOptLab

helm = helm_series(net; max_order=40, tol=1e-8)
helm.status       # :converged | :series_diverged | :max_order_reached
helm.residual     # full nonlinear current mismatch [A]
helm.pade_spread  # last-two Padé differences for every coefficient row
helm.coefficient_tail_ratios
helm.singularity_estimate # heuristic; validate before calling it a margin
```

The final test is physical: returned voltages are substituted into the full nonlinear current mismatch. `solve_pf_helm` returns the standard result shape.

### Why this is more than avoiding an initial guess

Iterative solvers can fail because of scaling, a local basin, line search, or a poor start. HELM fixes the studied branch by continuation from no load and therefore supplies different evidence. Persistent finite-order coefficient growth at ``s=1`` says this series evaluation did not resolve the requested state; another study is required before concluding that the solvable branch ends there.

## 3. Interpret status and margin precisely

| Status | Meaning | Response |
|---|---|---|
| `HELM_CONVERGED` | operational branch reached; physical residual met tolerance | inspect voltage/diagnostics and compare with nonlinear PF |
| `HELM_SERIES_DIVERGED` | mismatch failed and the exposed coefficient tail grows | inspect Padé spread/ratios and run a loading continuation or analytic check |
| `HELM_MAX_ORDER` | mismatch failed without the growing-tail classification | increase `max_order`; the result is inconclusive |

Coefficient-tail Domb-Sykes extrapolation estimates a coefficient-dominating singularity ``\lambda``. On the repository's analytic two-node saddle-node, a value near 2 means the present load is about half the known collapse loading, and a value below 1 agrees with the known branch ending before ``s=1``. General networks may have other real or complex singularities, so this interpretation does not transfer automatically. `NaN` means the finite tail could not produce an estimate.

### Pitfall: calling every failed solve voltage collapse

Both `HELM_SERIES_DIVERGED` and `HELM_MAX_ORDER` are inconclusive about non-existence. Likewise, failure of nonlinear `solve_pf` is not proof. For a defensible collapse claim, bracket the limit with continuation power flow or use an analytically known fixture, and report the common model and tolerances.

## 4. Four-wire and ideal-element strengths

HELM uses an augmented four-wire nodal system. Floating neutrals remain nodes, so neutral displacement is solved rather than silently grounded. Zero-impedance switches and ideal transformers are bordered constraints, not arbitrary huge admittances:

```math
\begin{bmatrix}Y&A^T\\A&0\end{bmatrix}
\begin{bmatrix}V\\w\end{bmatrix}=
\begin{bmatrix}I\\0\end{bmatrix}.
```

Their voltage identities hold term-by-term. With `switches=:constrain`, HELM also returns physical switch-conductor currents; `:alias` gives identical voltages via node fusion.

```julia
hr = helm_series(net; switches=:constrain, ideal_xfmrs=:constrain)
hr.couplings
```

### Pitfall: using a huge admittance for an ideal element

An arbitrary `1/z` changes conditioning and adds a fictitious voltage drop. Bordered constraints retain exact identities. Note that a transformer promoted to an ideal coupling is purely ideal: its magnetizing shunt and neutral-ground branch are absent. That modelling choice must match the comparison case.

## 5. HELM as an independent oracle

For a case inside both solvers' scope, compare conductor voltages and currents:

```julia
using BMOPFTools: solve_pf

h = solve_pf_helm(net)
n = solve_pf(net; per_unit=false)
@assert h["termination_status"] == "HELM_CONVERGED"
Δv = h["bus"]["loadbus"]["1"]["vm"] - n["bus"]["loadbus"]["1"]["vm"]
```

Repository validation combines closed-form feeders, a separate `ybus_linearized` current-mismatch oracle, OpenDSS decks, and three-way HELM-Ipopt-OpenDSS comparisons. A good research workflow is:

1. Specify common source, grounding, switch/transformer, and load physics.
2. Compare phase-to-earth voltages, neutral rise, and relevant currents.
3. Inspect residuals, Padé spreads, coefficient-tail ratios, and the singularity estimate.
4. Investigate disagreement before changing tolerances.
5. Add operational physics only to nonlinear PF, clearly marking HELM's scope.

### Pitfall: comparing different physics

HELM cannot validate nonlinear PF when the latter includes constant-current ZIP fractions, IBR controls, or operational limits absent from HELM. Reduce to the common supported model first; otherwise the discrepancy is expected model mismatch, not solver evidence.

## 6. The holomorphic load restriction is fundamental

The v1 embedding accepts constant-P and constant-Z components, including ZIP P/Z fractions. It rejects constant-current ZIP fractions and non-integer exponential loads because they involve ``|\Delta V|``, which is not holomorphic:

```julia
solve_pf_helm(net_with_constant_current_zip) # throws an informative ArgumentError
```

It also needs an energized WYE/SINGLE_PHASE source and a linear path from every node to a reference at the no-load germ. Generators and IBRs are not injected in v1. Use nonlinear OPF PF or the `ybus_linearized` fixed-point map when those models are required.

### Pitfall: treating rejection as a numerical weakness

The rejection protects the embedding's mathematical scope. Forcing a non-holomorphic load into the recursion makes the resulting diagnostics hard to interpret. An explicit scope is more scientifically useful than silently solving a different model.

## 7. Publication checklist

Report the load representation, source/reference and neutral treatment, ideal-element model, HELM order/tolerance/status/residual, Padé/coefficient diagnostics, singularity estimate, and cross-solver comparison on common physics. State every excluded control and unsupported load component. HELM's value is not replacing nonlinear PF: it provides a deterministic analytic-continuation reference whose claims can be checked independently.

For API details see [HELM power flow](../algorithms/helm.md).
