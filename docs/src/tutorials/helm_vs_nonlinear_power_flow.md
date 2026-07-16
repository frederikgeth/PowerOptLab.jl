# HELM versus nonlinear power flow

> **Audience:** power-system researchers · **Scope:** four-wire distribution power flow, solvability diagnostics, and cross-method validation.

The Holomorphic Embedding Load-flow Method (HELM) is not merely another nonlinear iteration. It constructs the branch connected to an energized no-load network, evaluates it by analytic continuation, and can distinguish a certified collapse from insufficient series order. This gives it unusual value beside the Ipopt-based nonlinear `solve_pf`.

The trade-off is explicit: the present embedding supports constant-power and constant-impedance load components only. It cannot replace the nonlinear engine for arbitrary ZIP/exponential loads, IBRs, generators, controls, or OPF limits.

## 1. Two solvers, different evidence

| Question | HELM `solve_pf_helm` | Nonlinear `solve_pf` |
|---|---|---|
| Operating point | no-load germ, power series, Padé continuation | iterative nonlinear feasibility/optimization solve |
| Initialization | none | numerical initialization and local convergence matter |
| Branch | continuously connected to no-load state | feasible local point reached by the solve |
| No solution | can certify collapse in its supported embedding | failure is not a nonexistence proof |
| Collapse distance | `load_margin` on every solve | needs a separate loading study |
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
helm.status       # :converged | :diverged_no_solution | :max_order_reached
helm.residual     # full nonlinear current mismatch [A]
helm.load_margin  # estimated collapse loading multiplier
```

The final test is physical: returned voltages are substituted into the full nonlinear current mismatch. `solve_pf_helm` returns the standard result shape.

### Why this is more than avoiding an initial guess

Iterative solvers can fail because of scaling, a local basin, line search, or a poor start. HELM fixes the operational branch by continuation from no load. For the supported holomorphic embedding, persistent coefficient growth at ``s=1`` means the solvable branch ends before the requested loading; it is not an instruction to try a different initialization.

## 3. Interpret status and margin precisely

| Status | Meaning | Response |
|---|---|---|
| `HELM_CONVERGED` | operational branch reached; physical residual met tolerance | inspect voltage/margin and compare with nonlinear PF |
| `HELM_NO_SOLUTION` | certified collapse for the supported embedding at ``s=1`` | redesign/reduce loading; do not hunt for a start point |
| `HELM_MAX_ORDER` | finite series order was insufficient | increase `max_order`; it is not collapse |

Coefficient-tail Domb-Sykes extrapolation estimates the collapse multiplier ``\lambda^*``. A margin near 2 means the present supported load is about half of its collapse loading; a margin below 1 is consistent with collapse before ``s=1``. `NaN` means the tail cannot estimate a reliable radius, often because collapse is too remote for light or constant-impedance-only loading.

### Pitfall: calling every failed solve voltage collapse

Only `HELM_NO_SOLUTION` has the collapse meaning. `HELM_MAX_ORDER` is inconclusive. Likewise, failure of nonlinear `solve_pf` is not proof of nonexistence without a solvability analysis such as HELM on a common model.

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
3. Inspect residuals and HELM margin.
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

The rejection protects the analytic claim. Forcing a non-holomorphic load into the recursion makes the certificate and load margin meaningless. An explicit scope is more scientifically useful than silently solving a different model.

## 7. Publication checklist

Report the load representation, source/reference and neutral treatment, ideal-element model, HELM order/tolerance/status/residual/margin, and cross-solver comparison on common physics. State every excluded control and unsupported load component. HELM's value is not replacing nonlinear PF: it provides deterministic analytic reference, branch-aware no-solution evidence, and a collapse-distance indicator where iterative PF alone is most ambiguous.

For API details see [HELM power flow](../algorithms/helm.md).
