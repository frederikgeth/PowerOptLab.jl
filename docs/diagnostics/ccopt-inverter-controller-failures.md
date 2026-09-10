# Full inverter controller comparison for CCOpt developers

This is a PR attachment for the smooth/complementarity controller implementation
on `codex/advanced-ibr-reliability`. It compares **CCOpt v0.1.0** with the same
physical controller represented by the established smooth graph and solved by
Ipopt or MadNLP. These are controlled AC power-flow problems, not isolated PWL
curve tests. The controller encoding changes; the plant, network, request,
measurement convention, and physical protection policy are held fixed.

No failed or nearly feasible solver status is promoted to success. In particular,
`LOCALLY_INFEASIBLE` is the solver's local termination report, not a proof that
the exact controller has no feasible equilibrium. A smooth solution alone also
does not prove existence of an exact equilibrium.

## Observed results (2026-09-07)

Julia 1.12.6 on macOS/aarch64; CCOpt 0.1.0, MadNLP 0.10.1, Ipopt 1.15.0,
JuMP 1.30.1, MOI 1.53.0, MathOptComplements 0.1.1 and NLPModelsJuMP 0.13.6.
The [compact result bundle](ccopt-inverter-controller-results.json) includes all
package versions, exact settings and source hashes.

| Case | Smooth / Ipopt | Smooth / MadNLP | CCOpt default | CCOpt tighter relaxation |
|:--|:--|:--|:--|:--|
| `balanced_constant` | Solved | Locally infeasible | Almost solved | Locally infeasible |
| `balanced_curve` | Solved | Solved | Almost solved | Solved |
| `unbalanced_droop` | Solved | Solved | Locally infeasible | Numerical error |
| `apparent_power_saturation` | Solved | Locally infeasible | Almost solved | Solved |
| `lcl_grid_target` | Solved | Solved | Locally infeasible | Locally infeasible |

“Solved” means `LOCALLY_SOLVED`; “almost solved” means `ALMOST_LOCALLY_SOLVED`
and is **not** a strict success. The tighter method strictly solves two of five
exact-controller cases, with full current-replay errors below `4e-11 A`. Both
also pass the public snapshot API publication audit (10 API assertions passed).

### Paired `bound_relax_factor=0.0` comparison

The final controller implementation was rerun for all 35 method/case combinations
in the same environment. The original 20 outcomes reproduced exactly. In each
new pair, **only `bound_relax_factor` changes**; `respect_comp_bounds` remains at
its default. MadNLP 0.10.1's default bound relaxation is `1e-8`. The tighter
CCOpt comparison changes the explicitly configured `1e-12` to `0.0`.

| Method | Original strict successes | With `bound_relax_factor=0.0` |
|:--|--:|--:|
| Smooth / MadNLP | 3/5 | 2/5 |
| CCOpt, default relaxation | 0/5 | 1/5 |
| CCOpt, tighter relaxation | 2/5 | 0/5 |

| Case | Smooth / MadNLP, zero bounds | CCOpt default relaxation, zero bounds | CCOpt tighter relaxation, zero bounds |
|:--|:--|:--|:--|
| `balanced_constant` | Locally infeasible | Locally infeasible | Locally infeasible |
| `balanced_curve` | Locally infeasible | Solved | Locally infeasible |
| `unbalanced_droop` | Locally infeasible | Locally infeasible | Iteration limit |
| `apparent_power_saturation` | Solved | Slow progress | Locally infeasible |
| `lcl_grid_target` | Solved | Locally infeasible | Locally infeasible |

- **Smooth MadNLP:** zero bound relaxation resolves apparent-power saturation,
  but the previously solved balanced volt-watt and unbalanced-droop cases now
  return `LOCALLY_INFEASIBLE`. Their largest scalar equality residuals remain
  about `1.78e-15` and `1.17e-15`, respectively, in raw row units. The curve-free
  case remains locally infeasible; the LCL case remains solved. Their strict
  termination classifications are sensitive to this option even when primal
  residuals are extremely small.
- **CCOpt with default relaxation:** zero bounds resolves the balanced volt-watt
  case (`LOCALLY_SOLVED`, `4.27e-6 A` replay error, `9.74e-9` complementarity
  error). The other four cases do not strictly solve. Curve-free replay error
  increases to `1.33e-2 A`, while the LCL replay error decreases to `9.94e-7 A`
  without a successful termination status.
- **CCOpt with tighter relaxation:** replacing `1e-12` by zero loses both prior
  successes (balanced volt-watt and apparent-power saturation). Unbalanced droop
  reaches the 500-iteration limit. The LCL candidate still has very small replay
  error (`4.88e-11 A`) and complementarity error (`6.28e-13`) despite reporting
  `LOCALLY_INFEASIBLE`.

Zero bound relaxation is therefore **not a general fix for these cases**. This
comparison does not establish why the statuses differ; the complete candidate
residuals and exact options are retained for the solver developers. No solver
default or publication rule was changed in PowerOptLab. The expanded comparison
completed successfully with 38 assertions checking the records and baseline.

### Failures to include in the PR

1. **Unbalanced droop — baseline smooth backends succeed; no tested CCOpt configuration strictly succeeds.** Default
   CCOpt returns `LOCALLY_INFEASIBLE`, with `7.93e-7 A` current-replay error and
   `1.33e-8` normalized complementarity error. Tighter relaxation returns
   `NUMERICAL_ERROR` despite `1.30e-8 A` replay error and `1.24e-10`
   complementarity error. Inspect restoration/linear-system termination and
   primal/dual stopping conditions; these residuals alone do not establish
   stationarity. With zero bounds, smooth MadNLP also loses strict success and
   tighter CCOpt reaches its iteration limit.
2. **LCL grid-current target — both smooth backends succeed; CCOpt does not.**
   Default CCOpt returns `LOCALLY_INFEASIBLE`, with `1.89e-4 A` replay error and
   `1.49e-6` complementarity error. Tighter relaxation still returns
   `LOCALLY_INFEASIBLE`, although replay error falls to `8.05e-11 A` and
   complementarity error to `9.83e-13`. The largest scalar equality residual is
   about `1.04e-9` in its stamped units.
3. **Balanced volt-watt — both smooth backends succeed; default CCOpt reaches
   only acceptable convergence.** The result is `ALMOST_LOCALLY_SOLVED`, with
   `1.25e-5 A` replay error. The explicit tighter-relaxation method resolves this
   case (`LOCALLY_SOLVED`, `3.00e-11 A` replay error). Default-relaxation CCOpt
   also solves it with zero bounds; tighter-relaxation CCOpt and smooth MadNLP
   lose their prior successes when their bound relaxation is set to zero.
4. **Curve-free balanced control — shared MadNLP/CCOpt status issue.** Ipopt
   succeeds, but smooth MadNLP returns `LOCALLY_INFEASIBLE` with maximum equality
   residual `1.78e-15` in raw row units. Default CCOpt returns
   `ALMOST_LOCALLY_SOLVED`; tighter CCOpt returns `LOCALLY_INFEASIBLE` with
   `2.80e-11 A` replay error and `5.67e-13` complementarity error. This should
   not be attributed solely to PWL curves or to the complementarity encoding.
5. **Apparent-power saturation — shared MadNLP/default-CCOpt difficulty.**
   Smooth Ipopt succeeds. Smooth MadNLP returns `LOCALLY_INFEASIBLE` with tiny
   equality residuals; default CCOpt returns `ALMOST_LOCALLY_SOLVED`. The
   tighter CCOpt method succeeds (`8.38e-12 A` replay error). Setting bound
   relaxation to zero resolves the smooth MadNLP case, but loses the tighter
   CCOpt success.
6. **CCOpt MOI diagnostic getter fails in every MPCC run.** Requesting
   `MOI.RawStatusString()` raises `FieldError: type CCOpt.RelaxationSolver has
   no field opt`. The v0.1.0 getter reads `optimizer.solver.opt`, while the
   relaxation solver has `opts` and `ipm`. The experiment runner records this
   secondary error separately without discarding the returned candidate.

These observations are reproducible counterexamples for investigation, not a
claim that every failed case is a CCOpt defect. In particular, equalities and
complementarity can be very accurate while dual convergence or stationarity is
unresolved. The full controller contains nonnegative norm equalities, redundant
inactive protection constraints and nonunique selector weights at ties.

## Reproduce

From the repository root:

```sh
julia --project=. scripts/instantiate_pinned.jl
julia scripts/formulations/setup.jl /tmp/pol-controller-encodings
POL_FORMULATION_RESULTS=/tmp/controller-encoding-comparison.toml \
  julia --project=/tmp/pol-controller-encodings \
  scripts/formulations/inverter_controller_comparison.jl
```

The optional environment preserves the production BMOPFTools source revision.
The comparison writes raw outcomes, candidate residuals, exact solver options,
package versions, and source/manifest hashes to TOML. CI runs the same script on
Julia 1.10 and the current stable Julia and uploads the result bundle. Every row
uses a fresh model and one solve; the tighter-relaxation method is an explicitly
separate experiment, not an automatic retry or replacement result.

## Cases and settings

All cases use `s_base=1e6`, three-leg inverters, RMS phase phasors, and the loss
selection objective. Each case's exact parameters are in the comparison script.

- `balanced_constant`: balanced feeder, average-voltage policy without curves,
  12 kW available, 20 kVA / 40 A ratings. This exercises exact limiting even with
  no PWL curve to blame.
- `balanced_curve`: the same feeder and ratings, with a flat-tail volt-watt curve.
- `unbalanced_droop`: unbalanced feeder, worst-phase volt-watt/volt-var, net
  conflict resolution, negative-sequence admittance droop and 50% ripple blend.
- `apparent_power_saturation`: 12 kW availability exceeds the 8 kVA rating;
  the controller and final converter capability allocator must saturate.
- `lcl_grid_target`: a damped LCL filter with 40 A converter / 35 A grid ratings,
  grid-current target, and 9 kW availability.

Smooth Ipopt, smooth MadNLP and the default CCOpt method use `tol=1e-8` and
`max_iter=500`. The separate tighter CCOpt method uses `tol=1e-9`,
`bound_relax_factor=1e-12`, and
`ProportionalRelaxationUpdate(sigma_mu_ratio=1e-4, sigma_min=1e-14)`.
Three additional paired methods set only `bound_relax_factor=0.0`: smooth
MadNLP, CCOpt with its default relaxation, and CCOpt with the tighter relaxation.
The last pair changes this setting from `1e-12` to `0.0`; all other options,
including `respect_comp_bounds` (left at its default), are unchanged.
The exact encoding uses MOI complementarity, bridged by MathOptComplements;
NLPModelsJuMP supplies CCOpt's nonlinear model interface.

## What to inspect

The TOML candidate metrics include:

- full exact-firmware current replay error in amperes, evaluated at the candidate
  voltage and filter state and including the final plant-capability backoff;
- normalized complementarity products and negative-component violations;
- maximum scalar equality and inequality violations, plus the worst equality;
- converter current and apparent-power margins in A and VA.

Scalar row residuals retain their **stamped row units**; they are not a common
physical norm and should not be compared across unrelated constraints without
considering scaling. Candidate diagnostics remain available in the research
runner even when the solver fails. They do not make that candidate publishable.
The snapshot API additionally checks exact current replay and complementarity
before publishing a strictly successful result.

The exact controller uses complementary simplex weights and nonnegative slacks
for max/min. With exact defining rows and nonnegativity, product tolerance τ
bounds selector output error by `2τ*scale`, including at ties. Weights can be
nonunique at ties. Exact nonnegative norm roots can be degenerate at zero;
sequence components are materialized before squaring to avoid cancellation near
balanced operation. These properties may matter to regularization/restoration.
