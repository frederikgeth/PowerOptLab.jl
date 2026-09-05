# Controller convergence investigation

The baseline investigation below describes `92e42a3`. Its diagnostic script was
committed at `6bb8dfd`, before the production changes. Run from a separate
checkout of `6bb8dfd` to reproduce that baseline; running the script on the
current branch evaluates the corrected formulation. Both snapshots are retained.

## Implementation follow-up

At `5827acb`, production removes provably zero roots, uses direct smooth
selectors and nondegenerate normalized expression definitions, and uses upper
norms for capability denominators. Current-linear losses and negative-sequence
voltage smoothing use BMOPFTools helpers. Dimensionless scale-selector widths
are independent of the OPF power base, and curve-free voltage initialization
uses the local nominal voltage rather than a fixed 230 V anchor.

The wrappers still perform one solve with fixed smoothing. No continuation,
retry, exact-law network re-solve, or relaxed result-publication rule was added.

`controller_convergence_after.toml` records 18 fresh single-solve comparisons
on Julia 1.12.6 with the same solver/package versions as the baseline.

| Fixture | Ipopt default | MadNLP default |
|---|---|---|
| Ripple 0.7 V | LOCALLY_SOLVED, 36 iterations | LOCALLY_INFEASIBLE, 106 iterations |
| Ripple 0.1 V | LOCALLY_SOLVED, 42 iterations | LOCALLY_SOLVED, 41 iterations |
| LCL grid target | LOCALLY_SOLVED, 6 iterations | LOCALLY_SOLVED, 23 iterations |

All nine Ipopt configurations converge. MadNLP converges in six of nine runs,
including an unscaled 0.7 V run, default and unscaled 0.1 V runs, and all three
LCL runs. Three ripple runs still return local-infeasibility status despite tiny
original-model residuals. Its initialization and configuration sensitivity
therefore remains an open issue. The corrected LCL controller has about 0.51 mA
exact-versus-smooth command discrepancy, within the 1 mA regression budget.
This is larger than the old unscaled Ipopt discrepancy: the changes prioritize
reliable convergence within a stated physical accuracy budget, rather than
promising that every numerical discrepancy decreases.

The regression suite checks all 36 combinations of four difficult cases, three
independent starts and three power bases (0.1, 1 and 10 MVA). It checks original
constraint residuals and physical current/power/ripple limits, and retains the
original 0.1 V ripple case. Ripple and LCL tests now use ordinary Ipopt settings
without fixture-specific scaling or barrier overrides.

A separate removal of the unused plant loss epigraphs regressed eight assertions
in the sequence/ripple test, so that removal was reverted. Their retention is
an empirical solver-path workaround, not a physically necessary loss term or
a mathematically justified regularization. Wider plant conditioning and a
reduced upstream MadNLP restoration reproducer remain separate work.

## Baseline investigation

The skipped ripple regression is numerically solvable. The existing controller
formulation nevertheless contains a demonstrable constraint-qualification defect;
changing solver options alone does not remove it. MadNLP reproduces failures, so
switching solvers is not a sufficient repair.

## Reproduce

From the repository root, first instantiate the pinned ordinary environment if
needed, then create a separate, empty diagnostic environment:

```sh
julia --project=. scripts/instantiate_pinned.jl
julia scripts/diagnostics/setup_controller_convergence.jl /tmp/pol-controller-diagnostics
POL_DIAG_OUTPUT=/tmp/controller-results.toml julia --project=/tmp/pol-controller-diagnostics scripts/diagnostics/controller_convergence.jl > /tmp/controller-solvers.log 2>&1
```

MadNLP is an optional diagnostic dependency, not a new package or CI dependency.
The setup preserves the versions in the local manifest and verifies the BMOPFTools
commit. The experiment uses CPU solvers, exact equality enforcement, `tol=1e-8`,
`max_iter=1000`, and one BLAS thread. Each run builds a fresh model with the same
physical data, formulation, loss objective and prescribed starts. The two solver
implementations still have different internal initialization and algorithms.

The default matrix has 48 runs: three fixtures, two solvers, four option sets, and
two formulations. Restrict it with comma-separated `POL_DIAG_CASES`
(`ripple_0.7,ripple_0.1,lcl_grid`), `POL_DIAG_SOLVERS` (`ipopt,madnlp`), or
`POL_DIAG_VARIANTS` (`default,unscaled,unscaled_adaptive,unscaled_unrelaxed`).
`POL_DIAG_ABLATIONS=no` runs only the original formulation. The script uses internal
handles to inspect failed candidates; these are diagnostic code, not a public API.

`unscaled` disables solver NLP scaling, while retaining the model's per-unit and
explicit row normalization. `unscaled_adaptive` additionally selects Ipopt's
adaptive barrier or MadNLP's quality-function update; these are analogous options,
not identical implementations. `unscaled_unrelaxed` instead disables bound
relaxation. Other solver defaults are retained. See the
[MadNLP options](https://madsuite.org/MadNLP.jl/stable/options/).

## Measured findings

The adjacent `controller_convergence_results.toml` records the baseline local matrix
on Julia 1.12.6 / aarch64 macOS, BMOPFTools
`5b51d2f361dab91bd7c16711019584407da79ed8`, Ipopt.jl 1.15.0 and MadNLP 0.10.1.
This is a platform-specific diagnostic snapshot, not a cross-platform guarantee
or a statistical estimate of failure probability.

Selected default-option runs (`tol` and iteration limit as above):

| Fixture | Solver | Original | Exact zero-root elimination |
|---|---|---|---|
| Ripple limit 0.7 V | Ipopt | LOCALLY_SOLVED, 111 iterations | LOCALLY_SOLVED, 29 iterations |
| Ripple limit 0.1 V | Ipopt | ITERATION_LIMIT, 1000 iterations | LOCALLY_SOLVED, 38 iterations |
| Ripple limit 0.7 V | MadNLP | NUMERICAL_ERROR, 116 iterations | LOCALLY_INFEASIBLE, 90 iterations |
| Ripple limit 0.1 V | MadNLP | LOCALLY_INFEASIBLE, 452 iterations | LOCALLY_INFEASIBLE, 308 iterations |
| LCL grid-current target | Ipopt | ALMOST_LOCALLY_SOLVED, 35 iterations | Identical: no zero-root rows |
| LCL grid-current target | MadNLP | LOCALLY_INFEASIBLE, 153 iterations | Identical: no zero-root rows |

The original 0.7 V MadNLP failure was also reproduced through the public
`solve_controlled_inverter(...; optimizer=MadNLP.Optimizer)` API. Its result remains
non-publishable and its published plant values are NaN, as required by the outcome
contract.

### 1. Exact zero-root lifting is a proven formulation defect

`_safe_direction_scale_implicit!` checks whether an offset is a *literal* numeric
zero before lifting its magnitude. Its later headroom check also recognizes
identically-zero JuMP expressions. That inconsistency introduces two redundant
roots in each ripple model: apparent-power offset and ripple-power offset.

They have precisely the form `a*y^2 == 0`, with `y >= 0`. At every feasible point,
`y=0` and the equality gradient is zero. Consequently LICQ and MFCQ fail for the
full formulation. This is an algebraic conclusion, not an inference from a solver
status. Tightening a squared-row residual tolerance also controls the root only
as its square root: `abs(y) <= sqrt(tolerance/a)`.

The experiment recognizes these rows using exact coefficients, deletes them and
fixes their variables to zero. This preserves the exact feasible set and objective;
it adds no smoothing or tolerance relaxation. Both ripple cases contain two such
rows; LCL contains none. The improvement in Ipopt, especially the recovered 0.1 V
case, demonstrates a material effect on this solver's convergence. It does not
establish that these are the only causes of failure.

### 2. MadNLP restoration status needs separate investigation

With zero roots eliminated and NLP scaling disabled, MadNLP returns
`LOCALLY_INFEASIBLE` for the 0.7 V fixture even though the maximum original-row
violation is about `3.9e-12` and exact-versus-smooth current discrepancy about
`2.3e-10 A`. Disabling bound relaxation does not remove this behavior. Other runs
fail with large violations; those are quite different outcomes.

In MadNLP 0.10.1, `robust!` in `src/IPM/solver.jl` returns
`INFEASIBLE_PROBLEM_DETECTED` when the restoration subproblem's residuals meet its
tolerance. Returning to the regular phase separately requires filter acceptance
and relative infeasibility reduction. The observed tiny-residual local-infeasibility
exits warrant a reduced upstream reproducer of that termination path. They do not
prove physical infeasibility or global nonexistence, and we do not override the
returned status. An upstream solver fix has not been established here.

`max_original_row_violation` evaluates the JuMP constraints and variable bounds,
without the solver's additional scaling. It mixes the formulation's working units
and is not a universal physical error norm. The separate candidate fields report
DC ripple in volts, conductor current in amperes, command-to-terminal disagreement,
and exact-to-smooth controller disagreement. Failed candidates are never extracted
by fabricating a successful `SolveStatus`. Low primal residuals alone are not a
stationarity or complementarity certificate; MadNLP's reported residuals are also
saved separately.

### 3. Shifted smoothing is present, but only in part of the implementation

PowerOptLab's inverter current-linear loss term uses a local
`_smooth_magnitude(radicand, epsilon) = sqrt(radicand + epsilon^2) - epsilon`.
The controller still uses `_implicit_sqrt!` and lifted min/max expressions. Thus
the earlier loss-norm repair never removed all lifted roots from the controller.
When `a_loss == 0`, the plant additionally retains unused current-magnitude
epigraphs. Their removal has not been studied in this matrix.

The pinned BMOPFTools exports `smooth_norm(ctx, x, y; scale, eps_rel, ...)` and a
vector overload. These return shifted expressions and record differentiability
annotations. PowerOptLab does not currently call them. Consolidating the loss
implementation onto that helper is sensible if its existing physical epsilon and
one-sided loss-error budget are preserved.

However, the shifted norm underestimates the exact norm by up to epsilon. Putting
it directly into an upper capability limit or a limiter denominator can enlarge
the accepted physical region or current command. The controller needs a derived
error allowance or conservative correction, plus physical replay. Further, the
shifted norm's curvature at the origin scales as `1/epsilon`; smoothness removes
the nondifferentiability but does not guarantee good conditioning for arbitrarily
small epsilon. For an identically-zero offset, using literal zero is exact and
requires no approximation at all.

## What remains

1. Remove provably zero roots in production, then validate the wider controller
   suite on supported Julia versions and Linux. Prior removal changed convergence
   elsewhere, so this diagnostic ablation is not yet a validated production fix.
2. Audit the remaining lifted selectors and unused loss epigraphs. Use the shared
   smooth helper where mathematically appropriate, with declared physical error
   budgets and exact-law checks. Do not replace every square root blindly.
3. Reduce the MadNLP restoration exit to a small standalone upstream case and
   examine its termination logic against original-problem residuals.
4. Fix the remaining Linux regression before calling PR #44 green. At head
   `92e42a3`, Julia 1.10 CI passed the restored ripple assertions but failed five
   assertions in the 255 V balanced Volt-watt/OpenDSS case after Ipopt returned
   `NUMERICAL_ERROR`. That is a distinct fixture and is outside this matrix.
   The Julia 1 job and documentation build passed at that head.

The evidence supports a formulation cleanup and better diagnostics, not relaxed
acceptance of failed solves or a claim that changing the barrier strategy has
solved numerical robustness generally.
