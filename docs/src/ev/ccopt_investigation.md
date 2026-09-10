# Numerical investigation: a full battery that should do nothing

This is a small reproducible sensitivity study, not a CCOpt benchmark or a
universal configuration recommendation. It records unsuccessful candidates as
well as successful references. Run `scripts/ev/zero_dispatch.jl` in the optional
solver environment described in [the model guide](model.md).

## An analytical solution before solver tuning

A single connected interval starts and ends at 40 kWh. Both AC power ratings
are 10 kW and both one-way efficiencies are 0.9. There is no trip or self-discharge.
Energy conservation and exact mode exclusion require

```math
0=0.9p^c-p^d/0.9,\qquad p^cp^d=0,\qquad p^c,p^d\ge0.
```

Hence ``p^d=0.81p^c`` and the only possible powers are ``p^c=p^d=0``.
The port is attached to a resistive AC feeder with a 10 kW baseline load and a
negative source price. The negative price rewards loss-making simultaneous
cycling whenever the numerical representation permits it.

The reference experiment restricts discharge to zero and solves the remaining
AC NLP with Ipopt. Energy conservation then forces charging to zero as well.
This reduction is equivalent **for this one-interval, fixed-inventory case**.
It is not valid for a multi-interval V2G cycle, where charging and discharging
may legitimately occur at different times. The reference is labeled separately;
it is not evidence that CCOpt solved the unreduced MPCC.

## What was varied

The experiment uses the installed CCOpt 0.1.0/MadNLP 0.10.1 stack, runs in both
SI and per-unit coordinates, and changes one configuration at a time relative
to the declared baseline. A few cases intentionally combine tighter tolerances
with a different strategy; they are comparisons, not a factorial causal study.

| Configuration | Observed maximum absolute raw power, approximately | Accepted MPCC result? |
|:--|:--|:--|
| Proportional continuation, NLP tolerance `1e-8`, continuation floor `1e-12` | 1.1605 W | No |
| Explicit zero device initialization | 1.1605 W | No |
| Network power base changed from 1 MVA to 10 kVA | 1.1605 W | No |
| NLP tolerance `1e-10`, continuation floor `1e-16` | 1.1116 W | No |
| Rolloff continuation with the tighter settings | 1.110–1.111 W | No |
| `respect_comp_bounds=true` with the tighter settings | Strongly coordinate-dependent failed candidates | No |
| Tighter settings and `bound_relax_factor=0` | 0.0274 W (PU), 0.0326 W (SI) | No |

The table records the original local macOS run. Linux CI subsequently returned
`LOCALLY_SOLVED` for the no-bound-relaxation case at approximately 0.0335 W:
per-unit on Julia 1.10, and both coordinate systems in the Julia `1` job. This still
exceeds the experiment's 0.01 W budget ([CI evidence](https://github.com/frederikgeth/PowerOptLab.jl/actions/runs/34478727627)). Thus the local failure status is not a
portable expectation, and a success status is not an accuracy certificate.

The saved TOML includes exact statuses, settings, raw powers, energy residuals,
package versions and source fingerprints. `s_base` is inactive in SI mode, so
its SI repeat is a control, not another scaling transformation. Zero device
starts retain the engine's network starts; they are not an exact feasible AC
initial point. Local failure/infeasibility statuses do not disprove the known
zero-dispatch solution.

## Why lowering the continuation floor was not enough

With equal 10 kW ratings, energy conservation gives the normalized product

```math
ab=0.81\left(p^c/10000\right)^2.
```

If the effective allowed product slack is ``\delta``, the negative-price
objective can select approximately

```math
p^c=10000\sqrt{\delta/0.81}.
```

For ``\delta=10^{-8}``, this predicts **1.1111 W**. That matches the plateau
in the tighter-continuation cases. In the tested implementation, the Scholtes
relation has an upper bound of zero after subtracting the continuation parameter,
and CCOpt passes MadNLP's `bound_relax_factor` into initialization. The observed
plateau is consistent with this additional bound relaxation. Removing it
substantially reduces leakage, but does not establish convergence of the
biactive MPCC. See the [CCOpt source](https://github.com/madsuite-org/CCOpt.jl),
particularly `src/Models/scholtes_relaxation.jl` and
`src/Solvers/relaxation/solver.jl` for the installed release.

This is why three quantities must remain distinct: the continuation parameter,
the backend's bound relaxation, and the stopping criteria. A small energy
residual can coexist with nonzero artificial cycling; check the individual
powers and complementarity as well.

## Decision for the library

The production model and publication contract are unchanged by this experiment.
There is no automatic acceptance of an “almost solved” point, no clipping of raw
powers, and no automatic branch selection. The experiment uses a 0.01 W power
budget and 0.01 Wh energy budget. `publishable` records the library's solver-status
contract; `study_accepted` additionally requires both physical budgets. Only
study-accepted powers populate `study_charge_w` and `study_discharge_w`; other
study powers are NaN while raw candidates and status-published dispatch remain
available for diagnosis. Neither budget is loosened for platform differences.
Other studies must choose their own justified budgets.

CI asserts that the analytic reference meets both requirements and that the
study gate rejects the Linux status-success/physical-failure counterexample.
The remaining configurations are characterization cases, not a promise that
every local solver success meets an independently chosen watt-level tolerance.

The analytic branch reference verifies the AC network separately. The existing
replenished-cycle and nonconcave-acceptance regressions continue to exercise
successful CCOpt compositions. None of these small examples establishes global
optimality or an independently checked MPCC stationarity certificate.
A general active-set/presolve integration would need its own model-equivalence
conditions and validation across multi-period cases before adoption.
