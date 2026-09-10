# Tutorial: taper, smoothing and the limits of aggregation

## Experiment 1: charging power is an upper bound

**Misconception:** every nonsmooth minimum needs smoothing or an integer model.

For an illustrative 40 kWh battery, let the AC acceptance cap be

```math
f(z)=\min\{7000,14000(1-z)\}\text{ W},\qquad z=E/40000.
```

This cap is concave. Imposing ``p^c\le7000`` and ``p^c\le14000(1-z)`` is exact.
It leaves the scheduler free to charge below the cap. An equality to the minimum
would instead prescribe a maximum-charging response and would be a different
model. `acceptance_formulation=:auto` performs this relation-aware lowering.

```@example ev_approx
using PowerOptLab
include(joinpath(dirname(pathof(PowerOptLab)), "..", "examples", "ev_charging.jl"))
r = EVChargingExamples.taper()
x = r.dispatch["taper"]
@assert solve_status(r).publishable
@assert isapprox(x.p_charge[1], 7000/1.35; atol=0.1)
(x.p_charge, x.energy_wh, x.diagnostics.acceptance_violation_w)
```

The example begins at ``z=0.5`` and uses a negative price to select the largest
feasible charge. With a one-hour constant dispatch and unity efficiency,
``z_2=0.5+p^c/40000``. Checking the terminal cap gives

```math
p^c\le7000-0.35p^c\quad\Rightarrow\quad p^c\le5185.185\ldots\text{ W}.
```

Checking only the initial cap would allow 7000 W for the entire hour, violating
the falling acceptance curve before the interval ended.

## Experiment 2: separate time discretization from smoothing

The model checks both energy boundaries. With a nonincreasing cap and affine
within-period energy this conservatively bounds constant dispatch throughout
the interval. It is not the continuous maximum-charge solution.

For this illustrative falling branch, the latter satisfies
``\dot z=0.35(1-z)`` and ``z(t)=1-0.5e^{-0.35t}``. Refining the interval grid
allows the constant dispatches to approach this continuous envelope from below.

```@example ev_approx
fine = EVChargingExamples.taper(dt_h=0.25, periods=4)
@assert solve_status(fine).publishable
continuous_energy = 40000*(1-0.5exp(-0.35))
@assert x.energy_wh[end] < fine.dispatch["taper"].energy_wh[end] < continuous_energy
(x.energy_wh[end], fine.dispatch["taper"].energy_wh[end], continuous_energy)
```

The analytical solution is a reference for this hand-specified curve, not
measured battery behavior. Power taper can depend on temperature, SoH, voltage,
and proprietary BMS decisions that this example omits.

## Experiment 3: a conservative smooth bound

For a nonconcave PWL cap, the intersection of all segment halfspaces is generally
too restrictive. `:auto` rejects an unresolved relation rather than silently
changing it. An explicit smooth approximation or `ComplementarityGraph()` stays
within the NLP/MPCC scope.

The optional `scripts/ev/ccopt.jl` regression uses breakpoints
``(0,0.5,1)`` with caps ``(7000,2000,0)`` W. Starting at 10 kWh in a 40 kWh
battery, a one-hour unity-efficiency interval remains on the first segment.
Its end-point bound is ``p^c\le4500-0.25p^c``, giving a 3600 W maximum.
The regression checks that analytical answer and the original physical cap,
and records its numerical tolerances separately from the charge/discharge
mode cases. It is a small solver experiment, not a scalability guarantee.

Although the concave example needs no smoothing, it is useful for comparing
representations. The library supports `LocalC2Formulation(width)` and
`SoftplusFormulation(width)`; width is in **energy-fraction units** here.

```@example ev_approx
smooth = EVChargingExamples.taper(LocalC2Formulation(0.01))
@assert solve_status(smooth).publishable
y = smooth.dispatch["taper"]
@assert maximum(y.diagnostics.acceptance_violation_w) < 0.1
(y.p_charge, y.energy_wh)
```

If the relation planner establishes
``\widetilde f(z)-f(z)\le\epsilon^+_P``, the corrected bound
``p^c\le\widetilde f(z)-\epsilon^+_P`` is conservative for the original scalar
upper relation in real arithmetic. Solver residuals are additional. Do not
infer conservativeness from a smooth curve's appearance.

A correction can make the cap negative near an original zero. Clipping it with
an unexamined `max(0,...)` changes the model and may lose smoothness. Use the
exact affine/concave representation where possible, restrict the operating
domain with scientific justification, or design and verify a zero-preserving
approximation. Disconnected intervals are removed from acceptance constraints,
so they remain exactly off even at a zero-cap state.

For fixed-efficiency PE accounting and a power-error bound in each interval,

```math
|\delta E_T|\le\sum_t\eta_c\Delta t_t\epsilon_{P,t}.
```

This is a trajectory-accounting bound for bounded power discrepancies. It is
not an objective-gap guarantee or a sensitivity theorem for optimized,
state-dependent trajectories. Keep four error sources distinct: curve fit,
smoothing, time discretization, and numerical solution. See
[physical error budgets](../formulations/error_budgets.md).

## Experiment 4: a fleet is not one battery

**Misconception:** summing energy capacities and kW ratings preserves charging
flexibility.

Consider two unity-efficiency vehicles with 1 kW chargers. A is connected only
in hour 1 and requires 1 kWh by that hour's end. B is connected only in hour 2
and requires 1 kWh by the horizon end. An aggregate description retaining only
2 kWh total demand and 2 kW nameplate power admits the schedule ``(0,2)`` kW.
The individual sessions cannot realize it: A misses departure, and B cannot
draw 2 kW. Time-varying availability caps repair part of this example; general
heterogeneous deadlines require the relevant cumulative/subset constraints or
explicit disaggregation.

For exact results under specified linear charging assumptions see
[Panda and Tindemans](https://arxiv.org/abs/2310.02729); robust aggregation with
uncertain requirements is developed by
[Mukhi et al.](https://arxiv.org/abs/2405.08232). Those guarantees do not
transfer automatically to nonconvex AC constraints, losses, or acceptance
curves. The implemented session model is a disaggregated deterministic building
block, not a certified aggregate-flexibility API.

## Experiment 5: perfect forecasts are not an operational guarantee

**Misconception:** an optimal day-ahead schedule ensures the driver leaves ready.

A one-kW vehicle needing one kWh may optimally postpone all charging to the
second of two hours. Departure after the first hour defeats the plan. Merely
re-solving at that departure cannot recover the lost opportunity. Useful
extensions include earlier reserve requirements, scenarios with common
pre-observation decisions, and explicit shortfall costs. Evaluate them against
unseen arrival/departure realizations, not only the optimization forecast.

`ChargingSession` currently treats these inputs as known. A future MPC driver
must transfer measured state between solves and distinguish predictions from
realized events; its optimization subproblems can remain NLP/MPCC. This tutorial
is a design experiment for that extension, not a claim that stochastic control
is already implemented.

## Minimum pilots and assignment choices

A continuous energy schedule is not necessarily an admissible pilot sequence.
ACN's [equipment API](https://acnportal.readthedocs.io/en/latest/acnsim/models.html)
includes continuous, deadband and finite-rate pilot models. A pilot is a current
permission, not a guarantee that a vehicle draws that current.

One possible exact MPCC encoding for an on/off permission is

```math
0\le y\perp(1-y)\ge0,\qquad I_{min}y\le u\le I_{max}y.
```

This encodes discrete choice despite using continuous variables. CCOpt does not
make the resulting combinatorial search globally easy. Finite-rate pilots and
vehicle assignment require further choices and replay. They are not currently
encoded by `EVSE`; no fractional pilot is advertised as field-ready.
