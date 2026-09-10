# Tutorial: artificial cycling and genuine energy arbitrage

## Experiment 1: the full battery at a negative price

**Misconception:** round-trip loss automatically excludes simultaneous charging
and discharging.

Take a full 40 kWh battery, one hour, 10 kW charge/discharge ratings,
``\eta_c=\eta_d=0.9``, and a requirement to remain full. Energy conservation gives
``p^d=0.81p^c``. The independent split can therefore draw 10 kW while exporting
8.1 kW: its energy is unchanged but it consumes 1.9 kW. A negative price rewards
that artificial dissipative mode.

```@example ev_modes
using PowerOptLab
include(joinpath(dirname(pathof(PowerOptLab)), "..", "examples", "ev_charging.jl"))
outer = EVChargingExamples.negative_price(:independent)
@assert solve_status(outer).publishable
(outer.dispatch["full"].p_charge, outer.dispatch["full"].p_discharge)
```

This is not justified by saying the vehicle alternates modes within the hour:
the implied duty fractions at the two 10 kW ratings sum to 1.81. A valid
averaged switching model would need time-sharing constraints and must also
respect intra-interval energy and electrical limits.

## Experiment 2: quantify the relaxed NLP leakage

With ``a=p^c/10000``, ``b=p^d/10000`` and ``ab\le\tau``, conservation implies
``b=0.81a``. Hence the largest charge power is

```math
p^c_{max}=10000\sqrt{\tau/0.81}\text{ W}.
```

At ``\tau=10^{-6}``, this is approximately 11.11 W, with 9 W discharge and
2.11 W net consumption. It is small, but it is not zero.

```@example ev_modes
relaxed = EVChargingExamples.negative_price(:relaxed; tolerance=1e-6)
x = relaxed.dispatch["full"]
@assert solve_status(relaxed).publishable
@assert isapprox(x.p_charge[1], 10000sqrt(1e-6/0.81); atol=0.2)
(x.p_charge, x.p_discharge, x.diagnostics.complementarity_product)
```

The script tightens Ipopt's feasibility settings and disables bound relaxation
for this analytical comparison. Solver residuals are separate from ``\tau``;
a default solver tolerance is not automatically a W-level accuracy promise.
Repeat the experiment in SI and per-unit coordinates and report both normalized
mode residuals and physical energy conservation.

For exact mathematical mode exclusion use `operation=:complementarity` with
[the CCOpt configuration](model.md). The executable `scripts/ev/ccopt.jl` checks
this full-battery example and a two-interval cycle on the complete AC network.
The replenished-cycle case requires a successful solve and physical residual
checks. The full-battery case is also a **solver characterization**: the tested
backend can report a local failure at this biactive zero-dispatch solution,
including in SI coordinates. Such a status cannot prove infeasibility, since
``p^c=p^d=0`` is analytically feasible. The runner records raw candidates
separately and checks that failed solves do not publish schedules. It does not
assert exact floating-point zero or global MPCC optimality.

## Experiment 3: distinguish initial inventory from arbitrage

**Misconception:** exporting in an expensive interval proves profitable repeated
charging/discharging operation.

A vehicle beginning at 20 kWh and allowed to leave empty can sell initial
inventory. That says little about the cost of replenishing it. Compare two
runs with identical ratings, 90% one-way efficiencies, and prices 0.05 and
0.25 currency/kWh:

```@example ev_modes
cycle = EVChargingExamples.arbitrage(replenish=true)
empty = EVChargingExamples.arbitrage(replenish=false)
@assert solve_status(cycle).publishable && solve_status(empty).publishable
(cycle.dispatch["cycle"].energy_wh,
 empty.dispatch["cycle"].energy_wh)
```

In the replenished run, 10 kWh is bought cheaply, 9 kWh enters storage, and
8.1 kWh can be sold later while returning to the original energy. The ideal-grid
energy margin is ``8.1(0.25)-10(0.05)=1.525`` currency, before degradation,
auxiliaries, network losses and other charges. These are dispatch economics,
not a lifetime profitability claim.

A terminal equality is one transparent boundary condition. A defensible
terminal value or a longer rolling horizon is another; simply removing the
boundary assigns a particular economic value to the remaining energy.

## Why a throughput penalty is not a theorem

A positive penalty can discourage cycling, but its required strength depends on
the full objective and constraints. It is not a universal replacement for mode
exclusion. Published sufficient-condition results apply to specified economic
dispatch or household formulations; their hypotheses must be established before
using them in an AC network model. See [the literature table](literature.md).
