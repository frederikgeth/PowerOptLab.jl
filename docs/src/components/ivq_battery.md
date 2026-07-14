# Current–voltage (IVQ) battery

> **Kind:** Component model · **Maturity:** prototype · **Direction:** forward · **Temporal:** single-snapshot + multi-period

[`IVQBattery`](@ref) models grid storage in the **voltage–current–charge**
variable space instead of the power–energy space of [`StorageDevice`](@ref). It
is the Julia counterpart to the model of

> P. Aaslid, F. Geth, M. Korpås, M. M. Belsnes, O. B. Fosso,
> *Non-linear charge-based battery storage optimization model with bi-variate
> cubic spline constraints*, **Journal of Energy Storage 32 (2020) 101979**.
> [doi:10.1016/j.est.2020.101979](https://doi.org/10.1016/j.est.2020.101979)

Where the power–energy ("PE") model — [`StorageDevice`](@ref) plus the
state-of-charge balance in [`solve_multiperiod_opf`](@ref) — tracks energy with a
**fixed** round-trip efficiency, the IVQ model tracks **charge** and represents
the cell **terminal voltage explicitly**, so voltage and current limits are
enforced *individually* rather than folded into a conservative power limit. Two
consequences matter for dispatch:

- **Safe operation up to the cell boundaries.** A PE-model must pad its power and
  SoC limits to guarantee the true voltage/current limits are never crossed; the
  IVQ model carries those limits directly and can use the full envelope.
- **Round-trip efficiency is derived, not assumed.**
  ``\eta = f_v(\text{soc}, i)\,/\,f_v(\text{soc}, -i)`` falls out of the voltage
  curve and *worsens with current* — no `eff_charge`/`eff_discharge` parameters.

## Terminal-voltage model

The paper represents the cell voltage as a bivariate cubic spline
``v = f_v(\text{soc}, i)`` fitted to multi-C-rate discharge data. That surface
needs data published for almost no chemistry, so PowerOptLab instead uses the
equivalent-circuit (Thévenin / *Rint*) decomposition

```math
v(\text{soc}, i) = \mathrm{OCV}(\text{soc}) - i \cdot R(\text{soc}),
\qquad i > 0 \ \text{discharge},
```

which encodes the same physics — terminal voltage drops under discharge and rises
under charge, and the round-trip efficiency above emerges — while only requiring
an **open-circuit-voltage curve** (abundant and citable for every chemistry) and
an **internal resistance** (from HPPC pulse data). The pack scales the cell by
its series/parallel counts: ``v_\text{pack} = n_\text{series}\,v_\text{cell}``,
``i_\text{pack} = n_\text{parallel}\,i_\text{cell}``, and the DC power is
``p_\text{dc} = v_\text{cell}\,i_\text{cell}\,n_\text{series}\,n_\text{parallel}``.

Three fidelity levels, cheapest first — each a [`BatteryChemistry`](@ref):

| Constructor | OCV(soc) | R | Data needed |
|---|---|---|---|
| [`thevenin_chemistry`](@ref) | constant | fixed | a nominal voltage + a resistance |
| [`linear_chemistry`](@ref) | linear `v_empty → v_full` | fixed | two endpoint voltages + a resistance |
| [`tabulated_chemistry`](@ref) | monotone cubic (PCHIP) | fixed or `R(soc)` table | an OCV curve (± HPPC) |

The PE-model is the further special case of `thevenin_chemistry` with `R = 0`
(constant voltage ⇒ power ∝ current).

!!! note "OCV(soc) is smooth and monotone by construction"
    In multi-period, `OCV(soc)` is evaluated at the SoC *variable*, so it must be
    twice-differentiable for the interior-point solver. `linear`/`thevenin` OCV is
    affine (embedded directly); `tabulated` OCV uses a **monotone cubic (PCHIP)**
    interpolant — smooth (C¹) *and* shape-preserving, so it registers as a smooth
    JuMP operator without overshoot. A plain cubic spline is C² but can overshoot
    outside its knots, handing an optimiser free voltage (= free energy) — the
    very defect this model avoids. OCV is also held **flat** outside the fitted
    range, and every chemistry carries a `soc_min`/`soc_max` window the solve
    clamps to.

## Reusing the converter

The AC-side converter is **not** re-implemented. An [`AdvancedInverter`](@ref)
owns the point-of-connection coupling, output filter, converter losses, and
apparent-power / three-phase topology limits; the battery attaches to its **DC
port**. In this first version the coupling is at the DC power level,

```math
v_\text{cell}\,i_\text{cell}\,n_\text{series}\,n_\text{parallel}
   \;=\; P_\text{dc}^{\text{(inverter)}} \quad [\text{W}],
```

which is exact for the default single-phase inverter (whose `v_dc` only enters an
optional modulation cap). Making `v_dc` a shared decision variable equal to the
battery terminal voltage — so the three-phase switching-polytope sees the true,
SoC-dependent DC rail — is the natural next step; until then set the inverter's
`v_dc` to the nominal pack voltage.

## Usage

```julia
using PowerOptLab

chem = lfp_chemistry()                       # or nmc/nca/lead_acid/leaf, or build your own
inv  = AdvancedInverter(id="bat", bus="poc", s_max=5000.0)
bat  = IVQBattery(id="bat", bus="poc", chemistry=chem,
                  n_series=100, n_parallel=1, soc_init=0.6, inverter=inv)

r = solve_ivq_battery(net, bat; objective=:max_export)
r.p_poc      # power delivered to the grid (W)
r.v_cell     # cell terminal voltage under load (V) — below OCV when discharging
r.i_cell     # signed cell current (A, > 0 discharge)
r.p_dc       # DC-link power = p_conv + p_loss
```

`objective` is `:max_export` (discharge to a cell voltage/current or converter
limit), `:max_charge`, or `:min_loss` (with a `p_set` delivery target). SoC is
fixed at `soc_init` — this single operating point reproduces the paper's second
case study: the deliverable power set by the cell's own voltage and current
limits, which a power-only model misses.

## Multi-period arbitrage

[`solve_multiperiod_ivq`](@ref) co-optimises a chronological sequence of
snapshots with the SoC linking the periods — now `soc` is a *decision variable*,
so `OCV(soc)` enters the model as a function of it (embedded for affine
chemistries; a registered smooth operator for tabulated ones).

```julia
nets = [net_expensive, net_cheap]            # e.g. differing voltage_source `cost`
bat  = IVQBattery(id="bat", bus="poc", chemistry=lfp_chemistry(),
                  n_series=300, n_parallel=1, soc_init=0.5, inverter=inv, cyclic=true)
res  = solve_multiperiod_ivq(nets, [bat]; dt_h=1.0)
res.dispatch["bat"].soc      # SoC trajectory (length T+1)
res.dispatch["bat"].i_cell   # signed cell current per period
```

The charge balance ``q_{t+1} = q_t - i\,\Delta t`` uses the `integration` rule on
the device: `:trapezoidal` (default) averages consecutive period currents;
`:forward` — exact for piecewise-constant period currents — takes the period's
own current. Terminal state is set by `cyclic=true` (return to `soc_init`) or
`soc_final`.

!!! warning "Nonconvex — expect local optima"
    The coupled cell + inverter model over several periods is nonconvex
    (`p = OCV(soc)·i` is bilinear on top of the rectangular power flow), so
    `solve_multiperiod_ivq` finds a *local* optimum and does not converge for
    every configuration. A non-`LOCALLY_SOLVED`/`OPTIMAL` status returns `NaN`
    trajectories rather than an unconverged point — always check the status. A
    PE-model warm start (solve [`solve_multiperiod_opf`](@ref) first, seed the
    dispatch) is the planned robustness improvement.

### Scaling: per-unit and SI

The battery is modelled in **pack quantities, per-unit on its own DC bases**
(pack nominal voltage and max discharge current), so the current and voltage
variables are ≈ O(1) *independent of the engine's `s_base`*. This is what keeps
the coupled solve conditioned — a naïve SI formulation mixes cell volts (~3.5 V),
network kV, and kW in one problem, spanning five orders of magnitude, and the
solver then returns poor points. The DC↔AC coupling carries the `s_base` factor
explicitly, so `solve_ivq_battery` and `solve_multiperiod_ivq` work in both the
engine's SI (`per_unit=false`, the default) and per-unit (`per_unit=true`) modes;
results are returned in SI either way.

## Chemistry library and data sources

The preloaded chemistries — [`lfp_chemistry`](@ref), [`nmc_chemistry`](@ref),
[`nca_chemistry`](@ref), [`lead_acid_chemistry`](@ref), [`leaf_chemistry`](@ref)
— use **stylised OCV shapes anchored to published nominal / min / max cell
voltages**. They capture each chemistry's qualitative OCV(soc) form (LFP's flat
plateau, NMC/NCA's slope, lead-acid's near-linear fall) for sensible defaults;
they are **not** fits to a specific cell's test data. For a calibrated curve,
load real data through [`tabulated_chemistry`](@ref).

**LFP is a deliberate first-class default**: it dominates residential and
stationary storage (Tesla Powerwall, BYD, sonnen), and its flat plateau with
sharp end-knees is the honest stress case — the PE-model's constant-voltage
assumption is least wrong mid-SoC yet the voltage/current limits bite hardest at
the knees.

Open datasets to calibrate from:

| Source | Provides | Notes |
|---|---|---|
| [**PyBaMM** parameter sets](https://docs.pybamm.org/en/latest/source/api/parameters/parameter_sets.html) (BSD-3) | OCV(soc) as cited functions per chemistry | `Chen2020` (NMC811, LG M50), `Prada2013` (LFP), NCA, LCO — the most reusable OCV source |
| [**Battery Archive**](https://batteryarchive.org) | OCV + cycling, standardised | aggregates CALCE, Sandia, Oxford, HNEI, SNL |
| **Sandia National Labs** (via Battery Archive) | OCV + `R(soc, T)` for NMC/NCA/LFP | multiple temperatures/DoDs — for a future temperature hook |
| [**CALCE** (U. Maryland)](https://calce.umd.edu/battery-data) | HPPC → `R(soc)`, OCV | classic equivalent-circuit identification data |
| [**Oxford Battery Degradation Dataset 1**](https://ora.ox.ac.uk/objects/uuid:03ba4b01-cfed-46d3-9b1a-7d4a7bdf6fac) (CC-BY) | Kokam pouch OCV + aging | for a future degradation hook |
| **Nissan-Leaf 2013 cell** ([Zenodo 2580327](https://doi.org/10.5281/zenodo.2580327), CC-BY) | full multi-C-rate `v(soc, i)` surface | the source paper's cell; use to reproduce its results |

Each preloaded chemistry records its provenance in the `source` field:

```julia
lfp_chemistry().source   # "Stylised LFP OCV (datasheet-anchored); PyBaMM Prada2013"
```

## Critique of the source model (and how this differs)

Honest notes on the 2020 paper, and the choices made here:

- **Data availability.** The bivariate `f_v(soc, i)` surface needs multi-C-rate
  data at matched SoC — a dead-end for a chemistry *library*. The OCV + `R`
  decomposition recovers the same physics from data that actually exists, turning
  "one chemistry" into "any chemistry with a published OCV curve".
- **Optimiser vs. simulation.** The paper's spline smoothing factor is chosen
  visually, and cubic extrapolation can produce non-physical voltages (its
  Fig. 2). Harmless when *simulating*; a correctness bug when an *optimiser* can
  park the state where the fit hands out free energy. Hence the monotone,
  clamped OCV here.
- **Converter assumptions.** The paper assumes a stable AC voltage and DC-side-
  only converter losses — which break in an unbalanced four-wire feeder, exactly
  where BMOPFTools operates. Reusing [`AdvancedInverter`](@ref), with its
  explicit internal AC node and current-dependent losses, repairs that.
- **Non-convexity.** `p = v·i` with `v = f(soc)` is bilinear-nonconvex on top of
  the already-nonconvex rectangular power flow; only local optima. The PE-model
  (`StorageDevice`) is the natural convex warm-start — a planned option for the
  multi-period solve.
- **Excluded.** Temperature and degradation (both affect OCV and especially `R`)
  are out of scope for v1; `R(soc)` leaves room to add a `T` dependence later.

## Status and next steps

Delivered: the chemistry library, the single-snapshot device
([`solve_ivq_battery`](@ref)), and multi-period arbitrage
([`solve_multiperiod_ivq`](@ref)) with the trapezoidal/forward charge balance and
per-unit conditioning. Planned:

1. **PE warm-start** for the multi-period solve — seed from
   [`solve_multiperiod_opf`](@ref) to put Ipopt in the right basin on the harder
   (e.g. net-charge-against-price) configurations that currently fail to converge.
2. **Shared `v_dc`**: couple the inverter DC rail to the battery terminal voltage
   (so the three-phase switching-polytope sees the true SoC-dependent DC rail).
3. **Temperature / degradation** hooks on `R(soc)` and capacity.
4. Co-optimise IVQ batteries alongside PE [`StorageDevice`](@ref)s in one horizon.

## API

See the API reference for [`IVQBattery`](@ref), [`solve_ivq_battery`](@ref),
[`IVQBatteryResult`](@ref), [`solve_multiperiod_ivq`](@ref),
[`MultiperiodIVQResult`](@ref), [`BatteryChemistry`](@ref), and the chemistry
constructors [`thevenin_chemistry`](@ref), [`linear_chemistry`](@ref),
[`tabulated_chemistry`](@ref), and the preloaded [`lfp_chemistry`](@ref),
[`nmc_chemistry`](@ref), [`nca_chemistry`](@ref), [`lead_acid_chemistry`](@ref),
[`leaf_chemistry`](@ref).
