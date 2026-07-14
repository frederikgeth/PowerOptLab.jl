# Current–voltage (IVQ) battery

> **Kind:** Component model · **Maturity:** prototype · **Direction:** forward · **Temporal:** single-snapshot (multi-period planned)

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
| [`tabulated_chemistry`](@ref) | clamped piecewise-linear | fixed or `R(soc)` table | an OCV curve (± HPPC) |

The PE-model is the further special case of `thevenin_chemistry` with `R = 0`
(constant voltage ⇒ power ∝ current).

!!! note "Extrapolation is clamped on purpose"
    A free cubic spline overshoots outside its sample hull; in an *optimiser*
    that hands out free voltage — i.e. free energy. `tabulated_chemistry` holds
    OCV **flat** beyond the fitted range and requires a **non-decreasing** OCV
    curve, and every chemistry carries a `soc_min`/`soc_max` usable window the
    device clamps `soc_init` to.

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

This branch delivers the chemistry library and the **single-snapshot** device.
Planned:

1. **Multi-period** IVQ: `soc` as a variable with the charge balance
   ``q_{t+1} = q_t - \Delta t\,\tfrac{i_t + i_{t+1}}{2}`` (trapezoidal — the
   paper's recommended integrator, with forward/backward Euler options), linked
   like [`StorageDevice`](@ref) in [`solve_multiperiod_opf`](@ref), plus a PE
   warm-start. For a linear chemistry this stays polynomial (directly
   embeddable); a tabulated OCV(soc) is registered as a smooth JuMP operator.
2. **Shared `v_dc`**: couple the inverter DC rail to the battery terminal voltage.
3. **Temperature / degradation** hooks on `R(soc)` and capacity.

## API

See the API reference for [`IVQBattery`](@ref), [`solve_ivq_battery`](@ref),
[`IVQBatteryResult`](@ref), [`BatteryChemistry`](@ref), and the chemistry
constructors [`thevenin_chemistry`](@ref), [`linear_chemistry`](@ref),
[`tabulated_chemistry`](@ref), and the preloaded [`lfp_chemistry`](@ref),
[`nmc_chemistry`](@ref), [`nca_chemistry`](@ref), [`lead_acid_chemistry`](@ref),
[`leaf_chemistry`](@ref).
