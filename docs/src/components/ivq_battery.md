# Current–voltage–charge (IVQ) battery

> **Kind:** Component model · **Maturity:** prototype · **Direction:** forward · **Temporal:** single-snapshot + multi-period

[`IVQBattery`](@ref) models grid storage in the **voltage–current–charge**
variable space instead of the power–energy space of [`StorageDevice`](@ref). It
uses a charge-based formulation motivated by

> P. Aaslid, F. Geth, M. Korpås, M. M. Belsnes, O. B. Fosso,
> *Non-linear charge-based battery storage optimization model with bi-variate
> cubic spline constraints*, **Journal of Energy Storage 32 (2020) 101979**.
> [doi:10.1016/j.est.2020.101979](https://doi.org/10.1016/j.est.2020.101979)

Where the power–energy ("PE") model — [`StorageDevice`](@ref) plus the
state-of-charge balance in [`solve_multiperiod_opf`](@ref) — tracks energy with
fixed one-way charging and discharging efficiencies, the IVQ model tracks
**charge** and represents the cell **terminal voltage explicitly**. Voltage and
current limits are constrained individually at the modeled operating points.
Two consequences matter for dispatch:

- **Explicit modeled cell limits.** IVQ constrains the voltage and current of the
  supplied equivalent-circuit model. Whether those bounds describe real equipment
  depends on calibration and omitted physics. The illustrative presets do not
  establish physical safety or a validated operating envelope. PE ratings and
  energy bounds likewise need evidence for the study's AC-to-stored-energy boundary.
- **A current-dependent efficiency proxy is derived, not assumed.**
  ``\eta = f_v(\text{soc}, i)\,/\,f_v(\text{soc}, -i)`` falls out of the voltage
  curve for equal positive current magnitudes at fixed SoC. For the Rint law with
  positive resistance it decreases with current — no cell-level
  `eff_charge`/`eff_discharge` parameters.
  Note this is an *instantaneous, equal-current, fixed-SoC cell* efficiency, **not**
  a full round-trip *energy* efficiency over a cycle, and it excludes the converter
  losses (which live on the [`AdvancedInverter`](@ref)).

## Terminal-voltage model

The paper represents the cell voltage as a bivariate cubic spline
``v = f_v(\text{soc}, i)`` fitted to multi-C-rate discharge data. An empirical
surface requires matched state-of-charge and current data. PowerOptLab uses the
0th-order equivalent-circuit (Thévenin / *Rint*) decomposition

```math
v(\text{soc}, i) = \mathrm{OCV}(\text{soc}) - i \cdot R(\text{soc}),
\qquad i > 0 \ \text{discharge},
```

This is the implemented 0th-order equivalent-circuit approximation, rather than
the paper's empirical spline surface. With positive resistance, terminal voltage
sags under discharge and rises under charge, using the supplied OCV and resistance
curves. Rint omits polarization dynamics, relaxation, rate dependence and hysteresis (which
matter especially for LFP — a hysteresis-aware RC model is the usual next step).
The pack scales the cell by
its series/parallel counts: ``v_\text{pack} = n_\text{series}\,v_\text{cell}``,
``i_\text{pack} = n_\text{parallel}\,i_\text{cell}``, and the DC power is
``p_\text{dc} = v_\text{cell}\,i_\text{cell}\,n_\text{series}\,n_\text{parallel}``.

Three fidelity levels, cheapest first — each a [`BatteryChemistry`](@ref):

| Constructor | OCV(soc) | R | Data needed |
|---|---|---|---|
| [`thevenin_chemistry`](@ref) | constant | fixed | a nominal voltage + a resistance |
| [`linear_chemistry`](@ref) | linear `v_empty → v_full` | fixed | two endpoint voltages + a resistance |
| [`tabulated_chemistry`](@ref) | monotone cubic (PCHIP) | fixed or `R(soc)` table | an OCV curve (± HPPC) |

With constant OCV and `R = 0`, ideal cell power is proportional to current and
stored energy is proportional to charge after choosing a consistent voltage and
capacity scale. Equivalence with a PE energy account additionally requires
consistent AC/DC boundaries, efficiencies and converter constraints. Setting
`R = 0` removes cell resistive loss; the [`AdvancedInverter`](@ref) still has its
own loss and capability model. See [PE versus IVQ](../tutorials/battery_storage_models.md)
for a comparison of the complete storage models.

!!! note "OCV(soc) smoothness — and its limits"
    In multi-period, `OCV(soc)` is evaluated at the SoC *variable*.
    `linear`/`thevenin` OCV is affine (C^∞, embedded directly). `tabulated` OCV
    uses a **monotone cubic (PCHIP)** interpolant, which is **C¹ but not C²** at
    interior knots and only **C⁰** at the two outer knots (where it clamps flat to
    forbid extrapolation). Ipopt's theory assumes twice-continuously-differentiable
    functions, so PCHIP is a pragmatic choice, not a guaranteed-convergent one — in
    practice it works because evaluations rarely land exactly on a knot, and the
    `soc_min`/`soc_max` window is required to sit strictly inside the outer knots so
    the C⁰ points are never in range. A monotone C² alternative (shape-constrained
    smoothing) is the robust upgrade. PCHIP is preferred over a plain cubic spline
    (C² but overshoots, handing an optimiser free voltage) because monotonicity
    matters more than the extra derivative here.

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

chem = illustrative_lfp()                       # or nmc/nca/lead_acid/leaf, or build your own
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
fixed at `soc_init`. This operating-point study computes deliverable power under
the supplied cell and converter limits. Reproducing published numerical results
requires the corresponding cell data, fitted voltage model and network assumptions.

## Multi-period arbitrage

[`solve_multiperiod_ivq`](@ref) co-optimises a chronological sequence of
snapshots with the SoC linking the periods — now `soc` is a *decision variable*,
so `OCV(soc)` enters the model as a function of it (embedded for affine
chemistries; a registered smooth operator for tabulated ones).

```julia
nets = [net_expensive, net_cheap]            # e.g. differing voltage_source `cost`
bat  = IVQBattery(id="bat", bus="poc", chemistry=illustrative_lfp(),
                  n_series=300, n_parallel=1, soc_init=0.5, inverter=inv, cyclic=true)
res  = solve_multiperiod_ivq(nets, [bat]; dt_h=1.0)
res.dispatch["bat"].soc      # SoC trajectory (length T+1)
res.dispatch["bat"].i_cell   # signed cell current per period
```

The charge balance is a forward update ``q_{t+1} = q_t - i_t\,\Delta t``, which is
**exact** for the piecewise-constant current of each period (there is one current
variable per period), so it conserves charge exactly — for a cyclic horizon the
closure ``\sum_t i_t\,\Delta t = 0`` is enforced term-for-term. (A true
trapezoidal rule would need currents at all ``T{+}1`` time *nodes* with consistent
endpoint semantics; the earlier per-interval average silently broke coulomb
conservation and was removed.) Terminal state is set by `cyclic=true` (return to
`soc_init`) or `soc_final`. The objective is the horizon energy cost (each
snapshot's cost rate ``\times\,\Delta t``).

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
engine's per-unit (`per_unit=true`, **the default**) and raw-SI (`per_unit=false`)
modes; results are returned in SI either way. Per-unit additionally conditions
the *AC-network* side of the coupled solve, and the nonconvex multi-period
problem converges far more reliably under it across platforms/Ipopt builds — a
raw-SI solve of the same problem is knife-edge (small solver-option changes flip
it between `LOCALLY_SOLVED` and `LOCALLY_INFEASIBLE`), which is why per-unit is
the default.

## Chemistry presets and data sources

!!! warning "The presets are illustrative, not measured"
    The `illustrative_*` presets — [`illustrative_lfp`](@ref),
    [`illustrative_nmc`](@ref), [`illustrative_nca`](@ref),
    [`illustrative_lead_acid`](@ref), [`illustrative_leaf`](@ref) — are **hand-drawn
    OCV shapes** in a plausible voltage band for each chemistry. They are **not**
    fits to any cell's measured data; they carry no temperature, state-of-health,
    or rest-protocol metadata; and their capacity / current / resistance defaults
    are representative round numbers, **not** tied to the cell that inspired the
    voltage band. Use them for demos, tests and defaults — **not** for scientific
    comparison or operational studies. For a calibrated chemistry, fit real data
    (below) and pass it through [`tabulated_chemistry`](@ref).

LFP is provided first because it is common in residential/stationary storage and
its flat plateau with sharp end-knees is a good stress case — the PE-model's
constant-voltage assumption is least wrong mid-SoC yet the voltage/current limits
bite hardest at the knees.

Open datasets to calibrate a real chemistry from:

| Source | Provides | Notes |
|---|---|---|
| [**PyBaMM** parameter sets](https://docs.pybamm.org/en/latest/source/api/parameters/parameter_sets.html) (BSD-3) | OCV(soc) as cited functions per chemistry | `Chen2020` (NMC811, LG M50), `Prada2013` (LFP), NCA, LCO — the most reusable OCV source |
| [**Battery Archive**](https://batteryarchive.org) | OCV + cycling, standardised | aggregates CALCE, Sandia, Oxford, HNEI, SNL |
| **Sandia National Labs** (via Battery Archive) | OCV + `R(soc, T)` for NMC/NCA/LFP | multiple temperatures/DoDs — for a future temperature hook |
| [**CALCE** (U. Maryland)](https://calce.umd.edu/battery-data) | HPPC → `R(soc)`, OCV | classic equivalent-circuit identification data |
| [**Oxford Battery Degradation Dataset 1**](https://ora.ox.ac.uk/objects/uuid:03ba4b01-cfed-46d3-9b1a-7d4a7bdf6fac) (CC-BY) | Kokam pouch OCV + aging | for a future degradation hook |
| **Nissan-Leaf 2013 cell** ([Zenodo 2580327](https://doi.org/10.5281/zenodo.2580327), CC-BY) | full multi-C-rate `v(soc, i)` surface | the source paper's cell; needed to actually reproduce its results |

Each preset records its (illustrative) provenance in the `source` field:

```julia
illustrative_lfp().source   # "ILLUSTRATIVE hand-drawn LFP OCV shape (not a fit); …"
```

A proper chemistry library — named, versioned cell parameter sets with stored
raw/licensed data, deterministic fitting scripts, separate provenance for OCV /
resistance / capacity / limits, and voltage-error validation over SoC, current,
temperature and charge/discharge direction — is future work, not what these
presets provide.

## Critique of the source model (and how this differs)

Honest notes on the 2020 paper, and the choices made here:

- **Data availability.** The bivariate `f_v(soc, i)` surface needs multi-C-rate
  data at matched SoC — a dead-end for a chemistry *library*. The OCV + `R`
  decomposition uses OCV and resistance data in a simpler equivalent-circuit
  approximation. It does not reproduce the empirical surface's rate dependence,
  polarization or hysteresis; an OCV curve alone does not calibrate those effects.
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
  (`StorageDevice`) has a simpler affine energy balance, but its coupled AC
  scheduling problem remains nonconvex. Automatic transfer of a PE dispatch into
  an IVQ warm start is not implemented.
- **Excluded.** Temperature and degradation (both affect OCV and especially `R`)
  are out of scope for v1; `R(soc)` leaves room to add a `T` dependence later.

## Status and next steps

Delivered: the chemistry library, the single-snapshot device
([`solve_ivq_battery`](@ref)), and multi-period arbitrage
([`solve_multiperiod_ivq`](@ref)) with the exact forward charge balance and
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
[`tabulated_chemistry`](@ref), and the preloaded [`illustrative_lfp`](@ref),
[`illustrative_nmc`](@ref), [`illustrative_nca`](@ref), [`illustrative_lead_acid`](@ref),
[`illustrative_leaf`](@ref).
