# Battery storage models: PE versus IVQ

> **Audience:** power-system researchers · **Scope:** selecting battery models
> for OPF, scheduling, and cell-limited network studies.

Battery models answer different questions. The power-energy (PE) model is often
the right abstraction for a feeder scheduling study; IVQ
(voltage-current-charge) is needed when cell voltage, current, and their
state-dependence determine deliverable power. A richer model with uncalibrated
data is not automatically more defensible than a transparent PE model with
well-supported limits.

This tutorial compares PowerOptLab's PE storage model with the three IVQ
chemistry levels. It explains what each model constrains, what data it needs,
and the conclusions it cannot support.

## 1. The families at a glance

| Model | State | Main decisions | Loss representation | Direct cell limits? | Best use |
|---|---|---|---|---|---|
| PE `StorageDevice` | energy (Wh) | AC charge/discharge power | fixed one-way efficiency | no | scheduling, arbitrage, EV studies, large scenario sets |
| IVQ Thévenin | charge/SoC | cell current, terminal voltage, DC power | constant OCV and Rint | yes | first cell-limited study |
| IVQ linear | charge/SoC | as above | linear OCV(soc) + Rint | yes | endpoint-based SoC sensitivity |
| IVQ tabulated | charge/SoC | as above | OCV(soc), optional R(soc) | yes | calibrated chemistry-specific study |

PE is used through [`solve_multiperiod_opf`](@ref). IVQ provides a fixed-SoC
limit study ([`solve_ivq_battery`](@ref)) and a chronological model
([`solve_multiperiod_ivq`](@ref)). Both operate in the same AC network.

## 2. PE: a contractual power-energy abstraction

[`StorageDevice`](@ref) splits net AC injection into non-negative charge and
discharge powers:

```math
P^{inj}=p^d-p^c,\qquad
E_{t+1}=E_t+\left(\eta^c p^c-\frac{p^d}{\eta^d}\right)\Delta t.
```

```julia
using PowerOptLab

pe = StorageDevice(id="bat", bus="poc",
    p_charge_max=40e3, p_discharge_max=40e3,
    energy_max=100e3, energy_init=40e3,
    eff_charge=0.95, eff_discharge=0.95, cyclic=true)

schedule = solve_multiperiod_opf(nets, [pe]; dt_h=1.0)
```

An energy terminal target can replace `cyclic=true`; `q_min`/`q_max` permit AC
reactive support (unity power factor is the default). Fixed losses normally make
simultaneous charge/discharge suboptimal without a binary complementarity model,
but inspect dispatch under unusual objectives or negative prices.

PE assumes its AC power limits remain valid across the usable energy window.
Cell voltage, current, internal resistance, and load-dependent efficiency are
already aggregated into ratings and efficiencies.

### Pitfall: treating PE energy as electrochemical charge

PE tracks Wh using assumed efficiencies. It neither produces a cell voltage nor
proves compliance with cell current/voltage limits. Do not calibrate `energy_max`
from nominal amp-hour capacity and then claim cell-safe operation without an
external derating calculation that supports the PE limits.

## 3. IVQ: direct cell current, voltage, and charge

IVQ uses a 0th-order equivalent circuit (Rint):

```math
v_{cell}(soc,i)=OCV(soc)-iR(soc),\qquad i>0\ \text{on discharge},
```

with pack power and chronological charge balance

```math
P_{dc}=v_{cell}i\,n_s n_p,\qquad
soc_{t+1}=soc_t-\frac{i_t\Delta t}{q_{cell}n_p}.
```

Cell voltage, charge/discharge current, and SoC bounds are direct constraints.
On discharge voltage falls below OCV; on charge it rises. A cell boundary can
therefore bind before the inverter, or an inverter can bind while the cell has
unused electrochemical headroom.

```julia
chem = illustrative_lfp()  # demonstration only, not a calibrated cell
inv = AdvancedInverter(id="bat", bus="poc", s_max=5e3)
ivq = IVQBattery(id="bat", bus="poc", chemistry=chem,
                 n_series=300, n_parallel=1, soc_init=0.5,
                 inverter=inv, cyclic=true)

limit = solve_ivq_battery(net, ivq; objective=:max_export)
limit.p_poc, limit.v_cell, limit.i_cell
```

The [`AdvancedInverter`](@ref) still owns AC coupling, filter/current limits,
apparent-power topology, and converter losses. IVQ attaches at its DC power
port rather than recreating converter physics.

### Pitfall: calling IVQ an electrochemical model

Rint represents OCV and ohmic sag/rise. It omits polarization, relaxation, rate
capacity, hysteresis, temperature, aging, and cell imbalance. It is more
cell-boundary-aware than PE, not a substitute for an RC or electrochemical
dynamic model.

## 4. The IVQ chemistry ladder

Use the least complex chemistry that answers the question and can be supported
by data.

### Thévenin: constant OCV and R

```julia
chem = thevenin_chemistry(name="screening-cell", v_nominal=3.6,
    r_internal=0.02, q_cell=10.0,
    i_charge_max=30.0, i_discharge_max=60.0)
```

This isolates current limits and voltage sag without SoC-voltage shape. With
`R=0`, it approaches constant-voltage/current behaviour, closest to PE, but it
still tracks charge and directly enforces cell limits.

Use it for first sensitivity and sizing studies. Do not claim a realistic
near-empty/full power envelope or chemistry-specific efficiency.

### Linear OCV plus R

```julia
chem = linear_chemistry(name="endpoint-fit", v_empty=3.0, v_full=3.6,
    r_internal=0.015, q_cell=50.0, soc_min=0.05, soc_max=0.95)
```

Linear OCV gives transparent SoC sensitivity from only endpoint data. It is
better than a flat voltage for exploration, but it erases plateaux and knees.
It is particularly unsuitable for claiming LFP-specific low/high-SoC behaviour.

### Tabulated OCV and optional R(soc)

```julia
chem = tabulated_chemistry(name="measured-cell-v1",
    soc_points=[0.0, 0.1, 0.5, 0.9, 1.0],
    ocv_points=[3.0, 3.2, 3.3, 3.35, 3.6],
    r_points=[0.030, 0.020, 0.012, 0.016, 0.030],
    q_cell=50.0, soc_min=0.05, soc_max=0.95,
    i_charge_max=40.0, i_discharge_max=80.0,
    source="documented OCV/HPPC test, 25 C, fresh cell")
```

Tabulation uses monotone PCHIP interpolation, preventing ordinary cubic-spline
overshoot that could hand an optimizer nonphysical voltage. It is C1 at interior
knots, not C2, and clamps outside outer knots; keep the usable SoC range strictly
inside those knots.

Use it only with explicit OCV, resistance, capacity, limit, temperature, state-
of-health, and test-protocol provenance. The `illustrative_*` presets are
hand-drawn demo shapes, not measured parameter sets.

### Pitfall: increasing fidelity without increasing data quality

Tabulating a guessed curve does not validate IVQ. If only nominal voltage,
rating, and cycle efficiency are known, PE or Thévenin sensitivity bands may be
more honest than a chemistry-labelled tabulated curve.

## 5. Efficiency has different meanings

PE takes `eff_charge` and `eff_discharge` as input; their product is an assumed
round-trip-energy approximation. IVQ derives the Rint proxy

```math
\eta_{proxy}=\frac{OCV-iR}{OCV+iR}
```

at fixed SoC and equal-magnitude current. It worsens with current, but is not a
full-cycle energy efficiency and excludes converter loss, which remains in the
inverter model.

### Pitfall: matching these two efficiencies numerically

They are different validation quantities. Calibrate PE efficiency against the
target cycle/energy throughput. Calibrate IVQ OCV/R against voltage/current
tests, then account for converter loss separately.

## 6. Operating point, scheduling, and numerical implications

Use `solve_ivq_battery` for a fixed-SoC physical-limit question:

```julia
export_limit = solve_ivq_battery(net, ivq; objective=:max_export)
charge_limit = solve_ivq_battery(net, ivq; objective=:max_charge)
```

Use `solve_multiperiod_ivq` for chronological dispatch:

```julia
result = solve_multiperiod_ivq(nets, [ivq]; dt_h=1.0)
d = result.dispatch["bat"]
d.soc, d.i_cell, d.v_cell, d.p_poc
```

The IVQ update conserves charge exactly for piecewise-constant period current.
Match initial/terminal condition, time step, inverter rating/losses, AC Q
policy, usable SoC window, and objective before comparing it with PE.

IVQ introduces `P=v*i` on top of nonlinear AC power flow, so multi-period IVQ
is nonconvex and returns local optima. It uses pack-specific DC bases and
per-unit AC mode by default; raw SI can be poorly conditioned. Non-success
statuses return NaN trajectories and must not be reported as schedules.

The present coupling is at DC power level: set inverter `v_dc` to nominal pack
voltage. The inverter's three-phase switching polytope does not yet see the
instantaneous SoC-dependent IVQ terminal rail.

### Pitfall: comparing unmatched battery studies

Equal kW rating does not make PE and IVQ runs comparable. Align every boundary
condition, repeat nonconvex IVQ solves from credible starts, compare against a
PE baseline, and report binding limits and solver status. A discrepancy may be
an important cell effect, an uncalibrated chemistry assumption, or a local basin.

## 7. A defensible comparison protocol

1. Define the target: energy throughput, AC delivery, peak current, voltage, or economics.
2. Use the same feeder, inverter, Q policy, horizon, and state boundaries.
3. Derive PE limits from the same BMS/test data used for IVQ where possible.
4. Progress IVQ from Thévenin to linear to tabulated only when additional data justify it.
5. Report SoC/current/voltage trajectories, AC/DC losses, binding constraints, and solver status.
6. Test sensitivity to OCV, R, temperature, SoH, and omitted hysteresis/polarization dynamics.

PE is a clear, efficient contractual abstraction; IVQ is a cell-boundary-aware
model. Their comparison matters when it tests whether PE is a sound aggregation
of battery physics for the decision at hand.

For API details see [Storage & EV devices](../components/devices.md) and
[Current-voltage (IVQ) battery](../components/ivq_battery.md).
