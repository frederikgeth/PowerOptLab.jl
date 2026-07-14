# Battery chemistry models for the current–voltage (IVQ) storage device.
#
# A `BatteryChemistry` describes ONE cell in the voltage–current–charge variable
# space of Aaslid, Geth, Korpås, Belsnes & Fosso, "Non-linear charge-based
# battery storage optimization model with bi-variate cubic spline constraints",
# Journal of Energy Storage 32 (2020) 101979 (https://doi.org/10.1016/j.est.2020.101979).
# In that model the cell terminal voltage is an empirical function of state of
# charge and current, `v = f_v(soc, i)`, and the state is counted as charge (Ah)
# rather than energy (Wh) — so voltage and current limits are represented
# individually instead of being folded into a single power/energy limit.
#
# We deliberately do NOT reproduce the paper's raw bivariate cubic spline as the
# primary interface: a full `f_v(soc, i)` surface needs multi-C-rate discharge
# data at matched SoC, which is published for almost no chemistry. Instead the
# terminal voltage is decomposed into an equivalent-circuit (Thévenin / Rint)
# form
#
#     v(soc, i) = OCV(soc) − i · R(soc)            (i > 0 discharge)
#
# where OCV(soc) is the open-circuit voltage (abundant and citable for every
# chemistry) and R(soc) the internal resistance (from HPPC pulse data). This
# recovers the same physics the paper's surface encodes — including the
# round-trip efficiency emerging as η = f_v(soc, i)/f_v(soc, −i) rather than
# being a fixed parameter — while only requiring data that actually exists.
#
# Three levels of fidelity, cheapest first:
#   • `thevenin_chemistry` — constant OCV behind a fixed resistance (the paper's
#     PE-model is the further special case R → 0).
#   • `linear_chemistry`   — OCV drops linearly between a full and an empty
#     voltage as SoC falls (a "state-dependent voltage source"), fixed R.
#   • `tabulated_chemistry`— OCV(soc) (and optionally R(soc)) from data points,
#     e.g. a published OCV curve or a slow-rate discharge test.
# Preloaded shapes for common chemistries (LFP, NMC, NCA, lead-acid, and the
# paper's Nissan-Leaf cell) build on these; see their docstrings for provenance.
#
# All quantities are SI and per cell: volts, amps (signed, i > 0 = discharge),
# ohms, amp-hours. Cell → pack scaling (series/parallel counts) lives on the
# `IVQBattery` device, not here.

"""
    BatteryChemistry

One battery cell in the voltage–current–charge (IVQ) variable space. The cell
terminal voltage follows the Thévenin/Rint decomposition
`v(soc, i) = OCV(soc) − i · R(soc)` with `i > 0` on discharge, so charging
(`i < 0`) raises the terminal voltage above OCV and discharging lowers it.

Construct one with [`thevenin_chemistry`](@ref), [`linear_chemistry`](@ref) or
[`tabulated_chemistry`](@ref), or use a preloaded shape ([`lfp_chemistry`](@ref),
[`nmc_chemistry`](@ref), [`nca_chemistry`](@ref), [`lead_acid_chemistry`](@ref),
[`leaf_chemistry`](@ref)).

# Fields
- `name::String` — chemistry label.
- `ocv::Function` — `soc ∈ [0,1] → open-circuit voltage (V)`; non-decreasing and
  smooth (C¹). `linear`/`thevenin` OCV is affine; `tabulated` OCV is a monotone
  cubic (PCHIP) so it is safe to embed as a function of the SoC *variable*.
- `r_internal::Function` — `soc ∈ [0,1] → internal resistance (Ω) ≥ 0`.
- `ocv_affine::Union{Tuple{Float64,Float64},Nothing}` — `(intercept, slope)` when
  OCV is affine (so the multi-period model embeds it as a plain expression rather
  than a registered operator), else `nothing`.
- `r_constant::Union{Float64,Nothing}` — the resistance value when `R` is constant,
  else `nothing`.
- `q_cell::Float64` — cell capacity (Ah).
- `v_cell_min`, `v_cell_max::Float64` — terminal-voltage operating bounds (V).
- `i_charge_max`, `i_discharge_max::Float64` — current magnitude limits (A ≥ 0).
- `soc_min`, `soc_max::Float64` — usable SoC window; the device clamps to it so
  the optimiser cannot exploit extrapolation outside the fitted range.
- `source::String` — data provenance / reference for the OCV and R values.
"""
struct BatteryChemistry
    name::String
    ocv::Function
    r_internal::Function
    ocv_affine::Union{Tuple{Float64,Float64},Nothing}
    r_constant::Union{Float64,Nothing}
    q_cell::Float64
    v_cell_min::Float64
    v_cell_max::Float64
    i_charge_max::Float64
    i_discharge_max::Float64
    soc_min::Float64
    soc_max::Float64
    source::String
end

# Monotone cubic (Fritsch–Carlson / PCHIP) interpolant of `(xs, ys)`, returned as
# a closure. It is **smooth** (C¹, continuous first derivative — what an
# interior-point solver needs) and **shape-preserving** (no overshoot, so a
# monotone OCV curve stays monotone), and is held flat outside `[xs[1], xs[end]]`
# so the optimiser cannot extrapolate into non-physical voltage. The closure is
# generic in its argument (accepts `ForwardDiff.Dual`), so it registers directly
# as a JuMP nonlinear operator with automatic-differentiation derivatives.
#
# A plain cubic spline would be C² but can overshoot — the very defect (free
# voltage ⇒ free energy) we criticise in the source paper's raw splines — so
# monotonicity is preferred over the extra derivative.
function _monotone_cubic(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs)
    h = [xs[i+1] - xs[i] for i in 1:n-1]
    Δ = [(ys[i+1] - ys[i]) / h[i] for i in 1:n-1]
    d = Vector{Float64}(undef, n)
    d[1] = Δ[1]; d[n] = Δ[n-1]                       # shape-preserving endpoints
    for i in 2:n-1
        if Δ[i-1] * Δ[i] <= 0
            d[i] = 0.0                                # local extremum ⇒ flat tangent
        else
            w1 = 2h[i] + h[i-1]; w2 = h[i] + 2h[i-1]  # Fritsch–Carlson weights
            d[i] = (w1 + w2) / (w1 / Δ[i-1] + w2 / Δ[i])
        end
    end
    return function (x)
        x <= xs[1] && return ys[1] + zero(x)
        x >= xs[n] && return ys[n] + zero(x)
        k = 1
        @inbounds while k < n - 1 && x > xs[k+1]
            k += 1
        end
        @inbounds begin
            t = (x - xs[k]) / h[k]
            t2 = t * t; t3 = t2 * t
            h00 =  2t3 - 3t2 + 1
            h10 =  t3 - 2t2 + t
            h01 = -2t3 + 3t2
            h11 =  t3 - t2
            return h00 * ys[k] + h10 * h[k] * d[k] + h01 * ys[k+1] + h11 * h[k] * d[k+1]
        end
    end
end

"""
    thevenin_chemistry(; name, v_nominal, r_internal, q_cell, kwargs...)

A Thévenin cell: a fixed open-circuit voltage `v_nominal` behind a constant
internal resistance `r_internal` (Ω). `v(soc, i) = v_nominal − i·R`. The paper's
PE-model is the further special case `r_internal = 0` (constant voltage, so
power and current are proportional and the current limit reduces to a power
limit). Useful as a floor-fidelity model and as a warm-start generator.

# Keywords
- `name="Thevenin"`, `v_nominal`, `r_internal`, `q_cell` (Ah).
- `v_cell_min = 0.8·v_nominal`, `v_cell_max = 1.2·v_nominal` — voltage bounds (V).
- `i_charge_max = q_cell`, `i_discharge_max = q_cell` — current limits (A; default 1C).
- `soc_min = 0.0`, `soc_max = 1.0`.
- `source="Thévenin / Rint model"`.
"""
function thevenin_chemistry(; name::String="Thevenin",
                            v_nominal::Float64, r_internal::Float64,
                            q_cell::Float64,
                            v_cell_min::Float64 = 0.8v_nominal,
                            v_cell_max::Float64 = 1.2v_nominal,
                            i_charge_max::Float64 = q_cell,
                            i_discharge_max::Float64 = q_cell,
                            soc_min::Float64 = 0.0, soc_max::Float64 = 1.0,
                            source::String = "Thévenin / Rint model")
    r_internal >= 0 || throw(ArgumentError("r_internal must be ≥ 0"))
    return BatteryChemistry(name, _ -> v_nominal, _ -> r_internal,
                            (v_nominal, 0.0), r_internal, q_cell,
                            v_cell_min, v_cell_max, i_charge_max, i_discharge_max,
                            soc_min, soc_max, source)
end

"""
    linear_chemistry(; name, v_full, v_empty, r_internal, q_cell, kwargs...)

A state-dependent voltage source: the open-circuit voltage falls linearly from
`v_full` at `soc = 1` to `v_empty` at `soc = 0`, behind a constant internal
resistance. `OCV(soc) = v_empty + (v_full − v_empty)·soc`. A good first model
whenever only the charged/discharged voltage endpoints and a resistance are
known, and the natural discretisation of a real OCV curve.

# Keywords
- `name="Linear"`, `v_full`, `v_empty` (V, with `v_full ≥ v_empty`), `r_internal`
  (Ω), `q_cell` (Ah).
- `v_cell_min = 0.95·v_empty`, `v_cell_max = 1.05·v_full` — voltage bounds (V).
- `i_charge_max`, `i_discharge_max`, `soc_min`, `soc_max`, `source` — as
  [`thevenin_chemistry`](@ref).
"""
function linear_chemistry(; name::String="Linear",
                          v_full::Float64, v_empty::Float64,
                          r_internal::Float64, q_cell::Float64,
                          v_cell_min::Float64 = 0.95v_empty,
                          v_cell_max::Float64 = 1.05v_full,
                          i_charge_max::Float64 = q_cell,
                          i_discharge_max::Float64 = q_cell,
                          soc_min::Float64 = 0.0, soc_max::Float64 = 1.0,
                          source::String = "linear OCV(soc) + fixed R")
    v_full >= v_empty || throw(ArgumentError("v_full must be ≥ v_empty"))
    r_internal >= 0 || throw(ArgumentError("r_internal must be ≥ 0"))
    ocv = soc -> v_empty + (v_full - v_empty) * soc
    return BatteryChemistry(name, ocv, _ -> r_internal,
                            (v_empty, v_full - v_empty), r_internal, q_cell,
                            v_cell_min, v_cell_max, i_charge_max, i_discharge_max,
                            soc_min, soc_max, source)
end

"""
    tabulated_chemistry(; name, soc_points, ocv_points, r_internal, q_cell, kwargs...)

Open-circuit voltage (and optionally resistance) from data points — e.g. a
published OCV curve or a slow-rate (≈ C/20) discharge test. `OCV(soc)` is the
clamped piecewise-linear interpolant of `(soc_points, ocv_points)`, which must
be non-decreasing (a physical OCV curve); outside the range it is held flat so
the optimiser cannot extrapolate into non-physical voltage. `r_internal` may be
a scalar Ω or a matching vector `r_points` interpolated the same way.

# Keywords
- `name`, `soc_points::Vector`, `ocv_points::Vector` (non-decreasing, V).
- `r_internal` — scalar Ω, **or** pass `r_points::Vector` (Ω) for an R(soc) table.
- `q_cell` (Ah), `v_cell_min = minimum(ocv_points)`, `v_cell_max = maximum(ocv_points)`,
  `i_charge_max`, `i_discharge_max`, `soc_min = minimum(soc_points)`,
  `soc_max = maximum(soc_points)`, `source`.
"""
function tabulated_chemistry(; name::String="Tabulated",
                             soc_points::Vector{Float64},
                             ocv_points::Vector{Float64},
                             r_internal::Union{Float64,Nothing} = nothing,
                             r_points::Union{Vector{Float64},Nothing} = nothing,
                             q_cell::Float64,
                             v_cell_min::Float64 = minimum(ocv_points),
                             v_cell_max::Float64 = maximum(ocv_points),
                             i_charge_max::Float64 = q_cell,
                             i_discharge_max::Float64 = q_cell,
                             soc_min::Float64 = minimum(soc_points),
                             soc_max::Float64 = maximum(soc_points),
                             source::String = "tabulated OCV(soc)")
    length(soc_points) == length(ocv_points) ||
        throw(ArgumentError("soc_points and ocv_points must have equal length"))
    length(soc_points) >= 2 || throw(ArgumentError("need at least two data points"))
    issorted(soc_points) || throw(ArgumentError("soc_points must be sorted ascending"))
    issorted(ocv_points) ||
        throw(ArgumentError("ocv_points must be non-decreasing (a physical OCV curve)"))

    ocv = _monotone_cubic(soc_points, ocv_points)   # smooth (C¹), monotone

    r_const = nothing
    rfun = if r_points !== nothing
        length(r_points) == length(soc_points) ||
            throw(ArgumentError("r_points must match soc_points in length"))
        all(>=(0), r_points) || throw(ArgumentError("r_points must be ≥ 0"))
        _monotone_cubic(soc_points, r_points)
    elseif r_internal !== nothing
        r_internal >= 0 || throw(ArgumentError("r_internal must be ≥ 0"))
        r_const = r_internal
        _ -> r_internal
    else
        throw(ArgumentError("provide either r_internal (scalar Ω) or r_points (Vector Ω)"))
    end

    return BatteryChemistry(name, ocv, rfun, nothing, r_const, q_cell,
                            v_cell_min, v_cell_max, i_charge_max, i_discharge_max,
                            soc_min, soc_max, source)
end

# ── Preloaded chemistry shapes ───────────────────────────────────────────────
#
# Each returns a `BatteryChemistry` for a single cell. The OCV shapes below are
# STYLISED curves anchored to published nominal / min / max cell voltages — they
# capture the qualitative OCV(soc) form of each chemistry (LFP's flat plateau,
# NMC/NCA's slope, lead-acid's near-linear fall) for illustration and defaults.
# They are NOT fits to a specific cell's test data. For a calibrated curve, load
# real data through `tabulated_chemistry` — see the module docs for open
# datasets (PyBaMM parameter sets, Sandia/CALCE/Oxford via batteryarchive.org,
# and the Nissan-Leaf cell dataset used by the source paper).

"""
    lfp_chemistry(; q_cell=100.0, r_internal=0.006, kwargs...)

Lithium iron phosphate (LFP/LiFePO₄): a flat ~3.2–3.3 V plateau between sharp
knees, the dominant chemistry in stationary and residential storage. The flat
OCV is where a constant-voltage PE-model is least wrong in mid-SoC yet where the
terminal-voltage and current limits still bite hardest at the knees — a good
honest stress case for the IVQ model.

Stylised OCV anchored to typical LFP cell voltages (charged ≈ 3.65 V, nominal
≈ 3.2 V, knee ≈ 2.5 V). Source: LFP cell datasheets (e.g. A123 APR18650,
CATL/BYD prismatic); PyBaMM `Prada2013` LFP parameter set.
"""
function lfp_chemistry(; q_cell::Float64=100.0, r_internal::Float64=0.006,
                       i_charge_max::Float64=q_cell, i_discharge_max::Float64=q_cell,
                       soc_min::Float64=0.05, soc_max::Float64=0.98)
    soc = [0.0, 0.03, 0.08, 0.15, 0.30, 0.50, 0.70, 0.85, 0.92, 0.97, 1.0]
    ocv = [2.50, 3.00, 3.20, 3.25, 3.27, 3.29, 3.31, 3.33, 3.35, 3.45, 3.65]
    return tabulated_chemistry(; name="LFP", soc_points=soc, ocv_points=ocv,
        r_internal=r_internal, q_cell=q_cell,
        v_cell_min=2.5, v_cell_max=3.65,
        i_charge_max=i_charge_max, i_discharge_max=i_discharge_max,
        soc_min=soc_min, soc_max=soc_max,
        source="Stylised LFP OCV (datasheet-anchored); PyBaMM Prada2013")
end

"""
    nmc_chemistry(; q_cell=5.0, r_internal=0.03, kwargs...)

Nickel-manganese-cobalt (NMC, ~811): a monotonic slope from ≈ 3.0 V to 4.2 V,
the workhorse of EV cells. Stylised OCV anchored to typical NMC811 cell voltages.
Source: LG INR21700-M50 (NMC811/graphite-SiOx) datasheet; PyBaMM `Chen2020`
parameter set (Chen et al., J. Electrochem. Soc. 167 (2020) 080534).
"""
function nmc_chemistry(; q_cell::Float64=5.0, r_internal::Float64=0.03,
                       i_charge_max::Float64=q_cell, i_discharge_max::Float64=2q_cell,
                       soc_min::Float64=0.05, soc_max::Float64=0.98)
    soc = [0.0, 0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95, 1.0]
    ocv = [3.00, 3.35, 3.55, 3.68, 3.80, 3.95, 4.08, 4.15, 4.20]
    return tabulated_chemistry(; name="NMC811", soc_points=soc, ocv_points=ocv,
        r_internal=r_internal, q_cell=q_cell,
        v_cell_min=3.0, v_cell_max=4.2,
        i_charge_max=i_charge_max, i_discharge_max=i_discharge_max,
        soc_min=soc_min, soc_max=soc_max,
        source="Stylised NMC811 OCV (datasheet-anchored); PyBaMM Chen2020")
end

"""
    nca_chemistry(; q_cell=3.2, r_internal=0.035, kwargs...)

Nickel-cobalt-aluminium (NCA): similar sloped profile to NMC, ≈ 2.5–4.2 V, used
in Panasonic/Tesla cylindrical cells. Stylised OCV anchored to typical NCA cell
voltages. Source: Panasonic NCR18650B datasheet.
"""
function nca_chemistry(; q_cell::Float64=3.2, r_internal::Float64=0.035,
                       i_charge_max::Float64=0.5q_cell, i_discharge_max::Float64=2q_cell,
                       soc_min::Float64=0.05, soc_max::Float64=0.98)
    soc = [0.0, 0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95, 1.0]
    ocv = [2.50, 3.30, 3.52, 3.65, 3.78, 3.95, 4.08, 4.16, 4.20]
    return tabulated_chemistry(; name="NCA", soc_points=soc, ocv_points=ocv,
        r_internal=r_internal, q_cell=q_cell,
        v_cell_min=2.5, v_cell_max=4.2,
        i_charge_max=i_charge_max, i_discharge_max=i_discharge_max,
        soc_min=soc_min, soc_max=soc_max,
        source="Stylised NCA OCV (datasheet-anchored); Panasonic NCR18650B")
end

"""
    lead_acid_chemistry(; q_cell=100.0, r_internal=0.004, kwargs...)

Lead-acid: a roughly linear open-circuit voltage from ≈ 1.75 V (empty) to
≈ 2.15 V (full) per cell — modelled here with [`linear_chemistry`](@ref). Still
common in off-grid and backup storage. Source: standard flooded/AGM lead-acid
OCV-SoC relation (e.g. IEEE 1188 / manufacturer curves).
"""
function lead_acid_chemistry(; q_cell::Float64=100.0, r_internal::Float64=0.004,
                             i_charge_max::Float64=0.2q_cell, i_discharge_max::Float64=q_cell,
                             soc_min::Float64=0.2, soc_max::Float64=1.0)
    return linear_chemistry(; name="Lead-acid", v_full=2.15, v_empty=1.75,
        r_internal=r_internal, q_cell=q_cell,
        v_cell_min=1.70, v_cell_max=2.40,
        i_charge_max=i_charge_max, i_discharge_max=i_discharge_max,
        soc_min=soc_min, soc_max=soc_max,
        source="Linear lead-acid OCV (1.75–2.15 V/cell); IEEE 1188 / datasheets")
end

"""
    leaf_chemistry(; kwargs...)

The 2013 Nissan-Leaf cell used by the source paper (Aaslid et al., 2020) — an
LMO/NMC blend cycled 3.0–4.2 V. Terminal voltage bounds and the ≈ 3.9 V average
match the paper's Table 2 (`Vmin/Vmax/Vavg = 3.20/4.15/3.90 V`, `Qmax = 29 Ah`,
`Ib,ch/Ib,dch = 30/90 A`); the OCV shape is a stylised monotone curve through
that band. Provided for reproducing the paper; for the true `f_v(soc, i)`
surface load the cell dataset (Wiggins, Allu & Wang, ORNL, 2020,
https://doi.org/10.5281/zenodo.2580327) via [`tabulated_chemistry`](@ref).
"""
function leaf_chemistry(; q_cell::Float64=29.0, r_internal::Float64=0.0035,
                        i_charge_max::Float64=30.0, i_discharge_max::Float64=90.0,
                        soc_min::Float64=0.05, soc_max::Float64=0.98)
    soc = [0.0, 0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95, 1.0]
    ocv = [3.20, 3.45, 3.60, 3.72, 3.85, 3.98, 4.06, 4.12, 4.15]
    return tabulated_chemistry(; name="Nissan-Leaf-2013", soc_points=soc, ocv_points=ocv,
        r_internal=r_internal, q_cell=q_cell,
        v_cell_min=3.20, v_cell_max=4.15,
        i_charge_max=i_charge_max, i_discharge_max=i_discharge_max,
        soc_min=soc_min, soc_max=soc_max,
        source="Aaslid et al. 2020 (JES 32:101979), Table 2; cell data Zenodo 2580327")
end
