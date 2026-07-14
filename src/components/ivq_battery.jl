# Current–voltage (IVQ) battery storage device.
#
# The IVQ counterpart to the energy/power `StorageDevice` (`components/devices.jl`),
# implementing the voltage–current–charge storage model of Aaslid, Geth, Korpås,
# Belsnes & Fosso, "Non-linear charge-based battery storage optimization model
# with bi-variate cubic spline constraints", Journal of Energy Storage 32 (2020)
# 101979. Where the PE-model tracks energy with a fixed round-trip efficiency,
# the IVQ device tracks charge and models the cell terminal voltage explicitly:
#
#     v(soc, i) = OCV(soc) − i · R(soc)            (i > 0 discharge, per cell)
#
# so voltage and current limits are enforced individually — enabling safe
# operation right up to the cell boundaries that a power-only model can only
# approximate with conservative padding — and the round-trip efficiency
# η = f_v(soc, i)/f_v(soc, −i) emerges from the physics instead of being a
# parameter (see [`BatteryChemistry`](@ref)).
#
# The AC-side converter is NOT re-implemented here: an [`AdvancedInverter`](@ref)
# owns the point-of-connection coupling, output filter, converter losses and
# apparent-power/topology limits, exactly as for a standalone inverter. The
# battery attaches to that inverter's DC port. In this first version the two are
# coupled at the DC power level,
#
#     v_cell · i_cell · n_series · n_parallel  =  P_dc (inverter)      [W],
#
# which is exact for the (default) single-phase inverter, whose DC-link voltage
# `v_dc` only enters an optional modulation cap. Making `v_dc` a shared decision
# variable equal to the battery terminal voltage (so the switching-polytope of
# the three-phase topologies sees the true, SoC-dependent DC rail) is a planned
# follow-up; set the inverter's `v_dc` to the nominal pack voltage until then.
#
# Signs and units: cell current `i` is signed with i > 0 on discharge (matching
# the paper); `p_dc > 0` means the pack delivers power to the DC link. Device
# parameters are SI; the inverter handles per-unit conditioning through
# `ctx.bases`.

"""
    IVQBattery(; id, bus, chemistry, n_series, n_parallel, soc_init, inverter, kwargs...)

A current–voltage battery: `n_series × n_parallel` cells of a
[`BatteryChemistry`](@ref), coupled to the grid through an
[`AdvancedInverter`](@ref) at the DC port. The pack terminal voltage is
`n_series · v_cell` and the pack current `n_parallel · i_cell`.

# Required
- `id::String`, `bus::String` — identifier and point-of-connection bus (must
  match the `inverter`'s `bus`).
- `chemistry::BatteryChemistry` — the cell model.
- `n_series::Int`, `n_parallel::Int` — pack configuration.
- `soc_init::Float64` — state of charge (0–1) at the operating point / horizon start.
- `inverter::AdvancedInverter` — the AC↔DC converter (reused as-is).

# Optional (reserved for multi-period use; ignored in the single-snapshot solve)
- `cyclic::Bool=true` — require terminal SoC to return to `soc_init`.
- `soc_final::Union{Float64,Nothing}=nothing` — pin terminal SoC (overrides `cyclic`).
- `integration::Symbol=:trapezoidal` — charge-balance rule (`:forward`,
  `:backward`, `:trapezoidal`) for the SoC ODE `dq/dt = −i`.
"""
Base.@kwdef struct IVQBattery
    id::String
    bus::String
    chemistry::BatteryChemistry
    n_series::Int
    n_parallel::Int
    soc_init::Float64
    inverter::AdvancedInverter
    cyclic::Bool = true
    soc_final::Union{Float64,Nothing} = nothing
    integration::Symbol = :trapezoidal
end

# Handle published by `_stamp_battery!` for reporting. All SI.
struct _BatHandles
    i_cell        # signed cell current variable (A; > 0 discharge)
    v_cell        # cell terminal-voltage expression (V)
    p_dc          # pack DC power expression (W; > 0 discharge)
    soc           # state of charge at the operating point (Float64, fixed here)
    n_series::Int
    n_parallel::Int
end

"""
    _stamp_battery!(ctx, battery, inv_handles) -> _BatHandles

Stamp the battery cell physics into `ctx` and couple its DC power to the
inverter's DC link (`inv_handles` from `_stamp_inverter!`). SoC is fixed at
`battery.soc_init` (single snapshot): the terminal voltage is affine in the cell
current and the pack DC power is quadratic in it, so no nonlinear-operator
registration is needed at this fidelity.
"""
function _stamp_battery!(ctx, battery::IVQBattery, inv_handles)
    m = ctx.model
    sb = inv_handles.sb                       # VA base (1.0 in SI mode)
    chem = battery.chemistry
    ns = battery.n_series; np = battery.n_parallel
    soc = battery.soc_init

    chem.soc_min <= soc <= chem.soc_max || throw(ArgumentError(
        "IVQBattery '$(battery.id)': soc_init=$soc outside the chemistry's usable " *
        "window [$(chem.soc_min), $(chem.soc_max)]"))

    ocv = chem.ocv(soc)                        # V, evaluated at the fixed SoC
    r   = chem.r_internal(soc)                 # Ω

    # Signed cell current (A): i > 0 discharge, i < 0 charge.
    i = JuMP.@variable(m, base_name = "i_cell_$(battery.id)")
    JuMP.@constraint(m, i >= -chem.i_charge_max)
    JuMP.@constraint(m, i <=  chem.i_discharge_max)

    # Cell terminal voltage v = OCV − i·R (drops on discharge, rises on charge),
    # bounded by the cell's operating window.
    v = JuMP.@expression(m, ocv - i * r)
    JuMP.@constraint(m, v >= chem.v_cell_min)
    JuMP.@constraint(m, v <= chem.v_cell_max)

    # Pack DC power (SI) = v_pack·i_pack = v·i·n_series·n_parallel, coupled to the
    # inverter DC-link power (converted from model units back to SI).
    p_dc = JuMP.@expression(m, v * i * ns * np)
    JuMP.@constraint(m, p_dc == inv_handles.p_dc * sb)

    return _BatHandles(i, v, p_dc, soc, ns, np)
end

"""
    IVQBatteryResult

Result of [`solve_ivq_battery`](@ref). SI units throughout.

# Fields
- `termination_status::String`
- `p_poc`, `q_poc` — active/reactive power at the grid POC (from the inverter).
- `p_conv` — converter-side active power at the internal node (W).
- `p_dc` — DC-link power (W; = pack discharge power, `p_dc = p_conv + p_loss`).
- `p_loss` — converter loss (W).
- `soc` — state of charge at the operating point (0–1).
- `v_cell`, `i_cell` — cell terminal voltage (V) and signed current (A; > 0 discharge).
- `v_pack`, `i_pack` — pack terminal voltage (V) and current (A).
- `bus::Dict{String,Any}` — the BMOPFTools `result["bus"]`.
"""
struct IVQBatteryResult
    termination_status::String
    p_poc::Float64
    q_poc::Float64
    p_conv::Float64
    p_dc::Float64
    p_loss::Float64
    soc::Float64
    v_cell::Float64
    i_cell::Float64
    v_pack::Float64
    i_pack::Float64
    bus::Dict{String,Any}
end

"""
    solve_ivq_battery(net, battery; objective=:max_export, kwargs...) -> IVQBatteryResult

Stamp `battery` (and its [`AdvancedInverter`](@ref)) into `net` and solve at a
single operating point (SoC fixed at `battery.soc_init`), demonstrating the
coupled cell + converter feasible region. The network supplies the surrounding
grid (a voltage source, lines); the battery exchanges power at its POC bus
through the inverter.

# Objective
- `:max_export` — maximise active power delivered to the grid (battery discharges
  until a cell voltage/current limit or the converter rating binds).
- `:max_charge` — maximise power drawn from the grid into the battery.
- `:min_loss` — minimise converter loss subject to a `p_set` (W) POC delivery.

# Keywords
- `p_set=nothing` — required active-power target (W) for `:min_loss`.
- `q_set=nothing` — optional reactive-power constraint at the POC (var).
- `per_unit=false`, `s_base=1e6`, `optimizer=Ipopt.Optimizer`, `verbose=false`,
  `solver_options=()`. Results are returned in SI regardless of `per_unit`.
"""
function solve_ivq_battery(net::Dict{String,Any}, battery::IVQBattery;
                           objective::Symbol=:max_export,
                           p_set::Union{Float64,Nothing}=nothing,
                           q_set::Union{Float64,Nothing}=nothing,
                           per_unit::Bool=false,
                           s_base::Float64=1e6,
                           optimizer=Ipopt.Optimizer,
                           verbose::Bool=false,
                           solver_options=())
    objective in (:max_export, :max_charge, :min_loss) ||
        throw(ArgumentError("objective must be :max_export, :max_charge or :min_loss, got :$objective"))
    objective == :min_loss && p_set === nothing &&
        throw(ArgumentError(":min_loss requires a p_set (W) active-power target"))
    battery.inverter.bus == battery.bus ||
        throw(ArgumentError("battery bus '$(battery.bus)' must match its inverter bus '$(battery.inverter.bus)'"))
    battery.integration in (:forward, :backward, :trapezoidal) ||
        throw(ArgumentError("integration must be :forward, :backward or :trapezoidal"))

    inv_h = Ref{Any}(); bat_h = Ref{Any}()
    hook! = ctx -> begin
        ih = _stamp_inverter!(ctx, battery.inverter)
        bh = _stamp_battery!(ctx, battery, ih)
        inv_h[] = ih; bat_h[] = bh
        q_set === nothing || JuMP.@constraint(ctx.model, ih.q_poc == q_set / ih.sb)
        if objective == :max_export
            JuMP.@objective(ctx.model, Max, ih.p_poc)
        elseif objective == :max_charge
            JuMP.@objective(ctx.model, Min, ih.p_poc)
        else
            JuMP.@constraint(ctx.model, ih.p_poc == p_set / ih.sb)
            JuMP.@objective(ctx.model, Min, ih.p_loss)
        end
    end

    ctx = build_opf_model(net; per_unit=per_unit, s_base=s_base,
                          add_objective=false, model_hook! = hook!,
                          optimizer=optimizer, verbose=verbose)
    for (name, value) in solver_options
        JuMP.set_attribute(ctx.model, string(name), value)
    end
    enforce_kcl!(ctx)
    JuMP.optimize!(ctx.model)

    status = string(JuMP.termination_status(ctx.model))
    solved = JuMP.primal_status(ctx.model) == JuMP.MOI.FEASIBLE_POINT

    ih = inv_h[]; bh = bat_h[]
    sb = ih.sb
    val(e) = solved ? JuMP.value(e) : NaN
    i_cell = val(bh.i_cell)
    v_cell = val(bh.v_cell)

    result = extract_result(ctx)
    return IVQBatteryResult(status,
                            solved ? JuMP.value(ih.p_poc) * sb : NaN,
                            solved ? JuMP.value(ih.q_poc) * sb : NaN,
                            solved ? JuMP.value(ih.p_conv) * sb : NaN,
                            solved ? JuMP.value(bh.p_dc) : NaN,
                            solved ? JuMP.value(ih.p_loss) * sb : NaN,
                            bh.soc,
                            v_cell, i_cell,
                            v_cell * bh.n_series, i_cell * bh.n_parallel,
                            result["bus"])
end
