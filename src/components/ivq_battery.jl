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

# DC per-unit bases tied to the battery's OWN rating — pack nominal voltage and
# pack max discharge current — so the pack current/voltage variables are ≈ O(1)
# independent of the engine's system `s_base`. This is what keeps the coupled
# battery+inverter solve well conditioned (the battery's ~3.5 V cell scale and
# the network's ~kV/~kW scales otherwise span 5 orders of magnitude and make
# Ipopt return infeasible points). Works in both SI and per-unit engine modes:
# the DC↔AC coupling carries the s_base factor explicitly.
function _dc_bases(chem::BatteryChemistry, ns::Int, np::Int)
    vbase = ns * chem.ocv(0.5 * (chem.soc_min + chem.soc_max))   # pack nominal V
    ibase = np * chem.i_discharge_max                            # pack max discharge A
    return vbase, ibase
end

# Handle published by `_stamp_battery!` for reporting. `ipu`/`vpu` are per-unit
# on the battery's own DC bases (`vbase`, `ibase`); multiply back for SI.
struct _BatHandles
    ipu           # signed per-unit pack current variable (> 0 discharge)
    vpu           # per-unit pack terminal-voltage expression
    vbase::Float64
    ibase::Float64
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

    # The battery is modelled in PACK quantities, per-unit on its own DC bases
    # (`_dc_bases`), so the current/voltage variables are ≈ O(1) and the coupled
    # solve is well conditioned. Cell relations: V_pack = n_series·v_cell,
    # I_pack = n_parallel·i_cell; per-unit: Ip = ipu·ibase, Vp = vpu·vbase.
    vbase, ibase = _dc_bases(chem, ns, np); pbase = vbase * ibase
    ipu = JuMP.@variable(m, base_name = "ipu_$(battery.id)", start = 0.0)   # > 0 discharge
    JuMP.@constraint(m, ipu >= -np * chem.i_charge_max / ibase)
    JuMP.@constraint(m, ipu <=  np * chem.i_discharge_max / ibase)

    # Pack terminal voltage (per-unit): vpu = [n_series·(OCV − i_cell·R)] / vbase.
    vpu = JuMP.@expression(m, (ns * ocv - ipu * ibase * (ns * r / np)) / vbase)
    JuMP.@constraint(m, vpu >= ns * chem.v_cell_min / vbase)
    JuMP.@constraint(m, vpu <= ns * chem.v_cell_max / vbase)

    # Pack DC power = vpu·ipu·pbase (SI) coupled to the inverter DC-link power.
    JuMP.@constraint(m, vpu * ipu * pbase == inv_handles.p_dc * sb)

    return _BatHandles(ipu, vpu, vbase, ibase, soc, ns, np)
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
    i_pack = solved ? JuMP.value(bh.ipu) * bh.ibase : NaN   # per-unit → SI
    v_pack = solved ? JuMP.value(bh.vpu) * bh.vbase : NaN

    result = extract_result(ctx)
    return IVQBatteryResult(status,
                            solved ? JuMP.value(ih.p_poc) * sb : NaN,
                            solved ? JuMP.value(ih.q_poc) * sb : NaN,
                            solved ? JuMP.value(ih.p_conv) * sb : NaN,
                            solved ? v_pack * i_pack : NaN,        # p_dc (SI)
                            solved ? JuMP.value(ih.p_loss) * sb : NaN,
                            bh.soc,
                            v_pack / bh.n_series, i_pack / bh.n_parallel,   # cell = pack / n
                            v_pack, i_pack,
                            result["bus"])
end

# ── Multi-period ─────────────────────────────────────────────────────────────
#
# Several snapshots co-optimised in one model with each battery's state of charge
# linking the periods — the IVQ analogue of `solve_multiperiod_opf`. Now the SoC
# is a decision variable, so the terminal voltage `OCV(soc)` is a function of a
# variable: for a `linear`/`thevenin` chemistry (affine OCV) it embeds as a plain
# expression; for a `tabulated` chemistry the monotone-cubic OCV(soc) (and any
# R(soc)) is registered ONCE as a smooth JuMP nonlinear operator and reused
# across snapshots.
#
# Because each snapshot's current is piecewise-constant over its period, the
# charge balance q[t+1] = q[t] − i[t]·Δt is exact (`:forward`); a `:trapezoidal`
# option averages consecutive period currents for users who model nodal currents.

"""
    MultiperiodIVQResult

Result of [`solve_multiperiod_ivq`](@ref). SI units throughout.

# Fields
- `termination_status::String`, `objective::Float64`.
- `snapshots::Vector{Dict{String,Any}}` — the per-period BMOPFTools result dict.
- `dispatch::Dict{String,NamedTuple}` — per battery id: `i_cell`, `v_cell`,
  `i_pack`, `v_pack`, `p_poc`, `q_poc`, `p_dc` (each length `T`) and `soc`
  (length `T+1`, 0–1, with `soc[1]` the initial state). SI throughout.

The coupled cell + inverter model is nonconvex, so `solve_multiperiod_ivq` finds
a local optimum and may not converge for every configuration; a non-`LOCALLY_SOLVED`
/`OPTIMAL` status returns `NaN` trajectories rather than an unconverged point.
"""
struct MultiperiodIVQResult
    termination_status::String
    objective::Float64
    snapshots::Vector{Dict{String,Any}}
    dispatch::Dict{String,NamedTuple}
end

# Stamp one battery into snapshot `ctx` at SoC variable `soc_t`, coupling its DC
# power to the inverter `ih`. Modelled in PACK quantities per-unit on the
# battery's own DC bases (well-conditioned; see `_stamp_battery!`). `ocv_op`/`r_op`
# are the registered smooth operators (or `nothing` when OCV is affine / R
# constant). Returns (ipu, vpu, vbase, ibase).
function _stamp_battery_mp!(ctx, battery::IVQBattery, soc_t, ih, ocv_op, r_op)
    m = ctx.model
    chem = battery.chemistry
    ns = battery.n_series; np = battery.n_parallel
    vbase, ibase = _dc_bases(chem, ns, np); pbase = vbase * ibase

    ipu = JuMP.@variable(m, base_name = "ipu_$(battery.id)", start = 0.0)   # > 0 discharge
    JuMP.@constraint(m, ipu >= -np * chem.i_charge_max / ibase)
    JuMP.@constraint(m, ipu <=  np * chem.i_discharge_max / ibase)

    ocv_expr = if chem.ocv_affine !== nothing
        a, b = chem.ocv_affine
        JuMP.@expression(m, a + b * soc_t)
    else
        ocv_op(soc_t)
    end
    r_expr = chem.r_constant !== nothing ? chem.r_constant : r_op(soc_t)

    # Pack terminal voltage (per-unit): vpu = [n_series·(OCV(soc) − i_cell·R)] / vbase.
    vpu = JuMP.@expression(m, (ns * ocv_expr - ipu * ibase * (ns * r_expr / np)) / vbase)
    JuMP.@constraint(m, vpu >= ns * chem.v_cell_min / vbase)
    JuMP.@constraint(m, vpu <= ns * chem.v_cell_max / vbase)

    JuMP.@constraint(m, vpu * ipu * pbase == ih.p_dc * ih.sb)
    return ipu, vpu, vbase, ibase
end

"""
    solve_multiperiod_ivq(nets, batteries; kwargs...) -> MultiperiodIVQResult

Co-optimise a chronological sequence of network snapshots `nets` with a set of
[`IVQBattery`](@ref) devices whose state of charge couples the periods. Each
battery arbitrages across time subject to its cell voltage/current limits, SoC
window, and terminal/cyclic condition, exchanging power through its
[`AdvancedInverter`](@ref). The snapshots share one model and one objective (the
total generation cost across the horizon); period economics come from the
snapshots (e.g. a time-varying slack import price via each net's `voltage_source`
`cost`).

# Arguments
- `nets::Vector` — `T` network dicts (`parse_bmopf` output), one per period.
- `batteries::Vector{IVQBattery}` — each battery's `bus`/inverter must exist in
  every snapshot. `soc_init` sets the horizon start; `cyclic`/`soc_final` the
  terminal condition; `integration` (`:forward` default, or `:trapezoidal`) the
  charge-balance rule.

# Keywords
- `dt_h=1.0` — period duration in hours.
- `per_unit=false`, `s_base=1e6`, `optimizer=Ipopt.Optimizer`, `verbose=false`,
  `solver_options=()`. Results are SI regardless of `per_unit`.
"""
function solve_multiperiod_ivq(nets::AbstractVector, batteries::AbstractVector;
                               dt_h::Float64=1.0,
                               per_unit::Bool=false,
                               s_base::Float64=1e6,
                               optimizer=Ipopt.Optimizer,
                               verbose::Bool=false,
                               solver_options=())
    T = length(nets)
    T >= 1 || throw(ArgumentError("need at least one snapshot"))
    ids = [b.id for b in batteries]
    allunique(ids) || throw(ArgumentError("battery ids must be unique: $ids"))
    for b in batteries
        b.inverter.bus == b.bus || throw(ArgumentError(
            "battery '$(b.id)' bus '$(b.bus)' must match its inverter bus '$(b.inverter.bus)'"))
        b.integration in (:forward, :trapezoidal) || throw(ArgumentError(
            "battery '$(b.id)': integration must be :forward or :trapezoidal"))
        chem = b.chemistry
        chem.soc_min <= b.soc_init <= chem.soc_max || throw(ArgumentError(
            "battery '$(b.id)': soc_init=$(b.soc_init) outside the chemistry window " *
            "[$(chem.soc_min), $(chem.soc_max)]"))
    end

    model = JuMP.Model(optimizer)
    verbose || JuMP.set_silent(model)
    for (name, value) in solver_options
        JuMP.set_attribute(model, string(name), value)
    end

    # Per battery: register smooth OCV/R operators once (tabulated only), and the
    # SoC trajectory bounded to the chemistry's usable window.
    ocv_ops = Dict{String,Any}(); r_ops = Dict{String,Any}()
    socs = Dict{String,Any}()
    for b in batteries
        chem = b.chemistry
        if chem.ocv_affine === nothing
            ocv_ops[b.id] = JuMP.add_nonlinear_operator(model, 1, chem.ocv;
                                name = Symbol("ocv_$(b.id)"))
        end
        if chem.r_constant === nothing
            r_ops[b.id] = JuMP.add_nonlinear_operator(model, 1, chem.r_internal;
                              name = Symbol("rint_$(b.id)"))
        end
        soc = JuMP.@variable(model, [1:T+1], base_name = "soc_$(b.id)",
                             lower_bound = chem.soc_min, upper_bound = chem.soc_max,
                             start = b.soc_init)
        JuMP.@constraint(model, soc[1] == b.soc_init)
        socs[b.id] = soc
    end

    ipack = Dict(id => Vector{Any}(undef, T) for id in ids)   # per-unit pack current
    vpack = Dict(id => Vector{Any}(undef, T) for id in ids)   # per-unit pack voltage
    ppoc  = Dict(id => Vector{Any}(undef, T) for id in ids)
    qpoc  = Dict(id => Vector{Any}(undef, T) for id in ids)

    dcbases = Dict(b.id => _dc_bases(b.chemistry, b.n_series, b.n_parallel) for b in batteries)

    stamp_all(t) = ctx -> begin
        for b in batteries
            ih = _stamp_inverter!(ctx, b.inverter)
            ipu, vpu, _, _ = _stamp_battery_mp!(ctx, b, socs[b.id][t], ih,
                                                get(ocv_ops, b.id, nothing),
                                                get(r_ops, b.id, nothing))
            ipack[b.id][t] = ipu; vpack[b.id][t] = vpu
            ppoc[b.id][t] = ih.p_poc; qpoc[b.id][t] = ih.q_poc
        end
    end

    ctxs = [build_opf_model(nets[t]; model=model, per_unit=per_unit, s_base=s_base,
                            add_objective=false, model_hook! = stamp_all(t))
            for t in 1:T]
    sb = ctxs[1].bases === nothing ? 1.0 : ctxs[1].bases.s_base

    # Charge balance q[t+1] = q[t] − i·Δt (soc = q_pack/qp), then the terminal
    # state. The per-unit pack current ipu maps to SI amps as ipu·ibase.
    for b in batteries
        soc = socs[b.id]; ip = ipack[b.id]
        _, ibase = dcbases[b.id]
        qp = b.chemistry.q_cell * b.n_parallel   # pack charge capacity (Ah)
        for t in 1:T
            if b.integration == :forward || t == T
                JuMP.@constraint(model, soc[t+1] == soc[t] - dt_h * ip[t] * ibase / qp)
            else
                JuMP.@constraint(model, soc[t+1] == soc[t] - dt_h * (ip[t] + ip[t+1]) * ibase / (2qp))
            end
        end
        if b.soc_final !== nothing
            JuMP.@constraint(model, soc[T+1] == b.soc_final)
        elseif b.cyclic
            JuMP.@constraint(model, soc[T+1] == soc[1])
        end
    end

    JuMP.@objective(model, Min, sum(generation_cost(ctxs[t]) for t in 1:T))
    foreach(enforce_kcl!, ctxs)
    JuMP.optimize!(model)

    status = string(JuMP.termination_status(model))
    solved = JuMP.primal_status(model) == JuMP.MOI.FEASIBLE_POINT
    obj = solved ? JuMP.objective_value(model) : NaN
    snapshots = [extract_result(ctxs[t]) for t in 1:T]

    dispatch = Dict{String,NamedTuple}()
    for b in batteries
        id = b.id; ns = b.n_series; np = b.n_parallel
        vbase, ibase = dcbases[id]
        val(x) = solved ? JuMP.value(x) : NaN
        ipk = [val(ipack[id][t]) * ibase for t in 1:T]   # per-unit → SI
        vpk = [val(vpack[id][t]) * vbase for t in 1:T]
        dispatch[id] = (
            i_cell = ipk ./ np,                           # cell = pack / n
            v_cell = vpk ./ ns,
            i_pack = ipk,
            v_pack = vpk,
            p_poc  = [val(ppoc[id][t]) * sb for t in 1:T],
            q_poc  = [val(qpoc[id][t]) * sb for t in 1:T],
            p_dc   = ipk .* vpk,
            soc    = [val(socs[id][k]) for k in 1:T+1],
        )
    end

    return MultiperiodIVQResult(status, obj, snapshots, dispatch)
end
