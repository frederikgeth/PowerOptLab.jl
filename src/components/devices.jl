# Storage-class devices (battery, EV) and the low-level machinery that stamps
# their inverter port into one network snapshot via a BMOPFTools `model_hook!`.
#
# A device contributes, per snapshot:
#   - a current injection (cr, ci) on each phase terminal, added to the engine's
#     KCL accumulators so the device draws/supplies real physical current;
#   - nonnegative charge/discharge power split (pc, pd) with the AC injection
#     P = pd − pc (discharge positive), so an efficiency-aware state-of-charge
#     can be linked across snapshots (see `multiperiod.jl`);
#   - reactive-power box bounds on the aggregate device Q.
#
# All model quantities are per-unit when the engine builds per-unit (the default);
# device parameters are SI (watts, vars, watt-hours) and are scaled by the
# snapshot's `s_base` at stamp time. A device with `s_base = nothing` (SI solve)
# uses a scale of 1.

"""
    StorageDevice(; id, bus, kwargs...)

A grid-connected battery (or generic energy-storage) inverter with an
inter-temporal state of charge. Powers are SI watts, energies SI watt-hours.

# Required
- `id::String` — unique device id.
- `bus::String` — connection bus.
- `p_charge_max`, `p_discharge_max` — charge / discharge power limits (W ≥ 0).
- `energy_max` — usable energy capacity (Wh).
- `energy_init` — energy at the start of the horizon (Wh).

# Optional
- `phase_terminals=["1"]`, `neutral="n"` — the inverter's phase conductor(s) and
  the return terminal (`nothing` if referenced directly to ground).
- `energy_min=0.0` — lower energy bound (Wh).
- `q_min=0.0`, `q_max=0.0` — reactive-power box (var); default is unity power factor.
- `eff_charge=1.0`, `eff_discharge=1.0` — one-way efficiencies (0,1].
- `cyclic=true` — require the terminal state of charge to equal `energy_init`.
- `energy_final=nothing` — if set, pin the terminal energy to this value (Wh),
  overriding `cyclic`.
"""
Base.@kwdef struct StorageDevice
    id::String
    bus::String
    phase_terminals::Vector{String} = ["1"]
    neutral::Union{String,Nothing} = "n"
    p_charge_max::Float64
    p_discharge_max::Float64
    q_min::Float64 = 0.0
    q_max::Float64 = 0.0
    energy_max::Float64
    energy_init::Float64
    energy_min::Float64 = 0.0
    eff_charge::Float64 = 1.0
    eff_discharge::Float64 = 1.0
    cyclic::Bool = true
    energy_final::Union{Float64,Nothing} = nothing
end

"""
    EVDevice(; id, bus, available, departure_energy, kwargs...)

An electric-vehicle charger: a storage device that is only controllable while
plugged in and must reach a target energy by departure. Set `p_discharge_max > 0`
for bidirectional (V2G) charging; the default `0.0` gives unidirectional (V1G).

# Required
- `id`, `bus` — as [`StorageDevice`](@ref).
- `p_charge_max` — charge power limit (W).
- `energy_max`, `energy_init` — battery capacity and plug-in energy (Wh).
- `available::Vector{Bool}` — per period, whether the vehicle is plugged in.
  While unavailable the charger is idle and the state of charge is held.
- `departure_energy` — energy (Wh) required by `departure_period`.

# Optional
- `p_discharge_max=0.0` — V2G discharge limit (W); `0.0` ⇒ charge-only.
- `departure_period=nothing` — period index by whose end `departure_energy` must
  be met; `nothing` ⇒ the end of the horizon.
- `phase_terminals`, `neutral`, `energy_min`, `q_min`, `q_max`,
  `eff_charge`, `eff_discharge` — as [`StorageDevice`](@ref).
"""
Base.@kwdef struct EVDevice
    id::String
    bus::String
    phase_terminals::Vector{String} = ["1"]
    neutral::Union{String,Nothing} = "n"
    p_charge_max::Float64
    p_discharge_max::Float64 = 0.0
    q_min::Float64 = 0.0
    q_max::Float64 = 0.0
    energy_max::Float64
    energy_init::Float64
    energy_min::Float64 = 0.0
    eff_charge::Float64 = 1.0
    eff_discharge::Float64 = 1.0
    available::Vector{Bool}
    departure_energy::Float64
    departure_period::Union{Int,Nothing} = nothing
end

# Uniform view over the fields the port-stamping and SOC-linking code needs, so
# StorageDevice and EVDevice share one implementation.
_dev_id(d) = d.id
_dev_bus(d) = d.bus
_dev_phases(d) = d.phase_terminals
_dev_neutral(d) = d.neutral
_dev_pcmax(d) = d.p_charge_max
_dev_pdmax(d) = d.p_discharge_max
_dev_qmin(d) = d.q_min
_dev_qmax(d) = d.q_max
_dev_emax(d) = d.energy_max
_dev_emin(d) = d.energy_min
_dev_einit(d) = d.energy_init
_dev_effc(d) = d.eff_charge
_dev_effd(d) = d.eff_discharge

function _validate_connection(id, bus, phases, neutral, nets)
    isempty(id) && throw(ArgumentError("device id must not be empty"))
    isempty(bus) && throw(ArgumentError("device '$id' bus must not be empty"))
    isempty(phases) && throw(ArgumentError("device '$id' needs at least one phase terminal"))
    allunique(phases) || throw(ArgumentError("device '$id' phase terminals must be unique"))
    neutral in phases && throw(ArgumentError(
        "device '$id' neutral cannot also be a phase terminal"))
    for (t, net) in enumerate(nets)
        buses = get(net, "bus", Dict())
        haskey(buses, bus) || throw(ArgumentError(
            "snapshot $t: device '$id' bus '$bus' not found"))
        terminals = Set(String.(get(buses[bus], "terminal_names", String[])))
        for phase in phases
            phase in terminals || throw(ArgumentError(
                "snapshot $t: device '$id' phase terminal '$phase' not found at bus '$bus'"))
        end
        neutral === nothing || neutral in terminals || throw(ArgumentError(
            "snapshot $t: device '$id' neutral '$neutral' not found at bus '$bus'"))
    end
    return nothing
end

function _validate_storage_values(d)
    values = (("p_charge_max", d.p_charge_max),
              ("p_discharge_max", d.p_discharge_max),
              ("energy_min", d.energy_min),
              ("energy_max", d.energy_max),
              ("energy_init", d.energy_init),
              ("q_min", d.q_min), ("q_max", d.q_max),
              ("eff_charge", d.eff_charge), ("eff_discharge", d.eff_discharge))
    all(isfinite(value) for (_, value) in values) || throw(ArgumentError(
        "device '$(d.id)' power, energy, reactive-power, and efficiency values must be finite"))
    d.p_charge_max >= 0 || throw(ArgumentError(
        "device '$(d.id)' p_charge_max must be >= 0"))
    d.p_discharge_max >= 0 || throw(ArgumentError(
        "device '$(d.id)' p_discharge_max must be >= 0"))
    d.q_min <= d.q_max || throw(ArgumentError(
        "device '$(d.id)' requires q_min <= q_max"))
    0 <= d.energy_min <= d.energy_init <= d.energy_max || throw(ArgumentError(
        "device '$(d.id)' requires 0 <= energy_min <= energy_init <= energy_max"))
    0 < d.eff_charge <= 1 || throw(ArgumentError(
        "device '$(d.id)' eff_charge must lie in (0, 1]"))
    0 < d.eff_discharge <= 1 || throw(ArgumentError(
        "device '$(d.id)' eff_discharge must lie in (0, 1]"))
    return nothing
end

function _validate_device(d::StorageDevice, T, nets)
    _validate_connection(d.id, d.bus, d.phase_terminals, d.neutral, nets)
    _validate_storage_values(d)
    if d.energy_final !== nothing
        isfinite(d.energy_final) && d.energy_min <= d.energy_final <= d.energy_max ||
            throw(ArgumentError(
                "device '$(d.id)' energy_final must be finite and within its energy bounds"))
    end
    return nothing
end

function _validate_device(d::EVDevice, T, nets)
    _validate_connection(d.id, d.bus, d.phase_terminals, d.neutral, nets)
    _validate_storage_values(d)
    length(d.available) == T || throw(ArgumentError(
        "EV '$(d.id)' availability must contain exactly $T entries"))
    isfinite(d.departure_energy) && d.energy_min <= d.departure_energy <= d.energy_max ||
        throw(ArgumentError(
            "EV '$(d.id)' departure_energy must be finite and within its energy bounds"))
    dp = d.departure_period === nothing ? T : d.departure_period
    1 <= dp <= T || throw(ArgumentError(
        "EV '$(d.id)': departure_period=$(d.departure_period) out of range 1:$T"))
    return nothing
end

# Per-period availability: batteries are always controllable; EVs follow their mask.
_dev_available(d::StorageDevice, t::Int) = true
_dev_available(d::EVDevice, t::Int) = t <= length(d.available) ? d.available[t] : false

# Handle returned by `_stamp_port!`: the JuMP variables/expressions a snapshot
# exposes for state-of-charge linking and result extraction. Powers are in the
# model's units (per-unit unless the solve is SI).
struct PortHandle
    pc            # charge power variable (≥ 0)
    pd            # discharge power variable (≥ 0)
    p             # net AC injection expression = pd − pc (discharge positive)
    q             # net AC reactive injection expression
end

# s_base for a context (VA), or 1.0 for an SI solve.
_sbase(ctx) = ctx.bases === nothing ? 1.0 : ctx.bases.s_base

# Phase-to-neutral voltage difference (VariableRef or AffExpr) at (bus, phase).
function _dv(ctx, bus, ph, neutral)
    vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
    if neutral === nothing
        return vr[(bus, ph)], vi[(bus, ph)]
    end
    m = ctx.model
    return JuMP.@expression(m, vr[(bus, ph)] - vr[(bus, neutral)]),
           JuMP.@expression(m, vi[(bus, ph)] - vi[(bus, neutral)])
end

"""
    _stamp_port!(ctx, dev; active=true) -> PortHandle

Stamp `dev`'s inverter port into snapshot `ctx`: declare per-phase current
injections and add them to KCL, declare the charge/discharge power split, and
bound the aggregate device P (through the split) and Q. Returns the [`PortHandle`]
used by the SOC linker and result extraction. `active=false` pins the port to
zero (an unplugged EV): no current, no power.
"""
function _stamp_port!(ctx, dev; active::Bool=true)
    m   = ctx.model
    sb  = _sbase(ctx)
    bus = _dev_bus(dev)
    phases = _dev_phases(dev)
    neutral = _dev_neutral(dev)
    id  = _dev_id(dev)

    kcl_r = ctx.kcl_r; kcl_i = ctx.kcl_i

    # Per-phase current injections summed into the aggregate device power.
    P = zero(JuMP.QuadExpr)
    Q = zero(JuMP.QuadExpr)
    for (k, ph) in enumerate(phases)
        cr = JuMP.@variable(m, base_name = "cr_$(id)_$(ph)")
        ci = JuMP.@variable(m, base_name = "ci_$(id)_$(ph)")
        dvr, dvi = _dv(ctx, bus, ph, neutral)
        P += JuMP.@expression(m, dvr*cr + dvi*ci)      # injected into network
        Q += JuMP.@expression(m, dvi*cr - dvr*ci)
        JuMP.add_to_expression!(kcl_r[(bus, ph)],  cr)
        JuMP.add_to_expression!(kcl_i[(bus, ph)],  ci)
        if neutral !== nothing
            JuMP.add_to_expression!(kcl_r[(bus, neutral)], -cr)
            JuMP.add_to_expression!(kcl_i[(bus, neutral)], -ci)
        end
    end

    # Charge/discharge power split (per-unit). Round-trip loss (eff < 1) makes
    # simultaneous pc,pd > 0 suboptimal, so the split stays physical without a
    # complementarity constraint.
    pc = JuMP.@variable(m, base_name = "pc_$(id)", lower_bound = 0.0)
    pd = JuMP.@variable(m, base_name = "pd_$(id)", lower_bound = 0.0)
    if active
        JuMP.@constraint(m, pc <= _dev_pcmax(dev) / sb)
        JuMP.@constraint(m, pd <= _dev_pdmax(dev) / sb)
        JuMP.@constraint(m, P == pd - pc)
        JuMP.@constraint(m, Q >= _dev_qmin(dev) / sb)
        JuMP.@constraint(m, Q <= _dev_qmax(dev) / sb)
    else
        # Unplugged: no exchange with the grid this period.
        JuMP.@constraint(m, pc == 0.0)
        JuMP.@constraint(m, pd == 0.0)
        JuMP.@constraint(m, P == 0.0)
        JuMP.@constraint(m, Q == 0.0)
    end

    return PortHandle(pc, pd, P, Q)
end
