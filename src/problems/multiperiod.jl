# Multi-period optimal power flow: several network snapshots co-optimised in one
# JuMP model, with storage/EV devices whose state of charge links each step to
# the next. Built entirely on the BMOPFTools staged API — `build_opf_model`
# (one shared model, `add_objective=false`), `generation_cost`, `enforce_kcl!`,
# `extract_result` — plus this package's device port stamping and SOC linking.

"""
    MultiperiodResult

Result of [`solve_multiperiod_opf`](@ref).

# Fields
- `termination_status::String` — JuMP status of the single shared solve.
- `objective::Float64` — optimal objective value.
- `snapshots::Vector{Dict{String,Any}}` — the per-period BMOPFTools result dict
  (SI), one per input net, in order.
- `dispatch::Dict{String,NamedTuple}` — per device id, SI trajectories over the
  horizon: `p_charge`, `p_discharge`, `p_net` (discharge positive), `q`
  (each length `T`), and `soc` (length `T+1`, energy in Wh at each step boundary,
  `soc[1]` = initial).
"""
struct MultiperiodResult
    termination_status::String
    objective::Float64
    snapshots::Vector{Dict{String,Any}}
    dispatch::Dict{String,NamedTuple}
end

# Terminal / departure state-of-charge constraint, dispatched per device type.
# `soc` is the per-unit SOC vector (length T+1); `sb` the VA base; `T` the horizon.
function _finalize_soc!(model, d::StorageDevice, soc, sb, T)
    if d.energy_final !== nothing
        JuMP.@constraint(model, soc[T+1] == d.energy_final / sb)
    elseif d.cyclic
        JuMP.@constraint(model, soc[T+1] == soc[1])
    end
end

function _finalize_soc!(model, d::EVDevice, soc, sb, T)
    dp = d.departure_period === nothing ? T : d.departure_period
    1 <= dp <= T || throw(ArgumentError(
        "EV '$(d.id)': departure_period=$(d.departure_period) out of range 1:$T"))
    # Must hold at least the required energy by the end of the departure period.
    JuMP.@constraint(model, soc[dp+1] >= d.departure_energy / sb)
end

# Link a device's per-period ports through its state of charge (all per-unit).
#   soc[t+1] = soc[t] + (eff_c·pc[t] − pd[t]/eff_d)·Δt ,   e_min ≤ soc ≤ e_max
function _link_soc!(model, dev, ports::Vector{PortHandle}, sb, dt_h, T)
    soc = JuMP.@variable(model, [1:T+1], base_name = "soc_$(_dev_id(dev))")
    JuMP.@constraint(model, soc[1] == _dev_einit(dev) / sb)
    effc = _dev_effc(dev); effd = _dev_effd(dev)
    for t in 1:T
        JuMP.@constraint(model,
            soc[t+1] == soc[t] + (effc*ports[t].pc - ports[t].pd/effd) * dt_h)
        JuMP.@constraint(model, soc[t+1] >= _dev_emin(dev) / sb)
        JuMP.@constraint(model, soc[t+1] <= _dev_emax(dev) / sb)
    end
    _finalize_soc!(model, dev, soc, sb, T)
    return soc
end

"""
    solve_multiperiod_opf(nets, devices; kwargs...) -> MultiperiodResult

Co-optimise a sequence of network snapshots `nets` (one BMOPFTools net dict per
period, in chronological order) with a set of storage/EV `devices` whose state of
charge couples the periods. The snapshots share one JuMP model and one objective
(the sum of each snapshot's generation cost); the devices arbitrage across time
subject to their power, energy, efficiency, and terminal/departure constraints.

Per-period economics come from the snapshots themselves — e.g. a time-varying
slack import price set via each net's `voltage_source` `cost`, or differing loads.

# Arguments
- `nets::Vector` — `T` network dicts (`parse_bmopf` output), one per period.
- `devices::Vector` — [`StorageDevice`](@ref) / [`EVDevice`](@ref) instances. A
  device's `bus`/terminals must exist in every snapshot.

# Keywords
- `dt_h=1.0` — period duration in hours (SOC integrates power over this).
- `per_unit=true`, `s_base=1e6` — engine unit handling (results are SI regardless).
- `optimizer=Ipopt.Optimizer`, `verbose=false`, `solver_options=()` — solver control.

# Returns
A [`MultiperiodResult`](@ref) with the per-period solutions and each device's SI
charge/discharge/SOC trajectory.
"""
function solve_multiperiod_opf(nets::AbstractVector, devices::AbstractVector;
                               dt_h::Real=1.0,
                               per_unit::Bool=true,
                               s_base::Float64=1e6,
                               optimizer=Ipopt.Optimizer,
                               verbose::Bool=false,
                               solver_options=())
    T = length(nets)
    T >= 1 || throw(ArgumentError("need at least one snapshot"))
    isfinite(dt_h) && dt_h > 0 || throw(ArgumentError(
        "dt_h must be a positive finite number of hours"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be finite and > 0"))
    all(d -> d isa Union{StorageDevice,EVDevice}, devices) || throw(ArgumentError(
        "devices must contain only StorageDevice or EVDevice values"))
    foreach(d -> _validate_device(d, T, nets), devices)
    ids = [_dev_id(d) for d in devices]
    allunique(ids) || throw(ArgumentError("device ids must be unique: $ids"))

    model = JuMP.Model(optimizer)
    verbose || JuMP.set_silent(model)
    _set_solver_options!(model, solver_options)

    # ports[dev.id][t] :: PortHandle. Filled as each snapshot's hook runs.
    ports = Dict{String,Vector{PortHandle}}(id => Vector{PortHandle}(undef, T) for id in ids)

    stamp_all(t) = ctx -> begin
        for d in devices
            ports[_dev_id(d)][t] = _stamp_port!(ctx, d; active=_dev_available(d, t))
        end
    end

    # Build every snapshot into the shared model; the engine adds no objective.
    ctxs = [build_opf_model(nets[t]; model=model, per_unit=per_unit, s_base=s_base,
                            add_objective=false, model_hook! = stamp_all(t))
            for t in 1:T]

    sb = _sbase(ctxs[1])

    # State-of-charge linking per device.
    socs = Dict{String,Any}(_dev_id(d) => _link_soc!(model, d, ports[_dev_id(d)], sb, dt_h, T)
                            for d in devices)

    # One combined objective: total generation cost across the horizon.
    JuMP.@objective(model, Min, Float64(dt_h) *
        sum(generation_cost(ctxs[t]) for t in 1:T))

    foreach(enforce_kcl!, ctxs)
    JuMP.optimize!(model)

    outcome = _solve_outcome(model)
    status = string(outcome.termination_status)
    solved = _publishable(outcome)
    obj = solved ? JuMP.objective_value(model) : NaN

    snapshots = [_extract_result(ctxs[t], outcome) for t in 1:T]

    dispatch = Dict{String,NamedTuple}()
    for d in devices
        id = _dev_id(d); ph = ports[id]; soc = socs[id]
        val(x) = solved ? JuMP.value(x) : NaN
        dispatch[id] = (
            p_charge    = [val(ph[t].pc) * sb for t in 1:T],
            p_discharge = [val(ph[t].pd) * sb for t in 1:T],
            p_net       = [val(ph[t].p)  * sb for t in 1:T],
            q           = [val(ph[t].q)  * sb for t in 1:T],
            soc         = [val(soc[k])   * sb for k in 1:T+1],
        )
    end

    return MultiperiodResult(status, obj, snapshots, dispatch)
end
