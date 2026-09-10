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
  `soc[1]` = initial). `energy_wh` is the explicit energy name (`soc` is retained
  for compatibility). PE devices also report phase P/Q, conductor currents,
  operating-mode metadata and independent physical residuals in `diagnostics`.
- `solve::SolveStatus` — exact normalized status and publication decision.
"""
struct MultiperiodResult <: AbstractSolveResult
    termination_status::String
    objective::Float64
    snapshots::Vector{Dict{String,Any}}
    dispatch::Dict{String,NamedTuple}
    solve::SolveStatus
end

# Terminal / departure state-of-charge constraint, dispatched per device type.
# `soc` is the per-unit SOC vector (length T+1); `sb` the VA base; `T` the horizon.
function _finalize_soc!(model, d::StorageDevice, soc, sb, T)
    if d.energy_final !== nothing
        JuMP.fix(soc[T+1], d.energy_final / sb; force=true)
    elseif d.cyclic
        JuMP.@constraint(model, soc[T+1] == soc[1])
    end
end

function _finalize_soc!(model, d::EVDevice, soc, sb, T)
    dp = d.departure_period === nothing ? T : d.departure_period
    1 <= dp <= T || throw(ArgumentError(
        "EV '$(d.id)': departure_period=$(d.departure_period) out of range 1:$T"))
    if d.energy_final !== nothing
        JuMP.fix(soc[T+1], d.energy_final / sb; force=true)
    end
    # Avoid a redundant inequality on an already fixed terminal state. It can
    # otherwise make a full-battery, biactive MPCC unnecessarily degenerate.
    if !(dp == T && d.energy_final !== nothing && d.energy_final >= d.departure_energy)
        JuMP.@constraint(model, soc[dp+1] >= d.departure_energy / sb)
    end
end

# Link a device's per-period ports through its state of charge (all per-unit).
#   soc[t+1] = soc[t] + (eff_c·pc[t] − pd[t]/eff_d)·Δt ,   e_min ≤ soc ≤ e_max
function link_device!(model, dev::Union{StorageDevice,EVDevice},
                      ports::AbstractVector, sb, grid::TimeGrid)
    T = length(grid)
    soc = JuMP.@variable(model, [1:T+1], base_name = "soc_$(_dev_id(dev))",
        lower_bound = _dev_emin(dev)/sb, upper_bound = _dev_emax(dev)/sb,
        start = _dev_einit(dev)/sb)
    JuMP.fix(soc[1], _dev_einit(dev)/sb; force=true)
    effc = _dev_effc(dev); effd = _dev_effd(dev)
    for t in 1:T
        JuMP.@constraint(model,
            soc[t+1] == soc[t] +
                        (effc*ports[t].pc - ports[t].pd/effd) * grid[t])
    end
    _finalize_soc!(model, dev, soc, sb, T)
    # Retain physical interval lengths for independent post-solve accounting.
    return (energy=soc, durations_h=copy(grid.durations_h))
end

_link_soc!(model, dev, ports::AbstractVector, sb, dt_h, T) =
    link_device!(model, dev, ports, sb, TimeGrid(T, dt_h))

function extract_device(dev::Union{StorageDevice,EVDevice},
                        ports::AbstractVector, soc, sb,
                        status::SolveStatus)
    T = length(ports)
    grid = soc.durations_h
    soc = soc.energy
    val(x) = status.publishable ? JuMP.value(x) : NaN
    energy = [val(soc[k]) * sb for k in 1:T+1]
    pc = [val(h.pc) * sb for h in ports]
    pd = [val(h.pd) * sb for h in ports]
    a = dev.p_charge_max > 0 ? pc ./ dev.p_charge_max : zeros(T)
    b = dev.p_discharge_max > 0 ? pd ./ dev.p_discharge_max : zeros(T)
    ir = [[val(c[1])*h.ib for c in h.currents] for h in ports]
    ii = [[val(c[2])*h.ib for c in h.currents] for h in ports]
    im = [hypot.(r,i) for (r,i) in zip(ir,ii)]
    neutral = [hypot(sum(r),sum(i)) for (r,i) in zip(ir,ii)]
    phase_p = [[val(p)*sb for p in h.phase_p] for h in ports]
    phase_q = [[val(q)*sb for q in h.phase_q] for h in ports]
    snorm = [hypot(val(h.p)*sb,val(h.q)*sb) for h in ports]
    departure = dev isa EVDevice ? (dev.departure_period === nothing ? T : dev.departure_period) : T
    target = dev isa EVDevice ? dev.departure_energy : nothing
    terminal = dev.energy_final !== nothing ? dev.energy_final :
        dev isa StorageDevice && dev.cyclic ? dev.energy_init : nothing
    diagnostics = (
        terminal_energy_error_wh=terminal === nothing ? nothing : energy[end]-terminal,
        ac_power_balance_w=[sum(phase_p[t])-(pd[t]-pc[t]) for t in 1:T],
        reactive_bound_violation_var=[ports[t].active ? max(0.0,dev.q_min-sum(phase_q[t]),sum(phase_q[t])-dev.q_max) : abs(sum(phase_q[t])) for t in 1:T],
        energy_balance_wh=diff(energy) .- grid .* (dev.eff_charge .* pc .- pd ./ dev.eff_discharge),
        complementarity_product=a .* b,
        complementarity_minimum=min.(a,b),
        simultaneous_power_w=min.(pc,pd),
        energy_bound_violation_wh=max.(dev.energy_min .- energy, energy .- dev.energy_max, 0.0),
        departure_shortfall_wh=target === nothing ? 0.0 : max(0.0,target-energy[departure+1]),
        power_limit_violation_w=max.(pc .- dev.p_charge_max,pd .- dev.p_discharge_max,-pc,-pd,0.0),
        current_limit_violation_a=dev.i_max === nothing ? nothing : [max(0.0,maximum(i)-dev.i_max) for i in im],
        neutral_limit_violation_a=dev.neutral_current_max === nothing ? nothing : max.(0.0,neutral .- dev.neutral_current_max),
        apparent_power_violation_va=dev.s_max === nothing ? nothing : max.(0.0,snorm .- dev.s_max),
        disconnected_current_a=[ports[t].active ? 0.0 : maximum(im[t]) for t in 1:T],
        phase_sharing_residual_w=dev.phase_policy == :equal_power ? [maximum(abs.(p .- p[1])) for p in phase_p] : nothing,
        phase_sharing_residual_var=dev.phase_policy == :equal_power ? [maximum(abs.(q .- q[1])) for q in phase_q] : nothing,
    )
    return (
        energy_wh=energy,
        current_real_a=ir, current_imag_a=ii, current_magnitude_a=im,
        neutral_current_a=dev.neutral === nothing ? nothing : neutral,
        phase_power_w=phase_p, phase_reactive_var=phase_q,
        operation=dev.operation, complementarity_tolerance=dev.operation == :relaxed ? dev.complementarity_tolerance : nothing,
        diagnostics=diagnostics,
        p_charge    = [val(ports[t].pc) * sb for t in 1:T],
        p_discharge = [val(ports[t].pd) * sb for t in 1:T],
        p_net       = [val(ports[t].p)  * sb for t in 1:T],
        q           = [val(ports[t].q)  * sb for t in 1:T],
        soc         = [val(soc[k])      * sb for k in 1:T+1],
    )
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
- `devices::Vector` — [`AbstractDevice`](@ref) instances implementing the
  validation/stamp/link/extract lifecycle. Built-in storage and EV devices
  require their bus/terminals to exist in every snapshot.

# Keywords
- `dt_h=1.0` — uniform period duration in hours (compatibility shorthand).
- `time_grid=nothing` — pass `TimeGrid([Δt₁, Δt₂, ...])` for nonuniform
  durations. When supplied it takes precedence over `dt_h`.
- `per_unit=true`, `s_base=1e6` — engine unit handling (results are SI regardless).
- `optimizer=Ipopt.Optimizer`, `verbose=false`, `solver_options=()` — solver control.
- `configure! = identity` — callback on the empty JuMP model before optimizer
  attachment, e.g. `MathOptComplements.Bridges.add_all_bridges` for CCOpt.

Bidirectional PE devices default to an explicitly approximate normalized product
relaxation. Exact mode exclusion uses `operation=:complementarity` and a suitable
backend. Published status does not certify exact complementarity or global
optimality. Inspect `dispatch[id].diagnostics` in physical units.

# Returns
A [`MultiperiodResult`](@ref) with the per-period solutions and each device's SI
charge/discharge/SOC trajectory.
"""
function solve_multiperiod_opf(nets::AbstractVector, devices::AbstractVector;
                               dt_h::Real=1.0,
                               time_grid::Union{Nothing,TimeGrid}=nothing,
                               per_unit::Bool=true,
                               s_base::Float64=1e6,
                               optimizer=Ipopt.Optimizer,
                               verbose::Bool=false,
                               solver_options=(),
                               configure!::Function=identity)
    T = length(nets)
    T >= 1 || throw(ArgumentError("need at least one snapshot"))
    grid = _resolve_time_grid(T, dt_h, time_grid)
    all(d -> d isa AbstractDevice, devices) || throw(ArgumentError(
        "devices must contain only AbstractDevice values"))
    foreach(d -> validate_device(d, nets; periods=T), devices)
    _validate_charging_assignments(devices)
    ids = [device_id(d) for d in devices]
    allunique(ids) || throw(ArgumentError("device ids must be unique: $ids"))

    # Compose each session once rather than rebuilding its horizon mask in every
    # snapshot. Linking/extraction still dispatch on the original session.
    compiled = Dict(device_id(d) => (d isa ChargingSession ? _session_device(d) : d) for d in devices)

    # ports[dev.id][t] :: PortHandle. Filled as each snapshot's hook runs.
    ports = Dict{String,Vector{Any}}(id => Vector{Any}(undef, T) for id in ids)

    stamp_all(t) = ctx -> begin
        for d in devices
            ports[device_id(d)][t] = stamp_device!(ctx, compiled[device_id(d)]; period=t)
        end
    end

    # Build every snapshot into the shared model; the engine adds no objective.
    multi = build_multi_context(nets; hook_factory=stamp_all, per_unit, s_base,
                                optimizer, verbose, solver_options, configure!)
    model = multi.model
    ctxs = multi.contexts

    sb = _sbase(ctxs[1])

    # State-of-charge linking per device.
    socs = Dict{String,Any}(device_id(d) =>
                            link_device!(model, d, ports[device_id(d)], sb, grid)
                            for d in devices)

    # One combined objective: total generation cost across the horizon.
    JuMP.@objective(model, Min,
        sum(grid[t] * (generation_cost(ctxs[t]) +
            sum(_snapshot_device_cost(d,ports[device_id(d)][t]) for d in devices;init=0.0)) for t in 1:T))

    foreach(enforce_kcl!, ctxs)
    JuMP.optimize!(model)

    outcome = _solve_outcome(model)
    status_contract = SolveStatus(outcome)
    status = string(outcome.termination_status)
    solved = _publishable(outcome)
    obj = solved ? JuMP.objective_value(model) : NaN

    snapshots = [_extract_result(ctxs[t], outcome) for t in 1:T]
    for t in 1:T, d in devices
        injection=_snapshot_device_injection(d,ports[device_id(d)][t],status_contract)
        injection===nothing && continue
        ledger=get!(snapshots[t],"custom_injection",Dict("p"=>0.0,"q"=>0.0))
        ledger["p"]+=injection.p; ledger["q"]+=injection.q
    end

    dispatch = Dict{String,NamedTuple}()
    for d in devices
        id = device_id(d)
        dispatch[id] = extract_device(d, ports[id], socs[id], sb, status_contract)
    end

    return MultiperiodResult(status, obj, snapshots, dispatch, status_contract)
end

solve_status(result::MultiperiodResult) = result.solve

solve_diagnostics(result::MultiperiodResult) =
    (objective=result.objective, periods=length(result.snapshots),
     devices=length(result.dispatch),
     device_diagnostics=Dict(id=>d.diagnostics for (id,d) in result.dispatch if hasproperty(d,:diagnostics)))

function link_device!(model, s::ChargingSession, ports::AbstractVector, sb, grid::TimeGrid)
    state = link_device!(model, _session_device(s), ports, sb, grid)
    f = s.ev.charge_acceptance
    if f !== nothing
        for t in 1:length(grid)
            ports[t].active || continue  # exact disconnection, even at a zero cap
            # Piecewise-constant powers imply affine PE energy within a period.
            # Both boundaries bound a nonincreasing cap over that whole path.
            for k in (t,t+1)
                fraction = JuMP.@expression(model, state.energy[k]*sb/s.ev.energy_max)
                formulate_pwl_relation!(model,f,fraction,ports[t].pc;
                    domain=(s.ev.energy_min/s.ev.energy_max,1.0),relation=:upper,
                    formulation=s.acceptance_formulation,output_scale=sb,
                    conservative=s.acceptance_conservative && s.acceptance_formulation isa AbstractPWLSmoothing)
            end
        end
    end
    return state
end

function extract_device(s::ChargingSession, ports::AbstractVector, soc, sb, status::SolveStatus)
    d = extract_device(_session_device(s), ports, soc, sb, status)
    f = s.ev.charge_acceptance
    acceptance_error = if f === nothing
        nothing
    elseif !status.publishable
        fill(NaN,length(ports))
    else
        [ports[t].active ? max(0.0,d.p_charge[t] - min(
            _pwl_exact(f,d.energy_wh[t]/s.ev.energy_max),
            _pwl_exact(f,d.energy_wh[t+1]/s.ev.energy_max))) : 0.0 for t in eachindex(ports)]
    end
    merge(d,(vehicle_id=s.ev.id, evse_id=s.evse.id, occupied=copy(s.available),
        acceptance_formulation=s.acceptance_formulation,
        acceptance_conservative=s.acceptance_conservative,
        diagnostics=merge(d.diagnostics,(acceptance_violation_w=acceptance_error,))))
end

function _snapshot_device_injection(::Union{StorageDevice,EVDevice,ChargingSession}, h::PortHandle, status::SolveStatus)
    (p=status.publishable ? JuMP.value(h.p)*h.sb : NaN,
     q=status.publishable ? JuMP.value(h.q)*h.sb : NaN)
end
