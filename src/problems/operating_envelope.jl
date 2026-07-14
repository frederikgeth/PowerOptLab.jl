# Dynamic operating envelopes (DOEs).
#
# A DOE is the time-varying active-power export limit a network operator can
# allocate to each connection point (DER) such that, when every point exports at
# its allocated limit simultaneously, the network's operational constraints
# (voltage bounds, thermal limits) still hold. "Dynamic" refers to recomputing
# the envelope per interval as baseline conditions (load) change — the intervals
# are otherwise independent (no inter-temporal coupling), so each is a separate
# constrained allocation solve.
#
# This is a DIFFERENT problem specification over the same physics, and — unlike
# state estimation — it KEEPS the engine's operational bounds: those bounds are
# exactly what limit the envelope. It is built on `solve_opf` with a `model_hook!`
# that (a) stamps a free export injection at each connection point, and (b)
# replaces the generation-cost objective with a fairness objective over the
# allocated exports; a `solution_hook!` reads the allocation back out.

"""
    ConnectionPoint(; id, bus, export_max, kwargs...)

A DER connection point whose export limit a [`solve_operating_envelope`](@ref)
allocation computes. `export_max` is the inverter/connection active-power ceiling
(W); the allocated envelope never exceeds it and is usually smaller, bound by
network voltage/thermal constraints.

- `id::String`, `bus::String` — identifier and connection bus.
- `export_max::Float64` — maximum exportable active power (W ≥ 0).
- `phase_terminals=["1"]`, `neutral="n"` — phase conductor(s) and return terminal
  (`nothing` if referenced directly to ground). Export is at unity power factor.
"""
Base.@kwdef struct ConnectionPoint
    id::String
    bus::String
    phase_terminals::Vector{String} = ["1"]
    neutral::Union{String,Nothing} = "n"
    export_max::Float64
end

"""
    OperatingEnvelopeResult

Result of [`solve_operating_envelope`](@ref). All powers SI (W).

# Fields
- `termination_status::Vector{String}` — solver status per interval.
- `envelope::Dict{String,Vector{Float64}}` — per connection-point id, the
  allocated export limit (W) for each interval.
- `total_export::Vector{Float64}` — sum of allocated exports across points, per
  interval.
- `snapshots::Vector{Dict{String,Any}}` — the per-interval BMOPFTools result dict
  (bus voltages at the allocation, etc.).
"""
struct OperatingEnvelopeResult
    termination_status::Vector{String}
    envelope::Dict{String,Vector{Float64}}
    total_export::Vector{Float64}
    snapshots::Vector{Dict{String,Any}}
end

# Stamp a connection point's export port into `ctx`: a free active-power export
# `pe ∈ [0, export_max]` injected at unity power factor, added to KCL. Returns the
# per-unit export variable.
function _stamp_export_port!(ctx, cp::ConnectionPoint)
    m = ctx.model
    sb = _sbase(ctx)
    bus = cp.bus; phases = cp.phase_terminals; neutral = cp.neutral; id = cp.id

    pe = JuMP.@variable(m, base_name = "export_$(id)", lower_bound = 0.0)
    JuMP.@constraint(m, pe <= cp.export_max / sb)

    P = zero(JuMP.QuadExpr)
    Q = zero(JuMP.QuadExpr)
    for ph in phases
        cr = JuMP.@variable(m, base_name = "cr_ex_$(id)_$(ph)")
        ci = JuMP.@variable(m, base_name = "ci_ex_$(id)_$(ph)")
        dvr, dvi = _dv(ctx, bus, ph, neutral)
        P += JuMP.@expression(m, dvr*cr + dvi*ci)
        Q += JuMP.@expression(m, dvi*cr - dvr*ci)
        JuMP.add_to_expression!(ctx.kcl_r[(bus, ph)], cr)
        JuMP.add_to_expression!(ctx.kcl_i[(bus, ph)], ci)
        if neutral !== nothing
            JuMP.add_to_expression!(ctx.kcl_r[(bus, neutral)], -cr)
            JuMP.add_to_expression!(ctx.kcl_i[(bus, neutral)], -ci)
        end
    end
    JuMP.@constraint(m, P == pe)      # export the allocated active power …
    JuMP.@constraint(m, Q == 0.0)     # … at unity power factor
    return pe
end

const _FAIRNESS = (:equal, :sum, :proportional)

"""
    solve_operating_envelope(nets, connection_points; fairness=:equal, kwargs...)
        -> OperatingEnvelopeResult

Compute a dynamic operating envelope: for each interval (each net in `nets`),
allocate an active-power export limit to every connection point such that all
points can export simultaneously without violating the network's operational
constraints (the voltage/thermal bounds the net declares).

Each interval is an independent constrained-allocation solve on
[`BMOPFTools.solve_opf`](https://github.com/frederikgeth/BMOPFTools.jl): a
`model_hook!` stamps a free export at each connection point and sets the fairness
objective, so the engine's voltage and thermal limits bound the result.

# Arguments
- `nets::Vector` — one BMOPFTools net dict per interval, carrying that interval's
  baseline loads and its operational limits (`v_min`/`v_max`, line `i_max`, …).
- `connection_points::Vector{ConnectionPoint}` — the DERs to allocate to.

# Keywords
- `fairness=:equal` — allocation rule:
  - `:equal` maximises a common export level assigned to every point (equitable;
    the level is capped by the weakest point and the tightest constraint);
  - `:sum` maximises the total allocated export (efficient; points with more
    network headroom get more, and a weak point may receive ≈0);
  - `:proportional` maximises `sum(log(pₑ))` over the points (proportional /
    Nash–Kelly fairness) — a middle ground where no point is starved but a point
    with more headroom still gets more.
- `per_unit=true`, `s_base=1e6`, `optimizer=Ipopt.Optimizer`, `verbose=false`,
  `solver_options=()` — solver control, as in `solve_opf`.

# Returns
An [`OperatingEnvelopeResult`](@ref) with each connection point's allocated
export (W) per interval.

Reactive power is held at zero (unity-PF export); import envelopes and
power-factor-flexible allocation are natural extensions.
"""
function solve_operating_envelope(nets::AbstractVector,
                                  connection_points::AbstractVector{ConnectionPoint};
                                  fairness::Symbol=:equal,
                                  per_unit::Bool=true,
                                  s_base::Float64=1e6,
                                  optimizer=Ipopt.Optimizer,
                                  verbose::Bool=false,
                                  solver_options=())
    fairness in _FAIRNESS || throw(ArgumentError(
        "fairness must be one of $(_FAIRNESS), got :$fairness"))
    T = length(nets)
    T >= 1 || throw(ArgumentError("need at least one interval"))
    cps = collect(connection_points)
    ids = [cp.id for cp in cps]
    allunique(ids) || throw(ArgumentError("connection-point ids must be unique: $ids"))

    envelope = Dict{String,Vector{Float64}}(id => Vector{Float64}(undef, T) for id in ids)
    total_export = Vector{Float64}(undef, T)
    statuses = Vector{String}(undef, T)
    snapshots = Vector{Dict{String,Any}}(undef, T)

    for t in 1:T
        pe = Dict{String,Any}()
        doe_hook! = ctx -> begin
            m = ctx.model
            for cp in cps
                pe[cp.id] = _stamp_export_port!(ctx, cp)
            end
            # Replace the generation-cost objective with the fairness objective.
            if fairness == :equal
                level = JuMP.@variable(m, base_name = "doe_level", lower_bound = 0.0)
                for cp in cps
                    JuMP.@constraint(m, pe[cp.id] == level)
                end
                JuMP.@objective(m, Max, level)
            elseif fairness == :sum
                JuMP.@objective(m, Max, sum(pe[cp.id] for cp in cps))
            else # :proportional — maximise Σ log(pe): no point is starved, but a
                 # point with more headroom still receives more (Nash/Kelly fairness).
                for cp in cps
                    JuMP.set_start_value(pe[cp.id], 1.0)   # keep the log argument off 0
                end
                JuMP.@objective(m, Max, sum(log(pe[cp.id]) for cp in cps))
            end
        end
        env_hook! = (ctx, result) -> begin
            sb = _sbase(ctx)
            result["operating_envelope"] =
                Dict(cp.id => JuMP.value(pe[cp.id]) * sb for cp in cps)
        end

        res = solve_opf(nets[t]; model_hook! = doe_hook!, solution_hook! = env_hook!,
                        per_unit=per_unit, s_base=s_base, optimizer=optimizer,
                        verbose=verbose, solver_options=solver_options)

        snapshots[t] = res
        statuses[t] = res["termination_status"]
        alloc = res["operating_envelope"]
        for cp in cps
            envelope[cp.id][t] = alloc[cp.id]
        end
        total_export[t] = sum(alloc[cp.id] for cp in cps)
    end

    return OperatingEnvelopeResult(statuses, envelope, total_export, snapshots)
end

# Single-interval convenience: a lone net becomes a one-element horizon.
function solve_operating_envelope(net::Dict{String,Any},
                                  connection_points::AbstractVector{ConnectionPoint};
                                  kwargs...)
    solve_operating_envelope([net], connection_points; kwargs...)
end
