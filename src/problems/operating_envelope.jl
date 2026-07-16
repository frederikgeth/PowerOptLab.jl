# Dynamic operating envelopes (DOEs).
#
# The implementation deliberately distinguishes a feasible allocation at the
# simultaneous bound (`security=:bound_point`) from an envelope whose corners
# have all been represented in the AC model (`security=:corners`).  The latter
# is still a local nonlinear-OPF result, not a proof that a non-convex feasible
# region has no interior holes; this limitation is carried into the result
# diagnostics instead of being hidden behind the word "robust".

"""
    ConnectionPoint(; id, bus, export_max=0.0, import_max=0.0, kwargs...)

An active connection point whose positive operating-envelope capacity is
calculated by [`solve_operating_envelope`](@ref).

- `export_max` / `import_max` are the connection's active-power nameplate limits
  (W, both non-negative). Select which one is used with `direction` on the solve.
- `ibr_id=nothing` retains the lightweight legacy port: aggregate unity-PF power
  is stamped directly at `bus` using `phase_terminals` and `neutral`.
- `ibr_id="pv1"` binds the envelope to an existing BMOPFTools `ibr`. This is the
  recommended representation for PV and batteries: the IBR keeps its prescribed
  Volt-VAr/Volt-Watt or fixed-power-factor control law, apparent/current limits,
  phase topology, and DC coupling while the DOE controls active power only.
- `requested` is an optional requested/forecast capacity (W), used by
  `FairnessPolicy(normalization=:request)`.
- `normalization` is an optional custom fairness reference (W), used by
  `FairnessPolicy(normalization=:custom)`.
"""
Base.@kwdef struct ConnectionPoint
    id::String
    bus::String
    phase_terminals::Vector{String} = ["1"]
    neutral::Union{String,Nothing} = "n"
    export_max::Float64 = 0.0
    import_max::Float64 = 0.0
    ibr_id::Union{String,Nothing} = nothing
    requested::Union{Float64,Nothing} = nothing
    normalization::Union{Float64,Nothing} = nothing
end

"""
    FairnessPolicy(; kind=:equal, normalization=:none, weights=Dict(),
                     alpha=1.0, epsilon=1e-6)

Parameterized policy for allocating active-power envelope capacity.

`normalization` defines the reference in `xᵢ = capacityᵢ/referenceᵢ`:

- `:none` — absolute watts (the reference is 1 W);
- `:capacity` — the connection's export/import nameplate for this direction;
- `:request` — `ConnectionPoint.requested`;
- `:custom` — `ConnectionPoint.normalization`.

Supported `kind` values:

- `:equal` — require equal normalized allocations and maximize their level;
- `:max_total` — maximize the weighted sum of normalized allocations;
- `:proportional` — weighted proportional fairness, `Σwᵢ log(xᵢ+ε)`;
- `:alpha` — weighted alpha fairness (`alpha=0` is weighted sum and `alpha=1`
  is proportional fairness);
- `:max_min` — maximize the minimum normalized allocation, then maximize total
  allocation while retaining that locally optimal minimum;
- `:equal_curtailment` — require equal normalized curtailment from nameplate and
  minimize it.

`weights` is keyed by connection-point id and defaults to one. All weights must
be finite and strictly positive. `epsilon` regularizes logarithmic/negative-power
utilities at zero; it does not create physical capacity.
"""
Base.@kwdef struct FairnessPolicy
    kind::Symbol = :equal
    normalization::Symbol = :none
    weights::Dict{String,Float64} = Dict{String,Float64}()
    alpha::Float64 = 1.0
    epsilon::Float64 = 1e-6
end

"""
    OperatingEnvelopeResult

Result of [`solve_operating_envelope`](@ref). Capacities are positive SI watts
for the selected `direction`.

`snapshots[t]` is the first scenario at the all-active upper corner. When an
interval has no feasible primal solution it contains only status metadata and
all capacities for that interval are `NaN`; infeasible solver iterates are never
published as envelopes.

`total_export` is retained as a backward-compatible alias of `total_capacity`.
For a new import study use `total_capacity` and inspect `direction`.

`fairness_metrics` reports allocation and curtailment metrics for the published
capacity at each interval. `schedule` records the issue/validity metadata and
whether a non-optimised fallback was published.
"""
struct OperatingEnvelopeResult
    termination_status::Vector{String}
    envelope::Dict{String,Vector{Float64}}
    total_export::Vector{Float64}
    snapshots::Vector{Dict{String,Any}}
    direction::Symbol
    total_capacity::Vector{Float64}
    diagnostics::Vector{Dict{String,Any}}
    fairness_metrics::Vector{Dict{String,Any}}
    schedule::Vector{Dict{String,Any}}
end

# Source compatibility for callers that constructed the original four-field
# result directly. New solves always populate the richer fields above.
OperatingEnvelopeResult(status, envelope, total, snapshots) =
    OperatingEnvelopeResult(status, envelope, copy(total), snapshots, :export,
                            total, Dict{String,Any}[], Dict{String,Any}[],
                            Dict{String,Any}[])

"""
    OperatingEnvelopeVerification

Result of [`verify_operating_envelope`](@ref). Each interval reports whether a
fixed, already-issued allocation was locally feasible at every requested
utilisation point and forecast/model scenario. This is a verification result,
not a new allocation.
"""
struct OperatingEnvelopeVerification
    termination_status::Vector{String}
    feasible::Vector{Bool}
    snapshots::Vector{Dict{String,Any}}
    diagnostics::Vector{Dict{String,Any}}
end

const _FAIRNESS_KINDS =
    (:equal, :max_total, :proportional, :alpha, :max_min, :equal_curtailment)
const _NORMALIZATIONS = (:none, :capacity, :request, :custom)
const _SECURITY_MODES = (:bound_point, :corners)
const _DIRECTIONS = (:export, :import)

_power_base(ctx) = ctx.bases === nothing ? 1.0 : Float64(ctx.bases.s_base)
_capacity_limit(cp::ConnectionPoint, direction::Symbol) =
    direction == :export ? cp.export_max : cp.import_max

function _as_policy(fairness)
    fairness isa FairnessPolicy && return fairness
    fairness isa Symbol || throw(ArgumentError(
        "fairness must be a Symbol or FairnessPolicy, got $(typeof(fairness))"))
    fairness == :sum          && return FairnessPolicy(kind=:max_total)
    fairness == :proportional && return FairnessPolicy(kind=:proportional)
    fairness == :equal        && return FairnessPolicy(kind=:equal)
    fairness in _FAIRNESS_KINDS && return FairnessPolicy(kind=fairness)
    throw(ArgumentError("unknown fairness policy :$fairness"))
end

function _scenario_groups(nets)
    nets isa Dict{String,Any} && return [[nets]]
    nets isa AbstractVector || throw(ArgumentError(
        "nets must be a network Dict, a vector of interval networks, or a vector of scenario vectors"))
    isempty(nets) && throw(ArgumentError("need at least one interval"))
    if all(n -> n isa Dict{String,Any}, nets)
        return [[n] for n in nets]
    elseif all(g -> g isa AbstractVector && !isempty(g) &&
                    all(n -> n isa Dict{String,Any}, g), nets)
        return [collect(g) for g in nets]
    end
    throw(ArgumentError(
        "nets must contain only network Dicts or only non-empty vectors of network Dicts"))
end

function _validate_policy(policy::FairnessPolicy, cps, direction)
    policy.kind in _FAIRNESS_KINDS || throw(ArgumentError(
        "fairness kind must be one of $(_FAIRNESS_KINDS), got :$(policy.kind)"))
    policy.normalization in _NORMALIZATIONS || throw(ArgumentError(
        "normalization must be one of $(_NORMALIZATIONS), got :$(policy.normalization)"))
    isfinite(policy.alpha) || throw(ArgumentError("fairness alpha must be finite"))
    isfinite(policy.epsilon) && policy.epsilon > 0 || throw(ArgumentError(
        "fairness epsilon must be finite and > 0"))
    unknown = setdiff(Set(keys(policy.weights)), Set(cp.id for cp in cps))
    isempty(unknown) || throw(ArgumentError("fairness weights contain unknown ids: $(collect(unknown))"))
    for cp in cps
        w = get(policy.weights, cp.id, 1.0)
        isfinite(w) && w > 0 || throw(ArgumentError(
            "fairness weight for '$(cp.id)' must be finite and > 0"))
        _fairness_reference(cp, policy.normalization, direction)
    end
end

function _fairness_reference(cp::ConnectionPoint, normalization::Symbol,
                             direction::Symbol)
    ref = if normalization == :none
        1.0
    elseif normalization == :capacity
        _capacity_limit(cp, direction)
    elseif normalization == :request
        cp.requested === nothing && throw(ArgumentError(
            "connection '$(cp.id)' needs requested for normalization=:request"))
        cp.requested
    else
        cp.normalization === nothing && throw(ArgumentError(
            "connection '$(cp.id)' needs normalization for normalization=:custom"))
        cp.normalization
    end
    isfinite(ref) && ref > 0 || throw(ArgumentError(
        "fairness reference for '$(cp.id)' must be finite and > 0, got $ref"))
    return Float64(ref)
end

function _validate_connection_points(groups, cps, policy, direction, security,
                                     max_exact_corners)
    direction in _DIRECTIONS || throw(ArgumentError(
        "direction must be one of $(_DIRECTIONS), got :$direction"))
    security in _SECURITY_MODES || throw(ArgumentError(
        "security must be one of $(_SECURITY_MODES), got :$security"))
    isempty(cps) && throw(ArgumentError("need at least one connection point"))
    ids = [cp.id for cp in cps]
    allunique(ids) || throw(ArgumentError("connection-point ids must be unique: $ids"))
    for cp in cps
        isempty(cp.id) && throw(ArgumentError("connection-point id must not be empty"))
        isempty(cp.bus) && throw(ArgumentError("connection '$(cp.id)' bus must not be empty"))
        isempty(cp.phase_terminals) && throw(ArgumentError(
            "connection '$(cp.id)' needs at least one phase terminal"))
        allunique(cp.phase_terminals) || throw(ArgumentError(
            "connection '$(cp.id)' phase terminals must be unique"))
        cp.neutral in cp.phase_terminals && throw(ArgumentError(
            "connection '$(cp.id)' neutral cannot also be a phase terminal"))
        for (name, value) in (("export_max", cp.export_max), ("import_max", cp.import_max))
            isfinite(value) && value >= 0 || throw(ArgumentError(
                "connection '$(cp.id)' $name must be finite and >= 0, got $value"))
        end
        for (name, value) in (("requested", cp.requested),
                              ("normalization", cp.normalization))
            value === nothing || (isfinite(value) && value > 0) || throw(ArgumentError(
                "connection '$(cp.id)' $name must be finite and > 0"))
        end
    end
    all(_capacity_limit(cp, direction) == 0 for cp in cps) && throw(ArgumentError(
        "all connection points have zero $(direction) capacity"))
    max_exact_corners >= 1 || throw(ArgumentError("max_exact_corners must be >= 1"))
    security == :corners && length(cps) > max_exact_corners && throw(ArgumentError(
        "security=:corners needs 2^N AC contexts; got N=$(length(cps)) > " *
        "max_exact_corners=$max_exact_corners"))
    _validate_policy(policy, cps, direction)

    for (t, group) in enumerate(groups), (s, net) in enumerate(group)
        buses = get(net, "bus", Dict())
        for cp in cps
            haskey(buses, cp.bus) || throw(ArgumentError(
                "interval $t scenario $s: connection '$(cp.id)' bus '$(cp.bus)' not found"))
            terminals = Set(String.(get(buses[cp.bus], "terminal_names", String[])))
            if cp.ibr_id === nothing
                for term in cp.phase_terminals
                    term in terminals || throw(ArgumentError(
                        "interval $t scenario $s: terminal '$term' for '$(cp.id)' not found at bus '$(cp.bus)'"))
                end
                cp.neutral === nothing || cp.neutral in terminals || throw(ArgumentError(
                    "interval $t scenario $s: neutral '$(cp.neutral)' for '$(cp.id)' not found"))
            else
                invs = get(net, "ibr", Dict())
                haskey(invs, cp.ibr_id) || throw(ArgumentError(
                    "interval $t scenario $s: IBR '$(cp.ibr_id)' for '$(cp.id)' not found"))
                inv = invs[cp.ibr_id]
                get(inv, "bus", nothing) == cp.bus || throw(ArgumentError(
                    "interval $t scenario $s: IBR '$(cp.ibr_id)' is not at bus '$(cp.bus)'"))
                topo = uppercase(String(get(inv, "topology", "FOUR_LEG")))
                topo in ("SINGLE_PHASE", "FOUR_LEG") || throw(ArgumentError(
                    "connection-bound IBR '$(cp.ibr_id)' topology '$topo' is not supported; " *
                    "use SINGLE_PHASE or FOUR_LEG for prescribed Q-V control"))
            end
        end
    end
end

# Legacy connection port. It is intentionally aggregate unity-PF and exists for
# backward compatibility and simple teaching examples. Real PV/battery studies
# should bind `ConnectionPoint.ibr_id` to the engine's prescribed-control model.
function _stamp_legacy_port!(ctx, cp::ConnectionPoint)
    m = ctx.model
    bus = cp.bus
    P = zero(JuMP.QuadExpr)
    Q = zero(JuMP.QuadExpr)
    for ph in cp.phase_terminals
        cr = JuMP.@variable(m, base_name="cr_doe_$(cp.id)_$(ph)")
        ci = JuMP.@variable(m, base_name="ci_doe_$(cp.id)_$(ph)")
        dvr, dvi = _dv(ctx, bus, ph, cp.neutral)
        P += JuMP.@expression(m, dvr*cr + dvi*ci)
        Q += JuMP.@expression(m, dvi*cr - dvr*ci)
        JuMP.add_to_expression!(ctx.kcl_r[(bus, ph)], cr)
        JuMP.add_to_expression!(ctx.kcl_i[(bus, ph)], ci)
        if cp.neutral !== nothing
            JuMP.add_to_expression!(ctx.kcl_r[(bus, cp.neutral)], -cr)
            JuMP.add_to_expression!(ctx.kcl_i[(bus, cp.neutral)], -ci)
        end
    end
    JuMP.@constraint(m, Q == 0.0)
    return P
end

# Recover the active-power expression already stamped by a BMOPFTools IBR. This
# adds no reactive decision: the engine's own constant-PF / Volt-VAr equality is
# retained unchanged.
function _ibr_active_power(ctx, cp::ConnectionPoint)
    m = ctx.model
    inv = ctx.net["ibr"][cp.ibr_id]
    bus = String(inv["bus"])
    tm = String.(inv["terminal_map"])
    topo = uppercase(String(get(inv, "topology", "FOUR_LEG")))
    vr = ctx.vars[:vr]; vi = ctx.vars[:vi]
    cri = ctx.vars[:cri]; cii = ctx.vars[:cii]

    if topo == "SINGLE_PHASE"
        ph, ref = tm[1], tm[2]
        dvr = JuMP.@expression(m, vr[(bus,ph)] - vr[(bus,ref)])
        dvi = JuMP.@expression(m, vi[(bus,ph)] - vi[(bus,ref)])
        return JuMP.@expression(m,
            dvr*cri[(cp.ibr_id,1)] + dvi*cii[(cp.ibr_id,1)])
    end

    neutral = cp.neutral !== nothing && cp.neutral in tm ? cp.neutral : tm[end]
    phases = [term for term in tm if term != neutral]
    terms = JuMP.QuadExpr[]
    for (idx, ph) in enumerate(phases)
        dvr = JuMP.@expression(m, vr[(bus,ph)] - vr[(bus,neutral)])
        dvi = JuMP.@expression(m, vi[(bus,ph)] - vi[(bus,neutral)])
        push!(terms, JuMP.@expression(m,
            dvr*cri[(cp.ibr_id,idx)] + dvi*cii[(cp.ibr_id,idx)]))
    end
    return sum(terms)
end

_connection_active_power(ctx, cp) =
    cp.ibr_id === nothing ? _stamp_legacy_port!(ctx, cp) : _ibr_active_power(ctx, cp)

function _dispatch_patterns(n::Int, security::Symbol)
    security == :bound_point && return [ones(Float64, n)]
    return [[Float64((mask >> (i-1)) & 1) for i in 1:n]
            for mask in 0:(Int(2)^n - 1)]
end

function _set_fairness_objective!(model, cap, cps, policy, direction, power_base;
                                  temporal_history=nothing, temporal_dt_h=1.0)
    refs = Dict(cp.id => _fairness_reference(cp, policy.normalization, direction) /
                              power_base for cp in cps)
    x = Dict(cp.id => JuMP.@expression(model, cap[cp.id] / refs[cp.id]) for cp in cps)
    weights = Dict(cp.id => get(policy.weights, cp.id, 1.0) for cp in cps)
    limits = Dict(cp.id => _capacity_limit(cp, direction) / power_base for cp in cps)

    if temporal_history !== nothing
        level = JuMP.@variable(model, base_name="doe_cumulative_fairness", lower_bound=0.0)
        for cp in cps
            prior = get(temporal_history, cp.id, 0.0)
            JuMP.@constraint(model, prior + temporal_dt_h * x[cp.id] / weights[cp.id] >= level)
        end
        JuMP.@objective(model, Max, level)
        return (kind=:cumulative_max_min, level=level, x=x, weights=weights)
    elseif policy.kind == :equal
        level = JuMP.@variable(model, base_name="doe_equal_level", lower_bound=0.0)
        for cp in cps
            JuMP.@constraint(model, x[cp.id] == level)
        end
        JuMP.@objective(model, Max, level)
        return (kind=:single_stage, level=nothing)
    elseif policy.kind == :max_total
        JuMP.@objective(model, Max, sum(weights[cp.id] * x[cp.id] for cp in cps))
        return (kind=:single_stage, level=nothing)
    elseif policy.kind in (:proportional, :alpha)
        α = policy.kind == :proportional ? 1.0 : policy.alpha
        ε = policy.epsilon
        if isapprox(α, 1.0; atol=1e-12, rtol=0.0)
            JuMP.@objective(model, Max,
                sum(weights[cp.id] * log(x[cp.id] + ε) for cp in cps))
        elseif isapprox(α, 0.0; atol=1e-12, rtol=0.0)
            JuMP.@objective(model, Max,
                sum(weights[cp.id] * x[cp.id] for cp in cps))
        else
            JuMP.@objective(model, Max,
                sum(weights[cp.id] * (x[cp.id] + ε)^(1.0-α) / (1.0-α)
                    for cp in cps))
        end
        return (kind=:single_stage, level=nothing)
    elseif policy.kind == :max_min
        level = JuMP.@variable(model, base_name="doe_max_min_level", lower_bound=0.0)
        for cp in cps
            JuMP.@constraint(model, x[cp.id] >= level)
        end
        JuMP.@objective(model, Max, level)
        return (kind=:max_min, level=level, x=x, weights=weights)
    else
        level = JuMP.@variable(model, base_name="doe_curtailment", lower_bound=0.0)
        for cp in cps
            JuMP.@constraint(model, (limits[cp.id] - cap[cp.id]) / refs[cp.id] == level)
        end
        JuMP.@objective(model, Min, level)
        return (kind=:single_stage, level=nothing)
    end
end

_has_primal(model) = _publishable(_solve_outcome(model))

function _result_margins(result, net)
    best = Dict{String,Tuple{Float64,String}}()
    consider!(kind, margin, label) = begin
        isfinite(margin) || return
        if !haskey(best, kind) || margin < best[kind][1]
            best[kind] = (margin, label)
        end
    end

    # Phase-to-ground voltage bounds. More specialized BMOPFTools voltage
    # constraints remain enforced even when not reducible to one scalar margin.
    for (bus_id, bus) in get(net, "bus", Dict())
        rb = get(get(result, "bus", Dict()), bus_id, nothing)
        rb isa Dict || continue
        terminals = String.(get(bus, "terminal_names", String[]))
        grounded = Set(String.(get(bus, "perfectly_grounded_terminals", String[])))
        neutral_term = findfirst(t -> lowercase(t) in ("n", "neutral"), terminals)
        phases = [t for (idx, t) in enumerate(terminals)
                  if !(t in grounded) && idx != neutral_term]
        for (field, sense) in (("v_min", :lower), ("v_max", :upper))
            limits = get(bus, field, nothing)
            limits isa AbstractVector || continue
            for (idx, limit) in enumerate(limits)
                idx <= length(phases) || break
                term = phases[idx]
                haskey(rb, term) || continue
                vm = Float64(rb[term]["vm"])
                margin = sense == :lower ? vm - Float64(limit) : Float64(limit) - vm
                consider!("voltage", margin, "bus:$bus_id:$term:$field")
            end
        end

        vneg_max = get(bus, "vneg_max", nothing)
        if vneg_max isa Number && length(phases) == 3
            vn = if neutral_term === nothing || !haskey(rb, terminals[neutral_term])
                0.0 + 0.0im
            else
                rn = rb[terminals[neutral_term]]
                Float64(rn["vr"]) + im*Float64(rn["vi"])
            end
            V = ComplexF64[]
            for term in phases
                rt = rb[term]
                push!(V, Float64(rt["vr"]) + im*Float64(rt["vi"]) - vn)
            end
            a = cis(2pi/3)
            V2 = (V[1] + a^2*V[2] + a*V[3]) / 3
            consider!("negative_sequence", Float64(vneg_max) - abs(V2),
                      "bus:$bus_id:vneg_max")
        end
    end

    # Per-conductor line ampacity inherited from the referenced linecode.
    for (line_id, line) in get(net, "line", Dict())
        rl = get(get(result, "line", Dict()), line_id, nothing)
        rl isa Dict || continue
        lc = get(get(net, "linecode", Dict()), get(line, "linecode", ""), Dict())
        raw = get(line, "i_max", get(lc, "i_max", nothing))
        raw === nothing && continue
        limits = raw isa Number ? fill(Float64(raw), length(get(line, "terminal_map_from", []))) :
                                  Float64.(raw)
        terms = String.(get(line, "terminal_map_from", String[]))
        for (idx, limit) in enumerate(limits)
            idx <= length(terms) || break
            term = terms[idx]
            haskey(rl, term) || continue
            current = max(Float64(get(rl[term], "cm_fr", 0.0)),
                          Float64(get(rl[term], "cm_to", 0.0)))
            consider!("thermal", limit - current, "line:$line_id:$term:i_max")
        end
    end
    return best
end

function _merge_margins(results_and_nets)
    worst = Dict{String,Tuple{Float64,String}}()
    for (result, net) in results_and_nets
        for (kind, item) in _result_margins(result, net)
            if !haskey(worst, kind) || item[1] < worst[kind][1]
                worst[kind] = item
            end
        end
    end
    tolerances = Dict("voltage"=>0.05, "thermal"=>0.01,
                      "negative_sequence"=>0.01)
    binding = [item[2] for (kind, item) in worst
               if item[1] <= get(tolerances, kind, 0.0)]
    sort!(binding)
    return Dict{String,Any}(
        "minimum_margins" => Dict(kind=>item[1] for (kind, item) in worst),
        "minimum_margin_locations" => Dict(kind=>item[2] for (kind, item) in worst),
        "binding_constraints" => binding)
end

function _optimize_fairness!(model, stage; max_min_tolerance)
    JuMP.optimize!(model)
    stage.kind in (:max_min, :cumulative_max_min) || return
    _has_primal(model) || return
    best = JuMP.value(stage.level)
    JuMP.@constraint(model, stage.level >= best - max_min_tolerance)
    JuMP.@objective(model, Max,
        sum(stage.weights[id] * stage.x[id] for id in keys(stage.x)))
    JuMP.optimize!(model)
end

function _fairness_metrics(alloc, cps, policy, direction; cumulative=nothing)
    normalized = Dict{String,Float64}()
    curtailment = Dict{String,Float64}()
    for cp in cps
        value = Float64(alloc[cp.id])
        ref = _fairness_reference(cp, policy.normalization, direction)
        normalized[cp.id] = value / ref
        limit = _capacity_limit(cp, direction)
        curtailment[cp.id] = limit > 0 ? 1 - value / limit : 0.0
    end
    values_ = collect(values(normalized))
    denominator = length(values_) * sum(x^2 for x in values_)
    # Roundoff can push the mathematically bounded index marginally above one
    # (for example 1.0000000000000002 for equal allocations).
    jain = denominator > 0 ? clamp(sum(values_)^2 / denominator, 0.0, 1.0) : 1.0
    out = Dict{String,Any}(
        "total_capacity_W" => sum(values(alloc)),
        "normalized_allocations" => normalized,
        "curtailment_fraction" => curtailment,
        "jain_index" => jain,
        "min_normalized" => minimum(values_),
        "max_normalized" => maximum(values_),
        "mean_normalized" => sum(values_) / length(values_))
    cumulative === nothing || (out["cumulative_normalized"] = copy(cumulative))
    return out
end

function _validate_temporal_fairness(mode, history, dt_h, cps)
    mode in (:none, :cumulative_max_min) || throw(ArgumentError(
        "temporal_fairness must be :none or :cumulative_max_min"))
    isfinite(dt_h) && dt_h > 0 || throw(ArgumentError("temporal_dt_h must be finite and > 0"))
    unknown = setdiff(Set(keys(history)), Set(cp.id for cp in cps))
    isempty(unknown) || throw(ArgumentError("fairness_history contains unknown ids: $(collect(unknown))"))
    for (id, value) in history
        isfinite(value) && value >= 0 || throw(ArgumentError(
            "fairness_history for '$id' must be finite and >= 0"))
    end
end

function _capacity_trajectory(capacities, cps, T, direction)
    source = capacities isa OperatingEnvelopeResult ? capacities.envelope : capacities
    source isa AbstractDict || throw(ArgumentError(
        "capacities must be an OperatingEnvelopeResult or a dictionary keyed by connection-point id"))
    result = Vector{Dict{String,Float64}}(undef, T)
    for t in 1:T
        item = Dict{String,Float64}()
        for cp in cps
            haskey(source, cp.id) || throw(ArgumentError("capacities missing id '$(cp.id)'"))
            raw = source[cp.id]
            value = raw isa AbstractVector ? (length(raw) == T || throw(ArgumentError(
                "capacity vector for '$(cp.id)' must have $T entries")); raw[t]) : raw
            value isa Number && isfinite(value) && value >= 0 || throw(ArgumentError(
                "capacity for '$(cp.id)' at interval $t must be finite and >= 0"))
            value <= _capacity_limit(cp, direction) + 1e-8 || throw(ArgumentError(
                "capacity for '$(cp.id)' exceeds its declared $(direction) limit"))
            item[cp.id] = Float64(value)
        end
        result[t] = item
    end
    return result
end

function _verification_patterns(utilizations, n)
    utilizations == :bound_point && return _dispatch_patterns(n, :bound_point)
    utilizations == :corners && return _dispatch_patterns(n, :corners)
    utilizations isa AbstractVector || throw(ArgumentError(
        "utilizations must be :bound_point, :corners, or a vector of utilization vectors"))
    patterns = Vector{Vector{Float64}}()
    for pattern in utilizations
        pattern isa AbstractVector && length(pattern) == n || throw(ArgumentError(
            "each utilization vector must contain one value per connection point"))
        values_ = Float64.(pattern)
        all(x -> isfinite(x) && 0 <= x <= 1, values_) || throw(ArgumentError(
            "utilization values must be finite and lie in [0, 1]"))
        push!(patterns, values_)
    end
    isempty(patterns) && throw(ArgumentError("need at least one utilization point"))
    return patterns
end

function _solve_interval_group(group, cps, policy;
                               direction, security, per_unit, s_base,
                               optimizer, verbose, solver_options,
                               volt_var_watt_eps, max_min_tolerance,
                               temporal_history=nothing,
                               temporal_dt_h=1.0,
                               fixed_capacity=nothing,
                               patterns_override=nothing)
    model = JuMP.Model(optimizer)
    pb = per_unit ? s_base : 1.0
    cap = Dict{String,Any}()
    for cp in cps
        upper = _capacity_limit(cp, direction) / pb
        if fixed_capacity === nothing
            cap[cp.id] = JuMP.@variable(model, base_name="doe_capacity_$(cp.id)",
                lower_bound=0.0, upper_bound=upper)
        else
            value = fixed_capacity[cp.id] / pb
            cap[cp.id] = JuMP.@variable(model, base_name="doe_capacity_$(cp.id)",
                lower_bound=value, upper_bound=value)
        end
    end

    patterns = patterns_override === nothing ? _dispatch_patterns(length(cps), security) : patterns_override
    sign = direction == :export ? 1.0 : -1.0
    specifications = [(net=net, scenario=scenario_index,
                       pattern=pattern_index, fractions=fractions)
                      for (scenario_index, net) in enumerate(group)
                      for (pattern_index, fractions) in enumerate(patterns)]
    hook_factory = context_index -> begin
        fractions = specifications[context_index].fractions
        ctx -> begin
            for (i, cp) in enumerate(cps)
                p = _connection_active_power(ctx, cp)
                JuMP.@constraint(ctx.model,
                    p == sign * fractions[i] * cap[cp.id])
            end
        end
    end
    multi = build_multi_context([spec.net for spec in specifications]; model,
        hook_factory, per_unit, s_base, optimizer, verbose, solver_options,
        context_options=(volt_var_watt_eps=volt_var_watt_eps, verbose=verbose))
    records = [(ctx=multi.contexts[index], net=spec.net,
                scenario=spec.scenario, pattern=spec.pattern,
                fractions=spec.fractions)
               for (index, spec) in enumerate(specifications)]
    foreach(r -> enforce_kcl!(r.ctx), records)
    if fixed_capacity === nothing
        stage = _set_fairness_objective!(model, cap, cps, policy, direction, pb;
            temporal_history=temporal_history, temporal_dt_h=temporal_dt_h)
        _optimize_fairness!(model, stage; max_min_tolerance=max_min_tolerance)
    else
        JuMP.@objective(model, Min, 0.0)
        JuMP.optimize!(model)
    end

    outcome = _solve_outcome(model)
    status = string(outcome.termination_status)
    primal = string(outcome.primal_status)
    feasible = _publishable(outcome)
    total = NaN
    alloc = Dict(cp.id => NaN for cp in cps)
    snapshot = Dict{String,Any}("termination_status"=>status,
                                "primal_status"=>primal)
    objective = NaN
    margin_diagnostics = Dict{String,Any}(
        "minimum_margins"=>Dict{String,Float64}(),
        "minimum_margin_locations"=>Dict{String,String}(),
        "binding_constraints"=>String[])
    if feasible
        for cp in cps
            raw = JuMP.value(cap[cp.id]) * pb
            alloc[cp.id] = clamp(raw, 0.0, _capacity_limit(cp, direction))
        end
        total = sum(values(alloc))
        objective = JuMP.objective_value(model)
        representative_index = findfirst(r -> r.scenario == 1 && all(isone, r.fractions), records)
        representative = records[something(representative_index, 1)]
        snapshot = extract_result(representative.ctx)
        checked = Tuple{Dict{String,Any},Dict{String,Any}}[]
        for record in records
            result = record.scenario == 1 && all(isone, record.fractions) ?
                     snapshot : extract_result(record.ctx)
            push!(checked, (result, record.net))
        end
        margin_diagnostics = _merge_margins(checked)
    end

    diag = Dict{String,Any}(
        "feasible" => feasible,
        "primal_status" => primal,
        "objective" => objective,
        "direction" => direction,
        "security" => security,
        "security_scope" => security == :bound_point ?
            :simultaneous_upper_bound_only : :all_box_corners,
        "guarantee" => :local_ac_feasibility_at_tested_dispatches,
        "scenario_count" => length(group),
        "dispatch_points_per_scenario" => length(patterns),
        "fairness_kind" => policy.kind,
        "normalization" => policy.normalization,
        "verification" => fixed_capacity !== nothing,
        "temporal_fairness" => temporal_history === nothing ? :none : :cumulative_max_min)
    merge!(diag, margin_diagnostics)
    return (status=status, alloc=alloc, total=total,
            snapshot=snapshot, diagnostics=diag)
end

"""
    solve_operating_envelope(nets, connection_points; kwargs...)
        -> OperatingEnvelopeResult

Calculate active-power operating-envelope capacity for each interval.

`nets` accepts three shapes:

- one network `Dict` — one interval, one forecast scenario;
- `Vector{Dict}` — several intervals, one scenario each (backward compatible);
- `Vector{Vector{Dict}}` — several intervals, each containing one or more
  forecast/model scenarios. One capacity allocation is shared by every scenario
  in an interval.

Important keywords:

- `direction=:export` or `:import`; returned capacities are positive magnitudes;
- `fairness=:equal`, the legacy symbols `:sum` / `:proportional`, or a
  [`FairnessPolicy`](@ref);
- `security=:bound_point` enforces only simultaneous full utilisation;
- `security=:corners` embeds all `2^N` zero/full-utilisation corners for every
  scenario. It is deliberately capped by `max_exact_corners=10` and is reported
  as local AC feasibility at tested points, not a global robust certificate;
- `volt_var_watt_eps=2e-3` controls the engine's smooth approximation of
  mandatory IBR Volt-VAr/Volt-Watt curve corners.

Loads retain their known P/Q from each network snapshot. Connection-bound IBRs
retain their prescribed Q-V law. Other network devices, including BMOPFTools
STATCOM IBRs, remain available to the OPF; comparing otherwise identical nets
with and without a STATCOM therefore quantifies its impact on active-power DOEs.
Network voltage, phase-to-neutral, negative-sequence (`vneg_max`), branch thermal,
neutral-current, and device limits declared by BMOPFTools remain in force.
"""
function solve_operating_envelope(nets,
                                  connection_points::AbstractVector{ConnectionPoint};
                                  fairness=:equal,
                                  direction::Symbol=:export,
                                  security::Symbol=:bound_point,
                                  per_unit::Bool=true,
                                  s_base::Float64=1e6,
                                  optimizer=Ipopt.Optimizer,
                                  verbose::Bool=false,
                                  solver_options=(),
                                  volt_var_watt_eps::Float64=2e-3,
                                  max_exact_corners::Int=10,
                                  max_min_tolerance::Float64=1e-7,
                                  temporal_fairness::Symbol=:none,
                                  fairness_history::AbstractDict=Dict{String,Float64}(),
                                  temporal_dt_h::Float64=1.0,
                                  issued_at::Union{Nothing,DateTime}=nothing,
                                  interval_seconds::Float64=300.0,
                                  validity_seconds::Float64=interval_seconds,
                                  fallback::Symbol=:missing)
    groups = _scenario_groups(nets)
    cps = collect(connection_points)
    policy = _as_policy(fairness)
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be finite and > 0"))
    isfinite(volt_var_watt_eps) && volt_var_watt_eps > 0 || throw(ArgumentError(
        "volt_var_watt_eps must be finite and > 0"))
    isfinite(max_min_tolerance) && max_min_tolerance >= 0 || throw(ArgumentError(
        "max_min_tolerance must be finite and >= 0"))
    _validate_connection_points(groups, cps, policy, direction, security,
                                max_exact_corners)
    _validate_temporal_fairness(temporal_fairness, fairness_history, temporal_dt_h, cps)
    isfinite(interval_seconds) && interval_seconds > 0 || throw(ArgumentError(
        "interval_seconds must be finite and > 0"))
    isfinite(validity_seconds) && validity_seconds > 0 || throw(ArgumentError(
        "validity_seconds must be finite and > 0"))
    fallback in (:missing, :zero, :last_feasible) || throw(ArgumentError(
        "fallback must be :missing, :zero, or :last_feasible"))

    T = length(groups)
    ids = [cp.id for cp in cps]
    envelope = Dict(id => fill(NaN, T) for id in ids)
    total = fill(NaN, T)
    statuses = Vector{String}(undef, T)
    snapshots = Vector{Dict{String,Any}}(undef, T)
    diagnostics = Vector{Dict{String,Any}}(undef, T)
    metrics = Vector{Dict{String,Any}}(undef, T)
    schedule = Vector{Dict{String,Any}}(undef, T)
    history = Dict{String,Float64}(cp.id => Float64(get(fairness_history, cp.id, 0.0)) for cp in cps)
    last_feasible = nothing

    for t in 1:T
        solved = _solve_interval_group(groups[t], cps, policy;
            direction=direction, security=security, per_unit=per_unit,
            s_base=s_base, optimizer=optimizer, verbose=verbose,
            solver_options=solver_options, volt_var_watt_eps=volt_var_watt_eps,
            max_min_tolerance=max_min_tolerance,
            temporal_history=temporal_fairness == :cumulative_max_min ? history : nothing,
            temporal_dt_h=temporal_dt_h)
        statuses[t] = solved.status
        snapshots[t] = solved.snapshot
        diagnostics[t] = solved.diagnostics
        total[t] = solved.total
        publication_source = :optimized
        published = solved.alloc
        if !solved.diagnostics["feasible"]
            if fallback == :zero
                published = Dict(id => 0.0 for id in ids)
                total[t] = 0.0
                publication_source = :zero_fallback
            elseif fallback == :last_feasible && last_feasible !== nothing
                published = copy(last_feasible)
                total[t] = sum(values(published))
                publication_source = :last_feasible_fallback
            else
                publication_source = :missing
            end
            diagnostics[t]["fallback_network_safe"] = false
        else
            last_feasible = copy(published)
            if temporal_fairness == :cumulative_max_min
                for cp in cps
                    history[cp.id] += temporal_dt_h * published[cp.id] /
                        _fairness_reference(cp, policy.normalization, direction)
                end
            end
        end
        for id in ids
            envelope[id][t] = published[id]
        end
        metrics[t] = if publication_source == :missing
            Dict{String,Any}("available"=>false, "publication_source"=>publication_source)
        else
            outcome = _fairness_metrics(published, cps, policy, direction;
                cumulative=temporal_fairness == :cumulative_max_min ? history : nothing)
            outcome["available"] = true
            outcome["publication_source"] = publication_source
            outcome
        end
        valid_from = issued_at === nothing ? nothing :
            issued_at + Millisecond(round(Int, 1000 * interval_seconds * (t - 1)))
        valid_until = valid_from === nothing ? nothing :
            valid_from + Millisecond(round(Int, 1000 * validity_seconds))
        schedule[t] = Dict{String,Any}("interval_index"=>t, "issued_at"=>issued_at,
            "valid_from"=>valid_from, "valid_until"=>valid_until,
            "publication_source"=>publication_source)
    end

    OperatingEnvelopeResult(statuses, envelope, copy(total), snapshots,
                            direction, total, diagnostics, metrics, schedule)
end

"""
    compare_operating_envelope_policies(nets, connection_points, policies; kwargs...)

Solve the same DOE study under several fairness policies. `policies` is a
dictionary or vector of `label => fairness` pairs. The returned dictionary maps
each label to an [`OperatingEnvelopeResult`](@ref), whose `fairness_metrics`
make the capacity/fairness trade-off directly comparable.
"""
function compare_operating_envelope_policies(nets, connection_points, policies; kwargs...)
    entries = policies isa AbstractDict ? collect(pairs(policies)) : collect(policies)
    isempty(entries) && throw(ArgumentError("need at least one fairness policy"))
    result = Dict{String,OperatingEnvelopeResult}()
    for entry in entries
        entry isa Pair || throw(ArgumentError("policies must contain label => policy pairs"))
        label = string(first(entry))
        haskey(result, label) && throw(ArgumentError("duplicate policy label '$label'"))
        result[label] = solve_operating_envelope(nets, connection_points;
            fairness=last(entry), kwargs...)
    end
    return result
end

"""
    verify_operating_envelope(nets, connection_points, capacities; kwargs...)
        -> OperatingEnvelopeVerification

Check an already-issued capacity trajectory without optimising it. The active
power capacities are fixed and the normal network physics, including prescribed
Q-V IBR controls and any STATCOM model present in the network, is solved at
every requested scenario and utilisation point. `utilizations` is
`:bound_point`, `:corners`, or explicit vectors in `[0, 1]^N`.
"""
function verify_operating_envelope(nets,
                                   connection_points::AbstractVector{ConnectionPoint},
                                   capacities;
                                   direction::Symbol=:export,
                                   utilizations=:bound_point,
                                   per_unit::Bool=true,
                                   s_base::Float64=1e6,
                                   optimizer=Ipopt.Optimizer,
                                   verbose::Bool=false,
                                   solver_options=(),
                                   volt_var_watt_eps::Float64=2e-3,
                                   max_exact_corners::Int=10)
    groups = _scenario_groups(nets)
    cps = collect(connection_points)
    policy = FairnessPolicy(kind=:max_total)
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be finite and > 0"))
    isfinite(volt_var_watt_eps) && volt_var_watt_eps > 0 || throw(ArgumentError(
        "volt_var_watt_eps must be finite and > 0"))
    security = utilizations == :corners ? :corners : :bound_point
    _validate_connection_points(groups, cps, policy, direction, security, max_exact_corners)
    patterns = _verification_patterns(utilizations, length(cps))
    trajectory = _capacity_trajectory(capacities, cps, length(groups), direction)
    statuses = String[]
    feasible = Bool[]
    snapshots = Dict{String,Any}[]
    diagnostics = Dict{String,Any}[]
    for (t, group) in enumerate(groups)
        solved = _solve_interval_group(group, cps, policy;
            direction=direction, security=security, per_unit=per_unit, s_base=s_base,
            optimizer=optimizer, verbose=verbose, solver_options=solver_options,
            volt_var_watt_eps=volt_var_watt_eps, max_min_tolerance=1e-7,
            fixed_capacity=trajectory[t], patterns_override=patterns)
        push!(statuses, solved.status)
        push!(feasible, solved.diagnostics["feasible"])
        push!(snapshots, solved.snapshot)
        push!(diagnostics, solved.diagnostics)
    end
    return OperatingEnvelopeVerification(statuses, feasible, snapshots, diagnostics)
end
