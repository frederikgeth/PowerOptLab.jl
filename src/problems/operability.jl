# Post-OPF static voltage operability checks.
#
# The first slice deliberately audits the native static load-flow scope exposed
# by BMOPFTools' ybus_linearized seam. It does not pretend that generators,
# IBRs, or controller equations absent from that seam have been checked.

const _OPERABILITY_STATUSES = (:pass, :fail, :inconclusive, :not_applicable)
# Scaling-policy types/accessors are being staged in BMOPFTools. Keep the
# package loadable against the released BMOPFTools 0.1.0 used by the docs
# environment; the checker itself requires the richer seam at call time.
const _OPERABILITY_SCALING_POLICY = isdefined(BMOPFTools, :AbstractOpfScalingPolicy) ?
    getfield(BMOPFTools, :AbstractOpfScalingPolicy) : Any

"""One independently evaluated operability check."""
struct OperabilityCheck
    status::Symbol
    value::Any
    limit::Any
    message::String
    function OperabilityCheck(status::Symbol, value, limit, message::AbstractString)
        status in _OPERABILITY_STATUSES || throw(ArgumentError(
            "unknown operability status :$status"))
        new(status, value, limit, String(message))
    end
end

"""
    OperabilitySpec(; scaling_policy, ...)

Configuration for [`check_opf_operability`](@ref). The scaling policy is
required because the checker must use the coordinate contract of the audited
OPF, rather than silently choosing its own per-unit convention. A staged
BMOPFTools context can be supplied to `check_opf_operability` to obtain the
policy, coordinate bases, and research provenance automatically.

The first implementation supports native static loads represented by
`ybus_linearized`: constant-P, constant-I, constant-Z, ZIP, and exponential
loads, including wye, single-phase, and delta connections. Generator, IBR, and
endogenous-control equations are reported as outside this scope until a public

Set `compute_helm=true` to request an independent no-load-connected HELM
cross-check on its narrower constant-power/constant-impedance scope. A HELM
failure is retained as `:inconclusive`; unsupported physics is
`:not_applicable`.
"""
struct OperabilitySpec
    scaling_policy::Union{Nothing,_OPERABILITY_SCALING_POLICY}
    scaling_bases::Any
    provenance::Any
    residual_atol::Float64
    residual_rtol::Float64
    jacobian_step::Float64
    jacobian_rank_rtol::Float64
    sensitivity_step::Float64
    voltage_min::Float64
    voltage_max::Float64
    vuf_max::Float64
    compute_sensitivity::Bool
    compute_helm::Bool
    helm_max_order::Int
    helm_tol::Float64
    helm_endpoint_atol::Float64
    helm_endpoint_rtol::Float64
    function OperabilitySpec(;
            scaling_policy=nothing,
            scaling_bases=nothing,
            provenance=nothing,
            residual_atol::Real=1e-8,
            residual_rtol::Real=1e-6,
            jacobian_step::Real=1e-6,
            jacobian_rank_rtol::Real=1e-8,
            sensitivity_step::Real=1e-5,
            voltage_min::Real=0.0,
            voltage_max::Real=Inf,
            vuf_max::Real=Inf,
            compute_sensitivity::Bool=true,
            compute_helm::Bool=false,
            helm_max_order::Integer=40,
            helm_tol::Real=1e-8,
            helm_endpoint_atol::Real=1e-6,
            helm_endpoint_rtol::Real=1e-6)
        scaling_policy === nothing ||
            (_OPERABILITY_SCALING_POLICY !== Any && scaling_policy isa _OPERABILITY_SCALING_POLICY) ||
            throw(ArgumentError("scaling_policy must be an AbstractOpfScalingPolicy"))
        vals = (residual_atol, residual_rtol, jacobian_step,
                jacobian_rank_rtol, sensitivity_step)
        all(x -> isfinite(Float64(x)) && Float64(x) > 0, vals) ||
            throw(ArgumentError("operability tolerances and steps must be finite and > 0"))
        isfinite(Float64(voltage_min)) && voltage_min >= 0 ||
            throw(ArgumentError("voltage_min must be finite and >= 0"))
        (isfinite(Float64(voltage_max)) && voltage_max >= 0) || voltage_max == Inf ||
            throw(ArgumentError("voltage_max must be finite and >= 0 or Inf"))
        (isfinite(Float64(vuf_max)) && vuf_max >= 0) || vuf_max == Inf ||
            throw(ArgumentError("vuf_max must be finite and >= 0 or Inf"))
        helm_max_order >= 1 || throw(ArgumentError("helm_max_order must be >= 1"))
        helm_vals = (helm_tol, helm_endpoint_atol, helm_endpoint_rtol)
        all(x -> isfinite(Float64(x)) && Float64(x) > 0, helm_vals) ||
            throw(ArgumentError("HELM tolerances must be finite and > 0"))
        new(scaling_policy, scaling_bases, provenance,
            Float64(residual_atol), Float64(residual_rtol), Float64(jacobian_step),
            Float64(jacobian_rank_rtol), Float64(sensitivity_step),
            Float64(voltage_min), Float64(voltage_max), Float64(vuf_max),
            compute_sensitivity, compute_helm, Int(helm_max_order), Float64(helm_tol),
            Float64(helm_endpoint_atol), Float64(helm_endpoint_rtol))
    end
end

"""
    OperabilityResult

Structured first-slice post-OPF operability evidence. `jacobian` is the
coordinate-scaled rectangular current-balance Jacobian. The state and voltage
maps retain explicit node ordering so critical singular vectors and sensitivities
can be attributed back to terminals.
"""
struct OperabilityResult <: AbstractSolveResult
    status::Symbol
    endpoint_residual::Float64
    endpoint_residual_normalized::Float64
    state_nodes::Vector{_Node}
    state::Vector{Float64}
    jacobian::Matrix{Float64}
    singular_values::Vector{Float64}
    condition_number::Float64
    node_voltages::Dict{_Node,ComplexF64}
    load_connections::Dict{String,Any}
    sequences::Dict{String,Any}
    sensitivities::Dict{String,Any}
    branch_evidence::Dict{String,Any}
    checks::Dict{String,OperabilityCheck}
    provenance::Dict{String,Any}
    unsupported::Vector{String}
end

function solve_status(result::OperabilityResult)
    reportable = result.status in (:pass, :fail)
    SolveStatus(string(result.status),
                result.status === :pass ? "FEASIBLE_POINT" : "DIAGNOSTIC_POINT",
                true, result.status === :pass, false, reportable)
end

solve_diagnostics(result::OperabilityResult) = (
    endpoint_residual=result.endpoint_residual,
    endpoint_residual_normalized=result.endpoint_residual_normalized,
    smallest_singular_value=isempty(result.singular_values) ? NaN : last(result.singular_values),
    condition_number=result.condition_number,
    branch_evidence=result.branch_evidence,
    check_status=Dict(k => v.status for (k, v) in result.checks),
    unsupported=result.unsupported,
)

Base.show(io::IO, result::OperabilityResult) = print(io,
    "OperabilityResult(status=$(result.status), " *
    "residual=$(result.endpoint_residual), " *
    "sigma_min=$(isempty(result.singular_values) ? NaN : last(result.singular_values)))")

struct _OperabilityCoordinates
    policy::_OPERABILITY_SCALING_POLICY
    voltage_base::Dict{String,Float64}
    current_base::Dict{String,Float64}
    provenance::Dict{String,Any}
end

function _operability_provenance(spec::OperabilitySpec, policy,
                                 context, provenance)
    record = provenance === nothing ? spec.provenance : provenance
    record === nothing && return Dict{String,Any}(
        "scaling_policy" => BMOPFTools.opf_scaling_policy_data(policy))
    record isa AbstractDict || throw(ArgumentError("provenance must be a dictionary"))
    Dict{String,Any}(string(k) => v for (k, v) in record)
end

function _operability_coordinates(net::Dict{String,Any}, spec::OperabilitySpec;
                                  context=nothing)
    policy = context === nothing ? spec.scaling_policy : BMOPFTools.opf_scaling_policy(context)
    policy === nothing && throw(ArgumentError(
        "operability checking requires the audited OPF scaling_policy or a context"))

    supplied = spec.scaling_bases
    bases = Dict{String,Any}()
    if context !== nothing
        for (bus, _) in get(net, "bus", Dict())
            bases[String(bus)] = BMOPFTools.opf_ac_coordinate_bases(context, String(bus))
        end
    elseif supplied !== nothing
        supplied isa AbstractDict || throw(ArgumentError("scaling_bases must be a dictionary keyed by bus"))
        for (bus, base) in supplied
            hasproperty(base, :voltage) && hasproperty(base, :current) ||
                throw(ArgumentError("scaling_bases[$bus] needs voltage and current fields"))
            bases[String(bus)] = base
        end
    elseif !(policy isa BMOPFTools.SIUnitsScaling)
        throw(ArgumentError(
            "non-SI operability checking requires scaling_bases or an OPF context"))
    end

    vbase = Dict{String,Float64}(); ibase = Dict{String,Float64}()
    for (bus, _) in get(net, "bus", Dict())
        bid = String(bus)
        base = get(bases, bid, (voltage=1.0, current=1.0))
        vb = Float64(base.voltage); ib = Float64(base.current)
        isfinite(vb) && vb > 0 && isfinite(ib) && ib > 0 || throw(ArgumentError(
            "invalid voltage/current coordinate base for bus $bid"))
        vbase[bid] = vb; ibase[bid] = ib
    end
    provenance = context === nothing ? spec.provenance :
        BMOPFTools.opf_research_provenance(context)
    _OperabilityCoordinates(policy, vbase, ibase,
                            _operability_provenance(spec, policy, context, provenance))
end

function _operability_source_buses(net)
    Set{String}(String(get(vs, "bus", ""))
                for (_, vs) in get(net, "voltage_source", Dict()))
end

function _operability_voltage(solution, node::_Node)
    bus, terminal = node
    buses = get(solution, "bus", Dict())
    haskey(buses, bus) || throw(ArgumentError("solution is missing bus $(repr(bus))"))
    data = get(buses[bus], terminal, nothing)
    data === nothing && throw(ArgumentError(
        "solution is missing voltage for terminal $(repr(node))"))
    vr = Float64(get(data, "vr", NaN)); vi = Float64(get(data, "vi", NaN))
    isfinite(vr) && isfinite(vi) || throw(ArgumentError(
        "solution voltage for $(repr(node)) is not finite"))
    complex(vr, vi)
end

function _operability_solution_map(solution, nodes)
    Dict{_Node,ComplexF64}(node => ComplexF64(_operability_voltage(solution, node))
                           for node in nodes)
end

function _operability_preflight(net::Dict{String,Any})
    reasons = String[]
    isempty(get(net, "generator", Dict())) || push!(reasons,
        "generator injections are outside the ybus_linearized audit scope")
    isempty(get(net, "ibr", Dict())) || push!(reasons,
        "IBR injections and controls are outside the ybus_linearized audit scope")
    supported_models = Set(("constant_power", "constant_current", "constant_impedance",
                            "zip", "exponential"))
    for (id, load) in get(net, "load", Dict())
        model = string(get(load, "model", "constant_power"))
        model in supported_models || push!(reasons,
            "load $(repr(String(id))) uses unsupported model $(repr(model))")
        cfg = uppercase(string(get(load, "configuration", "WYE")))
        cfg in ("WYE", "SINGLE_PHASE", "DELTA") || push!(reasons,
            "load $(repr(String(id))) uses unsupported configuration $(repr(cfg))")
    end
    reasons
end

function _operability_scale_network(net::Dict{String,Any}, λ::Real)
    result = deepcopy(net)
    scale = Float64(λ)
    for (_, load) in get(result, "load", Dict())
        for key in ("p_nom", "q_nom")
            value = get(load, key, nothing)
            value === nothing && continue
            load[key] = value isa AbstractVector ? scale .* Float64.(value) : scale * Float64(value)
        end
    end
    result
end

function _operability_state_meta(lin, vmap, source_buses)
    fixed = [i for (i, node) in enumerate(lin.nodes) if node[1] in source_buses]
    free = [i for i in eachindex(lin.nodes) if !(i in fixed)]
    fixed_values = ComplexF64[vmap[lin.nodes[i]] for i in fixed]
    state = vcat(real.(ComplexF64[vmap[lin.nodes[i]] for i in free]),
                 imag.(ComplexF64[vmap[lin.nodes[i]] for i in free]))
    (fixed=fixed, free=free, fixed_values=fixed_values, state=state,
     state_nodes=lin.nodes[free])
end

function _operability_voltage_vector(lin, meta, x)
    n = length(lin.nodes); nf = length(meta.free)
    length(x) == 2nf || throw(DimensionMismatch("operability state has wrong length"))
    V = zeros(ComplexF64, n)
    V[meta.fixed] = meta.fixed_values
    V[meta.free] = ComplexF64.(x[1:nf] .+ im .* x[nf+1:end])
    V
end

function _operability_residual(lin, meta, x)
    V = _operability_voltage_vector(lin, meta, x)
    r = lin.Y * V - lin.i_comp(V)
    free = meta.free
    vcat(real.(r[free]), imag.(r[free]))
end

function _operability_complex_residual(lin, meta, x)
    V = _operability_voltage_vector(lin, meta, x)
    r = lin.Y * V - lin.i_comp(V)
    r[meta.free]
end

function _operability_node_map(lin, V)
    Dict{_Node,ComplexF64}(lin.nodes[i] => ComplexF64(V[i]) for i in eachindex(lin.nodes))
end

function _operability_node_derivative(lin, meta, dx)
    nf = length(meta.free); dV = zeros(ComplexF64, length(lin.nodes))
    dV[meta.free] = ComplexF64.(dx[1:nf] .+ im .* dx[nf+1:end])
    _operability_node_map(lin, dV)
end

function _operability_load_records(net, vmap; dvmap=nothing)
    records = Dict{String,Any}()
    for (id_raw, load) in sort!(collect(get(net, "load", Dict())); by=x -> String(x[1]))
        id = String(id_raw)
        for (k, sl) in enumerate(_load_subloads(load, net))
            vp = get(vmap, sl.pos, 0.0im)
            vn = sl.neg === nothing ? 0.0im : get(vmap, sl.neg, 0.0im)
            dv = dvmap === nothing ? 0.0im :
                (get(dvmap, sl.pos, 0.0im) -
                 (sl.neg === nothing ? 0.0im : get(dvmap, sl.neg, 0.0im)))
            u = vp - vn; mag = abs(u)
            dmag = dvmap === nothing || mag == 0 ? NaN : real(conj(u) * dv) / mag
            key = "$id/$k"
            records[key] = Dict{String,Any}(
                "load" => id, "connection_index" => k,
                "positive" => sl.pos, "negative" => sl.neg,
                "voltage" => ComplexF64(u), "magnitude" => mag,
                "voltage_derivative" => ComplexF64(dv),
                "magnitude_derivative" => dmag,
            )
        end
    end
    records
end

function _operability_sequences(net, vmap; dvmap=nothing)
    out = Dict{String,Any}(); α = cis(2pi / 3)
    for (bus_raw, data) in sort!(collect(get(net, "bus", Dict())); by=x -> String(x[1]))
        bus = String(bus_raw); tm = String.(get(data, "terminal_names", String[]))
        phases = _phase_positions(tm, _neutral_labels(net))
        length(phases) == 3 || continue
        np = _neutral_pos(tm, _neutral_labels(net))
        neutral = np === nothing ? 0.0im : get(vmap, (bus, tm[np]), 0.0im)
        dneutral = dvmap === nothing || np === nothing ? 0.0im :
            get(dvmap, (bus, tm[np]), 0.0im)
        phase = ComplexF64[get(vmap, (bus, tm[p]), 0.0im) - neutral for p in phases]
        dphase = dvmap === nothing ? ComplexF64[0.0im for _ in phases] :
            ComplexF64[get(dvmap, (bus, tm[p]), 0.0im) - dneutral for p in phases]
        v0 = sum(phase) / 3
        # PowerOptLab's phase ordering follows the convention used by
        # `_sequence_components`: positive sequence is a + α b + α² c.
        v1 = (phase[1] + α * phase[2] + α^2 * phase[3]) / 3
        v2 = (phase[1] + α^2 * phase[2] + α * phase[3]) / 3
        dv0 = sum(dphase) / 3
        dv1 = (dphase[1] + α * dphase[2] + α^2 * dphase[3]) / 3
        dv2 = (dphase[1] + α^2 * dphase[2] + α * dphase[3]) / 3
        vuf = abs(v1) == 0 ? (abs(v2) == 0 ? 0.0 : Inf) : abs(v2) / abs(v1)
        dvuf = if dvmap === nothing || abs(v1) == 0 || abs(v2) == 0
            NaN
        else
            vuf * (real(conj(v2) * dv2) / abs(v2)^2 -
                   real(conj(v1) * dv1) / abs(v1)^2)
        end
        out[bus] = Dict{String,Any}(
            "phase_voltage" => phase, "phase_voltage_derivative" => dphase,
            "v0" => ComplexF64(v0), "v1" => ComplexF64(v1), "v2" => ComplexF64(v2),
            "v0_derivative" => ComplexF64(dv0), "v1_derivative" => ComplexF64(dv1),
            "v2_derivative" => ComplexF64(dv2), "vuf" => vuf,
            "vuf_derivative" => dvuf)
    end
    out
end

function _operability_scales(lin, meta, coords::_OperabilityCoordinates)
    vs = Float64[coords.voltage_base[node[1]] for node in meta.state_nodes]
    is = Float64[coords.current_base[node[1]] for node in meta.state_nodes]
    vcat(vs, vs), vcat(is, is)
end

function _operability_scaled_jacobian(J, state_scale, residual_scale)
    result = similar(J)
    for i in axes(J, 1), j in axes(J, 2)
        result[i, j] = J[i, j] * state_scale[j] / residual_scale[i]
    end
    result
end

function _operability_node_key(node)
    "$(node[1])/$(node[2])"
end

function _operability_critical_mode(factorization, state_nodes)
    factorization === nothing && return Dict{String,Any}(
        "status" => :not_applicable,
        "message" => "no voltage Jacobian was available")
    u = Float64.(factorization.U[:, end])
    v = Float64.(factorization.V[:, end])
    nf = length(state_nodes)
    right = Dict{String,Float64}(); left = Dict{String,Float64}()
    for (k, node) in enumerate(state_nodes)
        key = _operability_node_key(node)
        right[key] = hypot(v[k], v[k + nf])
        left[key] = hypot(u[k], u[k + nf])
    end
    Dict{String,Any}(
        "status" => :pass,
        "left_vector" => u,
        "right_vector" => v,
        "left_node_participation" => left,
        "right_node_participation" => right,
        "state_nodes" => copy(state_nodes),
        "coordinate_order" => "[real(state_nodes); imag(state_nodes)]",
        "message" => "critical left/right singular vectors for the smallest scaled singular value")
end

function _operability_helm_preflight(net::Dict{String,Any})
    reasons = String[]
    isempty(get(net, "generator", Dict())) || push!(reasons,
        "generator injections are outside the HELM cross-check scope")
    isempty(get(net, "ibr", Dict())) || push!(reasons,
        "IBR injections and controls are outside the HELM cross-check scope")
    for (id_raw, load) in get(net, "load", Dict())
        id = String(id_raw)
        model = lowercase(string(get(load, "model", "constant_power")))
        model in ("constant_power", "constant_impedance") || push!(reasons,
            "load $(repr(id)) uses $(repr(model)); HELM currently supports only constant-power and constant-impedance parts")
        cfg = uppercase(string(get(load, "configuration", "WYE")))
        cfg in ("WYE", "SINGLE_PHASE", "DELTA") || push!(reasons,
            "load $(repr(id)) uses unsupported HELM configuration $(repr(cfg))")
    end
    for (id_raw, source) in get(net, "voltage_source", Dict())
        cfg = uppercase(string(get(source, "configuration", "WYE")))
        cfg in ("WYE", "SINGLE_PHASE") || push!(reasons,
            "voltage source $(repr(String(id_raw))) uses unsupported HELM configuration $(repr(cfg))")
    end
    reasons
end

function _operability_helm_reachability(net, node_voltages, spec::OperabilitySpec)
    reasons = _operability_helm_preflight(net)
    if !isempty(reasons)
        return Dict{String,Any}(
            "status" => :not_applicable,
            "homotopy" => "helm_load_scale_0_to_1",
            "base" => "energized_no_load_germ",
            "reasons" => reasons,
            "message" => "HELM cross-check is outside the supported physics scope")
    end
    hr = try
        helm_series(net; max_order=spec.helm_max_order, tol=spec.helm_tol)
    catch err
        # Expected unsupported-physics errors should have been caught above;
        # preserve unexpected failures as implementation/model errors.
        rethrow(err)
    end
    evidence = Dict{String,Any}(
        "homotopy" => "helm_load_scale_0_to_1",
        "base" => "energized_no_load_germ",
        "helm_status" => hr.status,
        "converged" => hr.converged,
        "residual" => hr.residual,
        "order" => hr.n_order,
        "singularity_estimate" => hr.singularity_estimate)
    if !hr.converged || hr.status !== :converged
        evidence["status"] = :inconclusive
        evidence["message"] = "HELM did not establish a converged endpoint; this is not a non-existence certificate"
        return evidence
    end
    common = intersect(keys(node_voltages), keys(hr.V))
    if isempty(common)
        evidence["status"] = :inconclusive
        evidence["message"] = "HELM and the audited solution have no comparable voltage nodes"
        return evidence
    end
    mismatch = maximum(abs(hr.V[node] - node_voltages[node]) for node in common)
    limit = spec.helm_endpoint_atol + spec.helm_endpoint_rtol * max(
        1.0, maximum(abs(node_voltages[node]) for node in common))
    evidence["common_nodes"] = length(common)
    evidence["endpoint_mismatch"] = mismatch
    evidence["endpoint_limit"] = limit
    evidence["status"] = mismatch <= limit ? :pass : :fail
    evidence["message"] = mismatch <= limit ?
        "HELM endpoint agrees with the audited voltage vector on the no-load-connected branch" :
        "HELM converged to an endpoint that does not match the audited voltage vector"
    evidence
end

function _operability_empty_result(coords, unsupported, message)
    checks = Dict("scope" => OperabilityCheck(:not_applicable, unsupported, nothing, message))
    OperabilityResult(:not_applicable, NaN, NaN, _Node[], Float64[], zeros(0, 0),
        Float64[], Inf, Dict{_Node,ComplexF64}(), Dict{String,Any}(),
        Dict{String,Any}(), Dict{String,Any}(), Dict{String,Any}(), checks,
        coords.provenance, unsupported)
end

function _operability_overall(checks::Dict{String,OperabilityCheck})
    any(c.status === :fail for c in values(checks)) && return :fail
    any(c.status === :inconclusive for c in values(checks)) && return :inconclusive
    haskey(checks, "endpoint") && checks["endpoint"].status === :not_applicable &&
        return :not_applicable
    any(c.status === :pass for c in values(checks)) ? :pass : :not_applicable
end

"""
    check_opf_operability(net, solution; spec, context=nothing)

Audit a solved static OPF point in the native load-flow scope. `solution` is the
SI-valued result dictionary returned by BMOPFTools. Pass the staged OPF context
when available so the effective scaling policy, coordinate bases, and research
provenance are inherited; otherwise construct `OperabilitySpec` with an
explicit `scaling_policy` (and `scaling_bases` for non-SI policies).
"""
function check_opf_operability(net::Dict{String,Any}, solution::AbstractDict;
                               spec::Union{Nothing,OperabilitySpec}=nothing,
                               context=nothing)
    spec === nothing && throw(ArgumentError(
        "pass OperabilitySpec(scaling_policy=...) or an OPF context"))
    coords = _operability_coordinates(net, spec; context)
    unsupported = _operability_preflight(net)
    !isempty(unsupported) && return _operability_empty_result(coords, unsupported,
        "candidate contains physics outside the first operability residual scope")
    source_buses = _operability_source_buses(net)
    isempty(source_buses) && throw(ArgumentError("operability checking requires a voltage source"))

    lin = BMOPFTools.ybus_linearized(net; fold=:constant_z)
    vmap = _operability_solution_map(solution, lin.nodes)
    meta = _operability_state_meta(lin, vmap, source_buses)
    x = meta.state
    residual_complex = _operability_complex_residual(lin, meta, x)
    endpoint_residual = isempty(residual_complex) ? 0.0 : maximum(abs, residual_complex)
    state_scale, residual_scale = _operability_scales(lin, meta, coords)
    endpoint_normalized = isempty(residual_complex) ? 0.0 : maximum(
        abs(residual_complex[i]) / residual_scale[i] for i in eachindex(residual_complex))
    endpoint_tol = spec.residual_atol + spec.residual_rtol * max(1.0,
        isempty(residual_complex) ? 0.0 : maximum(abs, residual_complex))
    checks = Dict{String,OperabilityCheck}()
    checks["endpoint"] = OperabilityCheck(
        endpoint_residual <= endpoint_tol ? :pass : :fail,
        endpoint_residual, endpoint_tol,
        "max current-balance mismatch on non-source nodes")

    V = _operability_voltage_vector(lin, meta, x)
    node_voltages = _operability_node_map(lin, V)
    load_records = _operability_load_records(net, node_voltages)
    sequence_records = _operability_sequences(net, node_voltages)
    magnitudes = [Float64(record["magnitude"]) for record in values(load_records)]
    if isempty(magnitudes) || (spec.voltage_min == 0.0 && spec.voltage_max == Inf)
        checks["terminal_voltage_bounds"] = OperabilityCheck(:not_applicable,
            isempty(magnitudes) ? nothing : extrema(magnitudes),
            (spec.voltage_min, spec.voltage_max), "no voltage magnitude limits requested")
    else
        ok = all(v -> spec.voltage_min <= v <= spec.voltage_max, magnitudes)
        checks["terminal_voltage_bounds"] = OperabilityCheck(ok ? :pass : :fail,
            extrema(magnitudes), (spec.voltage_min, spec.voltage_max),
            "all modeled load-terminal voltage magnitudes")
    end
    vufs = [Float64(r["vuf"]) for r in values(sequence_records)]
    if spec.vuf_max == Inf || isempty(vufs)
        checks["sequence_unbalance"] = OperabilityCheck(:not_applicable,
            isempty(vufs) ? nothing : maximum(vufs), spec.vuf_max,
            isempty(vufs) ? "no complete three-phase bus" : "no VUF limit requested")
    else
        ok = all(v -> v <= spec.vuf_max, vufs)
        checks["sequence_unbalance"] = OperabilityCheck(ok ? :pass : :fail,
            maximum(vufs), spec.vuf_max, "maximum |V₂|/|V₁| over complete three-phase buses")
    end

    Jphys = isempty(x) ? zeros(0, 0) : finite_difference_jacobian(
        y -> _operability_residual(lin, meta, y), x; step=spec.jacobian_step)
    J = _operability_scaled_jacobian(Jphys, state_scale, residual_scale)
    jacobian_factorization = isempty(J) ? nothing : svd(J)
    singular_values = jacobian_factorization === nothing ? Float64[] :
        Float64.(jacobian_factorization.S)
    condition_number = isempty(singular_values) || last(singular_values) == 0 ?
        Inf : first(singular_values) / last(singular_values)
    rank_tol = isempty(singular_values) ? Inf :
        spec.jacobian_rank_rtol * max(first(singular_values), 1.0)
    if isempty(x)
        checks["jacobian_regular"] = OperabilityCheck(:not_applicable, nothing, nothing,
            "no non-source voltage state remains after source/reference elimination")
    else
        checks["jacobian_regular"] = OperabilityCheck(last(singular_values) > rank_tol ?
            :pass : :inconclusive, last(singular_values), rank_tol,
            "scaled rectangular current-balance Jacobian smallest singular value")
    end

    sensitivities = Dict{String,Any}()
    if spec.compute_sensitivity && !isempty(x) && last(singular_values) > rank_tol
        h = spec.sensitivity_step
        plus = _operability_scale_network(net, 1.0 + h)
        minus = _operability_scale_network(net, 1.0 - h)
        lp = BMOPFTools.ybus_linearized(plus; fold=:constant_z)
        lm = BMOPFTools.ybus_linearized(minus; fold=:constant_z)
        fp = _operability_residual(lp, meta, x); fm = _operability_residual(lm, meta, x)
        dF = (fp - fm) / (2h)
        dF_scaled = dF ./ residual_scale
        dx_scaled = -(J \ dF_scaled)
        dx = dx_scaled .* state_scale
        dvmap = _operability_node_derivative(lin, meta, dx)
        sensitivities["load_scale"] = Dict{String,Any}(
            "state_derivative" => dx,
            "node_voltage_derivative" => dvmap,
            "load_connections" => _operability_load_records(net, node_voltages; dvmap),
            "sequences" => _operability_sequences(net, node_voltages; dvmap))
        checks["load_scale_sensitivity"] = OperabilityCheck(:pass,
            sensitivities["load_scale"]["state_derivative"], nothing,
            "sensitivity to uniform scaling of load nameplate P/Q terms")
    else
        checks["load_scale_sensitivity"] = OperabilityCheck(:not_applicable, nothing, nothing,
            isempty(x) ? "no state to differentiate" : "Jacobian is too close to singular")
    end

    branch_evidence = Dict{String,Any}(
        "critical_mode" => _operability_critical_mode(jacobian_factorization, meta.state_nodes))
    if spec.compute_helm
        reachability = _operability_helm_reachability(net, node_voltages, spec)
        branch_evidence["reachability"] = reachability
        reach_status = Symbol(reachability["status"])
        checks["helm_reachability"] = OperabilityCheck(reach_status,
            get(reachability, "endpoint_mismatch", nothing),
            get(reachability, "endpoint_limit", nothing),
            String(reachability["message"]))
    else
        branch_evidence["reachability"] = Dict{String,Any}(
            "status" => :not_applicable,
            "message" => "HELM cross-check was not requested")
        checks["helm_reachability"] = OperabilityCheck(:not_applicable, nothing, nothing,
            "HELM cross-check was not requested")
    end

    provenance = deepcopy(coords.provenance)
    provenance["operability"] = Dict("scope" => "static_ybus_linearized",
        "source_buses" => sort!(collect(source_buses)), "state_nodes" => meta.state_nodes,
        "coordinate_policy" => BMOPFTools.opf_scaling_policy_data(coords.policy))
    OperabilityResult(_operability_overall(checks), endpoint_residual,
        endpoint_normalized, meta.state_nodes, x, J, singular_values, condition_number,
        node_voltages, load_records, sequence_records, sensitivities, branch_evidence, checks,
        provenance, unsupported)
end
