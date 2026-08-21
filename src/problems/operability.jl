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
equilibrium seam is available. The current closure is frozen-dispatch native
static equilibrium and is recorded in result provenance.

Set `compute_helm=true` to request an independent no-load-connected HELM
cross-check on its narrower constant-power/constant-impedance scope. A HELM
failure is retained as `:inconclusive`; unsupported physics is
`:not_applicable`. Set `compute_sensitivity_validation=true` to re-solve small
uniform-load and named P/Q perturbations and compare the implicit derivatives
with independent finite differences. Set
`compute_fixed_point_certificate=true` to evaluate a conservative
Bernstein-style Z-bus contraction condition. The certificate is intentionally
restricted to constant-power and constant-impedance sub-loads; a failed
sufficient condition is `:inconclusive`, never a collapse or multiplicity
claim.
"""
struct OperabilitySpec
    scaling_policy::Union{Nothing,_OPERABILITY_SCALING_POLICY}
    scaling_bases::Any
    provenance::Any
    closure::Symbol
    residual_atol::Float64
    residual_rtol::Float64
    jacobian_step::Float64
    jacobian_rank_rtol::Float64
    sensitivity_step::Float64
    voltage_min::Float64
    voltage_max::Float64
    vuf_max::Float64
    compute_sensitivity::Bool
    compute_sensitivity_validation::Bool
    sensitivity_validation_atol::Float64
    sensitivity_validation_rtol::Float64
    compute_helm::Bool
    helm_max_order::Int
    helm_tol::Float64
    helm_endpoint_atol::Float64
    helm_endpoint_rtol::Float64
    compute_fixed_point_certificate::Bool
    function OperabilitySpec(;
            scaling_policy=nothing,
            scaling_bases=nothing,
            provenance=nothing,
            closure::Symbol=:frozen_dispatch,
            residual_atol::Real=1e-8,
            residual_rtol::Real=1e-6,
            jacobian_step::Real=1e-6,
            jacobian_rank_rtol::Real=1e-8,
            sensitivity_step::Real=1e-5,
            voltage_min::Real=0.0,
            voltage_max::Real=Inf,
            vuf_max::Real=Inf,
            compute_sensitivity::Bool=true,
            compute_sensitivity_validation::Bool=false,
            sensitivity_validation_atol::Real=1e-6,
            sensitivity_validation_rtol::Real=1e-3,
            compute_helm::Bool=false,
            helm_max_order::Integer=40,
            helm_tol::Real=1e-8,
            helm_endpoint_atol::Real=1e-6,
            helm_endpoint_rtol::Real=1e-6,
            compute_fixed_point_certificate::Bool=false)
        scaling_policy === nothing ||
            (_OPERABILITY_SCALING_POLICY !== Any && scaling_policy isa _OPERABILITY_SCALING_POLICY) ||
            throw(ArgumentError("scaling_policy must be an AbstractOpfScalingPolicy"))
        closure === :frozen_dispatch || throw(ArgumentError(
            "only closure=:frozen_dispatch is supported in the native static scope"))
        vals = (residual_atol, residual_rtol, jacobian_step,
                jacobian_rank_rtol, sensitivity_step,
                sensitivity_validation_atol, sensitivity_validation_rtol)
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
        new(scaling_policy, scaling_bases, provenance, closure,
            Float64(residual_atol), Float64(residual_rtol), Float64(jacobian_step),
            Float64(jacobian_rank_rtol), Float64(sensitivity_step),
            Float64(voltage_min), Float64(voltage_max), Float64(vuf_max),
            compute_sensitivity, compute_sensitivity_validation,
            Float64(sensitivity_validation_atol), Float64(sensitivity_validation_rtol),
            compute_helm, Int(helm_max_order), Float64(helm_tol),
            Float64(helm_endpoint_atol), Float64(helm_endpoint_rtol),
            compute_fixed_point_certificate)
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

function _operability_scope_inventory(net::Dict{String,Any})
    models = Set{String}(); configurations = Set{String}()
    for (_, load) in get(net, "load", Dict())
        push!(models, lowercase(string(get(load, "model", "constant_power"))))
        push!(configurations, uppercase(string(get(load, "configuration", "WYE"))))
    end
    Dict{String,Any}(
        "load_models" => sort!(collect(models)),
        "load_configurations" => sort!(collect(configurations)),
        "equilibrium_scope" => "native_ybus_linearized")
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

function _operability_load_connection_count(load)
    count = 0
    for field in ("p_nom", "q_nom")
        value = get(load, field, nothing)
        count = max(count, value === nothing ? 0 :
            value isa AbstractVector ? length(value) : 1)
    end
    count
end

function _operability_perturb_load(net::Dict{String,Any}, id::String, k::Int,
                                   field::String, delta::Real)
    result = deepcopy(net)
    loads = get(result, "load", Dict{String,Any}())
    haskey(loads, id) || throw(ArgumentError("unknown load $(repr(id))"))
    load = loads[id]
    values = get(load, field, nothing)
    vector = values === nothing ? Float64[] :
        (values isa AbstractVector ? Float64.(values) : [Float64(values)])
    if length(vector) < k
        append!(vector, zeros(k - length(vector)))
    end
    1 <= k <= length(vector) || throw(BoundsError(vector, k))
    vector[k] += Float64(delta)
    load[field] = vector
    result
end

function _operability_directional_sensitivity(net, lin, meta, x, J,
                                              residual_scale, state_scale, spec,
                                              id::String, k::Int, field::String)
    load = net["load"][id]
    values = get(load, field, nothing)
    base = values === nothing || (values isa AbstractVector && k > length(values)) ? 0.0 :
        Float64(values isa AbstractVector ? values[k] : values)
    h = spec.sensitivity_step * max(abs(base), 1.0)
    plus = _operability_perturb_load(net, id, k, field, h)
    minus = _operability_perturb_load(net, id, k, field, -h)
    lp = BMOPFTools.ybus_linearized(plus; fold=:constant_z)
    lm = BMOPFTools.ybus_linearized(minus; fold=:constant_z)
    dF = (_operability_residual(lp, meta, x) - _operability_residual(lm, meta, x)) / (2h)
    dF_scaled = dF ./ residual_scale
    dx_scaled = -(J \ dF_scaled)
    dx = dx_scaled .* state_scale
    dvmap = _operability_node_derivative(lin, meta, dx)
    Dict{String,Any}(
        "parameter" => field, "connection" => "$id/$k", "step" => h,
        "units" => field == "p_nom" ? "W" : "var",
        "state_derivative" => dx,
        "node_voltage_derivative" => dvmap,
        "load_connections" => _operability_load_records(net,
            _operability_node_map(lin, _operability_voltage_vector(lin, meta, x)); dvmap),
        "sequences" => _operability_sequences(net,
            _operability_node_map(lin, _operability_voltage_vector(lin, meta, x)); dvmap))
end

function _operability_sensitivity_corrector(lin, meta, x0, coords, spec)
    x = copy(x0); last_residual = Inf
    for iteration in 1:20
        residual = _operability_residual(lin, meta, x)
        _, residual_scale = _operability_scales(lin, meta, coords)
        scaled_residual = residual ./ residual_scale
        last_residual = isempty(scaled_residual) ? 0.0 : maximum(abs, scaled_residual)
        tolerance = spec.residual_atol + spec.residual_rtol
        last_residual <= tolerance && return x, true, iteration, last_residual
        state_scale, residual_scale = _operability_scales(lin, meta, coords)
        Jphys = finite_difference_jacobian(
            y -> _operability_residual(lin, meta, y), x; step=spec.jacobian_step)
        J = _operability_scaled_jacobian(Jphys, state_scale, residual_scale)
        step_scaled = try
            -(J \ (residual ./ residual_scale))
        catch
            return x, false, iteration, last_residual
        end
        all(isfinite, step_scaled) || return x, false, iteration, last_residual
        step = step_scaled .* state_scale
        α = 1.0; accepted = false
        while α >= 1 / 64
            candidate = x .+ α .* step
            candidate_residual = _operability_residual(lin, meta, candidate)
            candidate_scaled = candidate_residual ./ residual_scale
            candidate_norm = isempty(candidate_scaled) ? 0.0 : maximum(abs, candidate_scaled)
            if all(isfinite, candidate_scaled) && candidate_norm < last_residual
                x = candidate; accepted = true; break
            end
            α /= 2
        end
        accepted || return x, false, iteration, last_residual
    end
    x, false, 20, last_residual
end

function _operability_validate_load_scale(net, lin, meta, x, coords, spec, analytic_dx)
    # The implicit derivative uses the small configured residual perturbation;
    # an independent equilibrium re-solve needs a larger parameter step so the
    # voltage change is not lost beneath the current-balance corrector tolerance.
    h = max(100.0 * spec.sensitivity_step, sqrt(eps(Float64)))
    plus = BMOPFTools.ybus_linearized(_operability_scale_network(net, 1.0 + h);
                                      fold=:constant_z)
    minus = BMOPFTools.ybus_linearized(_operability_scale_network(net, 1.0 - h);
                                       fold=:constant_z)
    xp, okp, ip, rp = _operability_sensitivity_corrector(plus, meta, x, coords, spec)
    xm, okm, im, rm = _operability_sensitivity_corrector(minus, meta, x, coords, spec)
    if !(okp && okm)
        return Dict{String,Any}(
            "status" => :inconclusive, "step" => h,
            "plus_iterations" => ip, "minus_iterations" => im,
            "plus_residual" => rp, "minus_residual" => rm,
            "message" => "finite-difference perturbed equilibrium corrector failed")
    end
    finite_difference_dx = (xp - xm) / (2h)
    absolute_error = maximum(abs.(analytic_dx .- finite_difference_dx))
    scale = max(1.0, maximum(abs, analytic_dx), maximum(abs, finite_difference_dx))
    relative_error = absolute_error / scale
    tolerance = spec.sensitivity_validation_atol + spec.sensitivity_validation_rtol * scale
    Dict{String,Any}(
        "status" => absolute_error <= tolerance ? :pass : :fail,
        "step" => h, "analytic_state_derivative" => analytic_dx,
        "finite_difference_state_derivative" => finite_difference_dx,
        "absolute_error" => absolute_error, "relative_error" => relative_error,
        "tolerance" => tolerance, "plus_iterations" => ip, "minus_iterations" => im,
        "plus_residual" => rp, "minus_residual" => rm,
        "message" => "implicit load-scale sensitivity versus re-solved finite difference")
end

function _operability_validate_direction(net, lin, meta, x, coords, spec,
                                         analytic_dx, id::String, k::Int, field::String)
    load = net["load"][id]
    values = get(load, field, nothing)
    base = values === nothing || (values isa AbstractVector && k > length(values)) ? 0.0 :
        Float64(values isa AbstractVector ? values[k] : values)
    h = max(100.0 * spec.sensitivity_step, 1e-3) * max(abs(base), 1.0)
    plus = BMOPFTools.ybus_linearized(_operability_perturb_load(net, id, k, field, h);
                                      fold=:constant_z)
    minus = BMOPFTools.ybus_linearized(_operability_perturb_load(net, id, k, field, -h);
                                       fold=:constant_z)
    xp, okp, ip, rp = _operability_sensitivity_corrector(plus, meta, x, coords, spec)
    xm, okm, im, rm = _operability_sensitivity_corrector(minus, meta, x, coords, spec)
    if !(okp && okm)
        return Dict{String,Any}(
            "status" => :inconclusive, "parameter" => field,
            "connection" => "$id/$k", "step" => h,
            "plus_iterations" => ip, "minus_iterations" => im,
            "plus_residual" => rp, "minus_residual" => rm,
            "message" => "finite-difference directional corrector failed")
    end
    finite_difference_dx = (xp - xm) / (2h)
    absolute_error = maximum(abs.(analytic_dx .- finite_difference_dx))
    scale = max(1.0, maximum(abs, analytic_dx), maximum(abs, finite_difference_dx))
    relative_error = absolute_error / scale
    tolerance = spec.sensitivity_validation_atol + spec.sensitivity_validation_rtol * scale
    Dict{String,Any}(
        "status" => absolute_error <= tolerance ? :pass : :fail,
        "parameter" => field, "connection" => "$id/$k", "step" => h,
        "analytic_state_derivative" => analytic_dx,
        "finite_difference_state_derivative" => finite_difference_dx,
        "absolute_error" => absolute_error, "relative_error" => relative_error,
        "tolerance" => tolerance, "plus_iterations" => ip, "minus_iterations" => im,
        "plus_residual" => rp, "minus_residual" => rm,
        "message" => "implicit directional sensitivity versus re-solved finite difference")
end

function _operability_path_dpdv_evidence(sensitivities)
    records = get(get(sensitivities, "load_scale", Dict{String,Any}()),
                  "load_connections", Dict{String,Any}())
    connections = Dict{String,Any}(); finite_slopes = Float64[]
    for (key_raw, record) in records
        key = String(key_raw)
        slope = Float64(get(record, "path_dP_dV", NaN))
        classification = if !isfinite(slope)
            "not_available"
        elseif slope < -1e-9
            push!(finite_slopes, slope)
            "negative_high_side_indicator"
        elseif slope > 1e-9
            push!(finite_slopes, slope)
            "positive_low_side_indicator"
        else
            push!(finite_slopes, slope)
            "near_nose_indicator"
        end
        connections[key] = Dict{String,Any}(
            "path_dP_dV" => slope, "classification" => classification,
            "sign_convention" => "realized consumption P versus terminal magnitude along uniform load scale")
    end
    Dict{String,Any}(
        "status" => isempty(finite_slopes) ? :not_applicable : :diagnostic,
        "connections" => connections, "finite_count" => length(finite_slopes),
        "interpretation" => "path-qualified branch indicator; not a universal voltage-stability certificate")
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

function _operability_connection_power_slope(sl, magnitude)
    magnitude == 0.0 && return 0.0im
    h = 1e-6 * max(abs(magnitude), 1.0)
    (_subload_S(sl, magnitude + h) - _subload_S(sl, max(magnitude - h, 0.0))) /
        (magnitude + h - max(magnitude - h, 0.0))
end

function _operability_load_records(net, vmap; dvmap=nothing, uniform_scale=false)
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
            requested = ComplexF64(sl.p0 + im * sl.q0)
            realized = ComplexF64(_subload_S(sl, mag))
            local_dS_dmag = ComplexF64(_operability_connection_power_slope(sl, mag))
            path_dS = if uniform_scale && dvmap !== nothing
                realized + local_dS_dmag * dmag
            else
                ComplexF64(NaN + im * NaN)
            end
            path_dP_dV = uniform_scale && dvmap !== nothing && isfinite(dmag) &&
                abs(dmag) > eps(Float64) ? real(path_dS) / dmag : NaN
            key = "$id/$k"
            records[key] = Dict{String,Any}(
                "load" => id, "connection_index" => k,
                "positive" => sl.pos, "negative" => sl.neg,
                "voltage" => ComplexF64(u), "magnitude" => mag,
                "requested_power" => requested, "realized_power" => realized,
                "power_error" => realized - requested,
                "realized_power_local_derivative" => local_dS_dmag,
                "voltage_derivative" => ComplexF64(dv),
                "magnitude_derivative" => dmag,
                "realized_power_derivative" => path_dS,
                "path_dP_dV" => path_dP_dV,
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

function _operability_sequence_sensitivity_evidence(sensitivities)
    records = get(get(sensitivities, "load_scale", Dict{String,Any}()),
                  "sequences", Dict{String,Any}())
    buses = Dict{String,Any}()
    for (bus_raw, record) in records
        bus = String(bus_raw)
        v1 = ComplexF64(get(record, "v1", NaN + im * NaN))
        v2 = ComplexF64(get(record, "v2", NaN + im * NaN))
        dv1 = ComplexF64(get(record, "v1_derivative", NaN + im * NaN))
        dv2 = ComplexF64(get(record, "v2_derivative", NaN + im * NaN))
        dmag1 = isfinite(abs(v1)) && abs(v1) > eps(Float64) ?
            real(conj(v1) * dv1) / abs(v1) : NaN
        dmag2 = isfinite(abs(v2)) && abs(v2) > eps(Float64) ?
            real(conj(v2) * dv2) / abs(v2) : NaN
        buses[bus] = Dict{String,Any}(
            "v1" => v1, "v2" => v2,
            "positive_sequence_magnitude_derivative" => dmag1,
            "negative_sequence_magnitude_derivative" => dmag2,
            "vuf" => get(record, "vuf", NaN),
            "vuf_derivative" => get(record, "vuf_derivative", NaN),
            "interpretation" => "uniform-load path derivative; positive sequence can hide phase-level weakness")
    end
    Dict{String,Any}(
        "status" => isempty(buses) ? :not_applicable : :available,
        "buses" => buses,
        "interpretation" => "positive/negative sequence summaries complement terminal and phase evidence")
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

function _operability_fixed_point_certificate(net::Dict{String,Any}, lin, meta,
                                              spec::OperabilitySpec)
    reasons = String[]
    for (id_raw, load) in get(net, "load", Dict())
        model = lowercase(string(get(load, "model", "constant_power")))
        model in ("constant_power", "constant_impedance") || push!(reasons,
            "load $(repr(String(id_raw))) uses $(repr(model)); the certificate " *
            "currently supports only constant-power and constant-impedance parts")
    end
    assumptions = [
        "frozen-dispatch native ybus_linearized equilibrium",
        "source terminals are fixed and all non-source voltage nodes are included",
        "constant-power compensation currents are represented by connection incidence",
        "constant-impedance parts are folded into the linear Ybus",
        "the contraction region is a uniform complex-voltage polydisc around the no-load solution",
    ]
    base = Dict{String,Any}(
        "status" => :not_applicable,
        "method" => "bernstein_style_zbus_contraction",
        "theorem" => "sufficient existence, uniqueness, and Jacobian nonsingularity condition",
        "assumptions" => assumptions,
        "scope" => "constant_power_and_constant_impedance",
    )
    if !isempty(reasons)
        base["reasons"] = reasons
        base["message"] = "fixed-point certificate is outside the supported load-model scope"
        return base
    end
    nf = length(meta.free)
    if nf == 0
        base["message"] = "no non-source voltage state is available for a Z-bus certificate"
        return base
    end
    Yll = Matrix{ComplexF64}(lin.Y[meta.free, meta.free])
    Yls = isempty(meta.fixed) ? zeros(ComplexF64, nf, 0) :
        Matrix{ComplexF64}(lin.Y[meta.free, meta.fixed])
    Z = try
        Yll \ Matrix{ComplexF64}(I, nf, nf)
    catch err
        base["status"] = :inconclusive
        base["message"] = "source-eliminated Ybus is singular; contraction region is unavailable"
        base["error"] = sprint(showerror, err)
        return base
    end
    w = isempty(meta.fixed) ? zeros(ComplexF64, nf) :
        -(Yll \ (Yls * meta.fixed_values))
    wfull = zeros(ComplexF64, length(lin.nodes))
    wfull[meta.fixed] = meta.fixed_values
    wfull[meta.free] = w
    local_index = Dict{Int,Int}(gi => li for
        (li, gi) in enumerate(meta.free))
    edges = NamedTuple[]
    for (id_raw, load) in sort!(collect(get(net, "load", Dict())); by=x -> String(x[1]))
        id = String(id_raw)
        for (k, sl) in enumerate(_load_subloads(load, net))
            # The constant-Z component has already been absorbed into Y.
            s = ComplexF64(sl.pt.cc + im * sl.qt.cc)
            abs(s) == 0.0 && continue
            pi = get(lin.index, sl.pos, 0)
            ni = sl.neg === nothing ? 0 : get(lin.index, sl.neg, 0)
            a = zeros(Float64, nf)
            haskey(local_index, pi) && (a[local_index[pi]] += 1.0)
            haskey(local_index, ni) && (a[local_index[ni]] -= 1.0)
            anorm = sum(abs, a)
            anorm == 0.0 && continue
            vp = pi == 0 ? 0.0im : wfull[pi]
            vn = ni == 0 ? 0.0im : wfull[ni]
            dv0 = vp - vn
            push!(edges, (key="$id/$k", pos=sl.pos, neg=sl.neg,
                sabs=abs(s), dv0=dv0, a=a, anorm=anorm,
                zalpha=abs.(Z * a)))
        end
    end
    base["state_nodes"] = copy(meta.state_nodes)
    base["coordinate_order"] = "complex free nodes in LinearizedYbus order"
    base["no_load_voltage"] = copy(w)
    base["edge_count"] = length(edges)
    if isempty(edges)
        candidate = _operability_voltage_vector(lin, meta, meta.state)
        distance = maximum(abs.(candidate[meta.free] .- w))
        base["radius"] = 0.0
        base["candidate_distance"] = distance
        base["contraction_factor"] = 0.0
        base["invariance_bound"] = 0.0
        base["status"] = distance <= spec.residual_atol ? :pass : :inconclusive
        base["message"] = base["status"] === :pass ?
            "zero constant-power injection gives the exact no-load certificate" :
            "candidate is not the no-load solution"
        return base
    end
    radius_limit = minimum(abs(e.dv0) / e.anorm for e in edges)
    if !(isfinite(radius_limit) && radius_limit > 0.0)
        base["message"] = "the no-load connection voltage is zero or has no positive safety radius"
        return base
    end
    contribution(rho) = begin
        lower = [abs(e.dv0) - e.anorm * rho for e in edges]
        any(v -> v <= 0.0, lower) && return (Inf, Inf, lower)
        b = maximum(sum(e.zalpha[i] * e.sabs / lower[j]
                        for (j, e) in enumerate(edges)) for i in 1:nf)
        q = maximum(sum(e.zalpha[i] * e.anorm * e.sabs / lower[j]^2
                        for (j, e) in enumerate(edges)) for i in 1:nf)
        b, q, lower
    end
    selected = nothing
    for rho in range(0.0, stop=0.999radius_limit, length=2001)
        b, q, lower = contribution(rho)
        if isfinite(b) && isfinite(q) && b <= rho && q < 1.0
            selected = (rho=rho, b=b, q=q, lower=lower)
            break
        end
    end
    # The map is independently iterated from the energized no-load germ. This
    # is evidence recorded alongside, but not substituted for, the theorem.
    map_complex(vfree) = begin
        vfull = copy(wfull); vfull[meta.free] = vfree
        w + Z * lin.i_comp(vfull)[meta.free]
    end
    pack(v) = vcat(real.(v), imag.(v))
    unpack(x) = ComplexF64.(x[1:nf] .+ im .* x[nf+1:end])
    oracle = try
        fixed_point_oracle(x -> pack(map_complex(unpack(x))), pack(w);
            max_iterations=200, atol=spec.residual_atol, rtol=spec.residual_rtol,
            store_trajectory=false)
    catch err
        nothing
    end
    oracle !== nothing && (base["oracle"] = Dict(
        "converged" => oracle.converged, "cycled" => oracle.cycled,
        "iterations" => oracle.iterations, "residual_norm" => oracle.residual_norm))
    candidate = _operability_voltage_vector(lin, meta, meta.state)
    distance = maximum(abs.(candidate[meta.free] .- w))
    base["candidate_distance"] = distance
    if selected === nothing
        base["message"] = "no invariant contraction polydisc was found; sufficient condition is inconclusive"
        return base
    end
    base["radius"] = selected.rho
    base["invariance_bound"] = selected.b
    base["contraction_factor"] = selected.q
    base["minimum_connection_voltage"] = minimum(selected.lower)
    base["condition_margin"] = 1.0 - selected.q
    inside = distance <= selected.rho + spec.residual_atol
    base["candidate_inside_region"] = inside
    if inside
        base["status"] = :pass
        base["message"] = "Bernstein-style Z-bus contraction certifies a unique no-load-connected solution containing the candidate"
    else
        base["status"] = :inconclusive
        base["message"] = "the contraction condition holds around the no-load germ, but the candidate lies outside its certified region"
    end
    base
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
    claim_keys = ("terminal_voltage_bounds", "sequence_unbalance", "helm_reachability",
                  "fixed_point_certificate")
    claim_statuses = [checks[key].status for key in claim_keys if haskey(checks, key)]
    any(status === :pass for status in claim_statuses) ? :pass : :not_applicable
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
    residual_real = _operability_residual(lin, meta, x)
    endpoint_normalized = isempty(residual_real) ? 0.0 : maximum(
        abs(residual_real[i]) / residual_scale[i] for i in eachindex(residual_real))
    endpoint_tol = spec.residual_atol + spec.residual_rtol
    checks = Dict{String,OperabilityCheck}()
    checks["endpoint"] = OperabilityCheck(
        endpoint_normalized <= endpoint_tol ? :pass : :fail,
        endpoint_normalized, endpoint_tol,
        "max audited-policy-normalized current-balance mismatch on non-source nodes")

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
    load_scale_dx = nothing
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
        load_scale_dx = dx
        dvmap = _operability_node_derivative(lin, meta, dx)
        sensitivities["load_scale"] = Dict{String,Any}(
            "state_derivative" => dx,
            "node_voltage_derivative" => dvmap,
            "load_connections" => _operability_load_records(net, node_voltages;
                dvmap, uniform_scale=true),
            "sequences" => _operability_sequences(net, node_voltages; dvmap))
        checks["load_scale_sensitivity"] = OperabilityCheck(:pass,
            sensitivities["load_scale"]["state_derivative"], nothing,
            "sensitivity to uniform scaling of load nameplate P/Q terms")
    else
        checks["load_scale_sensitivity"] = OperabilityCheck(:not_applicable, nothing, nothing,
            isempty(x) ? "no state to differentiate" : "Jacobian is too close to singular")
    end

    if spec.compute_sensitivity_validation && load_scale_dx !== nothing
        validation = _operability_validate_load_scale(
            net, lin, meta, x, coords, spec, load_scale_dx)
        sensitivities["validation"] = Dict("load_scale" => validation)
        checks["load_scale_sensitivity_validation"] = OperabilityCheck(
            Symbol(validation["status"]), get(validation, "absolute_error", nothing),
            get(validation, "tolerance", nothing), String(validation["message"]))
    elseif spec.compute_sensitivity_validation
        checks["load_scale_sensitivity_validation"] = OperabilityCheck(
            :not_applicable, nothing, nothing,
            "load-scale sensitivity is unavailable for independent validation")
    end

    directions = Dict{String,Any}()
    if spec.compute_sensitivity && !isempty(x) && last(singular_values) > rank_tol
        for (id_raw, load) in sort!(collect(get(net, "load", Dict())); by=x -> String(x[1]))
            id = String(id_raw)
            nconn = _operability_load_connection_count(load)
            for field in ("p_nom", "q_nom")
                family = field == "p_nom" ? "P" : "Q"
                family_records = get!(directions, family, Dict{String,Any}())
                for k in 1:nconn
                    family_records["$id/$k"] = _operability_directional_sensitivity(
                        net, lin, meta, x, J, residual_scale, state_scale, spec, id, k, field)
                end
            end
        end
        sensitivities["directions"] = directions
    end

    if spec.compute_sensitivity_validation && !isempty(directions)
        direction_validation = Dict{String,Any}()
        statuses = Symbol[]
        for (id_raw, load) in sort!(collect(get(net, "load", Dict())); by=x -> String(x[1]))
            id = String(id_raw)
            nconn = _operability_load_connection_count(load)
            for (field, family) in (("p_nom", "P"), ("q_nom", "Q"))
                family_validation = get!(direction_validation, family, Dict{String,Any}())
                for k in 1:nconn
                    key = "$id/$k"
                    record = directions[family][key]
                    validation = _operability_validate_direction(
                        net, lin, meta, x, coords, spec,
                        record["state_derivative"], id, k, field)
                    family_validation[key] = validation
                    push!(statuses, Symbol(validation["status"]))
                end
            end
        end
        sensitivities["validation"] = get(sensitivities, "validation", Dict{String,Any}())
        sensitivities["validation"]["directions"] = direction_validation
        aggregate = any(status -> status === :fail, statuses) ? :fail :
            (any(status -> status === :inconclusive, statuses) ? :inconclusive : :pass)
        checks["directional_sensitivity_validation"] = OperabilityCheck(
            aggregate, direction_validation, nothing,
            "implicit P/Q sensitivities versus re-solved finite differences")
    elseif spec.compute_sensitivity_validation
        checks["directional_sensitivity_validation"] = OperabilityCheck(
            :not_applicable, nothing, nothing,
            "directional sensitivities are unavailable for independent validation")
    end

    branch_evidence = Dict{String,Any}(
        "critical_mode" => _operability_critical_mode(jacobian_factorization, meta.state_nodes),
        "dP_dV" => _operability_path_dpdv_evidence(sensitivities),
        "sequence_sensitivity" => _operability_sequence_sensitivity_evidence(sensitivities))
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
    if spec.compute_fixed_point_certificate
        certificate = _operability_fixed_point_certificate(net, lin, meta, spec)
        branch_evidence["fixed_point_certificate"] = certificate
        certificate_status = Symbol(certificate["status"])
        checks["fixed_point_certificate"] = OperabilityCheck(
            certificate_status, get(certificate, "contraction_factor", nothing), 1.0,
            String(certificate["message"]))
    else
        branch_evidence["fixed_point_certificate"] = Dict{String,Any}(
            "status" => :not_applicable,
            "message" => "fixed-point certificate was not requested")
        checks["fixed_point_certificate"] = OperabilityCheck(:not_applicable,
            nothing, nothing, "fixed-point certificate was not requested")
    end

    provenance = deepcopy(coords.provenance)
    provenance["operability"] = Dict("scope" => "static_ybus_linearized",
        "closure" => String(spec.closure),
        "source_buses" => sort!(collect(source_buses)), "state_nodes" => meta.state_nodes,
        "coordinate_policy" => BMOPFTools.opf_scaling_policy_data(coords.policy),
        "model_inventory" => _operability_scope_inventory(net))
    OperabilityResult(_operability_overall(checks), endpoint_residual,
        endpoint_normalized, meta.state_nodes, x, J, singular_values, condition_number,
        node_voltages, load_records, sequence_records, sensitivities, branch_evidence, checks,
        provenance, unsupported)
end
