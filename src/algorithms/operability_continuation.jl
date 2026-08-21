# Static load-scale continuation for the post-OPF operability processor.
#
# The natural trace follows loading parameter λ from the energized no-load state
# to the audited endpoint. The opt-in pseudo-arclength trace below follows the
# same static scope with an augmented corrector; it is a first fold-capable
# diagnostic slice, not a global branch-discovery or fold-certification engine.

"""Validated options for a static operability continuation trace."""
struct OperabilityContinuationSpec
    initial_step::Float64
    min_step::Float64
    max_step::Float64
    corrector_tol::Float64
    max_corrector_iterations::Int
    jacobian_step::Float64
    fold_sigma_tol::Float64
    endpoint_atol::Float64
    endpoint_rtol::Float64
    function OperabilityContinuationSpec(;
            initial_step::Real=0.1,
            min_step::Real=1e-4,
            max_step::Real=0.25,
            corrector_tol::Real=1e-7,
            max_corrector_iterations::Integer=20,
            jacobian_step::Real=1e-6,
            fold_sigma_tol::Real=1e-6,
            endpoint_atol::Real=1e-5,
            endpoint_rtol::Real=1e-5)
        vals = (initial_step, min_step, max_step, corrector_tol, jacobian_step,
                fold_sigma_tol, endpoint_atol, endpoint_rtol)
        all(x -> isfinite(Float64(x)) && Float64(x) > 0, vals) ||
            throw(ArgumentError("continuation steps, tolerances, and thresholds must be finite and > 0"))
        min_step <= initial_step <= max_step || throw(ArgumentError(
            "continuation steps must satisfy min_step <= initial_step <= max_step"))
        max_corrector_iterations >= 1 || throw(ArgumentError(
            "max_corrector_iterations must be >= 1"))
        new(Float64(initial_step), Float64(min_step), Float64(max_step),
            Float64(corrector_tol), Int(max_corrector_iterations), Float64(jacobian_step),
            Float64(fold_sigma_tol), Float64(endpoint_atol), Float64(endpoint_rtol))
    end
end

"""Evidence recorded by [`continue_opf_operability`](@ref)."""
struct OperabilityContinuationResult <: AbstractSolveResult
    status::Symbol
    message::String
    lambdas::Vector{Float64}
    states::Vector{Vector{Float64}}
    node_voltages::Vector{Dict{_Node,ComplexF64}}
    residuals::Vector{Float64}
    singular_values::Vector{Float64}
    condition_numbers::Vector{Float64}
    corrector_iterations::Vector{Int}
    events::Vector{Dict{String,Any}}
    state_nodes::Vector{_Node}
    endpoint_match::Union{Nothing,Bool}
    endpoint_distance::Float64
    provenance::Dict{String,Any}
end

function solve_status(result::OperabilityContinuationResult)
    publishable = result.status === :pass
    _result_solve_status(string(result.status), publishable;
        primal_status=publishable ? "FEASIBLE_POINT" : "DIAGNOSTIC_POINT")
end

solve_diagnostics(result::OperabilityContinuationResult) = (
    status=result.status,
    lambdas=result.lambdas,
    residuals=result.residuals,
    smallest_singular_value=isempty(result.singular_values) ? NaN : last(result.singular_values),
    endpoint_match=result.endpoint_match,
    endpoint_distance=result.endpoint_distance,
    margin=operability_continuation_margin(result),
    event_count=length(result.events),
)

Base.show(io::IO, result::OperabilityContinuationResult) = print(io,
    "OperabilityContinuationResult(status=$(result.status), " *
    "points=$(length(result.lambdas)), endpoint_distance=$(result.endpoint_distance))")

function _continuation_empty_result(status, message, provenance, events=Dict{String,Any}[])
    OperabilityContinuationResult(status, message, Float64[], Vector{Float64}[],
        Dict{_Node,ComplexF64}[], Float64[], Float64[], Float64[], Int[], events,
        _Node[], nothing, NaN, provenance)
end

"""
    operability_continuation_margin(result; reference_lambda=1.0)

Summarize the first declared path boundary observed in a continuation trace.
The summary considers voltage-limit and fold-candidate events; a localized
bordered fold uses its refined λ when available. The returned margin is
relative to `reference_lambda` (normally the audited endpoint λ=1), so it is a
path-specific diagnostic rather than a global transfer-capacity certificate.
If the trace ended before observing one of these mechanisms, the result has
`status == :not_observed` rather than implying unlimited margin.
"""
function _continuation_margin(events; reference_lambda::Real=1.0)
    reference = Float64(reference_lambda)
    isfinite(reference) && reference > 0 || throw(ArgumentError(
        "reference_lambda must be finite and > 0"))
    candidates = NamedTuple{(:lambda, :kind, :event),
                            Tuple{Float64,Symbol,Dict{String,Any}}}[]
    for event in events
        kind_string = get(event, "kind", "")
        kind = Symbol(kind_string)
        kind in (:voltage_limit, :fold_candidate) || continue
        λ = get(event, "lambda", NaN)
        if kind === :fold_candidate
            localized = get(event, "fold_localization", nothing)
            if localized isa AbstractDict && get(localized, "status", nothing) === :pass
                λ = get(localized, "lambda", λ)
            end
        end
        isfinite(Float64(λ)) || continue
        push!(candidates, (lambda=Float64(λ), kind=kind, event=deepcopy(event)))
    end
    isempty(candidates) && return Dict{String,Any}(
        "status" => :not_observed,
        "mechanism" => nothing,
        "lambda" => NaN,
        "reference_lambda" => reference,
        "parameter_margin" => NaN,
        "relative_margin" => NaN,
        "pre_reference" => false,
        "message" => "no declared voltage-limit or fold boundary was observed")
    first_boundary = findmin(getfield.(candidates, :lambda))[2]
    candidate = candidates[first_boundary]
    parameter_margin = candidate.lambda - reference
    Dict{String,Any}(
        "status" => :observed,
        "mechanism" => candidate.kind,
        "lambda" => candidate.lambda,
        "reference_lambda" => reference,
        "parameter_margin" => parameter_margin,
        "relative_margin" => parameter_margin / reference,
        "pre_reference" => candidate.lambda < reference,
        "event" => candidate.event,
        "message" => "first declared path boundary on the recorded trace")
end

"""
    operability_continuation_margin(result; reference_lambda=1.0)

Summarize the first declared path boundary observed in a continuation result.
The summary considers voltage-limit and fold-candidate events; a localized
bordered fold uses its refined λ when available. The returned margin is
relative to `reference_lambda` (normally the audited endpoint λ=1), so it is a
path-specific diagnostic rather than a global transfer-capacity certificate.
If the trace ended before observing one of these mechanisms, the result has
`status == :not_observed` rather than implying unlimited margin.
"""
operability_continuation_margin(
    result::OperabilityContinuationResult; reference_lambda::Real=1.0) =
    _continuation_margin(result.events; reference_lambda)

"""
    operability_continuation_rows(result)

Return deterministic, one-row-per-accepted-point records for a continuation
trace. Each row preserves the index into `result.lambdas`, `result.states`, and
`result.node_voltages`, and includes residual, conditioning, corrector, and
curvature evidence plus event kinds observed at that λ. Curvature and
arclength are `NaN` for the no-load base and any fixed-λ endpoint refinement
that has no predictor tangent.
"""
function operability_continuation_rows(result::OperabilityContinuationResult)
    continuation = get(result.provenance, "continuation", Dict{String,Any}())
    curvature_history = get(continuation, "curvature_history", Float64[])
    arclength_steps = get(continuation, "arclength_steps", Float64[])
    rows = NamedTuple[]
    for index in eachindex(result.lambdas)
        λ = result.lambdas[index]
        event_kinds = Symbol[]
        for event in result.events
            event_lambda = get(event, "lambda", NaN)
            event_lambda isa Real || continue
            isfinite(Float64(event_lambda)) || continue
            tolerance = max(1e-10, 1e-8 * max(1.0, abs(λ)))
            abs(Float64(event_lambda) - λ) <= tolerance || continue
            kind = Symbol(get(event, "kind", ""))
            kind in event_kinds || push!(event_kinds, kind)
        end
        tangent_index = index - 1
        curvature = tangent_index >= 1 && tangent_index <= length(curvature_history) ?
            Float64(curvature_history[tangent_index]) : NaN
        arclength_step = tangent_index >= 1 && tangent_index <= length(arclength_steps) ?
            Float64(arclength_steps[tangent_index]) : NaN
        push!(rows, (
            index=index,
            lambda=Float64(λ),
            residual_norm=Float64(result.residuals[index]),
            sigma_min=Float64(result.singular_values[index]),
            condition_number=Float64(result.condition_numbers[index]),
            corrector_iterations=Int(result.corrector_iterations[index]),
            curvature=curvature,
            arclength_step=arclength_step,
            event_kinds=event_kinds,
        ))
    end
    rows
end

function _continuation_scaled_jacobian(lin, meta, x, coords, cfg)
    state_scale, residual_scale = _operability_scales(lin, meta, coords)
    Jphys = finite_difference_jacobian(
        y -> _operability_residual(lin, meta, y), x; step=cfg.jacobian_step)
    J = _operability_scaled_jacobian(Jphys, state_scale, residual_scale)
    J, state_scale, residual_scale
end

function _continuation_corrector(lin, meta, x0, coords, cfg)
    x = copy(x0)
    last_norm = Inf
    for iteration in 1:cfg.max_corrector_iterations
        residual = _operability_residual(lin, meta, x)
        _, residual_scale = _operability_scales(lin, meta, coords)
        scaled_residual = residual ./ residual_scale
        norm_residual = isempty(residual) ? 0.0 : maximum(abs, residual)
        last_norm = norm_residual
        norm_residual <= cfg.corrector_tol && return x, true, iteration, last_norm
        J, state_scale, _ = _continuation_scaled_jacobian(lin, meta, x, coords, cfg)
        step_scaled = try
            -(J \ scaled_residual)
        catch
            return x, false, iteration, last_norm
        end
        all(isfinite, step_scaled) || return x, false, iteration, last_norm
        step = step_scaled .* state_scale
        α = 1.0
        accepted = false
        while α >= 1 / 64
            candidate = x .+ α .* step
            candidate_residual = _operability_residual(lin, meta, candidate)
            candidate_norm = isempty(candidate_residual) ? 0.0 : maximum(abs, candidate_residual)
            if isfinite(candidate_norm) && candidate_norm < norm_residual
                x = candidate
                accepted = true
                break
            end
            α /= 2
        end
        accepted || return x, false, iteration, last_norm
    end
    x, false, cfg.max_corrector_iterations, last_norm
end

function _continuation_point(lin, meta, x, coords, cfg)
    residual = _operability_complex_residual(lin, meta, x)
    residual_norm = isempty(residual) ? 0.0 : maximum(abs, residual)
    J, _, _ = _continuation_scaled_jacobian(lin, meta, x, coords, cfg)
    singular_values = isempty(J) ? Float64[] : Float64.(svdvals(J))
    condition_number = isempty(singular_values) || last(singular_values) == 0 ?
        Inf : first(singular_values) / last(singular_values)
    node_voltages = _operability_node_map(lin, _operability_voltage_vector(lin, meta, x))
    (residual_norm=residual_norm, singular_values=singular_values,
     condition_number=condition_number, node_voltages=node_voltages)
end

function _continuation_voltage_events(net, node_voltages, λ, spec)
    events = Dict{String,Any}[]
    (spec.voltage_min == 0.0 && spec.voltage_max == Inf) && return events
    records = _operability_load_records(net, node_voltages)
    for (key, record) in records
        magnitude = Float64(record["magnitude"])
        if magnitude < spec.voltage_min || magnitude > spec.voltage_max
            push!(events, Dict{String,Any}(
                "kind" => "voltage_limit", "lambda" => λ, "connection" => key,
                "magnitude" => magnitude,
                "limits" => (spec.voltage_min, spec.voltage_max)))
        end
    end
    events
end

"""
    continue_opf_operability(net, solution; spec, continuation, context=nothing)

Trace the static native-load equilibrium from the energized no-load state
(`λ=0`) toward the supplied solution (`λ=1`) by scaling load nameplates and
using damped Newton correctors. The result records every accepted point,
residual, smallest singular value, corrector iteration count, and explicit
corrector/near-singular/voltage-limit events.

This is a natural-parameter predictor-corrector trace. It is useful branch
evidence and a high-voltage preflight, but it does not cross folds and is not a
pseudo-arclength reachability certificate. Pass `stop_on_voltage_limit=true`
to terminate with `:fail` at the first declared terminal-voltage violation;
the default keeps collecting the trace through such events.
"""
function continue_opf_operability(net::Dict{String,Any}, solution::AbstractDict;
                                  spec::OperabilitySpec,
                                  continuation::OperabilityContinuationSpec=
                                      OperabilityContinuationSpec(),
                                  context=nothing,
                                  stop_on_voltage_limit::Bool=false)
    coords = _operability_coordinates(net, spec; context)
    unsupported = _operability_preflight(net)
    !isempty(unsupported) && return _continuation_empty_result(
        :not_applicable, "candidate contains physics outside the continuation scope",
        coords.provenance, [Dict{String,Any}("kind" => "unsupported_physics",
                                             "reasons" => unsupported)])
    source_buses = _operability_source_buses(net)
    isempty(source_buses) && throw(ArgumentError("continuation requires a voltage source"))
    lin_target = BMOPFTools.ybus_linearized(net; fold=:constant_z)
    target_vmap = _operability_solution_map(solution, lin_target.nodes)
    meta = _operability_state_meta(lin_target, target_vmap, source_buses)
    target_x = meta.state
    isempty(target_x) && return _continuation_empty_result(
        :pass, "network has no free non-source voltage state", coords.provenance)

    # At λ=0 the load compensation is zero, so one Newton linearization gives
    # the energized no-load germ without starting at a zero-voltage singularity.
    net_zero = _operability_scale_network(net, 0.0)
    lin_zero = BMOPFTools.ybus_linearized(net_zero; fold=:constant_z)
    Jzero, state_scale, residual_scale = _continuation_scaled_jacobian(
        lin_zero, meta, target_x, coords, continuation)
    base_x = target_x - (Jzero \ (_operability_residual(lin_zero, meta, target_x) ./
        residual_scale)) .* state_scale
    base_x, base_ok, base_iterations, _ = _continuation_corrector(
        lin_zero, meta, base_x, coords, continuation)
    base_ok || return _continuation_empty_result(
        :inconclusive, "energized no-load corrector failed", coords.provenance,
        [Dict{String,Any}("kind" => "base_corrector_failure",
                           "iterations" => base_iterations)])

    lambdas = [0.0]; states = [copy(base_x)]
    points = [_continuation_point(lin_zero, meta, base_x, coords, continuation)]
    residuals = [points[1].residual_norm]
    singular_values = [isempty(points[1].singular_values) ? NaN : last(points[1].singular_values)]
    condition_numbers = [points[1].condition_number]
    corrector_iterations = [base_iterations]
    node_voltages = [points[1].node_voltages]
    events = _continuation_voltage_events(net_zero, node_voltages[1], 0.0, spec)
    if !isempty(singular_values) && singular_values[1] <= continuation.fold_sigma_tol
        push!(events, Dict{String,Any}("kind" => "near_singular", "lambda" => 0.0,
                                       "sigma_min" => singular_values[1]))
    end

    λ = 0.0; x = base_x; step = continuation.initial_step
    previous_step = nothing; previous_delta = nothing
    voltage_limit_stop = false
    reached = false; endpoint_match = nothing; endpoint_distance = NaN
    message = "continuation did not reach λ=1"
    while λ < 1.0 - eps(Float64)
        Δλ = min(step, 1.0 - λ)
        λtrial = λ + Δλ
        xpred = previous_delta === nothing ? x : x .+ previous_delta .* (Δλ / previous_step)
        net_trial = _operability_scale_network(net, λtrial)
        lin_trial = BMOPFTools.ybus_linearized(net_trial; fold=:constant_z)
        xnew, ok, iterations, _ = _continuation_corrector(
            lin_trial, meta, xpred, coords, continuation)
        if !ok
            push!(events, Dict{String,Any}("kind" => "corrector_failure",
                "lambda" => λtrial, "step" => Δλ, "iterations" => iterations))
            step /= 2
            if step < continuation.min_step
                message = "corrector step fell below min_step before λ=1"
                break
            end
            continue
        end
        point = _continuation_point(lin_trial, meta, xnew, coords, continuation)
        push!(lambdas, λtrial); push!(states, copy(xnew)); push!(points, point)
        push!(residuals, point.residual_norm)
        push!(singular_values, isempty(point.singular_values) ? NaN : last(point.singular_values))
        push!(condition_numbers, point.condition_number)
        push!(corrector_iterations, iterations); push!(node_voltages, point.node_voltages)
        voltage_events = _continuation_voltage_events(
            net_trial, point.node_voltages, λtrial, spec)
        append!(events, voltage_events)
        if stop_on_voltage_limit && !isempty(voltage_events)
            voltage_limit_stop = true
            message = "continuation stopped at a declared voltage limit"
            break
        end
        if !isempty(point.singular_values) && last(point.singular_values) <= continuation.fold_sigma_tol
            push!(events, Dict{String,Any}("kind" => "near_singular", "lambda" => λtrial,
                                           "sigma_min" => last(point.singular_values)))
        end
        previous_delta = xnew - x; previous_step = Δλ
        x = xnew; λ = λtrial
        if iterations <= 4
            step = min(continuation.max_step, step * 1.4)
        elseif iterations >= continuation.max_corrector_iterations ÷ 2
            step = max(continuation.min_step, step / 1.5)
        end
    end

    reached = !isempty(lambdas) && last(lambdas) >= 1.0 - eps(Float64)
    if reached && voltage_limit_stop
        status = :fail
        message = "natural load-scale trace reached a declared voltage limit"
    elseif reached
        endpoint_distance = maximum(abs(node_voltages[end][node] - target_vmap[node])
                                    for node in keys(target_vmap))
        endpoint_limit = continuation.endpoint_atol + continuation.endpoint_rtol * max(
            1.0, maximum(abs, values(target_vmap)))
        endpoint_match = endpoint_distance <= endpoint_limit
        pre_endpoint_fold = any(get(event, "kind", "") == "near_singular" &&
                                get(event, "lambda", 1.0) < 1.0 - eps(Float64)
                                for event in events)
        if endpoint_match && !pre_endpoint_fold
            message = "natural load-scale trace reached and matched λ=1 endpoint"
            status = :pass
        elseif endpoint_match
            message = "endpoint matched, but a near-singular point was recorded before λ=1"
            status = :inconclusive
        else
            message = "natural load-scale endpoint does not match the audited solution"
            status = :fail
        end
    elseif voltage_limit_stop
        status = :fail
    else
        status = :inconclusive
    end
    provenance = deepcopy(coords.provenance)
    provenance["continuation"] = Dict("homotopy" => "uniform_load_scale_0_to_1",
        "source_buses" => sort!(collect(source_buses)), "state_nodes" => meta.state_nodes,
        "natural_parameter" => true, "pseudo_arclength" => false,
        "stop_on_voltage_limit" => stop_on_voltage_limit,
        "margin" => _continuation_margin(events))
    OperabilityContinuationResult(status, message, lambdas, states, node_voltages,
        residuals, singular_values, condition_numbers, corrector_iterations, events,
        meta.state_nodes, endpoint_match, endpoint_distance, provenance)
end

"""Validated options for the pseudo-arclength operability trace."""
struct OperabilityPseudoArclengthSpec
    initial_step::Float64
    min_step::Float64
    max_step::Float64
    max_steps::Int
    corrector_tol::Float64
    max_corrector_iterations::Int
    jacobian_step::Float64
    fold_sigma_tol::Float64
    tangent_lambda_tol::Float64
    target_lambda_tol::Float64
    endpoint_atol::Float64
    endpoint_rtol::Float64
    curvature_low::Float64
    curvature_high::Float64
    function OperabilityPseudoArclengthSpec(;
            initial_step::Real=0.05,
            min_step::Real=1e-4,
            max_step::Real=0.1,
            max_steps::Integer=200,
            corrector_tol::Real=1e-7,
            max_corrector_iterations::Integer=20,
            jacobian_step::Real=1e-6,
            fold_sigma_tol::Real=1e-6,
            tangent_lambda_tol::Real=1e-3,
            target_lambda_tol::Real=5e-3,
            endpoint_atol::Real=1e-5,
            endpoint_rtol::Real=1e-5,
            curvature_low::Real=2.0,
            curvature_high::Real=5.0)
        vals = (initial_step, min_step, max_step, corrector_tol, jacobian_step,
                fold_sigma_tol, tangent_lambda_tol, target_lambda_tol,
                endpoint_atol, endpoint_rtol, curvature_low, curvature_high)
        all(x -> isfinite(Float64(x)) && Float64(x) > 0, vals) ||
            throw(ArgumentError("pseudo-arclength steps and tolerances must be finite and > 0"))
        min_step <= initial_step <= max_step || throw(ArgumentError(
            "pseudo-arclength steps must satisfy min_step <= initial_step <= max_step"))
        max_steps >= 1 || throw(ArgumentError("max_steps must be >= 1"))
        max_corrector_iterations >= 1 || throw(ArgumentError(
            "max_corrector_iterations must be >= 1"))
        curvature_low <= curvature_high || throw(ArgumentError(
            "curvature_low must be <= curvature_high"))
        new(Float64(initial_step), Float64(min_step), Float64(max_step), Int(max_steps),
            Float64(corrector_tol), Int(max_corrector_iterations), Float64(jacobian_step),
            Float64(fold_sigma_tol), Float64(tangent_lambda_tol),
            Float64(target_lambda_tol), Float64(endpoint_atol), Float64(endpoint_rtol),
            Float64(curvature_low), Float64(curvature_high))
    end
end

"""Result of an opt-in bordered-equation fold localization."""
struct OperabilityFoldResult <: AbstractSolveResult
    status::Symbol
    message::String
    lambda::Float64
    state::Vector{Float64}
    node_voltages::Dict{_Node,ComplexF64}
    residual_norm::Float64
    sigma_min::Float64
    iterations::Int
    critical_mode::Dict{String,Any}
    provenance::Dict{String,Any}
end

function solve_status(result::OperabilityFoldResult)
    publishable = result.status === :pass
    _result_solve_status(string(result.status), publishable;
        primal_status=publishable ? "FEASIBLE_POINT" : "DIAGNOSTIC_POINT")
end

solve_diagnostics(result::OperabilityFoldResult) = (
    status=result.status, lambda=result.lambda, residual_norm=result.residual_norm,
    sigma_min=result.sigma_min, iterations=result.iterations)

Base.show(io::IO, result::OperabilityFoldResult) = print(io,
    "OperabilityFoldResult(status=$(result.status), lambda=$(result.lambda), " *
    "residual=$(result.residual_norm))")

function _operability_fold_empty(status, message, provenance)
    OperabilityFoldResult(status, message, NaN, Float64[],
        Dict{_Node,ComplexF64}(), NaN, NaN, 0,
        Dict{String,Any}(), provenance)
end

function _fold_equations(net, u, meta, coords, cfg)
    n = length(meta.state)
    x_scaled = u[1:n]; λ = u[n + 1]; v = u[n + 2:end]
    lin = BMOPFTools.ybus_linearized(_operability_scale_network(net, λ);
                                    fold=:constant_z)
    state_scale, residual_scale = _operability_scales(lin, meta, coords)
    x = x_scaled .* state_scale
    f = _operability_residual(lin, meta, x) ./ residual_scale
    Jphys = finite_difference_jacobian(
        y -> _operability_residual(lin, meta, y), x; step=cfg.jacobian_step)
    J = _operability_scaled_jacobian(Jphys, state_scale, residual_scale)
    vcat(f, J * v, dot(v, v) - 1.0)
end

function _operability_solution_from_state(lin, meta, x)
    V = _operability_voltage_vector(lin, meta, x)
    buses = Dict{String,Any}()
    for (i, (bus, terminal)) in enumerate(lin.nodes)
        bus_data = get!(buses, bus, Dict{String,Any}())
        value = V[i]
        bus_data[terminal] = Dict{String,Any}("vr" => real(value), "vi" => imag(value))
    end
    Dict{String,Any}("bus" => buses)
end

"""
    locate_opf_operability_fold(net, solution; spec, lambda=1.0,
        max_iterations=30, tol=1e-8, jacobian_step=1e-6, context=nothing)

Refine a supplied fold candidate with the bordered equations ``F=0``, ``Jv=0``,
and ``‖v‖₂=1``. `solution` supplies the initial voltage state and `lambda`
specifies the load-scale at that candidate. This is a local diagnostic
localization, not a global fold or branch-discovery certificate.
"""
function locate_opf_operability_fold(
        net::Dict{String,Any}, solution::AbstractDict;
        spec::OperabilitySpec, lambda::Real=1.0, max_iterations::Integer=30,
        tol::Real=1e-8, jacobian_step::Real=1e-6, context=nothing)
    isfinite(Float64(lambda)) && lambda > 0 ||
        throw(ArgumentError("lambda must be finite and > 0"))
    max_iterations >= 1 || throw(ArgumentError("max_iterations must be >= 1"))
    isfinite(Float64(tol)) && tol > 0 || throw(ArgumentError("tol must be finite and > 0"))
    isfinite(Float64(jacobian_step)) && jacobian_step > 0 ||
        throw(ArgumentError("jacobian_step must be finite and > 0"))
    coords = _operability_coordinates(net, spec; context)
    unsupported = _operability_preflight(net)
    !isempty(unsupported) && return _operability_fold_empty(
        :not_applicable, "candidate contains physics outside the continuation scope",
        coords.provenance)
    source_buses = _operability_source_buses(net)
    isempty(source_buses) && throw(ArgumentError("fold localization requires a voltage source"))
    lin_nominal = BMOPFTools.ybus_linearized(net; fold=:constant_z)
    vmap = _operability_solution_map(solution, lin_nominal.nodes)
    meta = _operability_state_meta(lin_nominal, vmap, source_buses)
    isempty(meta.state) && return _operability_fold_empty(
        :not_applicable, "network has no free non-source voltage state", coords.provenance)
    state_scale, _ = _operability_scales(lin_nominal, meta, coords)
    J0, _, residual_scale = _continuation_scaled_jacobian(
        BMOPFTools.ybus_linearized(_operability_scale_network(net, lambda);
                                   fold=:constant_z), meta, meta.state, coords,
        OperabilityContinuationSpec(jacobian_step=jacobian_step))
    factorization = svd(J0)
    v0 = Float64.(factorization.V[:, end])
    u = vcat(meta.state ./ state_scale, Float64(lambda), v0)
    cfg = (jacobian_step=Float64(jacobian_step),)
    converged = false; last_norm = Inf; iterations = 0
    for iteration in 1:max_iterations
        iterations = iteration
        g = _fold_equations(net, u, meta, coords, cfg)
        last_norm = isempty(g) ? 0.0 : maximum(abs, g)
        last_norm <= tol && (converged = true; break)
        G = finite_difference_jacobian(
            w -> _fold_equations(net, w, meta, coords, cfg), u;
            step=Float64(jacobian_step))
        du = try
            -(G \ g)
        catch
            break
        end
        all(isfinite, du) || break
        α = 1.0; accepted = false
        while α >= 1 / 64
            candidate = u .+ α .* du
            gc = _fold_equations(net, candidate, meta, coords, cfg)
            nc = isempty(gc) ? 0.0 : maximum(abs, gc)
            if all(isfinite, gc) && nc < last_norm
                u = candidate; accepted = true; break
            end
            α /= 2
        end
        accepted || break
    end
    n = length(meta.state); λ = u[n + 1]; x = u[1:n] .* state_scale
    lin = BMOPFTools.ybus_linearized(_operability_scale_network(net, λ);
                                     fold=:constant_z)
    J, _, _ = _continuation_scaled_jacobian(
        lin, meta, x, coords, OperabilityContinuationSpec(jacobian_step=jacobian_step))
    sigma = isempty(J) ? NaN : last(svdvals(J))
    residual = _operability_complex_residual(lin, meta, x)
    residual_norm = isempty(residual) ? 0.0 : maximum(abs, residual)
    mode = isempty(J) ? Dict{String,Any}() :
        _operability_critical_mode(svd(J), meta.state_nodes)
    voltages = _operability_node_map(lin, _operability_voltage_vector(lin, meta, x))
    provenance = deepcopy(coords.provenance)
    provenance["fold_localization"] = Dict("equations" => "F=0,Jv=0,norm(v)=1",
        "source_buses" => sort!(collect(source_buses)), "state_nodes" => meta.state_nodes,
        "initial_lambda" => Float64(lambda), "iterations" => iterations)
    status = converged ? :pass : :inconclusive
    message = converged ? "bordered fold equations converged" :
        "bordered fold equations did not converge"
    OperabilityFoldResult(status, message, Float64(λ), x, voltages, residual_norm,
        Float64(sigma), iterations, mode, provenance)
end

function _pseudo_lambda_residual(net, λ, meta, x, residual_scale, h)
    if λ <= h
        lp = BMOPFTools.ybus_linearized(_operability_scale_network(net, λ + h);
                                        fold=:constant_z)
        l0 = BMOPFTools.ybus_linearized(_operability_scale_network(net, λ);
                                        fold=:constant_z)
        return (_operability_residual(lp, meta, x) - _operability_residual(l0, meta, x)) /
            h ./ residual_scale
    end
    lp = BMOPFTools.ybus_linearized(_operability_scale_network(net, λ + h);
                                    fold=:constant_z)
    lm = BMOPFTools.ybus_linearized(_operability_scale_network(net, λ - h);
                                    fold=:constant_z)
    (_operability_residual(lp, meta, x) - _operability_residual(lm, meta, x)) /
        (2h) ./ residual_scale
end

function _pseudo_arclength_jacobian(lin, meta, x, coords, cfg, arc_state_scale)
    J, state_scale, residual_scale = _continuation_scaled_jacobian(
        lin, meta, x, coords, cfg)
    # The audited Jacobian uses the caller's research scaling.  Pseudo-arclength
    # needs a separate metric so that SI volts and the dimensionless λ coordinate
    # have comparable magnitudes; transform the Jacobian columns accordingly.
    J_arc = copy(J)
    for j in axes(J_arc, 2)
        J_arc[:, j] .*= arc_state_scale[j] / state_scale[j]
    end
    J_arc, residual_scale
end

function _pseudo_corrector(net, ypred, tangent, meta, arc_state_scale, residual_scale,
                           coords, cfg::OperabilityPseudoArclengthSpec)
    # The continuation state contains real and imaginary components for every
    # free complex node.  `state_nodes` counts complex nodes, whereas `state`
    # (and the scaled Jacobian) has twice as many real entries.
    n = length(meta.state)
    y = copy(ypred)
    natural_cfg = OperabilityContinuationSpec(
        initial_step=cfg.initial_step, min_step=cfg.min_step, max_step=cfg.max_step,
        corrector_tol=cfg.corrector_tol,
        max_corrector_iterations=cfg.max_corrector_iterations,
        jacobian_step=cfg.jacobian_step, fold_sigma_tol=cfg.fold_sigma_tol,
        endpoint_atol=cfg.endpoint_atol, endpoint_rtol=cfg.endpoint_rtol)
    for iteration in 1:cfg.max_corrector_iterations
        x = y[1:n] .* arc_state_scale
        λ = y[end]
        lin = BMOPFTools.ybus_linearized(_operability_scale_network(net, λ);
                                         fold=:constant_z)
        f = _operability_residual(lin, meta, x) ./ residual_scale
        g = dot(tangent, y - ypred)
        G = vcat(f, g)
        norm_g = isempty(G) ? 0.0 : maximum(abs, G)
        norm_g <= cfg.corrector_tol && return y, true, iteration
        J, _ = _pseudo_arclength_jacobian(
            lin, meta, x, coords, natural_cfg, arc_state_scale)
        fλ = _pseudo_lambda_residual(net, λ, meta, x, residual_scale, 1e-5)
        A = [J fλ; reshape(tangent, 1, :)]
        dy = try
            -(A \ G)
        catch
            return y, false, iteration
        end
        all(isfinite, dy) || return y, false, iteration
        α = 1.0; accepted = false
        while α >= 1 / 64
            candidate = y .+ α .* dy
            xc = candidate[1:n] .* arc_state_scale; λc = candidate[end]
            linc = BMOPFTools.ybus_linearized(_operability_scale_network(net, λc);
                                               fold=:constant_z)
            fc = _operability_residual(linc, meta, xc) ./ residual_scale
            Gc = vcat(fc, dot(tangent, candidate - ypred))
            if all(isfinite, Gc) && maximum(abs, Gc) < norm_g
                y = candidate; accepted = true; break
            end
            α /= 2
        end
        accepted || return y, false, iteration
    end
    y, false, cfg.max_corrector_iterations
end

"""
    continue_opf_operability_pseudo_arclength(net, solution; spec,
        continuation=OperabilityPseudoArclengthSpec(), context=nothing)

Trace the same static load-scale equilibrium with a pseudo-arclength
predictor/corrector. The path may continue through a fold because λ is treated
as an unknown; fold candidates are recorded when the tangent reverses λ
direction or the scaled Jacobian becomes nearly singular, with critical
left/right singular-mode participation attached to those events. Endpoint crossings
are reported with their bracketing λ values and refined at fixed λ when the
crossing step does not land within `target_lambda_tol`. The result retains the
separate voltage-normalized arclength metric alongside the audited scaling
provenance, accepted arclength steps, and tangent-turning curvature history.
The step controller reduces the next step when curvature exceeds
`continuation.curvature_high`; when the corrector is cheap it uses faster
growth below `curvature_low` and conservative growth otherwise. Pass
`stop_at_target=false` to continue the stress path after λ=1;
the default remains endpoint matching. Pass `stop_on_voltage_limit=true` to
terminate with `:fail` at the first declared terminal-voltage violation.
"""
function continue_opf_operability_pseudo_arclength(
        net::Dict{String,Any}, solution::AbstractDict;
        spec::OperabilitySpec,
        continuation::OperabilityPseudoArclengthSpec=
            OperabilityPseudoArclengthSpec(), context=nothing,
        stop_at_target::Bool=true,
        stop_on_voltage_limit::Bool=false)
    coords = _operability_coordinates(net, spec; context)
    unsupported = _operability_preflight(net)
    !isempty(unsupported) && return _continuation_empty_result(
        :not_applicable, "candidate contains physics outside the continuation scope",
        coords.provenance, [Dict{String,Any}("kind" => "unsupported_physics",
                                             "reasons" => unsupported)])
    source_buses = _operability_source_buses(net)
    isempty(source_buses) && throw(ArgumentError("continuation requires a voltage source"))
    lin_target = BMOPFTools.ybus_linearized(net; fold=:constant_z)
    target_vmap = _operability_solution_map(solution, lin_target.nodes)
    meta = _operability_state_meta(lin_target, target_vmap, source_buses)
    target_x = meta.state
    isempty(target_x) && return _continuation_empty_result(
        :pass, "network has no free non-source voltage state", coords.provenance)

    natural_cfg = OperabilityContinuationSpec(
        initial_step=continuation.initial_step, min_step=continuation.min_step,
        max_step=continuation.max_step, corrector_tol=continuation.corrector_tol,
        max_corrector_iterations=continuation.max_corrector_iterations,
        jacobian_step=continuation.jacobian_step,
        fold_sigma_tol=continuation.fold_sigma_tol,
        endpoint_atol=continuation.endpoint_atol,
        endpoint_rtol=continuation.endpoint_rtol)
    net_zero = _operability_scale_network(net, 0.0)
    lin_zero = BMOPFTools.ybus_linearized(net_zero; fold=:constant_z)
    Jzero, state_scale, residual_scale = _continuation_scaled_jacobian(
        lin_zero, meta, target_x, coords, natural_cfg)
    base_x = target_x - (Jzero \ (_operability_residual(lin_zero, meta, target_x) ./
        residual_scale)) .* state_scale
    base_x, base_ok, base_iterations, _ = _continuation_corrector(
        lin_zero, meta, base_x, coords, natural_cfg)
    base_ok || return _continuation_empty_result(
        :inconclusive, "energized no-load corrector failed", coords.provenance,
        [Dict{String,Any}("kind" => "base_corrector_failure",
                           "iterations" => base_iterations)])

    base_point = _continuation_point(lin_zero, meta, base_x, coords, natural_cfg)
    Fλ = _pseudo_lambda_residual(net, 0.0, meta, base_x, residual_scale, 1e-5)
    # Keep the audited scaling for residual/Jacobian evidence, but use a
    # dimensionless voltage metric for the arclength coordinate.  This avoids
    # making SI-unit traces effectively voltage-only near a fold.
    arc_voltage_scale = max(1.0, maximum(abs, target_x))
    arc_state_scale = max.(state_scale, arc_voltage_scale)
    Jzero_arc = copy(Jzero)
    for j in axes(Jzero_arc, 2)
        Jzero_arc[:, j] .*= arc_state_scale[j] / state_scale[j]
    end
    dxλ = -(Jzero_arc \ Fλ)
    tangent = vcat(dxλ, 1.0); tangent ./= norm(tangent)
    y = vcat(base_x ./ arc_state_scale, 0.0)
    lambdas = [0.0]; states = [copy(base_x)]
    node_voltages = [base_point.node_voltages]
    residuals = [base_point.residual_norm]
    singular_values = [isempty(base_point.singular_values) ? NaN : last(base_point.singular_values)]
    condition_numbers = [base_point.condition_number]
    corrector_iterations = [base_iterations]
    events = _continuation_voltage_events(net_zero, base_point.node_voltages, 0.0, spec)
    step = continuation.initial_step
    curvature_history = Float64[]
    arclength_steps = Float64[]
    voltage_limit_stop = false
    endpoint_match = nothing; endpoint_distance = NaN; status = :inconclusive
    message = stop_at_target ?
        "pseudo-arclength trace did not produce a corrected λ=1 point" :
        "pseudo-arclength stress trace reached its step horizon"
    for _ in 1:continuation.max_steps
        step_used = step
        ypred = y .+ step .* tangent
        ynew, ok, iterations = _pseudo_corrector(
            net, ypred, tangent, meta, arc_state_scale, residual_scale, coords, continuation)
        if !ok
            push!(events, Dict{String,Any}("kind" => "corrector_failure",
                "lambda" => ypred[end], "step" => step, "iterations" => iterations))
            step /= 2
            if step < continuation.min_step
                message = "pseudo-arclength corrector step fell below min_step"
                break
            end
            continue
        end
        n = length(meta.state)
        xnew = ynew[1:n] .* arc_state_scale; λnew = ynew[end]
        lin_new = BMOPFTools.ybus_linearized(_operability_scale_network(net, λnew);
                                             fold=:constant_z)
        point = _continuation_point(lin_new, meta, xnew, coords, natural_cfg)
        λold = y[end]
        if (λold - 1.0) * (λnew - 1.0) <= 0 && λold != λnew
            push!(events, Dict{String,Any}("kind" => "target_crossing",
                "lambda_before" => λold, "lambda_after" => λnew))
            # A finite arclength step can bracket λ=1 without landing within
            # the endpoint tolerance.  Refine the crossing at fixed λ so the
            # endpoint comparison is made on an actual equilibrium, rather
            # than on a secant interpolation (which is especially important
            # when the crossing is close to a fold).
            if stop_at_target && abs(λnew - 1.0) > continuation.target_lambda_tol
                lin_endpoint = BMOPFTools.ybus_linearized(net; fold=:constant_z)
                xendpoint, endpoint_ok, endpoint_iterations, _ =
                    _continuation_corrector(lin_endpoint, meta, xnew, coords, natural_cfg)
                push!(events, Dict{String,Any}("kind" => "target_refinement",
                    "lambda" => 1.0, "iterations" => endpoint_iterations,
                    "status" => endpoint_ok ? :pass : :inconclusive))
                if endpoint_ok
                    endpoint_point = _continuation_point(
                        lin_endpoint, meta, xendpoint, coords, natural_cfg)
                    push!(lambdas, 1.0); push!(states, copy(xendpoint))
                    push!(node_voltages, endpoint_point.node_voltages)
                    push!(residuals, endpoint_point.residual_norm)
                    push!(singular_values, isempty(endpoint_point.singular_values) ?
                        NaN : last(endpoint_point.singular_values))
                    push!(condition_numbers, endpoint_point.condition_number)
                    push!(corrector_iterations, endpoint_iterations)
                    endpoint_voltage_events = _continuation_voltage_events(
                        net, endpoint_point.node_voltages, 1.0, spec)
                    append!(events, endpoint_voltage_events)
                    voltage_limit_stop = stop_on_voltage_limit &&
                        !isempty(endpoint_voltage_events)
                    endpoint_distance = maximum(
                        abs(endpoint_point.node_voltages[node] - target_vmap[node])
                        for node in keys(target_vmap))
                    endpoint_limit = continuation.endpoint_atol +
                        continuation.endpoint_rtol * max(
                            1.0, maximum(abs, values(target_vmap)))
                    endpoint_match = endpoint_distance <= endpoint_limit
                    status = endpoint_match ? :pass : :fail
                    message = endpoint_match ?
                        "pseudo-arclength trace refined and matched the λ=1 endpoint" :
                        "pseudo-arclength trace refined λ=1 but did not match the audited endpoint"
                    break
                end
            end
        end
        push!(lambdas, λnew); push!(states, copy(xnew)); push!(node_voltages, point.node_voltages)
        push!(residuals, point.residual_norm)
        push!(singular_values, isempty(point.singular_values) ? NaN : last(point.singular_values))
        push!(condition_numbers, point.condition_number); push!(corrector_iterations, iterations)
        voltage_events = _continuation_voltage_events(
            _operability_scale_network(net, λnew), point.node_voltages, λnew, spec)
        append!(events, voltage_events)
        if stop_on_voltage_limit && !isempty(voltage_events)
            voltage_limit_stop = true
            message = "pseudo-arclength trace stopped at a declared voltage limit"
        end
        sigma = isempty(point.singular_values) ? NaN : last(point.singular_values)
        fλ_new = _pseudo_lambda_residual(net, λnew, meta, xnew, residual_scale, 1e-5)
        Jnew, _ = _pseudo_arclength_jacobian(
            lin_new, meta, xnew, coords, natural_cfg, arc_state_scale)
        factorization_new = isempty(Jnew) ? nothing : svd(Jnew)
        near_singular = isfinite(sigma) && sigma <= continuation.fold_sigma_tol
        near_singular && push!(events, Dict{String,Any}(
            "kind" => "near_singular", "lambda" => λnew, "sigma_min" => sigma,
            "critical_mode" => _operability_critical_mode(
                factorization_new, meta.state_nodes)))
        dxλ_new = try -(Jnew \ fλ_new) catch; zeros(length(xnew)) end
        tangent_new = vcat(dxλ_new, 1.0); tangent_new ./= max(norm(tangent_new), eps())
        dot(tangent_new, tangent) < 0 && (tangent_new .*= -1)
        fold_event = nothing
        if abs(tangent[end]) > continuation.tangent_lambda_tol &&
            abs(tangent_new[end]) > continuation.tangent_lambda_tol &&
            signbit(tangent[end]) != signbit(tangent_new[end])
            fold_event = Dict{String,Any}("kind" => "fold_candidate",
                "lambda" => λnew, "tangent_lambda_before" => tangent[end],
                "tangent_lambda_after" => tangent_new[end])
        elseif abs(tangent_new[end]) <= continuation.tangent_lambda_tol
            fold_event = Dict{String,Any}("kind" => "fold_candidate",
                "lambda" => λnew, "tangent_lambda" => tangent_new[end])
        end
        if fold_event !== nothing
            fold_event["sigma_min"] = sigma
            fold_event["critical_mode"] = _operability_critical_mode(
                factorization_new, meta.state_nodes)
            fold_solution = _operability_solution_from_state(lin_new, meta, xnew)
            localized = locate_opf_operability_fold(net, fold_solution;
                spec=spec, lambda=λnew, context=context)
            fold_event["fold_localization"] = Dict{String,Any}(
                "status" => localized.status, "lambda" => localized.lambda,
                "residual_norm" => localized.residual_norm,
                "sigma_min" => localized.sigma_min,
                "iterations" => localized.iterations)
            push!(events, fold_event)
        end
        turning_angle = acos(clamp(real(dot(tangent, tangent_new)), -1.0, 1.0))
        curvature = turning_angle / max(step_used, eps(Float64))
        push!(curvature_history, curvature)
        push!(arclength_steps, step_used)
        y = ynew; tangent = tangent_new
        if curvature >= continuation.curvature_high
            step = max(continuation.min_step, step / 1.5)
        elseif iterations <= 4
            growth = curvature <= continuation.curvature_low ? 1.3 : 1.15
            step = min(continuation.max_step, step * growth)
        elseif iterations >= continuation.max_corrector_iterations ÷ 2
            step = max(continuation.min_step, step / 1.5)
        end
        if stop_at_target && abs(λnew - 1.0) <= continuation.target_lambda_tol
            endpoint_distance = maximum(abs(point.node_voltages[node] - target_vmap[node])
                                        for node in keys(target_vmap))
            endpoint_limit = continuation.endpoint_atol + continuation.endpoint_rtol * max(
                1.0, maximum(abs, values(target_vmap)))
            endpoint_match = endpoint_distance <= endpoint_limit
            status = voltage_limit_stop ? :fail : (endpoint_match ? :pass : :fail)
            message = voltage_limit_stop ?
                "pseudo-arclength trace reached a declared voltage limit" :
                (endpoint_match ?
                "pseudo-arclength trace reached and matched a λ=1 endpoint" :
                "pseudo-arclength trace reached λ=1 but did not match the audited endpoint")
            break
        end
        voltage_limit_stop && break
    end
    voltage_limit_stop && status === :inconclusive && (status = :fail)
    provenance = deepcopy(coords.provenance)
    provenance["continuation"] = Dict("homotopy" => "uniform_load_scale_0_to_1",
        "source_buses" => sort!(collect(source_buses)), "state_nodes" => meta.state_nodes,
        "natural_parameter" => false, "pseudo_arclength" => true,
        "arclength_state_scale" => arc_state_scale,
        "arclength_steps" => arclength_steps,
        "curvature_history" => curvature_history,
        "curvature_control" => Dict("low" => continuation.curvature_low,
                                     "high" => continuation.curvature_high),
        "stop_at_target" => stop_at_target,
        "stop_on_voltage_limit" => stop_on_voltage_limit,
        "margin" => _continuation_margin(events))
    OperabilityContinuationResult(status, message, lambdas, states, node_voltages,
        residuals, singular_values, condition_numbers, corrector_iterations, events,
        meta.state_nodes, endpoint_match, endpoint_distance, provenance)
end
