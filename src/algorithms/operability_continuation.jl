# Static load-scale continuation for the post-OPF operability processor.
#
# This first slice follows the natural loading parameter λ from the energized
# no-load state to λ=1. It records corrector failures and near-singular points,
# but deliberately does not call a natural-parameter trace a pseudo-arclength
# certificate. A full fold-crossing continuation engine remains a follow-up.

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
pseudo-arclength reachability certificate.
"""
function continue_opf_operability(net::Dict{String,Any}, solution::AbstractDict;
                                  spec::OperabilitySpec,
                                  continuation::OperabilityContinuationSpec=
                                      OperabilityContinuationSpec(),
                                  context=nothing)
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
        append!(events, _continuation_voltage_events(net_trial, point.node_voltages, λtrial, spec))
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
    if reached
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
    else
        status = :inconclusive
    end
    provenance = deepcopy(coords.provenance)
    provenance["continuation"] = Dict("homotopy" => "uniform_load_scale_0_to_1",
        "source_buses" => sort!(collect(source_buses)), "state_nodes" => meta.state_nodes,
        "natural_parameter" => true, "pseudo_arclength" => false)
    OperabilityContinuationResult(status, message, lambdas, states, node_voltages,
        residuals, singular_values, condition_numbers, corrector_iterations, events,
        meta.state_nodes, endpoint_match, endpoint_distance, provenance)
end
