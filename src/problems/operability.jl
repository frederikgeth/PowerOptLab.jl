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

"""Expected network/solution data failure while auditing one snapshot."""
struct OperabilityModelError <: Exception
    message::String
end

Base.showerror(io::IO, err::OperabilityModelError) = print(io, err.message)

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
Bernstein-style Z-bus contraction condition. The certificate covers the
native constant-P, constant-I, ZIP, exponential, and constant-Z load laws
when their declared voltage-domain parameters are finite; a failed sufficient
condition is `:inconclusive`, never a collapse or multiplicity claim.
The same report includes `fixed_point_local_region`, a separate
candidate-centered contraction result for local uniqueness evidence; it does
not establish reachability from the energized no-load germ and does not
downgrade the aggregate snapshot status. `fixed_point_euclidean_region`
provides an alternative no-load-connected sufficient condition using a
conservative Euclidean ball in the complex free-voltage state.
For large snapshots, `jacobian_spectrum=:extremes` estimates only the largest
and smallest singular values with deterministic iterative solves; the default
`:full` mode retains the complete SVD and critical-mode participation.
Set `jacobian_storage=:sparse` together with `jacobian_spectrum=:extremes` to
retain the finite-difference Jacobian as a sparse matrix and use sparse LU
right-hand-side solves; this reduced path omits the full SVD by design.
Set `record_jacobian_pattern=true` when sparse row/column provenance is needed;
the default `false` avoids retaining two explicit integer arrays in the sparse
large-network path.
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
    jacobian_spectrum::Symbol
    jacobian_storage::Symbol
    record_jacobian_pattern::Bool
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
            jacobian_spectrum::Symbol=:full,
            jacobian_storage::Symbol=:dense,
            record_jacobian_pattern::Bool=false,
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
        jacobian_spectrum in (:full, :extremes) || throw(ArgumentError(
            "jacobian_spectrum must be :full or :extremes"))
        jacobian_storage in (:dense, :sparse) || throw(ArgumentError(
            "jacobian_storage must be :dense or :sparse"))
        jacobian_storage === :sparse && jacobian_spectrum !== :extremes &&
            throw(ArgumentError("jacobian_storage=:sparse requires jacobian_spectrum=:extremes"))
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
            Float64(jacobian_rank_rtol), jacobian_spectrum, jacobian_storage,
            record_jacobian_pattern,
            Float64(sensitivity_step),
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
    jacobian::AbstractMatrix{Float64}
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

    # A non-SI audit must not silently fall back to unit bases for buses that
    # were omitted from the caller's coordinate contract.  SI audits may use
    # the unit defaults, but every bus needs an explicit base under any other
    # policy (or a context that supplies one for each bus).
    if !(policy isa BMOPFTools.SIUnitsScaling)
        for (bus, _) in get(net, "bus", Dict())
            haskey(bases, String(bus)) || throw(ArgumentError(
                "scaling_bases is missing voltage/current bases for bus $(repr(String(bus)))"))
        end
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
    haskey(buses, bus) || throw(OperabilityModelError(
        "solution is missing bus $(repr(bus))"))
    data = get(buses[bus], terminal, nothing)
    data === nothing && throw(OperabilityModelError(
        "solution is missing voltage for terminal $(repr(node))"))
    vr = Float64(get(data, "vr", NaN)); vi = Float64(get(data, "vi", NaN))
    isfinite(vr) && isfinite(vi) || throw(OperabilityModelError(
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

function _operability_scope_source_inventory(net::Dict{String,Any})
    bus_ids = Set{String}(String(id) for id in keys(get(net, "bus", Dict())))
    source_buses = sort!(collect(_operability_source_buses(net)))
    missing_source_buses = sort!(filter(bus -> !(bus in bus_ids), source_buses))
    has_voltage_source = !isempty(source_buses) && isempty(
        missing_source_buses)
    Dict{String,Any}(
        "bus_count" => length(bus_ids),
        "voltage_source_count" => length(get(net, "voltage_source", Dict())),
        "source_buses" => source_buses,
        "missing_source_buses" => missing_source_buses,
        "has_voltage_source" => has_voltage_source,
    )
end

"""
    operability_scope_audit(net)

Return a serialisable preflight record for the native frozen-dispatch
operability checker. This audit does not inspect a solved voltage state and does
not claim that an out-of-scope network is infeasible; it identifies the exact
model seam that would need to be extended before a snapshot can be checked.
Scope-audit `status` deliberately uses its own two-value readiness vocabulary,
`:supported` or `:not_applicable`; it is not an `OperabilityCheck` verdict.
"""
function operability_scope_audit(net::Dict{String,Any})
    reasons = _operability_preflight(net)
    inventory = _operability_scope_inventory(net)
    source_inventory = _operability_scope_source_inventory(net)
    if source_inventory["voltage_source_count"] == 0
        push!(reasons, "network has no voltage source for a native ybus_linearized equilibrium")
    elseif !isempty(source_inventory["missing_source_buses"])
        push!(reasons, "voltage source references missing bus(es) " *
            repr(source_inventory["missing_source_buses"]))
    end
    generator_count = length(get(net, "generator", Dict()))
    ibr_count = length(get(net, "ibr", Dict()))
    control_closure = generator_count == 0 && ibr_count == 0 ?
        "frozen_dispatch_native_static" : "outside_native_static_seam"
    Dict{String,Any}(
        "status" => isempty(reasons) ? :supported : :not_applicable,
        "scope" => "native_ybus_linearized",
        "closure" => :frozen_dispatch,
        "control_closure" => control_closure,
        "topology" => source_inventory,
        "generator_count" => generator_count,
        "ibr_count" => ibr_count,
        "model_inventory" => inventory,
        "unsupported_reasons" => reasons,
        "message" => isempty(reasons) ?
            "network is eligible for the native frozen-dispatch snapshot checker" :
            "network requires an extended equilibrium seam before operability claims can be evaluated")
end

"""Named finite stress direction for native static-load campaigns."""
struct OperabilityStressDirection
    name::Symbol
    p_scale::Float64
    q_scale::Float64
    load_weights::Dict{String,Float64}
    connection_weights::Dict{String,Vector{Float64}}
    function OperabilityStressDirection(name::Symbol=:uniform_load;
            p_scale::Real=1.0, q_scale::Real=1.0,
            load_weights=Dict{String,Float64}(), connection_weights=Dict())
        isfinite(Float64(p_scale)) && isfinite(Float64(q_scale)) ||
            throw(ArgumentError("stress-direction P/Q scales must be finite"))
        lw = Dict{String,Float64}()
        for (id, weight) in load_weights
            isfinite(Float64(weight)) || throw(ArgumentError(
                "stress-direction load weights must be finite"))
            lw[String(id)] = Float64(weight)
        end
        cw = Dict{String,Vector{Float64}}()
        for (id, weights) in connection_weights
            values = Float64.(collect(weights))
            all(isfinite, values) || throw(ArgumentError(
                "stress-direction connection weights must be finite"))
            isempty(values) && throw(ArgumentError(
                "stress-direction connection weights cannot be empty"))
            cw[String(id)] = values
        end
        new(name, Float64(p_scale), Float64(q_scale), lw, cw)
    end
end

function _operability_stress_network(net::Dict{String,Any}, λ::Real,
                                     direction::OperabilityStressDirection)
    result = deepcopy(net)
    scale = Float64(λ)
    isfinite(scale) && scale >= 0.0 || throw(ArgumentError(
        "stress λ must be finite and >= 0"))
    for (id_raw, load) in get(result, "load", Dict())
        id = String(id_raw)
        load_weight = get(direction.load_weights, id, 1.0)
        connection_weight = get(direction.connection_weights, id, nothing)
        nconn = _operability_load_connection_count(load)
        if connection_weight !== nothing && length(connection_weight) < nconn
            throw(ArgumentError("stress connection weights for $(repr(id)) need at least $nconn entries"))
        end
        for (key, factor) in (("p_nom", direction.p_scale),
                              ("q_nom", direction.q_scale))
            value = get(load, key, nothing)
            value === nothing && continue
            values = value isa AbstractVector ? Float64.(value) : [Float64(value)]
            load[key] = [scale * factor * load_weight *
                         (connection_weight === nothing ? 1.0 :
                          (k <= length(connection_weight) ? connection_weight[k] : 0.0)) *
                         values[k] for k in eachindex(values)]
        end
    end
    result
end

function _operability_scale_network(net::Dict{String,Any}, λ::Real)
    result = deepcopy(net)
    scale = Float64(λ)
    isfinite(scale) || throw(ArgumentError("load-scale λ must be finite"))
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

function _operability_directional_sensitivity(net, lin, meta, x, J, linear_solver,
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
    dx_scaled = linear_solver === nothing ? -(J \ dF_scaled) :
        -(linear_solver \ dF_scaled)
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

function _operability_sparse_finite_difference_jacobian(map, point; step::Real)
    x = Float64.(point)
    h = Float64(step)
    isfinite(h) && h > 0.0 || throw(ArgumentError("step must be finite and > 0"))
    y = Float64.(map(x))
    n = length(x)
    length(y) == n || throw(DimensionMismatch(
        "finite-difference map output length must equal input length"))
    rows = Int[]; columns = Int[]; values = Float64[]
    for column in 1:n
        delta = h * max(abs(x[column]), 1.0)
        plus = copy(x); plus[column] += delta
        minus = copy(x); minus[column] -= delta
        difference = (Float64.(map(plus)) - Float64.(map(minus))) / (2delta)
        length(difference) == n || throw(DimensionMismatch(
            "finite-difference map output length changed"))
        for row in eachindex(difference)
            value = difference[row]
            value == 0.0 && continue
            push!(rows, row); push!(columns, column); push!(values, value)
        end
    end
    sparse(rows, columns, values, n, n)
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
    if issparse(J)
        left = spdiagm(0 => 1.0 ./ residual_scale)
        right = spdiagm(0 => state_scale)
        return left * J * right
    end
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

"""Estimate the extreme singular values without materializing a full SVD.

The largest value uses power iteration on ``J'J``. The smallest value uses
inverse iteration through a caller-supplied LU factorization. This is a
reduced diagnostic, not a replacement for the full singular spectrum or its
left/right mode participation.
"""
function _operability_extreme_singular_values(J, factorization;
                                             max_iterations::Int=80,
                                             rtol::Float64=1e-7)
    isempty(J) && return Float64[]
    n = size(J, 2)
    n == 0 && return Float64[]
    lanczos_extreme(apply) = begin
        q = fill(inv(sqrt(Float64(n))), n)
        limit = min(max_iterations, n)
        basis = zeros(Float64, n, limit)
        alpha = zeros(Float64, limit); beta = zeros(Float64, limit)
        iterations = 0
        breakdown = false
        for k in 1:limit
            iterations = k
            basis[:, k] = q
            z = apply(q)
            a = dot(q, z)
            z .-= a .* q
            k > 1 && (z .-= basis[:, 1:k-1] * (basis[:, 1:k-1]' * z))
            b = norm(z)
            alpha[k] = a
            if b <= eps(Float64) * max(1.0, norm(apply(q)))
                breakdown = true
                break
            end
            beta[k] = b
            q = z / b
        end
        iterations == 0 && return nothing
        alpha_view = alpha[1:iterations]
        beta_view = iterations > 1 ? beta[1:iterations-1] : Float64[]
        eig = eigen(SymTridiagonal(alpha_view, beta_view))
        index = argmax(eig.values)
        residual = breakdown ? 0.0 : abs(beta[iterations] * eig.vectors[end, index])
        scale = max(1.0, abs(eig.values[index]))
        residual <= rtol * scale ? (value=eig.values[index], residual=residual) : nothing
    end
    maximum_result = lanczos_extreme(x -> J' * (J * x))
    maximum_result === nothing && return nothing
    factorization === nothing && return nothing
    minimum_inverse = lanczos_extreme(x -> factorization \ (factorization' \ x))
    minimum_inverse === nothing && return nothing
    maximum_result.value > 0.0 && minimum_inverse.value > 0.0 ?
        [sqrt(maximum_result.value), inv(sqrt(minimum_inverse.value))] : nothing
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
    supported_models = Set(("constant_power", "constant_current", "constant_impedance",
                            "zip", "exponential"))
    for (id_raw, load) in get(net, "load", Dict())
        model = lowercase(string(get(load, "model", "constant_power")))
        model in supported_models || push!(reasons,
            "load $(repr(String(id_raw))) uses unsupported model $(repr(model))")
    end
    assumptions = [
        "frozen-dispatch native ybus_linearized equilibrium",
        "source terminals are fixed and all non-source voltage nodes are included",
        "non-constant-Z compensation currents are represented by connection incidence",
        "constant-impedance parts are folded into the linear Ybus",
        "power/current/exponential current laws use conservative real Lipschitz bounds",
        "the contraction evidence uses connection-aware componentwise polydisc and Euclidean-ball geometries around the no-load solution",
    ]
    base = Dict{String,Any}(
        "status" => :not_applicable,
        "method" => "bernstein_style_zbus_contraction",
        "theorem" => "sufficient existence, uniqueness, and Jacobian nonsingularity condition",
        "assumptions" => assumptions,
        "scope" => "native_static_load_laws",
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
    Yll = ComplexF64.(lin.Y[meta.free, meta.free])
    Yls = isempty(meta.fixed) ? zeros(ComplexF64, nf, 0) :
        ComplexF64.(lin.Y[meta.free, meta.fixed])
    Yll_factor = try
        lu(Yll)
    catch err
        base["status"] = :inconclusive
        base["message"] = "source-eliminated Ybus is singular; contraction region is unavailable"
        base["error"] = sprint(showerror, err)
        return base
    end
    zsolve(rhs) = Yll_factor \ ComplexF64.(rhs)
    base["zbus_storage_mode"] = issparse(Yll) ?
        "implicit_sparse_factorization" : "implicit_dense_factorization"
    w = isempty(meta.fixed) ? zeros(ComplexF64, nf) :
        -(Yll_factor \ (Yls * meta.fixed_values))
    wfull = zeros(ComplexF64, length(lin.nodes))
    wfull[meta.fixed] = meta.fixed_values
    wfull[meta.free] = w
    local_index = Dict{Int,Int}(gi => li for
        (li, gi) in enumerate(meta.free))
    edges = NamedTuple[]
    for (id_raw, load) in sort!(collect(get(net, "load", Dict())); by=x -> String(x[1]))
        id = String(id_raw)
        for (k, sl) in enumerate(_load_subloads(load, net))
            # The constant-Z component has already been absorbed into Y. The
            # remaining terms are represented as current-law Lipschitz terms.
            terms = NamedTuple[]
            s = ComplexF64(sl.pt.cc + im * sl.qt.cc)
            abs(s) > 0.0 && push!(terms,
                (kind=:constant_power, coefficient=abs(s), gamma=0.0))
            current_coefficient = abs(sl.pt.cs + im * sl.qt.cs)
            current_coefficient > 0.0 && push!(terms,
                (kind=:constant_current, coefficient=current_coefficient, gamma=1.0))
            for component in (sl.pt, sl.qt)
                component.nl === nothing && continue
                coefficient, gamma = component.nl
                if !(isfinite(Float64(coefficient)) && isfinite(Float64(gamma)))
                    push!(reasons, "load $(repr(id))/$k has non-finite exponential parameters")
                    continue
                end
                coefficient_abs = abs(Float64(coefficient)) / sl.Vnom^Float64(gamma)
                coefficient_abs > 0.0 && push!(terms,
                    (kind=:exponential, coefficient=coefficient_abs, gamma=Float64(gamma)))
            end
            isempty(terms) && continue
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
            push!(edges, (key="$id/$k", pos=sl.pos, neg=sl.neg, sl=sl,
                terms=terms, dv0=dv0, a=a, anorm=anorm,
                zalpha=abs.(zsolve(a))))
        end
    end
    if !isempty(reasons)
        base["status"] = :not_applicable
        base["reasons"] = reasons
        base["message"] = "fixed-point certificate parameters are outside the supported finite domain"
        return base
    end
    base["state_nodes"] = copy(meta.state_nodes)
    base["coordinate_order"] = "complex free nodes in LinearizedYbus order"
    base["no_load_voltage"] = copy(w)
    base["edge_count"] = length(edges)
    base["connection_law_terms"] = Dict(
        e.key => [String(term.kind) for term in e.terms] for e in edges)
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
    term_bounds(term, lower, upper) = begin
        coefficient = term.coefficient
        if term.kind === :constant_power
            (coefficient / lower, coefficient / lower^2)
        elseif term.kind === :constant_current
            (coefficient, coefficient / lower)
        else
            magnitude_exponent = term.gamma - 1.0
            derivative_exponent = term.gamma - 2.0
            magnitude_radius = magnitude_exponent >= 0.0 ? upper : lower
            derivative_radius = derivative_exponent >= 0.0 ? upper : lower
            # In the radial/tangential frame of an exponential current law,
            # the two real Jacobian singular values are |γ-1| and 1 times
            # |C|m^(γ-2); the operator norm is therefore max(1, |γ-1|),
            # including the constant-power (γ=0) and constant-current (γ=1)
            # limits.
            (coefficient * magnitude_radius^magnitude_exponent,
             coefficient * max(1.0, abs(term.gamma - 1.0)) *
                 derivative_radius^derivative_exponent)
        end
    end
    component_contribution(radii) = begin
        lower = [abs(e.dv0) - sum(abs.(e.a) .* radii)
                 for e in edges]
        upper = [abs(e.dv0) + sum(abs.(e.a) .* radii)
                 for e in edges]
        any(v -> v <= 0.0, lower) &&
            return (fill(Inf, nf), Inf, lower, upper)
        bounds = [term_bounds(term, lower[j], upper[j])
                  for (j, e) in enumerate(edges) for term in e.terms]
        offsets = cumsum(vcat(0, [length(e.terms) for e in edges]))
        b = [sum(e.zalpha[i] * bounds[offsets[j] + t][1]
                 for (j, e) in enumerate(edges) for t in eachindex(e.terms))
             for i in 1:nf]
        q = maximum(
            sum(e.zalpha[i] * bounds[offsets[j] + t][2]
                * sum(abs.(e.a) .* radii) / radii[i]
                for (j, e) in enumerate(edges) for t in eachindex(e.terms))
            for i in 1:nf)
        b, q, lower, upper
    end
    contribution(rho) = begin
        lower = [abs(e.dv0) - e.anorm * rho for e in edges]
        upper = [abs(e.dv0) + e.anorm * rho for e in edges]
        any(v -> v <= 0.0, lower) && return (Inf, Inf, lower)
        bounds = [term_bounds(term, lower[j], upper[j])
                  for (j, e) in enumerate(edges) for term in e.terms]
        offsets = cumsum(vcat(0, [length(e.terms) for e in edges]))
        b = maximum(sum(e.zalpha[i] * bounds[offsets[j] + t][1]
                        for (j, e) in enumerate(edges) for t in eachindex(e.terms))
                    for i in 1:nf)
        q = maximum(sum(e.zalpha[i] * e.anorm * bounds[offsets[j] + t][2]
                        for (j, e) in enumerate(edges) for t in eachindex(e.terms))
                    for i in 1:nf)
        b, q, lower
    end
    # A coarse scan identifies the last feasible interval; bisection then
    # refines its upper boundary. The full predicate need not be globally
    # monotone because voltage-domain and contraction terms can trade off, so
    # the scan remains the guard against disconnected feasible bands. The
    # refinement is only applied upward from the last feasible grid point and
    # costs O(coarse_points + refinement_iterations) evaluations instead of the
    # former fixed 2001-point sweep.
    largest_feasible_scalar = function(start, stop, evaluate, feasible;
                                        coarse_points::Int=33,
                                        score=nothing,
                                        refinement_iterations::Int=12)
        start < stop || return nothing
        grid = collect(range(start, stop=stop, length=coarse_points))
        values = [evaluate(value) for value in grid]
        score_value(value) = begin
            score === nothing || value === nothing ? -Inf : begin
                score_candidate = Float64(score(value))
                isfinite(score_candidate) ? score_candidate : -Inf
            end
        end
        max_score = maximum(score_value(value) for value in values)
        feasible_indices = findall(feasible, values)
        isempty(feasible_indices) && return nothing
        index = last(feasible_indices)
        best = values[index]
        lo_scalar = grid[index]
        hi_scalar = index == length(grid) ? stop : grid[index + 1]
        if index == length(grid)
            endpoint = evaluate(hi_scalar)
            max_score = max(max_score, score_value(endpoint))
            if feasible(endpoint)
                return score === nothing ? endpoint : (value=endpoint, max_score=max_score)
            end
        end
        for _ in 1:refinement_iterations
            midpoint = (lo_scalar + hi_scalar) / 2.0
            candidate = evaluate(midpoint)
            max_score = max(max_score, score_value(candidate))
            if feasible(candidate)
                lo_scalar = midpoint
                best = candidate
            else
                hi_scalar = midpoint
            end
        end
        score === nothing ? best : (value=best, max_score=max_score)
    end
    # A second sufficient geometry uses a Euclidean ball in the complex free
    # voltage state.  For a state perturbation of norm ρ, connection j moves
    # by at most ||aⱼ||₂ρ.  Combining those bounds with |Z| and the scalar
    # load-law derivative bounds gives a conservative output Lipschitz factor
    # without assuming a componentwise box shape.
    euclidean_region = function(center_free, center_full, center_dv)
        offset = abs.(w + zsolve(lin.i_comp(center_full)[meta.free]) - center_free)
        offset_norm = norm(offset)
        connection_norms = Float64[norm(e.a) for e in edges]
        radius_limit = minimum(abs(center_dv[j]) / connection_norms[j]
                               for j in eachindex(edges)
                               if connection_norms[j] > 0.0)
        isfinite(radius_limit) && radius_limit > 0.0 || return nothing
        evaluate = rho -> begin
            lower = [abs(center_dv[j]) - connection_norms[j] * rho
                     for j in eachindex(edges)]
            upper = [abs(center_dv[j]) + connection_norms[j] * rho
                     for j in eachindex(edges)]
            any(v -> v <= 0.0, lower) && return nothing
            derivative_bounds = [sum(term_bounds(term, lower[j], upper[j])[2]
                                     for term in e.terms)
                                 for (j, e) in enumerate(edges)]
            # Each connection contributes through the incidence vector aⱼ;
            # summing |Z aⱼ| avoids assuming one nonlinear edge per free node
            # (e.g. a floating-neutral WYE load has one edge and four states).
            output_bound = zeros(Float64, nf)
            for (j, e) in enumerate(edges)
                # `zalpha = |Z aⱼ|` is precomputed when the connection edge is
                # built; reusing it keeps the scan linear in the number of
                # edges and free states instead of repeating dense matvecs.
                # A connection perturbation is bounded by ||a_j||₂ times the
                # free-state perturbation.  This incidence factor is essential
                # for phase-to-neutral edges with a floating neutral (and for
                # phase-to-phase/delta edges), where ||a_j||₂ can exceed one.
                output_bound .+= e.zalpha .* derivative_bounds[j] * connection_norms[j]
            end
            q = norm(output_bound)
            invariance = offset_norm + q * rho
            isfinite(q) && q < 1.0 && invariance <= rho || return nothing
            (radius=rho, q=q, offset_norm=offset_norm,
             lower=lower, upper=upper, invariance_bound=invariance)
        end
        largest_feasible_scalar(max(1e-9, radius_limit / 2000.0),
            0.999radius_limit, evaluate,
            value -> value !== nothing,
            score=value -> 1.0 - value.q)
    end
    euclidean_scan = euclidean_region(w, wfull, ComplexF64[e.dv0 for e in edges])
    euclidean = euclidean_scan === nothing ? nothing : euclidean_scan.value
    euclidean_max_condition_margin = euclidean_scan === nothing ? NaN :
        euclidean_scan.max_score
    selected_scan = largest_feasible_scalar(0.0, 0.999radius_limit,
        rho -> begin
            b, q, lower = contribution(rho)
            isfinite(b) && isfinite(q) && b <= rho && q < 1.0 || return nothing
            (rho=rho, b=b, q=q, lower=lower)
        end,
        value -> value !== nothing,
        score=value -> 1.0 - value.q)
    selected = selected_scan === nothing ? nothing : selected_scan.value
    selected_max_condition_margin = selected_scan === nothing ? NaN :
        selected_scan.max_score
    # A one-parameter family of componentwise polydiscs preserves the
    # connection-aware voltage geometry while allowing weak phases/nodes to
    # receive larger radii than the scalar uniform search. The scalar alpha
    # search keeps the certificate deterministic and auditable.
    b0, _, _, _ = component_contribution(fill(0.0, nf) .+ eps(Float64))
    componentwise = nothing
    if all(isfinite, b0) && any(>(0.0), b0)
        direction = max.(b0, eps(Float64))
        alpha_limit = minimum(abs(e.dv0) /
            sum(abs.(e.a) .* direction) for e in edges if
            sum(abs.(e.a) .* direction) > 0.0)
        if isfinite(alpha_limit) && alpha_limit > 1.0
            componentwise = largest_feasible_scalar(1.0, 0.999alpha_limit,
                alpha -> begin
                    radii = alpha .* direction
                    b, q, lower, upper = component_contribution(radii)
                    all(isfinite, b) && isfinite(q) && all(b .<= radii) && q < 1.0 ||
                        return nothing
                    (radii=radii, b=b, q=q, lower=lower, upper=upper, alpha=alpha)
                end,
                value -> value !== nothing,
                score=value -> 1.0 - value.q)
        end
    end
    componentwise_max_condition_margin = componentwise === nothing ? NaN :
        componentwise.max_score
    componentwise = componentwise === nothing ? nothing : componentwise.value
    candidate = _operability_voltage_vector(lin, meta, meta.state)
    candidate_delta = abs.(candidate[meta.free] .- w)
    if euclidean === nothing
        base["euclidean_region"] = Dict{String,Any}(
            "status" => :inconclusive,
            "message" => "no invariant Euclidean contraction ball was found")
    else
        euclidean_ratio = norm(candidate_delta) / euclidean.radius
        base["euclidean_region"] = Dict{String,Any}(
            "status" => euclidean_ratio <= 1.0 + spec.residual_atol ? :pass : :inconclusive,
            "radius" => euclidean.radius,
            "contraction_factor" => euclidean.q,
            "condition_margin" => 1.0 - euclidean.q,
            "max_condition_margin" => euclidean_max_condition_margin,
            "invariance_bound" => euclidean.invariance_bound,
            "invariance_form" => "offset_norm + q * radius <= radius in the Euclidean state norm",
            "candidate_ratio" => euclidean_ratio,
            "candidate_inside_region" => euclidean_ratio <= 1.0 + spec.residual_atol,
            "lower_connection_voltages" => euclidean.lower,
            "upper_connection_voltages" => euclidean.upper,
            "interpretation" => "conservative Euclidean-ball contraction around the energized no-load solution")
    end
    # In addition to the no-load-connected region, evaluate a local
    # candidate-centered contraction region.  Its invariance bound includes
    # the fixed-point residual at the candidate, so it can certify local
    # uniqueness even when the candidate is outside the no-load-connected
    # region.  This is deliberately reported as separate evidence: it says
    # nothing about reachability from the energized germ.
    local_region = try
        candidate_free = candidate[meta.free]
        center_full = copy(wfull); center_full[meta.free] = candidate_free
        map_center = w + zsolve(lin.i_comp(center_full)[meta.free])
        offset = abs.(map_center .- candidate_free)
        local_dv0 = [begin
            pi = get(lin.index, e.pos, 0)
            ni = e.neg === nothing ? 0 : get(lin.index, e.neg, 0)
            vp = pi == 0 ? 0.0im : center_full[pi]
            vn = ni == 0 ? 0.0im : center_full[ni]
            vp - vn
        end for e in edges]
        local_limit = minimum(abs(local_dv0[j]) / e.anorm
                              for (j, e) in enumerate(edges)
                              if e.anorm > 0.0)
        if !(isfinite(local_limit) && local_limit > 0.0)
            nothing
        else
            local_contribution(radii) = begin
                lower = [abs(local_dv0[j]) - sum(abs.(e.a) .* radii)
                         for (j, e) in enumerate(edges)]
                upper = [abs(local_dv0[j]) + sum(abs.(e.a) .* radii)
                         for (j, e) in enumerate(edges)]
                any(v -> v <= 0.0, lower) && return (fill(Inf, nf), Inf, lower, upper)
                derivative_bounds = [term_bounds(term, lower[j], upper[j])[2]
                                     for (j, e) in enumerate(edges) for term in e.terms]
                offsets = cumsum(vcat(0, [length(e.terms) for e in edges]))
                edge_delta = [sum(abs.(e.a) .* radii) for e in edges]
                b = [offset[i] + sum(e.zalpha[i] *
                        derivative_bounds[offsets[j] + t] * edge_delta[j]
                        for (j, e) in enumerate(edges) for t in eachindex(e.terms))
                     for i in 1:nf]
                q = maximum(sum(e.zalpha[i] *
                        derivative_bounds[offsets[j] + t] * edge_delta[j] / radii[i]
                        for (j, e) in enumerate(edges) for t in eachindex(e.terms))
                            for i in 1:nf)
                b, q, lower, upper
            end
            selected_local_scan = largest_feasible_scalar(
                max(1e-9, local_limit / 2000.0), 0.999local_limit,
                rho -> begin
                    radii = fill(rho, nf)
                    b, q, lower, upper = local_contribution(radii)
                    isfinite(q) && all(isfinite, b) && all(b .<= radii) && q < 1.0 ||
                        return nothing
                    (radii=radii, b=b, q=q, lower=lower, upper=upper)
                end,
                value -> value !== nothing,
                score=value -> 1.0 - value.q)
            selected_local = selected_local_scan === nothing ? nothing :
                selected_local_scan.value
            selected_local === nothing ? nothing :
                (kind=:candidate_local, radii=selected_local.radii,
                 b=selected_local.b, q=selected_local.q,
                 max_condition_margin=selected_local_scan.max_score,
                 lower=selected_local.lower, upper=selected_local.upper,
                 offset=offset)
        end
    catch
        nothing
    end
    uniform_region = if selected === nothing
        nothing
    else
        radii = fill(selected.rho, nf)
        b, q, lower, upper = component_contribution(radii)
        (kind=:uniform, radii=radii, b=b, q=q, lower=lower, upper=upper,
         max_condition_margin=selected_max_condition_margin)
    end
    component_region = if componentwise === nothing
        nothing
    else
        (kind=:componentwise, radii=componentwise.radii, b=componentwise.b,
         q=componentwise.q, lower=componentwise.lower, upper=componentwise.upper,
         max_condition_margin=componentwise_max_condition_margin)
    end
    candidate_ratio(region) = maximum(candidate_delta ./ region.radii)
    candidate_regions = [(region, candidate_ratio(region)) for region in
                         (uniform_region, component_region) if region !== nothing]
    inside_regions = [(region, ratio) for (region, ratio) in candidate_regions
                      if ratio <= 1.0 + spec.residual_atol]
    # Prefer the anisotropic componentwise polydisc whenever it contains the
    # candidate.  It retains phase-specific coupling information and is a
    # strictly more informative report for unbalanced networks; use volume as
    # the deterministic tie-breaker among regions of the same geometry.
    chosen = if !isempty(inside_regions)
        first(sort(inside_regions; by=x ->
            (x[1].kind == :componentwise ? 0 : 1, -sum(log, x[1].radii))))
    elseif !isempty(candidate_regions)
        first(sort(candidate_regions; by=x ->
            (x[1].kind == :componentwise ? 0 : 1, -sum(log, x[1].radii))))
    else
        nothing
    end
    chosen_region = chosen === nothing ? nothing : chosen[1]
    if chosen !== nothing
        current_jacobian(sl, dv) = begin
            h = 1e-6 * max(abs(dv), 1.0)
            current(z) = conj(_subload_S_nz(sl, abs(z))) / conj(z)
            plus_r = (current(dv + h) - current(dv - h)) / (2h)
            plus_i = (current(dv + im * h) - current(dv - im * h)) / (2h)
            [real(plus_r) real(plus_i); imag(plus_r) imag(plus_i)]
        end
        validation_rows = Dict{String,Any}()
        validation_ok = true
        for (j, edge) in enumerate(edges)
            lower = chosen_region.lower[j]
            upper = chosen_region.upper[j]
            radii = (lower, (lower + upper) / 2.0, upper)
            ratios = Float64[]; numerical = Float64[]; bounds = Float64[]
            for radius in radii
                dv = radius * edge.dv0 / abs(edge.dv0)
                numerical_norm = opnorm(current_jacobian(edge.sl, dv), 2)
                bound = sum(term_bounds(term, radius, radius)[2]
                            for term in edge.terms)
                ratio = bound == 0.0 ? 0.0 : numerical_norm / bound
                push!(ratios, ratio); push!(numerical, numerical_norm); push!(bounds, bound)
            end
            row_status = all(isfinite, ratios) && maximum(ratios) <= 1.0 + 1e-3
            validation_ok &= row_status
            validation_rows[edge.key] = Dict{String,Any}(
                "radii" => collect(radii), "numerical_operator_norm" => numerical,
                "analytic_bound" => bounds, "ratio" => ratios,
                "status" => row_status ? :pass : :inconclusive)
        end
        base["law_bound_validation"] = Dict{String,Any}(
            "status" => validation_ok ? :pass : :inconclusive,
            "method" => "central_finite_difference_real_current_jacobian",
            "samples_per_connection" => 3,
            "connections" => validation_rows,
            "interpretation" => "diagnostic validation of the analytic Lipschitz bound; not the sufficient-condition proof")
    else
        base["law_bound_validation"] = Dict{String,Any}(
            "status" => :not_applicable,
            "message" => "no invariant contraction region was available for law-bound validation")
    end
    # The map is independently iterated from the energized no-load germ. This
    # is evidence recorded alongside, but not substituted for, the theorem.
    map_complex(vfree) = begin
        vfull = copy(wfull); vfull[meta.free] = vfree
        w + zsolve(lin.i_comp(vfull)[meta.free])
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
    distance = maximum(candidate_delta)
    base["candidate_distance"] = distance
    base["uniform_region"] = uniform_region === nothing ? nothing : Dict(
        "radii" => uniform_region.radii, "contraction_factor" => uniform_region.q,
        "max_condition_margin" => uniform_region.max_condition_margin,
        "invariance_bound" => maximum(uniform_region.b),
        "invariance_form" => "componentwise b_i <= radius_i in the weighted infinity norm",
        "candidate_ratio" => candidate_ratio(uniform_region))
    base["componentwise_region"] = component_region === nothing ? nothing : Dict(
        "radii" => component_region.radii, "contraction_factor" => component_region.q,
        "max_condition_margin" => component_region.max_condition_margin,
        "invariance_bound" => maximum(component_region.b),
        "invariance_form" => "componentwise b_i <= radius_i in the weighted infinity norm",
        "candidate_ratio" => candidate_ratio(component_region))
    base["local_candidate_region"] = local_region === nothing ? Dict(
        "status" => :inconclusive,
        "message" => "no candidate-centered invariant contraction region was found") : Dict(
        "status" => :pass,
        "radii" => local_region.radii,
        "contraction_factor" => local_region.q,
        "invariance_bound" => maximum(local_region.b),
        "invariance_form" => "candidate-centered offset + q * radius <= radius in the componentwise norm",
        "condition_margin" => 1.0 - local_region.q,
        "max_condition_margin" => local_region.max_condition_margin,
        "candidate_residual_offset" => local_region.offset,
        "interpretation" => "local uniqueness evidence around the candidate; does not establish no-load reachability")
    margins = Float64[
        uniform_region === nothing ? NaN : uniform_region.max_condition_margin,
        component_region === nothing ? NaN : component_region.max_condition_margin,
        euclidean_max_condition_margin,
    ]
    finite_margins = filter(isfinite, margins)
    base["max_condition_margin"] = isempty(finite_margins) ? NaN : maximum(finite_margins)
    if chosen === nothing
        base["message"] = "no invariant contraction polydisc was found; sufficient condition is inconclusive"
        return base
    end
    region, region_ratio = chosen
    base["selected_region"] = String(region.kind)
    base["region_radii"] = copy(region.radii)
    base["radius"] = maximum(region.radii)
    base["invariance_bound"] = maximum(region.b)
    base["invariance_form"] = "componentwise b_i <= radius_i in the weighted infinity norm"
    base["contraction_factor"] = region.q
    base["minimum_connection_voltage"] = minimum(region.lower)
    base["condition_margin"] = 1.0 - region.q
    inside = region_ratio <= 1.0 + spec.residual_atol
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

function _operability_empty_result(coords, unsupported, message; net=nothing)
    checks = Dict("scope" => OperabilityCheck(:not_applicable, unsupported, nothing, message))
    provenance = deepcopy(coords.provenance)
    if net !== nothing
        generator_count = length(get(net, "generator", Dict()))
        ibr_count = length(get(net, "ibr", Dict()))
        provenance["operability"] = Dict(
            "scope" => "static_ybus_linearized",
            "status" => :not_applicable,
            "closure" => "frozen_dispatch",
            "control_closure" => generator_count == 0 && ibr_count == 0 ?
                "frozen_dispatch_native_static" : "outside_native_static_seam",
            "unsupported_reasons" => copy(unsupported),
            "model_inventory" => _operability_scope_inventory(net),
            "topology" => _operability_scope_source_inventory(net),
            "coordinate_policy" => BMOPFTools.opf_scaling_policy_data(coords.policy),
        )
    end
    OperabilityResult(:not_applicable, NaN, NaN, _Node[], Float64[], zeros(0, 0),
        Float64[], Inf, Dict{_Node,ComplexF64}(), Dict{String,Any}(),
        Dict{String,Any}(), Dict{String,Any}(), Dict{String,Any}(), checks,
        provenance, unsupported)
end

function _operability_overall(checks::Dict{String,OperabilityCheck})
    # Candidate-centered local uniqueness is complementary evidence. It must
    # not downgrade the primary snapshot verdict when the no-load-connected
    # certificate or another requested operational claim has its own outcome.
    primary = (c for (key, c) in checks if key != "fixed_point_local_region")
    any(c.status === :fail for c in primary) && return :fail
    primary = (c for (key, c) in checks if key != "fixed_point_local_region")
    any(c.status === :inconclusive for c in primary) && return :inconclusive
    claim_keys = ("terminal_voltage_bounds", "sequence_unbalance", "helm_reachability",
                  "fixed_point_certificate", "fixed_point_euclidean_region")
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
    scope_audit = operability_scope_audit(net)
    unsupported = copy(scope_audit["unsupported_reasons"])
    !isempty(unsupported) && return _operability_empty_result(coords, unsupported,
        "candidate contains physics outside the first operability residual scope";
        net=net)
    source_buses = _operability_source_buses(net)
    isempty(source_buses) && throw(ArgumentError("operability checking requires a voltage source"))

    lin = try
        BMOPFTools.ybus_linearized(net; fold=:constant_z)
    catch err
        err isa ArgumentError ? throw(OperabilityModelError(sprint(showerror, err))) : rethrow()
    end
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

    Jphys = if isempty(x)
        zeros(0, 0)
    elseif spec.jacobian_storage === :sparse
        _operability_sparse_finite_difference_jacobian(
            y -> _operability_residual(lin, meta, y), x; step=spec.jacobian_step)
    else
        finite_difference_jacobian(
            y -> _operability_residual(lin, meta, y), x; step=spec.jacobian_step)
    end
    J = _operability_scaled_jacobian(Jphys, state_scale, residual_scale)
    linear_solver = if !isempty(J) &&
        (spec.compute_sensitivity || spec.jacobian_spectrum === :extremes)
        try
            lu(J)
        catch
            nothing
        end
    else
        nothing
    end
    jacobian_factorization = spec.jacobian_spectrum === :full && !isempty(J) ? svd(J) : nothing
    effective_spectrum_mode = String(spec.jacobian_spectrum)
    extreme_estimator_budget = max(80, min(size(J, 2), 256))
    singular_values = if spec.jacobian_spectrum === :full
        jacobian_factorization === nothing ? Float64[] : Float64.(jacobian_factorization.S)
    else
        extremes = _operability_extreme_singular_values(J, linear_solver;
            max_iterations=extreme_estimator_budget)
        if extremes === nothing && !issparse(J) && !isempty(J)
            # A dense report can preserve the regularity/sensitivity contract
            # by falling back to its reference SVD when the reduced estimator
            # does not converge. Sparse mode remains explicitly inconclusive
            # rather than silently materializing a dense matrix.
            jacobian_factorization = svd(J)
            effective_spectrum_mode = "full_fallback"
            Float64.(jacobian_factorization.S)
        else
            extremes === nothing ? Float64[] : extremes
        end
    end
    condition_number = isempty(singular_values) || last(singular_values) == 0 ?
        Inf : first(singular_values) / last(singular_values)
    rank_tol = if isempty(singular_values) || !(first(singular_values) > 0.0)
        Inf
    else
        spec.jacobian_rank_rtol * first(singular_values)
    end
    if isempty(x)
        checks["jacobian_regular"] = OperabilityCheck(:not_applicable, nothing, nothing,
            "no non-source voltage state remains after source/reference elimination")
    elseif isempty(singular_values)
        checks["jacobian_regular"] = OperabilityCheck(:inconclusive, nothing, rank_tol,
            "reduced Jacobian extreme-singular-value estimator did not converge")
    else
        checks["jacobian_regular"] = OperabilityCheck(last(singular_values) > rank_tol ?
            :pass : :inconclusive, last(singular_values), rank_tol,
            "scaled rectangular current-balance Jacobian smallest singular value")
    end
    jacobian_unavailable_message = isempty(singular_values) &&
        spec.jacobian_spectrum === :extremes ?
        "extreme-singular-value estimator did not converge" :
        "Jacobian is too close to singular"

    sensitivities = Dict{String,Any}()
    load_scale_dx = nothing
    if spec.compute_sensitivity && !isempty(x) && !isempty(singular_values) &&
        last(singular_values) > rank_tol
        h = spec.sensitivity_step
        plus = _operability_scale_network(net, 1.0 + h)
        minus = _operability_scale_network(net, 1.0 - h)
        lp = BMOPFTools.ybus_linearized(plus; fold=:constant_z)
        lm = BMOPFTools.ybus_linearized(minus; fold=:constant_z)
        fp = _operability_residual(lp, meta, x); fm = _operability_residual(lm, meta, x)
        dF = (fp - fm) / (2h)
        dF_scaled = dF ./ residual_scale
        dx_scaled = linear_solver === nothing ? -(J \ dF_scaled) :
            -(linear_solver \ dF_scaled)
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
            isempty(x) ? "no state to differentiate" : jacobian_unavailable_message)
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
    if spec.compute_sensitivity && !isempty(x) && !isempty(singular_values) &&
        last(singular_values) > rank_tol
        for (id_raw, load) in sort!(collect(get(net, "load", Dict())); by=x -> String(x[1]))
            id = String(id_raw)
            nconn = _operability_load_connection_count(load)
            for field in ("p_nom", "q_nom")
                family = field == "p_nom" ? "P" : "Q"
                family_records = get!(directions, family, Dict{String,Any}())
                for k in 1:nconn
                    family_records["$id/$k"] = _operability_directional_sensitivity(
                        net, lin, meta, x, J, linear_solver,
                        residual_scale, state_scale, spec, id, k, field)
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

    nfree = length(meta.free)
    nreal = length(x)
    jacobian_nonzero_count = issparse(J) ? nnz(J) : count(v -> v != 0.0, J)
    nconnections = sum((_operability_load_connection_count(load)
                        for load in values(get(net, "load", Dict()))); init=0)
    ndirections = spec.compute_sensitivity && !isempty(x) &&
        last(singular_values) > rank_tol ? 1 + 2nconnections : 0
    jacobian_pattern = if issparse(J) && spec.record_jacobian_pattern
        rows, columns, _ = findnz(J)
        Dict("rows" => rows, "columns" => columns,
             "interpretation" => "finite-difference nonzero pattern at this snapshot")
    else
        nothing
    end
    branch_complexity = Dict{String,Any}(
        "free_node_count" => nfree,
        "real_state_dimension" => nreal,
        "load_connection_count" => nconnections,
        "sensitivity_direction_count" => ndirections,
        "jacobian_spectrum_mode" => String(spec.jacobian_spectrum),
        "jacobian_spectrum_effective_mode" => effective_spectrum_mode,
        "jacobian_extreme_estimator_iteration_budget" =>
            spec.jacobian_spectrum === :extremes ? extreme_estimator_budget : 0,
        "jacobian_storage_mode" => String(spec.jacobian_storage),
        "jacobian_nonzero_count" => jacobian_nonzero_count,
        "jacobian_pattern_recorded" => spec.record_jacobian_pattern && issparse(J),
        "jacobian_pattern" => jacobian_pattern,
        "linear_solver_factorization_reused" => linear_solver !== nothing,
        "jacobian_storage_bytes_dense" => sizeof(Float64) * nreal * nreal,
        "jacobian_storage_bytes_estimate" => Base.summarysize(J),
        "zbus_storage_bytes_dense" => spec.compute_fixed_point_certificate ?
            sizeof(ComplexF64) * nfree * nfree : 0,
        "zbus_storage_mode" => !spec.compute_fixed_point_certificate ? "not_requested" :
            (issparse(lin.Y) ? "implicit_sparse_factorization" :
             "implicit_dense_factorization"),
        "fixed_point_scan_points_per_geometry" => spec.compute_fixed_point_certificate ? 33 : 0,
        "fixed_point_refinement_iterations" => spec.compute_fixed_point_certificate ? 12 : 0,
        "fixed_point_scan_strategy" => spec.compute_fixed_point_certificate ?
            "coarse_grid_bisection" : "not_requested",
        "fixed_point_geometry_count" => spec.compute_fixed_point_certificate ? 4 : 0,
        "interpretation" => "diagnostic size indicators; dense Jacobian/SVD and Z-bus certificate storage can dominate large-network runs")
    critical_mode = jacobian_factorization !== nothing ?
        _operability_critical_mode(jacobian_factorization, meta.state_nodes) :
        Dict{String,Any}(
            "status" => :not_applicable,
            "message" => "critical singular vectors are omitted in jacobian_spectrum=:extremes mode",
            "spectrum_mode" => String(spec.jacobian_spectrum))
    branch_evidence = Dict{String,Any}(
        "critical_mode" => critical_mode,
        "dP_dV" => _operability_path_dpdv_evidence(sensitivities),
        "sequence_sensitivity" => _operability_sequence_sensitivity_evidence(sensitivities),
        "complexity" => branch_complexity)
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
        local_region = get(certificate, "local_candidate_region", Dict{String,Any}())
        local_status = Symbol(get(local_region, "status",
            haskey(certificate, "local_candidate_region") ? :inconclusive : :not_applicable))
        checks["fixed_point_local_region"] = OperabilityCheck(
            local_status, get(local_region, "contraction_factor", nothing), 1.0,
            local_status === :pass ?
                "candidate-centered contraction gives local uniqueness evidence; reachability remains separate" :
                "no candidate-centered contraction region was certified")
        euclidean = get(certificate, "euclidean_region", Dict{String,Any}())
        euclidean_status = Symbol(get(euclidean, "status", :not_applicable))
        checks["fixed_point_euclidean_region"] = OperabilityCheck(
            euclidean_status, get(euclidean, "contraction_factor", nothing), 1.0,
            euclidean_status === :pass ?
                "Euclidean-ball contraction certifies the candidate on the no-load-connected branch" :
                "Euclidean-ball sufficient condition is inconclusive")
    else
        branch_evidence["fixed_point_certificate"] = Dict{String,Any}(
            "status" => :not_applicable,
            "message" => "fixed-point certificate was not requested")
        checks["fixed_point_certificate"] = OperabilityCheck(:not_applicable,
            nothing, nothing, "fixed-point certificate was not requested")
    end

    provenance = deepcopy(coords.provenance)
    provenance["operability"] = Dict("scope" => "static_ybus_linearized",
        "status" => scope_audit["status"],
        "closure" => String(spec.closure),
        "control_closure" => "frozen_dispatch_native_static",
        "jacobian_spectrum_mode" => String(spec.jacobian_spectrum),
        "jacobian_storage_mode" => String(spec.jacobian_storage),
        "source_buses" => sort!(collect(source_buses)), "state_nodes" => meta.state_nodes,
        "coordinate_policy" => BMOPFTools.opf_scaling_policy_data(coords.policy),
        "model_inventory" => _operability_scope_inventory(net),
        "topology" => _operability_scope_source_inventory(net))
    OperabilityResult(_operability_overall(checks), endpoint_residual,
        endpoint_normalized, meta.state_nodes, x, J, singular_values, condition_number,
        node_voltages, load_records, sequence_records, sensitivities, branch_evidence, checks,
        provenance, unsupported)
end

"""
    operability_stress_network(net, lambda, direction)

Return a copied native-load network with a named finite stress direction
applied. Every nameplate P/Q term is zero at `lambda=0`; at other values it is
multiplied by the direction's P/Q scale, optional load weight, and optional
connection weight. This builder does not solve the network or claim that the
resulting endpoint is the audited OPF point.
"""
function operability_stress_network(net::Dict{String,Any}, lambda::Real,
                                    direction::OperabilityStressDirection)
    _operability_stress_network(net, lambda, direction)
end

"""
    operability_stress_rows(net, solution; spec, directions, lambdas, solve)

Run a bounded, deterministic campaign of named native-static stress snapshots.
`solve` is an explicit caller-supplied callback accepting either
`(network, previous_solution, direction, lambda)`,
`(network, previous_solution, direction)`, `(network, previous_solution)`, or
`(network)`, and must return an SI-valued solution dictionary. Set
`solve_arity` to `1`, `2`, `3`, or `4` when the callback is variadic or its
intended arity cannot be inferred safely. Each row
retains the stress name, lambda, endpoint residual, certificate status, region
margin, and any solver/checking error. These rows are finite study evidence;
they do not establish a global operating envelope or a continuation branch.
"""
function operability_stress_rows(net::Dict{String,Any}, solution::AbstractDict;
                                 spec::OperabilitySpec,
                                 directions=(OperabilityStressDirection(),),
                                 lambdas=(0.0, 0.5, 1.0), solve=nothing,
                                 context=nothing, solve_arity=nothing)
    solve === nothing && throw(ArgumentError(
        "operability_stress_rows requires an explicit solve callback"))
    direction_list = collect(directions)
    isempty(direction_list) && throw(ArgumentError("directions cannot be empty"))
    all(direction -> direction isa OperabilityStressDirection, direction_list) ||
        throw(ArgumentError("directions must contain OperabilityStressDirection values"))
    names = Symbol[direction.name for direction in direction_list]
    length(unique(names)) == length(names) ||
        throw(ArgumentError("stress-direction names must be unique"))
    lambda_list = sort!(Float64.(collect(lambdas)))
    isempty(lambda_list) && throw(ArgumentError("lambdas cannot be empty"))
    all(λ -> isfinite(λ) && λ >= 0.0, lambda_list) ||
        throw(ArgumentError("lambdas must be finite and >= 0"))
    if solve_arity !== nothing
        solve_arity isa Integer && solve_arity in 1:4 || throw(ArgumentError(
            "solve_arity must be nothing or an integer in 1:4"))
        solve_arity = Int(solve_arity)
    end
    rows = NamedTuple[]
    for direction in direction_list
        previous = solution
        for λ in lambda_list
            stressed = operability_stress_network(net, λ, direction)
            try
                solved = if solve_arity == 4 ||
                    (solve_arity === nothing && applicable(solve, stressed, previous, direction, λ))
                    solve(stressed, previous, direction, λ)
                elseif solve_arity == 3 ||
                    (solve_arity === nothing && applicable(solve, stressed, previous, direction))
                    solve(stressed, previous, direction)
                elseif solve_arity == 2 ||
                    (solve_arity === nothing && applicable(solve, stressed, previous))
                    solve(stressed, previous)
                elseif solve_arity == 1 ||
                    (solve_arity === nothing && applicable(solve, stressed))
                    solve(stressed)
                else
                    throw(ArgumentError("solve callback has no supported arity"))
                end
                solved isa AbstractDict || throw(ArgumentError(
                    "solve callback must return an AbstractDict solution"))
                report = check_opf_operability(stressed, solved; spec, context)
                certificate = get(report.branch_evidence, "fixed_point_certificate", Dict())
                previous = solved
                push!(rows, (
                    direction=direction.name, lambda=λ, status=report.status,
                    p_scale=direction.p_scale, q_scale=direction.q_scale,
                    load_weights=deepcopy(direction.load_weights),
                    connection_weights=deepcopy(direction.connection_weights),
                    endpoint_status=get(report.checks, "endpoint", OperabilityCheck(
                        :not_applicable, nothing, nothing, "endpoint unavailable")).status,
                    certificate_status=get(report.checks, "fixed_point_certificate",
                        OperabilityCheck(:not_applicable, nothing, nothing,
                            "certificate unavailable")).status,
                    endpoint_residual_normalized=report.endpoint_residual_normalized,
                    selected_region=get(certificate, "selected_region", missing),
                    contraction_factor=Float64(get(certificate, "contraction_factor", NaN)),
                    max_condition_margin=Float64(get(certificate,
                        "max_condition_margin", NaN)),
                    minimum_connection_voltage=Float64(get(certificate,
                        "minimum_connection_voltage", NaN)),
                    candidate_inside_region=get(certificate,
                        "candidate_inside_region", missing),
                    message="stress snapshot checked", error=nothing))
            catch err
                if err isa OperabilityModelError
                    # Expected snapshot data/model failures remain row-level
                    # inconclusive evidence in a finite campaign.
                    nothing
                elseif err isa InterruptException || err isa MethodError ||
                    err isa UndefVarError || err isa ArgumentError ||
                    err isa TypeError || err isa DimensionMismatch ||
                    err isa BoundsError || err isa LoadError
                    rethrow()
                end
                push!(rows, (
                    direction=direction.name, lambda=λ, status=:inconclusive,
                    p_scale=direction.p_scale, q_scale=direction.q_scale,
                    load_weights=deepcopy(direction.load_weights),
                    connection_weights=deepcopy(direction.connection_weights),
                    endpoint_status=:inconclusive, certificate_status=:inconclusive,
                    endpoint_residual_normalized=NaN, selected_region=missing,
                    contraction_factor=NaN, max_condition_margin=NaN,
                    minimum_connection_voltage=NaN,
                    candidate_inside_region=missing, message="stress snapshot failed",
                    error=sprint(showerror, err)))
            end
        end
    end
    rows
end

function _operability_stress_status(rows)
    statuses = Symbol[get(row, :status, :inconclusive) for row in rows]
    any(==(:fail), statuses) && return :fail
    any(==(:inconclusive), statuses) && return :inconclusive
    any(==(:pass), statuses) && return :pass
    :not_applicable
end

"""
    operability_stress_summary(rows; reference_lambda=1.0)

Aggregate one finite stress campaign into deterministic direction-level rows.
Status precedence is `:fail` > `:inconclusive` > `:pass` >
`:not_applicable`. The first failing or inconclusive row in ascending lambda
order is a path-specific observed boundary; absence of one is reported as
`:not_observed`, never as unlimited margin.
"""
function operability_stress_summary(rows; reference_lambda::Real=1.0)
    reference = Float64(reference_lambda)
    isfinite(reference) && reference > 0.0 || throw(ArgumentError(
        "reference_lambda must be finite and > 0"))
    row_list = collect(rows)
    isempty(row_list) && return NamedTuple[]
    all(hasproperty(row, :direction) && hasproperty(row, :lambda) for row in row_list) ||
        throw(ArgumentError("stress rows need direction and lambda fields"))
    grouped = Dict{Symbol,Vector{Any}}()
    for row in row_list
        push!(get!(grouped, Symbol(row.direction), Any[]), row)
    end
    summaries = NamedTuple[]
    for direction in sort!(collect(keys(grouped)); by=string)
        direction_rows = sort(grouped[direction]; by=row -> Float64(row.lambda))
        status = _operability_stress_status(direction_rows)
        boundary = nothing
        for row in direction_rows
            get(row, :status, :inconclusive) in (:fail, :inconclusive) || continue
            boundary = row
            break
        end
        finite_margins = Float64[]
        for row in direction_rows
            margin = if hasproperty(row, :max_condition_margin)
                Float64(row.max_condition_margin)
            elseif hasproperty(row, :contraction_factor)
                1.0 - Float64(row.contraction_factor)
            else
                NaN
            end
            isfinite(margin) && push!(finite_margins, margin)
        end
        finite_voltages = Float64[Float64(row.minimum_connection_voltage)
                                  for row in direction_rows
                                  if hasproperty(row, :minimum_connection_voltage) &&
                                     isfinite(Float64(row.minimum_connection_voltage))]
        residuals = Float64[Float64(row.endpoint_residual_normalized)
                            for row in direction_rows
                            if hasproperty(row, :endpoint_residual_normalized) &&
                               isfinite(Float64(row.endpoint_residual_normalized))]
        errors = count(row -> hasproperty(row, :error) && row.error !== nothing,
                       direction_rows)
        push!(summaries, (
            direction=direction, status=status, row_count=length(direction_rows),
            pass_count=count(row -> get(row, :status, :inconclusive) === :pass,
                             direction_rows),
            fail_count=count(row -> get(row, :status, :inconclusive) === :fail,
                             direction_rows),
            inconclusive_count=count(row ->
                get(row, :status, :inconclusive) === :inconclusive, direction_rows),
            error_count=errors, boundary_status=boundary === nothing ? :not_observed :
                get(boundary, :status, :inconclusive),
            boundary_lambda=boundary === nothing ? NaN : Float64(boundary.lambda),
            reference_lambda=reference,
            parameter_margin=boundary === nothing ? NaN :
                Float64(boundary.lambda) - reference,
            relative_margin=boundary === nothing ? NaN :
                (Float64(boundary.lambda) - reference) / reference,
            minimum_condition_margin=isempty(finite_margins) ? NaN : minimum(finite_margins),
            minimum_connection_voltage=isempty(finite_voltages) ? NaN : minimum(finite_voltages),
            maximum_endpoint_residual_normalized=isempty(residuals) ? NaN : maximum(residuals),
            p_scale=Float64(get(first(direction_rows), :p_scale, NaN)),
            q_scale=Float64(get(first(direction_rows), :q_scale, NaN)),
            message=boundary === nothing ?
                "no non-pass boundary observed on the finite stress rows" :
                "first non-pass row in the finite stress campaign"))
    end
    summaries
end

"""
    operability_stress_ensemble_rows(campaigns; reference_lambda=1.0)

Summarize multiple finite stress campaigns, supplied as a dictionary or
iterable of `(model_id, rows)` pairs. The returned rows retain the model label
and direction-level path margins; differences are campaign comparisons, not a
robustness guarantee over an uncertainty set.
"""
function operability_stress_ensemble_rows(campaigns; reference_lambda::Real=1.0)
    entries = campaigns isa AbstractDict ? collect(campaigns) : collect(campaigns)
    rows = NamedTuple[]
    for (model_id, campaign_rows) in entries
        summaries = operability_stress_summary(campaign_rows;
            reference_lambda=reference_lambda)
        for summary in summaries
            push!(rows, merge((model=String(model_id),), summary))
        end
    end
    sort!(rows; by=row -> (row.model, String(row.direction)))
end

"""
    operability_snapshot_row(result; snapshot_id=nothing)

Return one compact, table-ready row for a single [`OperabilityResult`](@ref).
The row preserves the snapshot label, endpoint/regularity evidence, voltage
and VUF extrema, branch-indicator counts, no-load and candidate-local
certificate/HELM statuses, scaling/storage metadata, and the number of
unsupported-scope reasons. It also retains scope status, closure/control
closure, and source-topology readiness so pooled rows cannot hide an
out-of-scope snapshot. It is a reporting projection of one snapshot, not a
contingency or operating-envelope assessment.
"""
function operability_snapshot_row(result::OperabilityResult; snapshot_id=nothing)
    magnitudes = Float64[Float64(get(record, "magnitude", NaN))
                         for record in values(result.load_connections)
                         if isfinite(Float64(get(record, "magnitude", NaN)))]
    vufs = Float64[Float64(get(record, "vuf", NaN))
                   for record in values(result.sequences)
                   if isfinite(Float64(get(record, "vuf", NaN)))]
    dpdv = get(get(result.branch_evidence, "dP_dV", Dict{String,Any}()),
               "connections", Dict{String,Any}())
    classifications = String[get(record, "classification", "not_available")
                             for record in values(dpdv)]
    certificate = get(result.branch_evidence, "fixed_point_certificate", Dict{String,Any}())
    local_certificate = get(certificate, "local_candidate_region", Dict{String,Any}())
    euclidean_certificate = get(certificate, "euclidean_region", Dict{String,Any}())
    reachability = get(result.branch_evidence, "reachability", Dict{String,Any}())
    complexity = get(result.branch_evidence, "complexity", Dict{String,Any}())
    operability = get(result.provenance, "operability", Dict{String,Any}())
    topology = get(operability, "topology", Dict{String,Any}())
    (
        snapshot_id=snapshot_id,
        status=result.status,
        endpoint_residual=result.endpoint_residual,
        endpoint_residual_normalized=result.endpoint_residual_normalized,
        smallest_singular_value=isempty(result.singular_values) ? NaN : last(result.singular_values),
        condition_number=result.condition_number,
        minimum_terminal_voltage=isempty(magnitudes) ? NaN : minimum(magnitudes),
        maximum_terminal_voltage=isempty(magnitudes) ? NaN : maximum(magnitudes),
        maximum_vuf=isempty(vufs) ? missing : maximum(vufs),
        maximum_vuf_status=isempty(vufs) ? :not_applicable : :available,
        high_side_indicator_count=count(==("negative_high_side_indicator"), classifications),
        near_nose_indicator_count=count(==("near_nose_indicator"), classifications),
        low_side_indicator_count=count(==("positive_low_side_indicator"), classifications),
        fixed_point_certificate_status=Symbol(get(certificate, "status", :not_applicable)),
        fixed_point_condition_margin=Float64(get(certificate, "condition_margin", NaN)),
        fixed_point_max_condition_margin=Float64(get(certificate,
            "max_condition_margin", NaN)),
        fixed_point_selected_region=get(certificate, "selected_region", missing),
        fixed_point_local_region_status=Symbol(get(local_certificate, "status", :not_applicable)),
        fixed_point_local_condition_margin=Float64(get(local_certificate, "condition_margin", NaN)),
        fixed_point_euclidean_region_status=Symbol(get(euclidean_certificate, "status", :not_applicable)),
        fixed_point_euclidean_condition_margin=Float64(get(euclidean_certificate, "condition_margin", NaN)),
        helm_reachability_status=Symbol(get(reachability, "status", :not_applicable)),
        jacobian_spectrum_mode=String(get(complexity, "jacobian_spectrum_mode", "not_available")),
        jacobian_storage_mode=String(get(complexity, "jacobian_storage_mode", "not_available")),
        jacobian_nonzero_count=get(complexity, "jacobian_nonzero_count", missing),
        jacobian_storage_bytes_estimate=get(complexity, "jacobian_storage_bytes_estimate", missing),
        zbus_storage_mode=String(get(complexity, "zbus_storage_mode", "not_applicable")),
        unsupported_count=length(result.unsupported),
        scope_status=Symbol(get(operability, "status",
            result.status === :not_applicable ? :not_applicable : :supported)),
        equilibrium_scope=String(get(operability, "scope", "static_ybus_linearized")),
        closure=Symbol(get(operability, "closure", "frozen_dispatch")),
        control_closure=String(get(operability, "control_closure",
            "frozen_dispatch_native_static")),
        topology_has_voltage_source=Bool(get(topology, "has_voltage_source", false)),
        topology_missing_source_buses=get(topology, "missing_source_buses", String[]),
        scope="single_snapshot_static_ybus",
    )
end

"""
    operability_snapshot_rows(results)

Project a dictionary or iterable of `(snapshot_id, result)` pairs into a
deterministically ordered vector of table-ready snapshot rows. Ordering is
lexicographic on `String(snapshot_id)`, so callers using numeric labels should
zero-pad them when numeric order matters. This adapter is intentionally limited
to already-solved single snapshots; it does not infer contingency, temporal, or
operating-envelope semantics.
"""
function operability_snapshot_rows(results)
    entries = results isa AbstractDict ? collect(results) : collect(results)
    rows = NamedTuple[]
    for entry in entries
        entry isa Pair || throw(ArgumentError(
            "operability_snapshot_rows expects (snapshot_id, OperabilityResult) pairs"))
        snapshot_id, result = entry
        result isa OperabilityResult || throw(ArgumentError(
            "snapshot $(repr(snapshot_id)) is not an OperabilityResult"))
        push!(rows, operability_snapshot_row(result; snapshot_id=String(snapshot_id)))
    end
    sort!(rows; by=row -> String(row.snapshot_id))
end
