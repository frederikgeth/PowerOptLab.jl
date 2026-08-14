# Generic fixed-point and local loop-gain diagnostics.

"""
    FixedPointIterationResult

Result of [`fixed_point_oracle`](@ref). `trajectory` stores column vectors,
including the initial state, when trajectory storage is enabled.

This is an empirical equilibrium/convergence diagnostic. It is not a proof of
global existence, uniqueness, or closed-loop stability.
"""
struct FixedPointIterationResult
    converged::Bool
    cycled::Bool
    iterations::Int
    cycle_period::Int
    residual_norm::Float64
    final_state::Vector{Float64}
    trajectory::Matrix{Float64}
end

"""Result of a finite-difference local fixed-point loop-gain screen."""
struct FixedPointGainScreen
    jacobian::Matrix{Float64}
    eigenvalues::Vector{ComplexF64}
    spectral_radius::Float64
    induced_inf_norm::Float64
    local_contractive::Bool
    margin::Float64
    finite_difference_step::Float64
end

"""
    InverterControlScalingAudit

SI-valued audit of the auxiliary starts, scales, and smoothing widths used by
the smooth sequence-controller formulation. The report is intentionally pure
numeric: it can be recorded with a study manifest without constructing a
JuMP model.
"""
struct InverterControlScalingAudit
    voltage_start_V::Float64
    voltage_scale_V::Float64
    sequence_voltage_start_V::Float64
    power_start_VA::Float64
    power_scale_VA::Float64
    capacity_aux_start_VA::Float64
    current_start_A::Float64
    current_scale_A::Float64
    voltage_floor_V::Float64
    current_epsilon_A::Float64
    power_epsilon_VA::Float64
    priority_headroom_fraction::Float64
end

function _real_state(x, name::String)
    values = Float64.(collect(x))
    all(isfinite, values) || throw(ArgumentError("$name must be finite"))
    values
end

function _control_voltage_anchor(policy::AbstractPositiveSequencePolicy,
                                fallback::Real)
    anchor = Float64(fallback)
    isfinite(anchor) && anchor > 0 || throw(ArgumentError(
        "voltage_fallback must be finite and > 0"))
    points = Float64[]
    for law in _policy_laws(policy)
        law === nothing || append!(points, law.breakpoints)
    end
    isempty(points) ? anchor : (minimum(points) + maximum(points))/2
end

"""
    inverter_control_scaling_audit(controller, ratings; voltage_fallback=230)
        -> InverterControlScalingAudit

Return the physical scales used by the smooth controller auxiliaries. Curve
breakpoints determine the voltage anchor when present; otherwise
`voltage_fallback` is used. The capacity-auxiliary start is the small residual
capability at the declared P/Q-priority headroom, rather than the full
apparent-power rating. This makes heterogeneous fleets auditable and exposes
poorly conditioned starts before a network solve is attempted.
"""
function inverter_control_scaling_audit(
        controller::SequenceController,
        ratings::InverterControlRatings;
        voltage_fallback::Real=230.0)
    voltage_start = _control_voltage_anchor(
        controller.positive, voltage_fallback)
    headroom = controller.limiter.priority_headroom_fraction
    capacity_start = ratings.s_max * sqrt(2headroom - headroom^2)
    InverterControlScalingAudit(
        voltage_start,
        voltage_start,
        0.05voltage_start,
        ratings.s_max,
        ratings.s_max,
        capacity_start,
        max(0.5ratings.i_max,
            _current_epsilon(controller.limiter, ratings)),
        ratings.i_max,
        controller.power_voltage_floor,
        _current_epsilon(controller.limiter, ratings),
        _power_epsilon(controller.limiter, ratings),
        headroom)
end

"""
    fixed_point_oracle(map, initial; kwargs...) -> FixedPointIterationResult

Iterate `x ↦ map(x)` with optional under-relaxation. The callable must accept a
finite real vector and return a vector of the same length. Convergence uses
the update residual with an absolute-plus-relative tolerance. A repeated state
within `cycle_tolerance` is reported as a detected cycle when convergence has
not already been reached.

The oracle is intentionally independent of a network solver. A caller may use
a reduced feeder map, a full quasi-static re-solve, or a measured/discrete-time
map. No claim about a physical response time follows unless the supplied map
contains the relevant sampling, delay, filtering, and ramp dynamics.
"""
function fixed_point_oracle(
        map, initial;
        relaxation::Real=1.0,
        max_iterations::Integer=100,
        atol::Real=1e-8,
        rtol::Real=1e-8,
        cycle_window::Integer=8,
        cycle_tolerance::Real=1e-8,
        store_trajectory::Bool=true)
    x = _real_state(initial, "initial state")
    0 < relaxation <= 1 || throw(ArgumentError(
        "relaxation must lie in (0, 1]"))
    max_iterations >= 1 || throw(ArgumentError(
        "max_iterations must be >= 1"))
    atol >= 0 && isfinite(atol) || throw(ArgumentError(
        "atol must be finite and >= 0"))
    rtol >= 0 && isfinite(rtol) || throw(ArgumentError(
        "rtol must be finite and >= 0"))
    cycle_window >= 1 || throw(ArgumentError(
        "cycle_window must be >= 1"))
    cycle_tolerance >= 0 && isfinite(cycle_tolerance) || throw(ArgumentError(
        "cycle_tolerance must be finite and >= 0"))
    relaxation_value = Float64(relaxation)
    atol_value = Float64(atol)
    rtol_value = Float64(rtol)
    cycle_tolerance_value = Float64(cycle_tolerance)
    states = Vector{Vector{Float64}}()
    history = Vector{Tuple{Int,Vector{Float64}}}([(0, copy(x))])
    store_trajectory && push!(states, copy(x))
    converged = false
    cycled = false
    cycle_period = 0
    residual = Inf
    iterations = 0
    for iteration in 1:Int(max_iterations)
        mapped = _real_state(map(x), "fixed-point map output")
        length(mapped) == length(x) || throw(DimensionMismatch(
            "fixed-point map returned length $(length(mapped)); expected " *
            "$(length(x))"))
        next = _real_state(
            (1 - relaxation_value).*x + relaxation_value.*mapped,
            "fixed-point iterate")
        residual = norm(next - x)
        threshold = atol_value + rtol_value * max(norm(x), norm(next), 1.0)
        iterations = iteration
        push!(history, (iteration, copy(next)))
        store_trajectory && push!(states, copy(next))
        if residual <= threshold
            x = next
            converged = true
            break
        end
        if cycle_window > 0 && length(history) > 1 &&
                residual > threshold && residual > cycle_tolerance_value
            first_index = max(1, length(history) - 1 - cycle_window)
            for index in first_index:(length(history) - 1)
                if norm(next - history[index][2]) <= cycle_tolerance_value
                    cycle_period = iteration - history[index][1]
                    cycled = true
                    break
                end
            end
        end
        length(history) > cycle_window + 1 && deleteat!(history, 1)
        x = next
        cycled && break
    end
    trajectory = if store_trajectory && !isempty(states)
        reduce(hcat, states)
    else
        Matrix{Float64}(undef, length(x), 0)
    end
    FixedPointIterationResult(
        converged, cycled, iterations, cycle_period, Float64(residual),
        x, trajectory)
end

"""Compute a central finite-difference Jacobian of a real vector map."""
function finite_difference_jacobian(
        map, point;
        step::Real=1e-6,
        relative_step::Bool=true)
    x = _real_state(point, "finite-difference point")
    h = Float64(step)
    isfinite(h) && h > 0 || throw(ArgumentError(
        "step must be finite and > 0"))
    y = _real_state(map(x), "finite-difference map output")
    n = length(x)
    length(y) == n || throw(DimensionMismatch(
        "finite-difference map output length must equal input length"))
    jacobian = Matrix{Float64}(undef, n, n)
    for column in 1:n
        delta = relative_step ? h * max(abs(x[column]), 1.0) : h
        plus = copy(x); plus[column] += delta
        minus = copy(x); minus[column] -= delta
        y_plus = _real_state(map(plus), "finite-difference map output")
        y_minus = _real_state(map(minus), "finite-difference map output")
        length(y_plus) == n && length(y_minus) == n || throw(DimensionMismatch(
            "finite-difference map output length changed"))
        jacobian[:, column] = (y_plus - y_minus) / (2delta)
    end
    jacobian
end

"""
    screen_fixed_point_gain(map, point; step=1e-6, threshold=1.0)

Return the finite-difference Jacobian and local spectral-radius screen for a
fixed-point map. `local_contractive` means only that the computed linearization
has spectral radius below `threshold`; it is not a global or nonlinear
stability certificate.
"""
function screen_fixed_point_gain(
        map, point;
        step::Real=1e-6,
        relative_step::Bool=true,
        threshold::Real=1.0)
    threshold_value = Float64(threshold)
    isfinite(threshold_value) && threshold_value > 0 || throw(ArgumentError(
        "threshold must be finite and > 0"))
    jacobian = finite_difference_jacobian(
        map, point; step=step, relative_step=relative_step)
    eigenvalues = ComplexF64.(eigvals(jacobian))
    spectral_radius = isempty(eigenvalues) ? 0.0 : maximum(abs, eigenvalues)
    induced_inf_norm = opnorm(jacobian, Inf)
    FixedPointGainScreen(
        jacobian, eigenvalues, spectral_radius, induced_inf_norm,
        spectral_radius < threshold_value, threshold_value - spectral_radius,
        Float64(step))
end

function _phasor_state(values)
    vector = _real_state(values, "phasor state")
    length(vector) == 6 || throw(DimensionMismatch(
        "a three-phase phasor state needs six real coordinates"))
    ntuple(index -> complex(vector[2index - 1], vector[2index]), 3)
end

function _phasor_vector(values::NTuple{3,<:Complex})
    result = Vector{Float64}(undef, 6)
    for index in 1:3
        result[2index - 1] = real(values[index])
        result[2index] = imag(values[index])
    end
    result
end

"""Finite-difference current-command Jacobian for the exact inverter law."""
function inverter_control_current_jacobian(
        controller::SequenceController,
        measurement::InverterControlMeasurement,
        request::InverterControlRequest,
        ratings::InverterControlRatings;
        step::Real=1e-6,
        relative_step::Bool=true)
    point = _phasor_vector(measurement.phase_voltage)
    map = state -> _phasor_vector(evaluate_exact(
        controller, InverterControlMeasurement(_phasor_state(state)),
        request, ratings).phase_current)
    finite_difference_jacobian(
        map, point; step=step, relative_step=relative_step)
end

"""Result of the reduced affine-feeder exact controller fixed-point oracle."""
struct InverterControlFixedPointResult
    oracle::FixedPointIterationResult
    reference_measurement::InverterControlMeasurement
    final_measurement::InverterControlMeasurement
    final_control::InverterControlResult
    voltage_sensitivity::Matrix{Float64}
end

"""Result of a plant-backed exact-law equilibrium iteration."""
struct InverterControlNetworkFixedPointResult
    smooth::ControlledInverterResult
    converged::Bool
    cycled::Bool
    iterations::Int
    cycle_period::Int
    residual_voltage_inf_V::Float64
    residual_current_inf_A::Float64
    exact_control::InverterControlResult
    plant::Union{Nothing,NamedTuple}
    bus::Dict{String,Any}
    solve::SolveStatus
    voltage_trajectory::Matrix{Float64}
end

function _validated_voltage_sensitivity(voltage_sensitivity::AbstractMatrix)
    sensitivity = Matrix{Float64}(voltage_sensitivity)
    size(sensitivity) == (6, 6) || throw(DimensionMismatch(
        "voltage_sensitivity must be a 6×6 real matrix"))
    all(isfinite, sensitivity) || throw(ArgumentError(
        "voltage_sensitivity must be finite"))
    sensitivity
end

"""
    inverter_control_fixed_point_oracle(controller, initial, request, ratings,
        reference, voltage_sensitivity; kwargs...) -> InverterControlFixedPointResult

Run the exact controller law against a supplied affine feeder sensitivity. The
reduced map is

``v_next = v_ref + Z (i_exact(v) - i_exact(v_ref))``

in six real rectangular phase-voltage coordinates. `reference` is the voltage
at which the feeder operating point and the reference controller current are
defined; it is not inferred from the initial condition. This is a lightweight
equilibrium oracle for campaign screening. It does not replace a full network
re-solve, and its convergence has no physical time interpretation unless the
supplied sensitivity/map is itself dynamic.
"""
function inverter_control_fixed_point_oracle(
        controller::SequenceController,
        initial::InverterControlMeasurement,
        request::InverterControlRequest,
        ratings::InverterControlRatings,
        reference::InverterControlMeasurement,
        voltage_sensitivity::AbstractMatrix;
        kwargs...)
    sensitivity = _validated_voltage_sensitivity(voltage_sensitivity)
    reference_state = _phasor_vector(reference.phase_voltage)
    initial_state = _phasor_vector(initial.phase_voltage)
    reference_current = _phasor_vector(evaluate_exact(
        controller, reference, request, ratings).phase_current)
    current_map = state -> _phasor_vector(evaluate_exact(
        controller, InverterControlMeasurement(_phasor_state(state)),
        request, ratings).phase_current)
    map = state -> reference_state + sensitivity *
        (current_map(state) - reference_current)
    oracle = fixed_point_oracle(map, initial_state; kwargs...)
    final_measurement = InverterControlMeasurement(
        _phasor_state(oracle.final_state))
    final_control = evaluate_exact(
        controller, final_measurement, request, ratings)
    InverterControlFixedPointResult(
        oracle, reference, final_measurement, final_control, sensitivity)
end

function _solve_exact_current_plant(
        net::Dict{String,Any}, inverter::AdvancedInverter,
        current::NTuple{3,ComplexF64}, target::AbstractCurrentTarget;
        s_base::Float64=1e6, optimizer=Ipopt.Optimizer,
        verbose::Bool=false, solver_options=())
    handles = Ref{_InvHandles}()
    hook! = ctx -> begin
        h = stamp_device!(ctx, inverter)
        handles[] = h
        target_re, target_im = target isa ConverterCurrentTarget ?
            (h.cri, h.cii) : (h.gri, h.gii)
        for k in 1:3
            JuMP.@constraint(_opf_model(ctx),
                target_re[k] == real(current[k]) / h.ib)
            JuMP.@constraint(_opf_model(ctx),
                target_im[k] == imag(current[k]) / h.ib)
        end
        JuMP.@objective(_opf_model(ctx), Min, h.p_loss + h.p_cap_loss)
    end
    ctx = build_opf_model(
        net; per_unit=true, s_base=s_base, add_objective=false,
        model_hook! = hook!, optimizer=optimizer, verbose=verbose)
    _set_solver_options!(_opf_model(ctx), solver_options)
    enforce_kcl!(ctx)
    JuMP.optimize!(_opf_model(ctx))
    outcome = _solve_outcome(_opf_model(ctx))
    status = SolveStatus(outcome)
    h = handles[]
    plant = extract_device(inverter, h, status)
    network = _extract_result(ctx, outcome)
    (status=status, plant=plant, network=network)
end

"""Finite-difference feeder voltage-sensitivity result with solve provenance."""
struct InverterControlNetworkSensitivityResult
    sensitivity::Matrix{Float64}
    base_voltage::NTuple{3,ComplexF64}
    base_current::NTuple{3,ComplexF64}
    perturbation_step_A::Float64
    base_status::SolveStatus
    plus_status::Vector{SolveStatus}
    minus_status::Vector{SolveStatus}
end

"""Finite-difference fleet voltage-sensitivity result with solve provenance."""
struct ControlledInverterFleetNetworkSensitivityResult
    sensitivity::Matrix{Float64}
    device_ids::Vector{String}
    base_voltage::Dict{String,NTuple{3,ComplexF64}}
    base_currents::Dict{String,NTuple{3,ComplexF64}}
    perturbation_step_A::Float64
    base_status::SolveStatus
    plus_status::Vector{SolveStatus}
    minus_status::Vector{SolveStatus}
end

function _current_tuple(current)
    values = ComplexF64.(collect(current))
    length(values) == 3 || throw(DimensionMismatch(
        "operating_current must contain three phase phasors"))
    all(v -> isfinite(real(v)) && isfinite(imag(v)), values) ||
        throw(ArgumentError("operating_current must be finite"))
    Tuple(values)
end

"""
    inverter_control_network_voltage_sensitivity(net, controlled,
        operating_current; step=1e-3, kwargs...) -> InverterControlNetworkSensitivityResult

Estimate the 6×6 real rectangular POC voltage sensitivity to the selected
converter/grid current target by central finite differences of fresh physical
AdvancedInverter solves. Columns are ordered
`(Re Ia, Im Ia, Re Ib, Im Ib, Re Ic, Im Ic)` and rows use the same voltage
ordering. The perturbation is `step * max(abs(component), 1 A)` in amperes.

Failed plus/minus solves leave the corresponding column as `NaN` and retain
their full `SolveStatus`; the function does not turn a physically infeasible
perturbation into a numerical sensitivity. Use the returned matrix with
[`inverter_control_loop_gain`](@ref) only when all required columns are finite.
"""
function inverter_control_network_voltage_sensitivity(
        net::Dict{String,Any},
        controlled::ControlledDevice{<:AdvancedInverter},
        operating_current;
        step::Real=1e-3,
        s_base::Real=1e6,
        optimizer=Ipopt.Optimizer,
        verbose::Bool=false,
        solver_options=())
    h = Float64(step)
    isfinite(h) && h > 0 || throw(ArgumentError(
        "step must be finite and > 0"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError(
        "s_base must be finite and > 0"))
    base_current = _current_tuple(operating_current)
    inverter = inverter_spec(controlled)
    target = controlled.controller.current_target
    solve_current(current) = _solve_exact_current_plant(
        net, inverter, Tuple(current), target; s_base=Float64(s_base),
        optimizer=optimizer, verbose=verbose, solver_options=solver_options)
    base = solve_current(base_current)
    nan_voltage = ntuple(_ -> ComplexF64(NaN, NaN), 3)
    if !base.status.publishable
        return InverterControlNetworkSensitivityResult(
            fill(NaN, 6, 6), nan_voltage, base_current, h,
            base.status, SolveStatus[], SolveStatus[])
    end
    base_voltage = _network_phase_state(base.network["bus"], inverter)
    point = _phasor_vector(base_current)
    sensitivity = fill(NaN, 6, 6)
    plus_status = SolveStatus[]
    minus_status = SolveStatus[]
    for column in 1:6
        delta = h * max(abs(point[column]), 1.0)
        plus = copy(point); plus[column] += delta
        minus = copy(point); minus[column] -= delta
        plus_result = solve_current(_phasor_state(plus))
        minus_result = solve_current(_phasor_state(minus))
        push!(plus_status, plus_result.status)
        push!(minus_status, minus_result.status)
        if plus_result.status.publishable && minus_result.status.publishable
            plus_voltage = _network_phase_state(
                plus_result.network["bus"], inverter)
            minus_voltage = _network_phase_state(
                minus_result.network["bus"], inverter)
            sensitivity[:, column] = (_phasor_vector(plus_voltage) -
                _phasor_vector(minus_voltage)) / (2delta)
        end
    end
    InverterControlNetworkSensitivityResult(
        sensitivity, base_voltage, base_current, h, base.status,
        plus_status, minus_status)
end

function _validated_fleet_currents(
        spec::ControlledInverterFleetSpec, currents::AbstractDict)
    ids = sort!(collect(keys(spec.devices)))
    canonical = Dict{String,NTuple{3,ComplexF64}}()
    for (id, value) in currents
        canonical_id = String(id)
        haskey(canonical, canonical_id) && throw(ArgumentError(
            "duplicate fleet current identifier '$canonical_id' after string conversion"))
        canonical[canonical_id] = _current_tuple(value)
    end
    Set(keys(canonical)) == Set(ids) || throw(ArgumentError(
        "operating_currents keys must match the controlled fleet ids"))
    canonical, ids
end

function _fleet_current_state(
        currents::Dict{String,NTuple{3,ComplexF64}}, ids)
    reduce(vcat, [_phasor_vector(currents[id]) for id in ids])
end

function _fleet_current_dict(state::AbstractVector, ids)
    length(state) == 6 * length(ids) || throw(DimensionMismatch(
        "fleet current state has length $(length(state)); expected " *
        "$(6 * length(ids))"))
    Dict(id => _phasor_state(state[6 * (index - 1) + 1:6 * index])
         for (index, id) in enumerate(ids))
end

function _fleet_bus_voltage(
        bus::AbstractDict, spec::ControlledInverterFleetSpec, ids)
    Dict(id => _network_phase_state(
        bus, inverter_spec(spec.devices[id])) for id in ids)
end

"""
    controlled_inverter_fleet_network_voltage_sensitivity(net, spec,
        operating_currents; step=1e-3, kwargs...)
        -> ControlledInverterFleetNetworkSensitivityResult

Estimate the real rectangular feeder voltage sensitivity for a selected fleet
by central finite differences of simultaneous physical fleet solves. The
returned matrix has `6N` rows and columns for `N` devices, with blocks ordered
by sorted device id and each device using
`(Re Ia, Im Ia, Re Ib, Im Ib, Re Ic, Im Ic)`. A column is retained only when
both fresh solves are publishable; failed perturbations remain `NaN` while
their `SolveStatus` values are preserved for auditability.
"""
function controlled_inverter_fleet_network_voltage_sensitivity(
        net::Dict{String,Any}, spec::ControlledInverterFleetSpec,
        operating_currents::AbstractDict;
        step::Real=1e-3,
        s_base::Real=1e6,
        optimizer=Ipopt.Optimizer,
        verbose::Bool=false,
        solver_options=())
    h = Float64(step)
    isfinite(h) && h > 0 || throw(ArgumentError(
        "step must be finite and > 0"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError(
        "s_base must be finite and > 0"))
    _validate_controlled_fleet(net, spec)
    base_currents, ids = _validated_fleet_currents(spec, operating_currents)
    solve_currents(values) = _solve_exact_current_fleet_plant(
        net, spec, values; s_base=Float64(s_base), optimizer=optimizer,
        verbose=verbose, solver_options=solver_options)
    base = solve_currents(base_currents)
    n = 6 * length(ids)
    nan_voltage = Dict(id => ntuple(_ -> ComplexF64(NaN, NaN), 3) for id in ids)
    if !base.status.publishable
        return ControlledInverterFleetNetworkSensitivityResult(
            fill(NaN, n, n), ids, nan_voltage, base_currents, h,
            base.status, SolveStatus[], SolveStatus[])
    end
    base_voltage = _fleet_bus_voltage(base.network["bus"], spec, ids)
    point = _fleet_current_state(base_currents, ids)
    sensitivity = fill(NaN, n, n)
    plus_status = SolveStatus[]
    minus_status = SolveStatus[]
    for column in 1:n
        delta = h * max(abs(point[column]), 1.0)
        plus = copy(point); plus[column] += delta
        minus = copy(point); minus[column] -= delta
        plus_result = solve_currents(_fleet_current_dict(plus, ids))
        minus_result = solve_currents(_fleet_current_dict(minus, ids))
        push!(plus_status, plus_result.status)
        push!(minus_status, minus_result.status)
        if plus_result.status.publishable && minus_result.status.publishable
            plus_voltage = _fleet_bus_voltage(plus_result.network["bus"], spec, ids)
            minus_voltage = _fleet_bus_voltage(minus_result.network["bus"], spec, ids)
            sensitivity[:, column] = (_fleet_voltage_state(plus_voltage, ids) -
                _fleet_voltage_state(minus_voltage, ids)) / (2delta)
        end
    end
    ControlledInverterFleetNetworkSensitivityResult(
        sensitivity, ids, base_voltage, base_currents, h, base.status,
        plus_status, minus_status)
end

"""Return deterministic long-form rows for a fleet sensitivity matrix."""
function controlled_inverter_fleet_network_sensitivity_rows(
        result::ControlledInverterFleetNetworkSensitivityResult)
    coordinates = NamedTuple[
        (device_id=id, phase=phase, component=component)
        for id in result.device_ids for phase in 1:3
        for component in (:real, :imag)]
    rows = NamedTuple[]
    n = length(coordinates)
    for output in 1:n, input in 1:n
        plus_ok = input <= length(result.plus_status) &&
            result.plus_status[input].publishable
        minus_ok = input <= length(result.minus_status) &&
            result.minus_status[input].publishable
        output_coordinate = coordinates[output]
        input_coordinate = coordinates[input]
        push!(rows, (
            output_device_id=output_coordinate.device_id,
            output_phase=output_coordinate.phase,
            output_component=output_coordinate.component,
            input_device_id=input_coordinate.device_id,
            input_phase=input_coordinate.phase,
            input_component=input_coordinate.component,
            voltage_sensitivity_V_per_A=result.sensitivity[output, input],
            base_publishable=result.base_status.publishable,
            plus_publishable=plus_ok,
            minus_publishable=minus_ok,
        ))
    end
    rows
end

function _validated_fleet_measurements(
        spec::ControlledInverterFleetSpec, measurements::AbstractDict)
    ids = sort!(collect(keys(spec.devices)))
    canonical = Dict{String,InverterControlMeasurement}()
    for (id, value) in measurements
        canonical_id = String(id)
        haskey(canonical, canonical_id) && throw(ArgumentError(
            "duplicate fleet measurement identifier '$canonical_id' after string conversion"))
        canonical[canonical_id] = value isa InverterControlMeasurement ?
            value : InverterControlMeasurement(ComplexF64.(collect(value)))
    end
    Set(keys(canonical)) == Set(ids) || throw(ArgumentError(
        "measurements keys must match the controlled fleet ids"))
    canonical, ids
end

"""
    controlled_inverter_fleet_loop_gain(spec, measurements,
        voltage_sensitivity; kwargs...) -> FixedPointGainScreen

Screen the simultaneous fleet loop
`Δv ↦ Z⋅diag(∂iₖ/∂vₖ)⋅Δv` in real rectangular coordinates. Device blocks are
ordered by sorted fleet id and each block uses the declared a-b-c phase order.
The supplied `voltage_sensitivity` must therefore be the `6N × 6N` matrix
returned by [`controlled_inverter_fleet_network_voltage_sensitivity`](@ref),
or an equivalent reduced/network map. This is a local quasi-static diagnostic,
not a dynamic stability certificate.
"""
function controlled_inverter_fleet_loop_gain(
        spec::ControlledInverterFleetSpec,
        measurements::AbstractDict,
        voltage_sensitivity::AbstractMatrix;
        step::Real=1e-6,
        relative_step::Bool=true,
        threshold::Real=1.0)
    validated, ids = _validated_fleet_measurements(spec, measurements)
    n = 6 * length(ids)
    sensitivity = Matrix{Float64}(voltage_sensitivity)
    size(sensitivity) == (n, n) || throw(DimensionMismatch(
        "fleet voltage_sensitivity must be a $(n)×$(n) real matrix"))
    all(isfinite, sensitivity) || throw(ArgumentError(
        "fleet voltage_sensitivity must be finite"))
    point = _fleet_voltage_state(
        Dict(id => validated[id].phase_voltage for id in ids), ids)
    controller_jacobian = zeros(Float64, n, n)
    current0 = Vector{Float64}(undef, n)
    for (index, id) in enumerate(ids)
        controlled = spec.devices[id]
        ratings = InverterControlRatings(
            inverter_spec(controlled), controlled.controller.current_target)
        request = spec.requests[id]
        measurement = validated[id]
        local_jacobian = inverter_control_current_jacobian(
            controlled.controller, measurement, request, ratings;
            step=step, relative_step=relative_step)
        block = 6 * (index - 1) + 1:6 * index
        controller_jacobian[block, block] = local_jacobian
        current0[block] = _phasor_vector(evaluate_exact(
            controlled.controller, measurement, request, ratings).phase_current)
    end
    map = state -> begin
        currents = Vector{Float64}(undef, n)
        for (index, id) in enumerate(ids)
            controlled = spec.devices[id]
            ratings = InverterControlRatings(
                inverter_spec(controlled), controlled.controller.current_target)
            request = spec.requests[id]
            block = 6 * (index - 1) + 1:6 * index
            measurement = InverterControlMeasurement(
                _phasor_state(state[block]))
            currents[block] = _phasor_vector(evaluate_exact(
                controlled.controller, measurement, request, ratings).phase_current)
        end
        point + sensitivity * (currents - current0)
    end
    screen = screen_fixed_point_gain(
        map, point; step=step, relative_step=relative_step,
        threshold=threshold)
    isapprox(screen.jacobian, sensitivity * controller_jacobian;
             atol=1e-7, rtol=1e-6) ||
        throw(ArgumentError("fleet loop-gain finite-difference inconsistency"))
    screen
end

function _network_phase_state(
        bus::AbstractDict, inverter::AdvancedInverter)
    bus_data = get(bus, inverter.bus, nothing)
    bus_data isa AbstractDict || throw(ArgumentError(
        "network result does not contain inverter bus '$(inverter.bus)'"))
    ntuple(3) do k
        phase = inverter.phase_terminals[k]
        values = get(bus_data, phase, nothing)
        values isa AbstractDict || throw(ArgumentError(
            "network result does not contain phase '$phase' at bus " *
            "'$(inverter.bus)'"))
        vr = Float64(get(values, "vr", NaN))
        vi = Float64(get(values, "vi", NaN))
        isfinite(vr) && isfinite(vi) || throw(ArgumentError(
            "network result has non-finite voltage at '$phase'"))
        ComplexF64(vr, vi)
    end
end

"""
    solve_inverter_control_network_fixed_point(net, controlled, request; kwargs...)
        -> InverterControlNetworkFixedPointResult

Compare the smooth controlled solve with a plant-backed fixed point of the
exact controller law. Each iteration evaluates the exact local law, fixes the
selected converter- or grid-side phase current in a fresh AdvancedInverter
solve, and feeds the resulting POC voltage into the next iteration. Physical
limits therefore remain part of the oracle rather than being replaced by a
linear sensitivity.

The result retains non-publishable final solves and cycling diagnostics. This
is a quasi-static equilibrium comparison: it is not a dynamic stability or
response-time simulation.
"""
function solve_inverter_control_network_fixed_point(
        net::Dict{String,Any},
        controlled::ControlledDevice{<:AdvancedInverter},
        request::InverterControlRequest;
        smooth_result::Union{Nothing,ControlledInverterResult}=nothing,
        initial_voltage=nothing,
        max_iterations::Integer=20,
        atol_voltage::Real=1e-6,
        atol_current::Real=1e-6,
        rtol::Real=1e-8,
        cycle_window::Integer=8,
        cycle_tolerance::Real=1e-7,
        s_base::Real=1e6,
        optimizer=Ipopt.Optimizer,
        verbose::Bool=false,
        solver_options=())
    max_iterations >= 1 || throw(ArgumentError(
        "max_iterations must be >= 1"))
    0 <= atol_voltage && isfinite(atol_voltage) || throw(ArgumentError(
        "atol_voltage must be finite and >= 0"))
    0 <= atol_current && isfinite(atol_current) || throw(ArgumentError(
        "atol_current must be finite and >= 0"))
    0 <= rtol && isfinite(rtol) || throw(ArgumentError(
        "rtol must be finite and >= 0"))
    cycle_window >= 1 || throw(ArgumentError(
        "cycle_window must be >= 1"))
    0 <= cycle_tolerance && isfinite(cycle_tolerance) || throw(ArgumentError(
        "cycle_tolerance must be finite and >= 0"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError(
        "s_base must be finite and > 0"))

    smooth = smooth_result === nothing ? solve_controlled_inverter(
        net, controlled, request; per_unit=true, s_base=s_base,
        optimizer=optimizer, verbose=verbose, solver_options=solver_options) :
        smooth_result
    if !smooth.solve.publishable
        return InverterControlNetworkFixedPointResult(
            smooth, false, false, 0, 0, Inf, Inf, smooth.exact_control, nothing,
            Dict{String,Any}(), smooth.solve,
            Matrix{Float64}(undef, 6, 0))
    end
    inverter = inverter_spec(controlled)
    ratings = InverterControlRatings(inverter, controlled.controller.current_target)
    initial = if initial_voltage === nothing
        smooth.control.phase_voltage
    elseif initial_voltage isa InverterControlMeasurement
        initial_voltage.phase_voltage
    else
        Tuple(ComplexF64.(collect(initial_voltage)))
    end
    length(initial) == 3 || throw(DimensionMismatch(
        "initial_voltage must contain three phase phasors"))
    voltage = ComplexF64.(initial)
    history = Vector{Vector{Float64}}([_phasor_vector(voltage)])
    trajectory = Vector{Vector{Float64}}([_phasor_vector(voltage)])
    exact_control = evaluate_exact(
        controlled.controller, InverterControlMeasurement(voltage),
        request, ratings)
    last_plant = nothing
    last_bus = Dict{String,Any}()
    last_status = smooth.solve
    converged = false
    cycled = false
    iterations = 0
    cycle_period = 0
    residual_voltage = Inf
    residual_current = Inf
    for iteration in 1:Int(max_iterations)
        iterations = iteration
        measurement = InverterControlMeasurement(voltage)
        exact_control = evaluate_exact(
            controlled.controller, measurement, request, ratings)
        solved = _solve_exact_current_plant(
            net, inverter, exact_control.phase_current,
            controlled.controller.current_target;
            s_base=Float64(s_base), optimizer=optimizer, verbose=verbose,
            solver_options=solver_options)
        last_status = solved.status
        last_plant = solved.plant
        network = solved.network
        last_bus = Dict{String,Any}(
            String(k) => v for (k, v) in network["bus"])
        if !last_status.publishable
            break
        end
        next_voltage = _network_phase_state(last_bus, inverter)
        next_measurement = InverterControlMeasurement(next_voltage)
        next_control = evaluate_exact(
            controlled.controller, next_measurement, request, ratings)
        residual_voltage = maximum(abs, next_voltage .- voltage)
        residual_current = maximum(abs,
            next_control.phase_current .- exact_control.phase_current)
        push!(trajectory, _phasor_vector(next_voltage))
        if residual_voltage <= Float64(atol_voltage) + Float64(rtol) *
                max(maximum(abs, voltage), maximum(abs, next_voltage), 1.0) &&
                residual_current <= Float64(atol_current) + Float64(rtol) *
                max(maximum(abs, exact_control.phase_current),
                    maximum(abs, next_control.phase_current), 1.0)
            voltage = next_voltage
            exact_control = next_control
            converged = true
            break
        end
        next_state = _phasor_vector(next_voltage)
        history_start = max(1, length(history) - cycle_window + 1)
        for index in history_start:length(history)
            if norm(next_state - history[index]) <= Float64(cycle_tolerance)
                cycled = true
                cycle_period = length(history) - index + 1
                voltage = next_voltage
                exact_control = next_control
                break
            end
        end
        cycled && break
        push!(history, next_state)
        voltage = next_voltage
    end
    final_bus = last_bus
    final_plant = last_plant === nothing ? nothing : last_plant
    InverterControlNetworkFixedPointResult(
        smooth, converged, cycled, iterations, cycle_period,
        Float64(residual_voltage), Float64(residual_current), exact_control,
        final_plant, final_bus, last_status, reduce(hcat, trajectory))
end

"""Result of a simultaneous plant-backed exact-law fleet iteration."""
struct ControlledInverterFleetNetworkFixedPointResult
    smooth::ControlledInverterFleetResult
    converged::Bool
    cycled::Bool
    iterations::Int
    cycle_period::Int
    residual_voltage_inf_V::Float64
    residual_current_inf_A::Float64
    exact_controls::Dict{String,InverterControlResult}
    plants::Dict{String,Any}
    bus::Dict{String,Any}
    solve::SolveStatus
    voltage_trajectory::Matrix{Float64}
end

function _solve_exact_current_fleet_plant(
        net::Dict{String,Any}, spec::ControlledInverterFleetSpec,
        currents::Dict{String,NTuple{3,ComplexF64}};
        s_base::Float64=1e6, optimizer=Ipopt.Optimizer,
        verbose::Bool=false, solver_options=())
    handles = Dict{String,_InvHandles}()
    builder = BMOPFTools.OpfDeviceBuilder(:PowerOptLab, function (ctx, ids)
        for id in ids
            _disable_native_ibr_port!(ctx, id)
            inverter = inverter_spec(spec.devices[id])
            h = stamp_device!(ctx, inverter)
            handles[id] = h
            target = spec.devices[id].controller.current_target
            target_re, target_im = target isa ConverterCurrentTarget ?
                (h.cri, h.cii) : (h.gri, h.gii)
            current = currents[id]
            for k in 1:3
                JuMP.@constraint(_opf_model(ctx),
                    target_re[k] == real(current[k]) / h.ib)
                JuMP.@constraint(_opf_model(ctx),
                    target_im[k] == imag(current[k]) / h.ib)
            end
        end
        ctx
    end)
    build_spec = BMOPFTools.OpfBuildSpec(component_builders=Dict(
        (:ibr, id) => builder for id in keys(spec.devices)))
    hook! = ctx -> JuMP.@objective(_opf_model(ctx), Min,
        sum(h.p_loss + h.p_cap_loss
            for h in values(handles)))
    ctx = build_opf_model(
        net; per_unit=true, s_base=s_base, add_objective=false,
        build_spec=build_spec, model_hook! = hook!, optimizer=optimizer,
        verbose=verbose)
    _set_solver_options!(_opf_model(ctx), solver_options)
    enforce_kcl!(ctx)
    JuMP.optimize!(_opf_model(ctx))
    outcome = _solve_outcome(_opf_model(ctx))
    status = SolveStatus(outcome)
    plants = Dict{String,Any}()
    for id in sort!(collect(keys(spec.devices)))
        plants[id] = extract_device(
            inverter_spec(spec.devices[id]), handles[id], status)
    end
    (status=status, plants=plants, network=_extract_result(ctx, outcome))
end

function _fleet_voltage_state(voltage_by_id::Dict{String,NTuple{3,ComplexF64}}, ids)
    reduce(vcat, [_phasor_vector(voltage_by_id[id]) for id in ids])
end

"""
    solve_controlled_inverter_fleet_network_fixed_point(net, spec; kwargs...)
        -> ControlledInverterFleetNetworkFixedPointResult

Run a simultaneous plant-backed fixed point for every selected fleet device.
The smooth fleet solve is used as the initial state; each exact iteration fixes
all selected converter/grid current targets in one fresh network solve. This
preserves feeder coupling and native unselected IBRs, and is the campaign-scale
equilibrium comparison for [`solve_controlled_inverter_fleet`](@ref).
"""
function solve_controlled_inverter_fleet_network_fixed_point(
        net::Dict{String,Any}, spec::ControlledInverterFleetSpec;
        smooth_result::Union{Nothing,ControlledInverterFleetResult}=nothing,
        initial_voltage=nothing,
        max_iterations::Integer=20,
        atol_voltage::Real=1e-6,
        atol_current::Real=1e-6,
        rtol::Real=1e-8,
        cycle_window::Integer=8,
        cycle_tolerance::Real=1e-7,
        s_base::Real=1e6,
        optimizer=Ipopt.Optimizer,
        verbose::Bool=false,
        selection_objective::Symbol=:loss,
        solver_options=())
    max_iterations >= 1 || throw(ArgumentError(
        "max_iterations must be >= 1"))
    0 <= atol_voltage && isfinite(atol_voltage) || throw(ArgumentError(
        "atol_voltage must be finite and >= 0"))
    0 <= atol_current && isfinite(atol_current) || throw(ArgumentError(
        "atol_current must be finite and >= 0"))
    0 <= rtol && isfinite(rtol) || throw(ArgumentError(
        "rtol must be finite and >= 0"))
    cycle_window >= 1 || throw(ArgumentError(
        "cycle_window must be >= 1"))
    0 <= cycle_tolerance && isfinite(cycle_tolerance) || throw(ArgumentError(
        "cycle_tolerance must be finite and >= 0"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError(
        "s_base must be finite and > 0"))
    _validate_controlled_fleet(net, spec)
    smooth = smooth_result === nothing ? solve_controlled_inverter_fleet(
        net, spec; s_base=s_base, optimizer=optimizer, verbose=verbose,
        selection_objective=selection_objective, solver_options=solver_options) :
        smooth_result
    ids = sort!(collect(keys(spec.devices)))
    if !smooth.solve.publishable
        return ControlledInverterFleetNetworkFixedPointResult(
            smooth, false, false, 0, 0, Inf, Inf,
            Dict{String,InverterControlResult}(), Dict{String,Any}(),
            Dict{String,Any}(), smooth.solve,
            Matrix{Float64}(undef, 3*length(ids), 0))
    end
    voltage_by_id = Dict{String,NTuple{3,ComplexF64}}(
        id => smooth.devices[id].control.phase_voltage for id in ids)
    if initial_voltage !== nothing
        initial_voltage isa AbstractDict || throw(ArgumentError(
            "initial_voltage must map fleet device ids to three phase phasors"))
        Set(String.(collect(keys(initial_voltage)))) == Set(ids) || throw(ArgumentError(
            "initial_voltage keys must match the controlled fleet ids"))
        for id in ids
            value = initial_voltage[id]
            phase_voltage = value isa InverterControlMeasurement ?
                value.phase_voltage : Tuple(ComplexF64.(collect(value)))
            length(phase_voltage) == 3 || throw(DimensionMismatch(
                "initial_voltage[$id] must contain three phase phasors"))
            all(v -> isfinite(real(v)) && isfinite(imag(v)), phase_voltage) ||
                throw(ArgumentError("initial_voltage[$id] must be finite"))
            voltage_by_id[id] = phase_voltage
        end
    end
    history = Vector{Vector{Float64}}([_fleet_voltage_state(voltage_by_id, ids)])
    trajectory = Vector{Vector{Float64}}(copy(history))
    exact_controls = Dict{String,InverterControlResult}()
    plants = Dict{String,Any}()
    last_bus = Dict{String,Any}()
    last_status = smooth.solve
    converged = false
    cycled = false
    iterations = 0
    cycle_period = 0
    residual_voltage = Inf
    residual_current = Inf
    for iteration in 1:Int(max_iterations)
        iterations = iteration
        currents = Dict{String,NTuple{3,ComplexF64}}()
        for id in ids
            inverter = inverter_spec(spec.devices[id])
            ratings = InverterControlRatings(
                inverter, spec.devices[id].controller.current_target)
            control = evaluate_exact(
                spec.devices[id].controller,
                InverterControlMeasurement(voltage_by_id[id]),
                spec.requests[id], ratings)
            exact_controls[id] = control
            currents[id] = control.phase_current
        end
        solved = _solve_exact_current_fleet_plant(
            net, spec, currents; s_base=Float64(s_base), optimizer=optimizer,
            verbose=verbose, solver_options=solver_options)
        last_status = solved.status
        plants = solved.plants
        last_bus = Dict{String,Any}(
            String(k) => v for (k, v) in solved.network["bus"])
        if !last_status.publishable
            break
        end
        next_voltage_by_id = Dict{String,NTuple{3,ComplexF64}}(
            id => _network_phase_state(last_bus, inverter_spec(spec.devices[id]))
            for id in ids)
        next_controls = Dict{String,InverterControlResult}()
        for id in ids
            inverter = inverter_spec(spec.devices[id])
            ratings = InverterControlRatings(
                inverter, spec.devices[id].controller.current_target)
            next_controls[id] = evaluate_exact(
                spec.devices[id].controller,
                InverterControlMeasurement(next_voltage_by_id[id]),
                spec.requests[id], ratings)
        end
        residual_voltage = maximum(abs,
            _fleet_voltage_state(next_voltage_by_id, ids) -
            _fleet_voltage_state(voltage_by_id, ids))
        residual_current = maximum(abs,
            reduce(vcat, [collect(next_controls[id].phase_current .-
                exact_controls[id].phase_current) for id in ids]))
        next_state = _fleet_voltage_state(next_voltage_by_id, ids)
        push!(trajectory, next_state)
        voltage_threshold = Float64(atol_voltage) + Float64(rtol) *
            max(maximum(abs, _fleet_voltage_state(voltage_by_id, ids)),
                maximum(abs, next_state), 1.0)
        current_threshold = Float64(atol_current) + Float64(rtol) *
            max(maximum(abs, reduce(vcat, [collect(exact_controls[id].phase_current)
                for id in ids])),
                maximum(abs, reduce(vcat, [collect(next_controls[id].phase_current)
                for id in ids])), 1.0)
        if residual_voltage <= voltage_threshold &&
                residual_current <= current_threshold
            voltage_by_id = next_voltage_by_id
            exact_controls = next_controls
            converged = true
            break
        end
        history_start = max(1, length(history) - cycle_window + 1)
        for index in history_start:length(history)
            if norm(next_state - history[index]) <= Float64(cycle_tolerance)
                cycled = true
                cycle_period = length(history) - index + 1
                break
            end
        end
        voltage_by_id = next_voltage_by_id
        exact_controls = next_controls
        cycled && break
        push!(history, next_state)
    end
    ControlledInverterFleetNetworkFixedPointResult(
        smooth, converged, cycled, iterations, cycle_period,
        Float64(residual_voltage), Float64(residual_current), exact_controls,
        plants, last_bus, last_status, reduce(hcat, trajectory))
end

"""Deterministic per-device comparison rows for a fleet fixed-point result."""
function controlled_inverter_network_fixed_point_rows(
        result::ControlledInverterFleetNetworkFixedPointResult)
    rows = NamedTuple[]
    for id in sort!(collect(keys(result.smooth.devices)))
        smooth = result.smooth.devices[id]
        exact = get(result.exact_controls, id, nothing)
        plant = get(result.plants, id, nothing)
        push!(rows, (
            device_id=id,
            converged=result.converged,
            cycled=result.cycled,
            cycle_period=result.cycle_period,
            iterations=result.iterations,
            residual_voltage_inf_V=result.residual_voltage_inf_V,
            residual_current_inf_A=result.residual_current_inf_A,
            solve_publishable=result.solve.publishable,
            smooth_exact_current_residual_A=
                smooth.exact_smooth_current_residual,
            exact_smooth_voltage_inf_V=exact === nothing ? NaN :
                maximum(abs, exact.phase_voltage .-
                    smooth.control.phase_voltage),
            exact_smooth_current_inf_A=exact === nothing ? NaN :
                maximum(abs, exact.phase_current .-
                    smooth.control.phase_current),
            smooth_p_poc_W=smooth.plant.p_poc,
            exact_p_poc_W=plant === nothing ? NaN : plant.p_poc,
            smooth_q_poc_var=smooth.plant.q_poc,
            exact_q_poc_var=plant === nothing ? NaN : plant.q_poc,
        ))
    end
    rows
end

"""Deterministic collection of exact fleet equilibrium runs from multiple starts."""
struct ControlledInverterFleetMultiStartResult
    smooth::ControlledInverterFleetResult
    start_ids::Vector{String}
    results::Vector{ControlledInverterFleetNetworkFixedPointResult}
end

function _fleet_result_voltage_vector(
        result::ControlledInverterFleetNetworkFixedPointResult)
    ids = sort!(collect(keys(result.smooth.spec.devices)))
    isempty(result.bus) && return nothing
    try
        reduce(vcat, [
            _phasor_vector(_network_phase_state(
                result.bus, inverter_spec(result.smooth.spec.devices[id])))
            for id in ids])
    catch error
        error isa ArgumentError || rethrow()
        nothing
    end
end

"""
    solve_controlled_inverter_fleet_multistart(net, spec, starts; kwargs...)
        -> ControlledInverterFleetMultiStartResult

Run the simultaneous exact fleet equilibrium oracle from several named initial
voltage states. `starts` is a dictionary from a reproducible string label to a
dictionary keyed by controlled device id, with each value a three-phase phasor
vector or `InverterControlMeasurement`. The smooth fleet solve is constructed
once and reused for every start. This makes initial-condition dependence
visible without conflating repeated smooth solves with equilibrium evidence.
"""
function solve_controlled_inverter_fleet_multistart(
        net::Dict{String,Any}, spec::ControlledInverterFleetSpec,
        starts::AbstractDict;
        smooth_result::Union{Nothing,ControlledInverterFleetResult}=nothing,
        max_iterations::Integer=20,
        atol_voltage::Real=1e-6,
        atol_current::Real=1e-6,
        rtol::Real=1e-8,
        cycle_window::Integer=8,
        cycle_tolerance::Real=1e-7,
        s_base::Real=1e6,
        optimizer=Ipopt.Optimizer,
        verbose::Bool=false,
        selection_objective::Symbol=:loss,
        solver_options=())
    isempty(starts) && throw(ArgumentError(
        "a multi-start equilibrium audit needs at least one start"))
    smooth = smooth_result === nothing ? solve_controlled_inverter_fleet(
        net, spec; s_base=s_base, optimizer=optimizer, verbose=verbose,
        selection_objective=selection_objective, solver_options=solver_options) :
        smooth_result
    labels = String[]
    canonical = Dict{String,Any}()
    for (label, value) in starts
        id = String(label)
        isempty(strip(id)) && throw(ArgumentError(
            "multi-start labels must be non-empty"))
        haskey(canonical, id) && throw(ArgumentError(
            "duplicate multi-start label '$id' after string conversion"))
        canonical[id] = value
        push!(labels, id)
    end
    sort!(labels)
    results = ControlledInverterFleetNetworkFixedPointResult[]
    for label in labels
        initial = canonical[label]
        initial isa AbstractDict || throw(ArgumentError(
            "multi-start '$label' must map device ids to phasors"))
        result = solve_controlled_inverter_fleet_network_fixed_point(
            net, spec; smooth_result=smooth, initial_voltage=initial,
            max_iterations=max_iterations, atol_voltage=atol_voltage,
            atol_current=atol_current, rtol=rtol, cycle_window=cycle_window,
            cycle_tolerance=cycle_tolerance, s_base=s_base, optimizer=optimizer,
            verbose=verbose, selection_objective=selection_objective,
            solver_options=solver_options)
        push!(results, result)
    end
    ControlledInverterFleetMultiStartResult(smooth, labels, results)
end

"""Return one deterministic robustness row per multi-start initial condition."""
function controlled_inverter_fleet_multistart_rows(
        result::ControlledInverterFleetMultiStartResult)
    baseline = nothing
    rows = NamedTuple[]
    for (label, run) in zip(result.start_ids, result.results)
        state = _fleet_result_voltage_vector(run)
        spread = if state === nothing || baseline === nothing
            NaN
        else
            maximum(abs, state - baseline)
        end
        run.converged && baseline === nothing && (baseline = state)
        push!(rows, (
            start_id=label,
            converged=run.converged,
            cycled=run.cycled,
            cycle_period=run.cycle_period,
            iterations=run.iterations,
            publishable=run.solve.publishable,
            residual_voltage_inf_V=run.residual_voltage_inf_V,
            residual_current_inf_A=run.residual_current_inf_A,
            final_equilibrium_voltage_spread_inf_V=spread,
        ))
    end
    rows
end

"""
    inverter_control_loop_gain(controller, measurement, request, ratings,
        voltage_sensitivity; kwargs...) -> FixedPointGainScreen

Screen the local loop ``Δv ↦ Z⋅(∂i/∂v)⋅Δv`` in real rectangular coordinates.
`voltage_sensitivity` is a 6×6 real matrix mapping phase-current perturbations
to phase-voltage perturbations. The exact controller evaluator supplies the
current-command Jacobian. This is a local diagnostic for the supplied
quasi-static sensitivity; it does not include dynamics absent from `Z` or the
controller map.
"""
function inverter_control_loop_gain(
        controller::SequenceController,
        measurement::InverterControlMeasurement,
        request::InverterControlRequest,
        ratings::InverterControlRatings,
        voltage_sensitivity::AbstractMatrix;
        step::Real=1e-6,
        relative_step::Bool=true,
        threshold::Real=1.0)
    sensitivity = _validated_voltage_sensitivity(voltage_sensitivity)
    controller_jacobian = inverter_control_current_jacobian(
        controller, measurement, request, ratings;
        step=step, relative_step=relative_step)
    point = _phasor_vector(measurement.phase_voltage)
    current0 = _phasor_vector(evaluate_exact(
        controller, measurement, request, ratings).phase_current)
    map = state -> point + sensitivity * (_phasor_vector(evaluate_exact(
        controller, InverterControlMeasurement(_phasor_state(state)),
        request, ratings).phase_current) - current0)
    screen = screen_fixed_point_gain(
        map, point; step=step, relative_step=relative_step,
        threshold=threshold)
    # Keep the independently calculated controller Jacobian visible to callers
    # through the identity J_loop = Z * J_controller.
    isapprox(screen.jacobian, sensitivity * controller_jacobian;
             atol=1e-7, rtol=1e-6) ||
        throw(ArgumentError("loop-gain finite-difference inconsistency"))
    screen
end
