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
