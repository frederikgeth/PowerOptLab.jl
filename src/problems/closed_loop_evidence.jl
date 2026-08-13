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

function _real_state(x, name::String)
    values = Float64.(collect(x))
    all(isfinite, values) || throw(ArgumentError("$name must be finite"))
    values
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
        store_trajectory && push!(states, copy(next))
        if residual <= threshold
            x = next
            converged = true
            break
        end
        if cycle_window > 0 && length(states) > 1 &&
                residual > threshold && residual > cycle_tolerance_value
            first_index = max(1, length(states) - 1 - cycle_window)
            for index in first_index:(length(states) - 1)
                if norm(next - states[index]) <= cycle_tolerance_value
                    cycle_period = length(states) - index
                    cycled = true
                    break
                end
            end
        end
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
    sensitivity = Matrix{Float64}(voltage_sensitivity)
    size(sensitivity) == (6, 6) || throw(DimensionMismatch(
        "voltage_sensitivity must be a 6×6 real matrix"))
    all(isfinite, sensitivity) || throw(ArgumentError(
        "voltage_sensitivity must be finite"))
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
