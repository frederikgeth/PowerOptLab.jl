# Local, phase-aware inverter control laws.
#
# The numeric evaluator is the firmware oracle: exact PWL corners, exact
# min/max phase selection, and an algebraic common-scale limiter. The JuMP
# evaluator stamps the same fixed computation graph with BMOPFTools' cached
# softplus PWL expressions and smooth algebraic max/min surrogates. Neither path
# solves a controller-side optimization problem.

"""Supertype for local inverter controllers composed with a physical device."""
abstract type AbstractInverterControlLaw end

"""Supertype for policies that form a positive-sequence power request."""
abstract type AbstractPositiveSequencePolicy end

"""Supertype for policies that form a negative-sequence current request."""
abstract type AbstractUnbalancePolicy end

"""Supertype for algebraic controller-output limiters."""
abstract type AbstractLimiterPolicy end

"""Supertype for the physical current location commanded by a controller."""
abstract type AbstractCurrentTarget end

"""Command the phase-leg current at the converter terminals."""
struct ConverterCurrentTarget <: AbstractCurrentTarget end

"""Command the grid-side current after the output filter."""
struct GridCurrentTarget <: AbstractCurrentTarget end

"""
    PiecewiseLinearLaw(breakpoints, values; smoothing_epsilon)

A continuous, flat-clamped piecewise-linear scalar law. `smoothing_epsilon` is
an absolute width in the same units as `breakpoints`; it is used only by the
smooth JuMP representation. [`evaluate_exact`](@ref) retains exact corners.
"""
struct PiecewiseLinearLaw
    breakpoints::Vector{Float64}
    values::Vector{Float64}
    smoothing_epsilon::Float64
    function PiecewiseLinearLaw(breakpoints::AbstractVector{<:Real},
                                values::AbstractVector{<:Real};
                                smoothing_epsilon::Real)
        xs = Float64.(breakpoints)
        ys = Float64.(values)
        length(xs) == length(ys) || throw(ArgumentError(
            "piecewise-linear breakpoints and values must have equal length"))
        length(xs) >= 2 || throw(ArgumentError(
            "a piecewise-linear law needs at least two points"))
        all(isfinite, xs) || throw(ArgumentError("breakpoints must be finite"))
        all(isfinite, ys) || throw(ArgumentError("values must be finite"))
        all(diff(xs) .> 0) || throw(ArgumentError(
            "breakpoints must be strictly increasing"))
        eps = Float64(smoothing_epsilon)
        isfinite(eps) && eps > 0 || throw(ArgumentError(
            "smoothing_epsilon must be finite and > 0"))
        new(xs, ys, eps)
    end
end

"""
    WorstPhaseVoltVarWatt(; volt_watt=nothing, volt_var=nothing,
                           conflict_policy=:dominant, extrema_epsilon=0.05,
                           conflict_epsilon=0.01,
                           volt_watt_basis=:available)

Direction-aware positive-sequence policy for three-phase local voltage control.
Volt-watt observes the largest phase-voltage magnitude. Volt-var observes both
the smallest and largest phase: low-voltage injection and high-voltage
absorption are retained. When both occur, `conflict_policy=:net` adds the two
opposing requests, `:dominant` continuously blends toward the branch with the
larger normalized response, `:low_voltage` prioritizes injection, and
`:high_voltage` prioritizes absorption.

Volt-watt ordinates are fractions in `[0, 1]`; Volt-var ordinates are fractions
of `InverterControlRequest.q_scale`. Both curves must be non-increasing.
`extrema_epsilon` is an SI voltage width for smooth phase extrema;
`conflict_epsilon` is the dimensionless transition width in normalized branch
severity and is part of both the firmware and JuMP laws. `volt_watt_basis` is
`:available` or `:rated`.
"""
struct WorstPhaseVoltVarWatt <: AbstractPositiveSequencePolicy
    volt_watt::Union{Nothing,PiecewiseLinearLaw}
    volt_var::Union{Nothing,PiecewiseLinearLaw}
    conflict_policy::Symbol
    extrema_epsilon::Float64
    conflict_epsilon::Float64
    volt_watt_basis::Symbol
end

function WorstPhaseVoltVarWatt(; volt_watt=nothing, volt_var=nothing,
                               conflict_policy::Symbol=:dominant,
                               extrema_epsilon::Real=0.05,
                               conflict_epsilon::Real=0.01,
                               volt_watt_basis::Symbol=:available)
    volt_watt isa Union{Nothing,PiecewiseLinearLaw} || throw(ArgumentError(
        "volt_watt must be a PiecewiseLinearLaw or nothing"))
    volt_var isa Union{Nothing,PiecewiseLinearLaw} || throw(ArgumentError(
        "volt_var must be a PiecewiseLinearLaw or nothing"))
    conflict_policy in (:dominant, :net, :low_voltage, :high_voltage) ||
        throw(ArgumentError("conflict_policy must be :dominant, :net, " *
                            ":low_voltage, or :high_voltage"))
    volt_watt_basis in (:available, :rated) || throw(ArgumentError(
        "volt_watt_basis must be :available or :rated"))
    veps = Float64(extrema_epsilon)
    qeps = Float64(conflict_epsilon)
    isfinite(veps) && veps > 0 || throw(ArgumentError(
        "extrema_epsilon must be finite and > 0"))
    isfinite(qeps) && qeps > 0 || throw(ArgumentError(
        "conflict_epsilon must be finite and > 0"))
    for (name, law) in (("volt_watt", volt_watt), ("volt_var", volt_var))
        law === nothing && continue
        all(diff(law.values) .<= 0) || throw(ArgumentError(
            "$name values must be non-increasing for worst-phase selection"))
    end
    if volt_watt !== nothing
        all(v -> 0 <= v <= 1, volt_watt.values) || throw(ArgumentError(
            "volt_watt values must lie in [0, 1]"))
    end
    return WorstPhaseVoltVarWatt(
        volt_watt, volt_var, conflict_policy, veps, qeps, volt_watt_basis)
end

"""
    AverageVoltageVoltVarWatt(; volt_watt=nothing, volt_var=nothing,
                                extrema_epsilon=0.05,
                                volt_watt_basis=:available)

Legacy balanced-current comparator. Both PWL laws observe the arithmetic mean
of the three local phase-voltage magnitudes.
"""
struct AverageVoltageVoltVarWatt <: AbstractPositiveSequencePolicy
    volt_watt::Union{Nothing,PiecewiseLinearLaw}
    volt_var::Union{Nothing,PiecewiseLinearLaw}
    extrema_epsilon::Float64
    volt_watt_basis::Symbol
end

"""
    PositiveSequenceVoltVarWatt(; volt_watt=nothing, volt_var=nothing,
                                  worst_phase_watt_guard=true,
                                  guard_epsilon=0.05,
                                  volt_watt_basis=:available)

Balanced-current sequence comparator. Volt-var observes `|U1|`. Volt-watt
observes the largest phase magnitude when `worst_phase_watt_guard=true`, so an
overvoltage phase cannot be hidden by the positive-sequence magnitude.
"""
struct PositiveSequenceVoltVarWatt <: AbstractPositiveSequencePolicy
    volt_watt::Union{Nothing,PiecewiseLinearLaw}
    volt_var::Union{Nothing,PiecewiseLinearLaw}
    worst_phase_watt_guard::Bool
    guard_epsilon::Float64
    volt_watt_basis::Symbol
end

function _validate_positive_curves(volt_watt, volt_var)
    volt_watt isa Union{Nothing,PiecewiseLinearLaw} || throw(ArgumentError(
        "volt_watt must be a PiecewiseLinearLaw or nothing"))
    volt_var isa Union{Nothing,PiecewiseLinearLaw} || throw(ArgumentError(
        "volt_var must be a PiecewiseLinearLaw or nothing"))
    for (name, law) in (("volt_watt", volt_watt), ("volt_var", volt_var))
        law === nothing && continue
        all(diff(law.values) .<= 0) || throw(ArgumentError(
            "$name values must be non-increasing"))
    end
    if volt_watt !== nothing
        all(v -> 0 <= v <= 1, volt_watt.values) || throw(ArgumentError(
            "volt_watt values must lie in [0, 1]"))
    end
    return nothing
end

function _validate_volt_watt_basis(basis::Symbol)
    basis in (:available, :rated) || throw(ArgumentError(
        "volt_watt_basis must be :available or :rated"))
    return basis
end

function AverageVoltageVoltVarWatt(; volt_watt=nothing, volt_var=nothing,
        extrema_epsilon::Real=0.05,
        volt_watt_basis::Symbol=:available)
    _validate_positive_curves(volt_watt, volt_var)
    eps = Float64(extrema_epsilon)
    isfinite(eps) && eps > 0 || throw(ArgumentError(
        "extrema_epsilon must be finite and > 0"))
    AverageVoltageVoltVarWatt(
        volt_watt, volt_var, eps, _validate_volt_watt_basis(volt_watt_basis))
end

function PositiveSequenceVoltVarWatt(; volt_watt=nothing, volt_var=nothing,
        worst_phase_watt_guard::Bool=true, guard_epsilon::Real=0.05,
        volt_watt_basis::Symbol=:available)
    _validate_positive_curves(volt_watt, volt_var)
    eps = Float64(guard_epsilon)
    isfinite(eps) && eps > 0 || throw(ArgumentError(
        "guard_epsilon must be finite and > 0"))
    PositiveSequenceVoltVarWatt(
        volt_watt, volt_var, worst_phase_watt_guard, eps,
        _validate_volt_watt_basis(volt_watt_basis))
end

"""Balanced-current baseline: request no negative-sequence current."""
struct NoUnbalanceControl <: AbstractUnbalancePolicy end

"""
    NegativeSequenceAdmittanceDroop(gain; impedance_angle=0,
        ripple_blend=0, voltage_floor=1)

Request negative-sequence current from local voltage only. `gain` maps
`eta = |U2| / sqrt(|U1|^2 + voltage_floor^2)` to an admittance magnitude in
A/V. `impedance_angle` fixes the configured negative-sequence actuation angle,
and `ripple_blend` blends the voltage-oriented request with the analytical
double-frequency-ripple-cancelling target.
"""
struct NegativeSequenceAdmittanceDroop <: AbstractUnbalancePolicy
    gain::PiecewiseLinearLaw
    impedance_angle::Float64
    ripple_blend::Float64
    voltage_floor::Float64
end

function NegativeSequenceAdmittanceDroop(gain::PiecewiseLinearLaw;
        impedance_angle::Real=0.0, ripple_blend::Real=0.0,
        voltage_floor::Real=1.0)
    all(gain.values .>= 0) || throw(ArgumentError(
        "negative-sequence admittance gains must be non-negative"))
    phi = Float64(impedance_angle)
    blend = Float64(ripple_blend)
    floor = Float64(voltage_floor)
    isfinite(phi) || throw(ArgumentError("impedance_angle must be finite"))
    isfinite(blend) && 0 <= blend <= 1 || throw(ArgumentError(
        "ripple_blend must lie in [0, 1]"))
    isfinite(floor) && floor > 0 || throw(ArgumentError(
        "voltage_floor must be finite and > 0"))
    NegativeSequenceAdmittanceDroop(gain, phi, blend, floor)
end

"""
    CommonScaleLimiter(; current_epsilon=1e-3, power_epsilon=1e-3,
                         pq_priority=:proportional)

Apply a common scalar to positive- and negative-sequence current commands so
the reconstructed currents at the configured current target respect its
declared per-leg limit. The physical plant independently retains every
converter- and grid-side current constraint.
The positive-sequence `P,Q` request is first scaled to the converter apparent-
power rating. The exact evaluator uses hard maxima. The smooth model represents
every magnitude by an implicit nonnegative square root and uses the stated SI
epsilons only in smooth maximum selectors. `pq_priority` is `:proportional`,
`:watt`, or `:var` for the positive-sequence apparent-power allocation.
"""
struct CommonScaleLimiter <: AbstractLimiterPolicy
    current_epsilon::Float64
    power_epsilon::Float64
    pq_priority::Symbol
end

function CommonScaleLimiter(; current_epsilon::Real=1e-3,
                            power_epsilon::Real=1e-3,
                            pq_priority::Symbol=:proportional)
    ieps = Float64(current_epsilon)
    seps = Float64(power_epsilon)
    isfinite(ieps) && ieps > 0 || throw(ArgumentError(
        "current_epsilon must be finite and > 0"))
    isfinite(seps) && seps > 0 || throw(ArgumentError(
        "power_epsilon must be finite and > 0"))
    pq_priority in (:proportional, :watt, :var) || throw(ArgumentError(
        "pq_priority must be :proportional, :watt, or :var"))
    CommonScaleLimiter(ieps, seps, pq_priority)
end

"""
    SequenceController(positive, unbalance, limiter;
                       current_target=ConverterCurrentTarget(),
                       power_voltage_floor=1)

Closed-form three-leg controller combining a positive-sequence Volt-var/Watt
policy, a negative-sequence policy, and an algebraic feasibility limiter. The
curves always observe local phase-to-neutral POC voltage. `current_target`
selects converter-leg or post-filter grid current. `power_voltage_floor` is the
strictly positive SI voltage floor in the regularized `P,Q`-to-`I1` conversion;
it is not a magnitude-smoothing epsilon.
"""
struct SequenceController{P<:AbstractPositiveSequencePolicy,
                          U<:AbstractUnbalancePolicy,
                          L<:AbstractLimiterPolicy,
                          T<:AbstractCurrentTarget} <: AbstractInverterControlLaw
    positive::P
    unbalance::U
    limiter::L
    current_target::T
    power_voltage_floor::Float64
end

_current_target(target::Union{ConverterCurrentTarget,GridCurrentTarget}) = target
_current_target(::AbstractCurrentTarget) = throw(ArgumentError(
    "only ConverterCurrentTarget and GridCurrentTarget are currently supported"))
_current_target(target::Symbol) = target === :converter ? ConverterCurrentTarget() :
    target === :grid ? GridCurrentTarget() : throw(ArgumentError(
        "current_target must be :converter, :grid, or an AbstractCurrentTarget"))

function SequenceController(
        positive::AbstractPositiveSequencePolicy,
        unbalance::AbstractUnbalancePolicy,
        limiter::AbstractLimiterPolicy;
        current_target=ConverterCurrentTarget(),
        power_voltage_floor::Real=1.0)
    floor = Float64(power_voltage_floor)
    isfinite(floor) && floor > 0 || throw(ArgumentError(
        "power_voltage_floor must be finite and > 0"))
    SequenceController(
        positive, unbalance, limiter, _current_target(current_target), floor)
end

SequenceController(positive::AbstractPositiveSequencePolicy;
                   unbalance::AbstractUnbalancePolicy=NoUnbalanceControl(),
                   limiter::AbstractLimiterPolicy=CommonScaleLimiter(),
                   current_target=ConverterCurrentTarget(),
                   power_voltage_floor::Real=1.0) =
    SequenceController(
        positive, unbalance, limiter;
        current_target=current_target,
        power_voltage_floor=power_voltage_floor)

"""Local three-phase RMS voltage phasors in the declared a-b-c terminal order."""
struct InverterControlMeasurement
    phase_voltage::NTuple{3,ComplexF64}
    function InverterControlMeasurement(
            phase_voltage::NTuple{3,<:Complex})
        values = ComplexF64.(phase_voltage)
        all(v -> isfinite(real(v)) && isfinite(imag(v)), values) ||
            throw(ArgumentError("phase-voltage phasors must be finite"))
        new(values)
    end
end

function InverterControlMeasurement(phase_voltage::AbstractVector{<:Complex})
    length(phase_voltage) == 3 || throw(ArgumentError(
        "inverter control measurement requires exactly three phase phasors"))
    InverterControlMeasurement(Tuple(phase_voltage))
end

"""Per-operating-point active availability, rated active power, and reactive-power base, in SI."""
struct InverterControlRequest
    p_available::Float64
    p_rated::Float64
    q_scale::Float64
    function InverterControlRequest(p_available::Real, p_rated::Real,
                                    q_scale::Real)
        p = Float64(p_available)
        pr = Float64(p_rated)
        q = Float64(q_scale)
        isfinite(p) && p >= 0 || throw(ArgumentError(
            "p_available must be finite and non-negative"))
        isfinite(pr) && pr > 0 || throw(ArgumentError(
            "p_rated must be finite and positive"))
        isfinite(q) && q >= 0 || throw(ArgumentError(
            "q_scale must be finite and non-negative"))
        new(p, pr, q)
    end
end

InverterControlRequest(p_available::Real, q_scale::Real) =
    InverterControlRequest(p_available, max(Float64(p_available), eps()), q_scale)

function InverterControlRequest(; p_available::Real,
        p_rated::Real=max(Float64(p_available), eps()), q_scale::Real)
    InverterControlRequest(p_available, p_rated, q_scale)
end

"""Fundamental apparent-power and per-phase conductor-current ratings, in SI."""
struct InverterControlRatings
    s_max::Float64
    i_max::Float64
    function InverterControlRatings(s_max::Real, i_max::Real)
        s = Float64(s_max)
        i = Float64(i_max)
        isfinite(s) && s > 0 || throw(ArgumentError(
            "s_max must be finite and > 0"))
        isfinite(i) && i > 0 || throw(ArgumentError(
            "i_max must be finite and > 0"))
        new(s, i)
    end
end

function InverterControlRatings(; s_max::Real, i_max::Real)
    InverterControlRatings(s_max, i_max)
end

function InverterControlRatings(inverter::AdvancedInverter)
    inverter.i_max === nothing && throw(ArgumentError(
        "a controlled AdvancedInverter requires a finite i_max"))
    InverterControlRatings(s_max=inverter.s_max, i_max=inverter.i_max)
end

InverterControlRatings(
    inverter::AdvancedInverter, ::ConverterCurrentTarget) =
    InverterControlRatings(inverter)

function InverterControlRatings(
        inverter::AdvancedInverter, ::GridCurrentTarget)
    converter_limit = inverter.i_max
    converter_limit === nothing && throw(ArgumentError(
        "a grid-current-controlled AdvancedInverter requires i_max"))
    grid_limit = something(inverter.i_grid_max, converter_limit)
    grid_limit === nothing && throw(ArgumentError(
        "a grid-current-controlled AdvancedInverter requires i_grid_max or i_max"))
    InverterControlRatings(
        s_max=inverter.s_max, i_max=min(converter_limit, grid_limit))
end

"""
Numeric controller record shared by exact evaluation and smooth-model
extraction. `sequence_power`, `total_power`, and `ripple_power` are command-space
diagnostics formed from the POC voltage measurement and commanded current at the
configured current target. They are not converter-terminal powers when an
output filter separates those locations. Complex powers use the injection sign
convention.
"""
struct InverterControlResult
    phase_voltage::NTuple{3,ComplexF64}
    voltage_sequence::NTuple{3,ComplexF64}
    voltage_min::Float64
    voltage_max::Float64
    eta::Float64
    p_request::Float64
    q_request::Float64
    power_scale::Float64
    i1_request::ComplexF64
    i2_voltage::ComplexF64
    i2_ripple::ComplexF64
    i2_request::ComplexF64
    current_scale::Float64
    i1::ComplexF64
    i2::ComplexF64
    phase_current::NTuple{3,ComplexF64}
    sequence_power::ComplexF64
    total_power::ComplexF64
    ripple_power::ComplexF64
end

"""
Converter-terminal phasors and powers derived from the physical plant solution.
All phasors are RMS and powers follow the injection sign convention. The
authoritative total converter power is `sum(phase_power)`. Sequence powers use
`S_k = 3U_k*conj(I_k)`, and `ripple_power` is the complex twice-fundamental
coefficient `sum(U_phase*I_phase)`.
"""
struct ConverterTerminalResult
    phase_voltage::NTuple{3,ComplexF64}
    phase_current::NTuple{3,ComplexF64}
    voltage_sequence::NTuple{3,ComplexF64}
    current_sequence::NTuple{3,ComplexF64}
    phase_power::NTuple{3,ComplexF64}
    zero_sequence_power::ComplexF64
    positive_sequence_power::ComplexF64
    negative_sequence_power::ComplexF64
    total_power::ComplexF64
    ripple_power::ComplexF64
end

function _converter_terminal_result(
        phase_voltage::NTuple{3,ComplexF64},
        phase_current::NTuple{3,ComplexF64})
    voltage_sequence = _phasor_sequences(phase_voltage)
    current_sequence = _phasor_sequences(phase_current)
    phase_power = ntuple(
        k -> phase_voltage[k]*conj(phase_current[k]), 3)
    sequence_power = ntuple(
        k -> 3voltage_sequence[k]*conj(current_sequence[k]), 3)
    total_power = sum(phase_power)
    ripple_power = sum(
        phase_voltage[k]*phase_current[k] for k in 1:3)
    return ConverterTerminalResult(
        phase_voltage, phase_current, voltage_sequence, current_sequence,
        phase_power, sequence_power[1], sequence_power[2], sequence_power[3],
        total_power, ripple_power)
end

_pwl_exact(law::PiecewiseLinearLaw, input::Real) =
    BMOPFTools.piecewise_linear_value(input, law.breakpoints, law.values)

function _phasor_sequences(values::NTuple{3,ComplexF64})
    seq = _sequence_components(real.(values), imag.(values))
    return ntuple(k -> ComplexF64(seq[k][1], seq[k][2]), 3)
end

_volt_watt_basis(policy::AbstractPositiveSequencePolicy) =
    policy.volt_watt_basis

function _volt_watt_request(policy::AbstractPositiveSequencePolicy,
                            request::InverterControlRequest, fraction)
    return _volt_watt_basis(policy) === :available ?
        request.p_available*fraction :
        request.p_rated*fraction
end

function _phase_from_sequences(i1, i2)
    a = cis(2pi / 3)
    return (i1 + i2, a^2*i1 + a*i2, a*i1 + a^2*i2)
end

function _volt_var_branch_scales(policy::WorstPhaseVoltVarWatt)
    values = policy.volt_var.values
    return max(maximum(values), eps()), max(-minimum(values), eps())
end

function _dominant_blend_numeric(low, high, policy::WorstPhaseVoltVarWatt)
    positive_scale, negative_scale = _volt_var_branch_scales(policy)
    severity_difference = low/positive_scale + high/negative_scale
    epsilon = policy.conflict_epsilon
    weight = (1 + severity_difference /
              sqrt(severity_difference^2 + epsilon^2))/2
    return weight*low + (1-weight)*high
end

function _volt_var_conflict_exact(policy::WorstPhaseVoltVarWatt, low, high)
    qlo = max(low, 0.0)
    qhi = min(high, 0.0)
    return policy.conflict_policy === :low_voltage ? qlo :
           policy.conflict_policy === :high_voltage ? qhi :
           policy.conflict_policy === :dominant ?
               _dominant_blend_numeric(qlo, qhi, policy) : qlo + qhi
end

function _worst_phase_exact(policy::WorstPhaseVoltVarWatt,
                            magnitudes, request)
    vmin, vmax = extrema(magnitudes)
    p_fraction = policy.volt_watt === nothing ? 1.0 :
        _pwl_exact(policy.volt_watt, vmax)
    q_fraction = if policy.volt_var === nothing
        0.0
    else
        _volt_var_conflict_exact(
            policy, _pwl_exact(policy.volt_var, vmin),
            _pwl_exact(policy.volt_var, vmax))
    end
    return _volt_watt_request(policy, request, p_fraction),
           request.q_scale*q_fraction,
           vmin, vmax
end

function _average_voltage_exact(policy::AverageVoltageVoltVarWatt,
                                magnitudes, request)
    voltage = sum(magnitudes)/3
    p_fraction = policy.volt_watt === nothing ? 1.0 :
        _pwl_exact(policy.volt_watt, voltage)
    q_fraction = policy.volt_var === nothing ? 0.0 :
        _pwl_exact(policy.volt_var, voltage)
    return _volt_watt_request(policy, request, p_fraction),
           request.q_scale*q_fraction,
           minimum(magnitudes), maximum(magnitudes)
end

function _positive_sequence_exact(policy::PositiveSequenceVoltVarWatt,
                                  u1, magnitudes, request)
    voltage = abs(u1)
    watt_voltage = policy.worst_phase_watt_guard ? maximum(magnitudes) : voltage
    p_fraction = policy.volt_watt === nothing ? 1.0 :
        _pwl_exact(policy.volt_watt, watt_voltage)
    q_fraction = policy.volt_var === nothing ? 0.0 :
        _pwl_exact(policy.volt_var, voltage)
    return _volt_watt_request(policy, request, p_fraction),
           request.q_scale*q_fraction,
           minimum(magnitudes), maximum(magnitudes)
end

_positive_exact(policy::WorstPhaseVoltVarWatt, u1, magnitudes, request) =
    _worst_phase_exact(policy, magnitudes, request)
_positive_exact(policy::AverageVoltageVoltVarWatt, u1, magnitudes, request) =
    _average_voltage_exact(policy, magnitudes, request)
_positive_exact(policy::PositiveSequenceVoltVarWatt, u1, magnitudes, request) =
    _positive_sequence_exact(policy, u1, magnitudes, request)

function _unbalance_exact(::NoUnbalanceControl, u1, u2, i1)
    return 0.0im, 0.0im, 0.0im, 0.0
end

_pwl_numeric_smooth(law::PiecewiseLinearLaw, input::Real) =
    BMOPFTools.piecewise_linear_value(
        input, law.breakpoints, law.values;
        epsilon=law.smoothing_epsilon)

function _numeric_smooth_extrema(magnitudes, epsilon)
    vmin = _smooth_min(_smooth_min(magnitudes[1], magnitudes[2], epsilon),
                       magnitudes[3], epsilon)
    vmax = _smooth_max(_smooth_max(magnitudes[1], magnitudes[2], epsilon),
                       magnitudes[3], epsilon)
    return vmin, vmax
end

function _worst_phase_numeric_smooth(policy::WorstPhaseVoltVarWatt,
                                     magnitudes, request)
    eps_v = policy.extrema_epsilon
    vmin, vmax = _numeric_smooth_extrema(magnitudes, eps_v)
    p_fraction = policy.volt_watt === nothing ? 1.0 :
        _pwl_numeric_smooth(policy.volt_watt, vmax)
    q_fraction = if policy.volt_var === nothing
        0.0
    else
        qlo = _smooth_positive(
            _pwl_numeric_smooth(policy.volt_var, vmin),
            policy.conflict_epsilon)
        qhi = _smooth_negative(
            _pwl_numeric_smooth(policy.volt_var, vmax),
            policy.conflict_epsilon)
        policy.conflict_policy === :low_voltage ? qlo :
        policy.conflict_policy === :high_voltage ? qhi :
        policy.conflict_policy === :dominant ?
            _dominant_blend_numeric(qlo, qhi, policy) : qlo + qhi
    end
    return _volt_watt_request(policy, request, p_fraction),
           request.q_scale*q_fraction,
           vmin, vmax
end

function _average_voltage_numeric_smooth(
        policy::AverageVoltageVoltVarWatt, magnitudes, request)
    voltage = sum(magnitudes)/3
    vmin, vmax = _numeric_smooth_extrema(
        magnitudes, policy.extrema_epsilon)
    p_fraction = policy.volt_watt === nothing ? 1.0 :
        _pwl_numeric_smooth(policy.volt_watt, voltage)
    q_fraction = policy.volt_var === nothing ? 0.0 :
        _pwl_numeric_smooth(policy.volt_var, voltage)
    return _volt_watt_request(policy, request, p_fraction),
           request.q_scale*q_fraction,
           vmin, vmax
end

function _positive_sequence_numeric_smooth(
        policy::PositiveSequenceVoltVarWatt, u1, magnitudes, request)
    voltage = abs(u1)
    vmin, vmax = _numeric_smooth_extrema(magnitudes, policy.guard_epsilon)
    watt_voltage = policy.worst_phase_watt_guard ? vmax : voltage
    p_fraction = policy.volt_watt === nothing ? 1.0 :
        _pwl_numeric_smooth(policy.volt_watt, watt_voltage)
    q_fraction = policy.volt_var === nothing ? 0.0 :
        _pwl_numeric_smooth(policy.volt_var, voltage)
    return _volt_watt_request(policy, request, p_fraction),
           request.q_scale*q_fraction,
           vmin, vmax
end

_positive_numeric_smooth(
    policy::WorstPhaseVoltVarWatt, u1, magnitudes, request) =
    _worst_phase_numeric_smooth(policy, magnitudes, request)
_positive_numeric_smooth(
    policy::AverageVoltageVoltVarWatt, u1, magnitudes, request) =
    _average_voltage_numeric_smooth(policy, magnitudes, request)
_positive_numeric_smooth(
    policy::PositiveSequenceVoltVarWatt, u1, magnitudes, request) =
    _positive_sequence_numeric_smooth(policy, u1, magnitudes, request)

function _unbalance_numeric_smooth(::NoUnbalanceControl, u1, u2, i1)
    return 0.0im, 0.0im, 0.0im, 0.0
end


function _unbalance_numeric_smooth(policy::NegativeSequenceAdmittanceDroop,
                                   u1, u2, i1)
    denominator = abs2(u1) + policy.voltage_floor^2
    eta = abs(u2) / sqrt(denominator)
    gain = _pwl_numeric_smooth(policy.gain, eta)
    i2v = -gain * cis(-policy.impedance_angle) * u2
    i2r = -(u2 * conj(u1) / denominator) * i1
    blend = policy.ripple_blend
    return i2v, i2r, (1-blend)*i2v + blend*i2r, eta
end

function _limit_positive_power_exact(p, q, smax, priority)
    raw_magnitude = hypot(p, q)
    if priority === :proportional
        scale = raw_magnitude > smax ? smax/raw_magnitude : 1.0
        return scale*p, scale*q, scale
    elseif priority === :watt
        p_limited = min(p, smax)
        q_capacity = sqrt(max(smax^2-p_limited^2, 0.0))
        q_limited = clamp(q, -q_capacity, q_capacity)
    else
        q_limited = clamp(q, -smax, smax)
        p_capacity = sqrt(max(smax^2-q_limited^2, 0.0))
        p_limited = min(p, p_capacity)
    end
    scale = iszero(raw_magnitude) ? 1.0 :
        hypot(p_limited, q_limited)/raw_magnitude
    return p_limited, q_limited, scale
end

function _limit_positive_power_numeric_smooth(p, q, smax, epsilon, priority)
    raw_magnitude = hypot(p, q)
    if priority === :proportional
        scale = smax/_smooth_max(smax, raw_magnitude, epsilon)
        return scale*p, scale*q, scale
    elseif priority === :watt
        p_limited = _smooth_min(p, smax, epsilon)
        q_capacity = sqrt(max(smax^2-p_limited^2, 0.0))
        q_scale = q_capacity/_smooth_max(q_capacity, abs(q), epsilon)
        q_limited = q_scale*q
    else
        q_scale = smax/_smooth_max(smax, abs(q), epsilon)
        q_limited = q_scale*q
        p_capacity = sqrt(max(smax^2-q_limited^2, 0.0))
        p_scale = p_capacity/_smooth_max(p_capacity, p, epsilon)
        p_limited = p_scale*p
    end
    scale = iszero(raw_magnitude) ? 1.0 :
        hypot(p_limited, q_limited)/raw_magnitude
    return p_limited, q_limited, scale
end

"""
    evaluate_smooth(controller, measurement, request, ratings)

Numerically evaluate the same differentiable computation graph stamped by
[`stamp_smooth_control!`](@ref). This is intended for inexpensive extraction,
SI/per-unit regression, and comparison with [`evaluate_exact`](@ref); it does
not solve an optimization problem.
"""
function evaluate_smooth(controller::SequenceController,
                         measurement::InverterControlMeasurement,
                         request::InverterControlRequest,
                         ratings::InverterControlRatings)
    phase_voltage = measurement.phase_voltage
    u0, u1, u2 = _phasor_sequences(phase_voltage)
    magnitudes = abs.(phase_voltage)
    p_raw, q_raw, vmin, vmax = _positive_numeric_smooth(
        controller.positive, u1, magnitudes, request)
    if _volt_watt_basis(controller.positive) === :rated
        p_raw = _smooth_min(
            request.p_available, p_raw, controller.limiter.power_epsilon)
    end

    p, q, power_scale = _limit_positive_power_numeric_smooth(
        p_raw, q_raw, ratings.s_max, controller.limiter.power_epsilon,
        controller.limiter.pq_priority)

    floor = controller.power_voltage_floor
    denominator = 3(abs2(u1) + floor^2)
    i1_request = conj(complex(p, q)) * u1 / denominator
    i2v, i2r, i2_request, eta = _unbalance_numeric_smooth(
        controller.unbalance, u1, u2, i1_request)
    phase_request = _phase_from_sequences(i1_request, i2_request)

    ieps = controller.limiter.current_epsilon
    imags = abs.(phase_request)
    current_denominator = foldl(
        (a, b) -> _smooth_max(a, b, ieps), imags; init=ratings.i_max)
    current_scale = ratings.i_max / current_denominator
    i1 = current_scale*i1_request
    i2 = current_scale*i2_request
    phase_current = _phase_from_sequences(i1, i2)
    sequence_power = 3u1*conj(i1)
    total_power = 3(u1*conj(i1) + u2*conj(i2))
    ripple = 3(u1*i2 + u2*i1)
    return InverterControlResult(
        phase_voltage, (u0, u1, u2), vmin, vmax, eta,
        p, q, power_scale, i1_request, i2v, i2r, i2_request,
        current_scale, i1, i2, phase_current, sequence_power,
        total_power, ripple)
end

function _unbalance_exact(policy::NegativeSequenceAdmittanceDroop,
                          u1, u2, i1)
    denom = abs2(u1) + policy.voltage_floor^2
    eta = abs(u2) / sqrt(denom)
    gain = _pwl_exact(policy.gain, eta)
    i2v = -gain * cis(-policy.impedance_angle) * u2
    i2r = -(u2 * conj(u1) / denom) * i1
    blend = policy.ripple_blend
    return i2v, i2r, (1-blend)*i2v + blend*i2r, eta
end

"""
    evaluate_exact(controller, measurement, request, ratings)

Evaluate the exact local controller without JuMP. This is the firmware oracle
for regression tests and exact-versus-smooth residuals. It uses only local RMS
phasors and fixed controller/rating data.
"""
function evaluate_exact(controller::SequenceController,
                        measurement::InverterControlMeasurement,
                        request::InverterControlRequest,
                        ratings::InverterControlRatings)
    phase_voltage = measurement.phase_voltage
    u0, u1, u2 = _phasor_sequences(phase_voltage)
    magnitudes = abs.(phase_voltage)
    p_raw, q_raw, vmin, vmax = _positive_exact(
        controller.positive, u1, magnitudes, request)
    if _volt_watt_basis(controller.positive) === :rated
        p_raw = min(request.p_available, p_raw)
    end

    p, q, power_scale = _limit_positive_power_exact(
        p_raw, q_raw, ratings.s_max, controller.limiter.pq_priority)

    floor = controller.power_voltage_floor
    denominator = 3(abs2(u1) + floor^2)
    i1_request = conj(complex(p, q)) * u1 / denominator
    i2v, i2r, i2_request, eta = _unbalance_exact(
        controller.unbalance, u1, u2, i1_request)
    phase_request = _phase_from_sequences(i1_request, i2_request)
    peak = maximum(abs, phase_request)
    current_scale = peak > ratings.i_max ? ratings.i_max / peak : 1.0
    i1 = current_scale*i1_request
    i2 = current_scale*i2_request
    phase_current = _phase_from_sequences(i1, i2)
    sequence_power = 3u1*conj(i1)
    total_power = 3(u1*conj(i1) + u2*conj(i2))
    ripple = 3(u1*i2 + u2*i1)
    return InverterControlResult(
        phase_voltage, (u0, u1, u2), vmin, vmax, eta,
        p, q, power_scale, i1_request, i2v, i2r, i2_request,
        current_scale, i1, i2, phase_current, sequence_power,
        total_power, ripple)
end

"""
    ControlledDevice(device, controller)

Typed composition of a physical inverter plant and a local control law. The
first supported plant is [`AdvancedInverter`](@ref).
"""
struct ControlledDevice{D<:AbstractDevice,C<:AbstractInverterControlLaw} <:
       AbstractDevice
    device::D
    controller::C
end

device_id(device::ControlledDevice) = device_id(device.device)
"""Return the physical [`AdvancedInverter`](@ref) owned by a control wrapper."""
inverter_spec(device::ControlledDevice{<:AdvancedInverter}) = device.device

"""Return stamped advanced-inverter handles through a control wrapper."""
inverter_handles(::ControlledDevice{<:AdvancedInverter},
                 handles::_InvHandles) = handles

function validate_device(controlled::ControlledDevice{<:AdvancedInverter}, nets=();
                         periods::Int=length(nets))
    inv = controlled.device
    validate_device(inv, nets; periods=periods)
    inv.topology == :THREE_LEG || throw(ArgumentError(
        "the initial SequenceController implementation requires topology=:THREE_LEG"))
    length(inv.phase_terminals) == 3 || throw(ArgumentError(
        "a SequenceController requires exactly three ordered phase terminals"))
    InverterControlRatings(inv)
    return controlled
end

_smooth_max(a, b, epsilon) = (a + b + sqrt((a-b)^2 + epsilon^2)) / 2
_smooth_min(a, b, epsilon) = (a + b - sqrt((a-b)^2 + epsilon^2)) / 2
_smooth_positive(a, epsilon) = (a + sqrt(a^2 + epsilon^2)) / 2
_smooth_negative(a, epsilon) = (a - sqrt(a^2 + epsilon^2)) / 2

function _implicit_sqrt!(m, radicand; start::Real=1.0, scale::Real=start)
    initial = max(Float64(start), 0.0)
    normalization = max(abs(Float64(scale)), eps(Float64))
    root = JuMP.@variable(m, lower_bound=0.0, start=initial)
    JuMP.@constraint(
        m, (root/normalization)^2 == radicand/normalization^2)
    return root
end

function _smooth_max_implicit!(m, a, b, epsilon; start::Real=1.0)
    root = _implicit_sqrt!(
        m, (a-b)^2 + epsilon^2; start=max(Float64(start), epsilon))
    return (a + b + root) / 2
end

function _smooth_min_implicit!(m, a, b, epsilon; start::Real=1.0)
    root = _implicit_sqrt!(
        m, (a-b)^2 + epsilon^2; start=max(Float64(start), epsilon))
    return (a + b - root) / 2
end

function _smooth_positive_implicit!(m, a, epsilon; start::Real=epsilon)
    root = _implicit_sqrt!(
        m, a^2 + epsilon^2; start=max(Float64(start), epsilon))
    return (a + root) / 2
end

function _smooth_negative_implicit!(m, a, epsilon; start::Real=epsilon)
    root = _implicit_sqrt!(
        m, a^2 + epsilon^2; start=max(Float64(start), epsilon))
    return (a - root) / 2
end

function _smooth_dominant_implicit!(
        m, low, high, policy::WorstPhaseVoltVarWatt)
    positive_scale, negative_scale = _volt_var_branch_scales(policy)
    difference = low/positive_scale + high/negative_scale
    epsilon = policy.conflict_epsilon
    root = _implicit_sqrt!(
        m, difference^2 + epsilon^2; start=epsilon, scale=1.0)
    weight = (1 + difference/root)/2
    return weight*low + (1-weight)*high
end

function _pwl_smooth(ctx, law::PiecewiseLinearLaw, input, input_scale)
    BMOPFTools.opf_piecewise_linear_expression(
        ctx, input, law.breakpoints ./ input_scale, law.values;
        epsilon=law.smoothing_epsilon / input_scale)
end

function _sequence_pair(values_re, values_im)
    return _sequence_components(values_re, values_im)
end

_cmul(ar, ai, br, bi) = (ar*br - ai*bi, ar*bi + ai*br)

function _phase_pairs(i1, i2)
    i1r, i1i = i1; i2r, i2i = i2
    h = _HALF_SQRT3
    return ((i1r+i2r, i1i+i2i),
            (-0.5i1r + h*i1i - 0.5i2r - h*i2i,
             -h*i1r - 0.5i1i + h*i2r - 0.5i2i),
            (-0.5i1r - h*i1i - 0.5i2r + h*i2i,
              h*i1r - 0.5i1i - h*i2r - 0.5i2i))
end

struct _ControlHandles
    phase_voltage
    voltage_sequence
    voltage_min
    voltage_max
    eta
    p_request
    q_request
    power_scale
    i1_request
    i2_voltage
    i2_ripple
    i2_request
    current_scale
    i1
    i2
    phase_current
    sequence_power
    total_power
    ripple_power
    sb::Float64
    vb::Float64
    ib::Float64
    controller::SequenceController
    request::InverterControlRequest
    ratings::InverterControlRatings
end

_policy_laws(policy::WorstPhaseVoltVarWatt) =
    (policy.volt_watt, policy.volt_var)
_policy_laws(policy::AverageVoltageVoltVarWatt) =
    (policy.volt_watt, policy.volt_var)
_policy_laws(policy::PositiveSequenceVoltVarWatt) =
    (policy.volt_watt, policy.volt_var)

function _control_voltage_start(policy::AbstractPositiveSequencePolicy, vb)
    voltage_points = Float64[]
    for law in _policy_laws(policy)
        law === nothing || append!(voltage_points, law.breakpoints)
    end
    return isempty(voltage_points) ? 230.0/vb :
        (minimum(voltage_points) + maximum(voltage_points))/(2vb)
end

function _smooth_extrema_implicit!(
        m, magnitudes, epsilon; voltage_start)
    spread_start = max(0.05voltage_start, epsilon)
    vmin12 = _smooth_min_implicit!(
        m, magnitudes[1], magnitudes[2], epsilon; start=spread_start)
    vmax12 = _smooth_max_implicit!(
        m, magnitudes[1], magnitudes[2], epsilon; start=spread_start)
    vmin = _smooth_min_implicit!(
        m, vmin12, magnitudes[3], epsilon; start=spread_start)
    vmax = _smooth_max_implicit!(
        m, vmax12, magnitudes[3], epsilon; start=spread_start)
    return vmin, vmax
end

function _safe_direction_scale_implicit!(
        m, offset, direction, limit, epsilon; magnitude_start, scale)
    offset_magnitude = all(x -> x isa Real && iszero(x), offset) ? 0.0 :
        _implicit_sqrt!(
            m, offset[1]^2 + offset[2]^2;
            start=magnitude_start, scale=scale)
    direction_magnitude = _implicit_sqrt!(
        m, direction[1]^2 + direction[2]^2;
        start=magnitude_start, scale=scale)
    available = limit - offset_magnitude
    denominator = _smooth_max_implicit!(
        m, available, direction_magnitude, epsilon;
        start=max(Float64(limit), epsilon))
    return available / denominator
end

function _limit_positive_power_smooth!(
        m, p_raw, q_raw, smax, epsilon, priority)
    raw_magnitude = _implicit_sqrt!(
        m, p_raw^2 + q_raw^2; start=smax, scale=smax)
    if priority === :proportional
        denominator = _smooth_max_implicit!(
            m, smax, raw_magnitude, epsilon; start=smax)
        scale = smax/denominator
        return scale*p_raw, scale*q_raw, scale
    elseif priority === :watt
        p_limited = _smooth_min_implicit!(
            m, p_raw, smax, epsilon; start=smax)
        q_capacity = _implicit_sqrt!(
            m, smax^2-p_limited^2; start=0.5smax, scale=smax)
        q_magnitude = _implicit_sqrt!(
            m, q_raw^2; start=0.25smax, scale=smax)
        denominator = _smooth_max_implicit!(
            m, q_capacity, q_magnitude, epsilon; start=smax)
        q_limited = q_raw*q_capacity/denominator
    else
        q_magnitude = _implicit_sqrt!(
            m, q_raw^2; start=0.25smax, scale=smax)
        q_denominator = _smooth_max_implicit!(
            m, smax, q_magnitude, epsilon; start=smax)
        q_limited = q_raw*smax/q_denominator
        p_capacity = _implicit_sqrt!(
            m, smax^2-q_limited^2; start=0.5smax, scale=smax)
        p_denominator = _smooth_max_implicit!(
            m, p_capacity, p_raw, epsilon; start=smax)
        p_limited = p_raw*p_capacity/p_denominator
    end
    limited_magnitude = _implicit_sqrt!(
        m, p_limited^2 + q_limited^2; start=smax, scale=smax)
    raw_denominator = _smooth_max_implicit!(
        m, raw_magnitude, epsilon, epsilon; start=smax)
    scale = limited_magnitude / raw_denominator
    return p_limited, q_limited, scale
end

function _combine_safe_scale_implicit!(m, scale, candidate, epsilon)
    return _smooth_min_implicit!(
        m, scale, candidate, epsilon; start=1.0)
end

function _worst_phase_smooth(ctx, policy::WorstPhaseVoltVarWatt,
                             magnitudes, request, sb, vb)
    m = _opf_model(ctx)
    eps_v = policy.extrema_epsilon / vb
    voltage_start = _control_voltage_start(policy, vb)
    vmin, vmax = _smooth_extrema_implicit!(
        m, magnitudes, eps_v; voltage_start=voltage_start)
    p_fraction = policy.volt_watt === nothing ? 1.0 :
        _pwl_smooth(ctx, policy.volt_watt, vmax, vb)
    q_fraction = if policy.volt_var === nothing
        0.0
    else
        qlo = _smooth_positive_implicit!(
            m, _pwl_smooth(ctx, policy.volt_var, vmin, vb),
            policy.conflict_epsilon; start=1.0)
        qhi = _smooth_negative_implicit!(
            m, _pwl_smooth(ctx, policy.volt_var, vmax, vb),
            policy.conflict_epsilon; start=1.0)
        policy.conflict_policy === :low_voltage ? qlo :
        policy.conflict_policy === :high_voltage ? qhi :
        policy.conflict_policy === :dominant ? _smooth_dominant_implicit!(
            m, qlo, qhi, policy) : qlo + qhi
    end
    return _volt_watt_request(policy, request, p_fraction)/sb,
           request.q_scale/sb*q_fraction, vmin, vmax
end


function _average_voltage_smooth(
        ctx, policy::AverageVoltageVoltVarWatt,
        magnitudes, u1, request, sb, vb)
    voltage = sum(magnitudes)/3
    m = _opf_model(ctx)
    voltage_start = _control_voltage_start(policy, vb)
    vmin, vmax = _smooth_extrema_implicit!(
        m, magnitudes, policy.extrema_epsilon/vb;
        voltage_start=voltage_start)
    p_fraction = policy.volt_watt === nothing ? 1.0 :
        _pwl_smooth(ctx, policy.volt_watt, voltage, vb)
    q_fraction = policy.volt_var === nothing ? 0.0 :
        _pwl_smooth(ctx, policy.volt_var, voltage, vb)
    return _volt_watt_request(policy, request, p_fraction)/sb,
           request.q_scale/sb*q_fraction, vmin, vmax
end


function _positive_sequence_smooth(
        ctx, policy::PositiveSequenceVoltVarWatt,
        magnitudes, u1, request, sb, vb)
    m = _opf_model(ctx)
    voltage_start = _control_voltage_start(policy, vb)
    voltage = _implicit_sqrt!(
        m, u1[1]^2 + u1[2]^2;
        start=voltage_start, scale=voltage_start)
    eps_v = policy.guard_epsilon/vb
    vmin, vmax = _smooth_extrema_implicit!(
        m, magnitudes, eps_v; voltage_start=voltage_start)
    watt_voltage = policy.worst_phase_watt_guard ? vmax : voltage
    p_fraction = policy.volt_watt === nothing ? 1.0 :
        _pwl_smooth(ctx, policy.volt_watt, watt_voltage, vb)
    q_fraction = policy.volt_var === nothing ? 0.0 :
        _pwl_smooth(ctx, policy.volt_var, voltage, vb)
    return _volt_watt_request(policy, request, p_fraction)/sb,
           request.q_scale/sb*q_fraction, vmin, vmax
end


_positive_smooth(ctx, policy::WorstPhaseVoltVarWatt,
                 magnitudes, u1, request, sb, vb) =
    _worst_phase_smooth(ctx, policy, magnitudes, request, sb, vb)
_positive_smooth(ctx, policy::AverageVoltageVoltVarWatt,
                 magnitudes, u1, request, sb, vb) =
    _average_voltage_smooth(ctx, policy, magnitudes, u1, request, sb, vb)
_positive_smooth(ctx, policy::PositiveSequenceVoltVarWatt,
                 magnitudes, u1, request, sb, vb) =
    _positive_sequence_smooth(ctx, policy, magnitudes, u1, request, sb, vb)

function _unbalance_smooth(
        ctx, ::NoUnbalanceControl, u1, u2, i1, vb, ib, voltage_start)
    zero_pair = (0.0, 0.0)
    return zero_pair, zero_pair, zero_pair, 0.0
end

function _unbalance_smooth(ctx, policy::NegativeSequenceAdmittanceDroop,
                           u1, u2, i1, vb, ib, voltage_start)
    u1r, u1i = u1; u2r, u2i = u2
    m = _opf_model(ctx)
    floor = policy.voltage_floor / vb
    denom = u1r^2 + u1i^2 + floor^2
    u2mag = _implicit_sqrt!(
        m, u2r^2 + u2i^2; start=0.05voltage_start,
        scale=voltage_start)
    u1denom = _implicit_sqrt!(
        m, denom; start=voltage_start, scale=voltage_start)
    eta = u2mag / u1denom
    gain = _pwl_smooth(ctx, policy.gain, eta, 1.0) * vb / ib
    c = cos(policy.impedance_angle); s = sin(policy.impedance_angle)
    rotated = _cmul(c, -s, u2r, u2i)
    i2v = (-gain*rotated[1], -gain*rotated[2])
    ratio = _cmul(u2r, u2i, u1r, -u1i)
    ratio = (ratio[1]/denom, ratio[2]/denom)
    product = _cmul(ratio[1], ratio[2], i1[1], i1[2])
    i2r = (-product[1], -product[2])
    blend = policy.ripple_blend
    i2 = ((1-blend)*i2v[1] + blend*i2r[1],
          (1-blend)*i2v[2] + blend*i2r[2])
    return i2v, i2r, i2, eta
end

"""
    stamp_smooth_control!(ctx, controller, request, inverter, handles)

Constrain an already-stamped [`AdvancedInverter`](@ref) to the smooth,
fixed-structure counterpart of [`evaluate_exact`](@ref). The controller uses
local phase-to-neutral POC phasors and commands the configured converter- or
grid-side phase currents.
"""
function stamp_smooth_control!(ctx, controller::SequenceController,
                               request::InverterControlRequest,
                               inverter::AdvancedInverter,
                               handles::_InvHandles)
    inverter.topology == :THREE_LEG || throw(ArgumentError(
        "smooth SequenceController stamping currently supports THREE_LEG only"))
    ratings = InverterControlRatings(inverter, controller.current_target)
    m = _opf_model(ctx)
    sb, vb, ib = handles.sb, handles.vb, handles.ib
    phase_voltage = ntuple(3) do k
        ph = inverter.phase_terminals[k]
        _dv(ctx, inverter.bus, ph, inverter.neutral)
    end
    vre = [phase_voltage[k][1] for k in 1:3]
    vim = [phase_voltage[k][2] for k in 1:3]
    voltage_start = _control_voltage_start(controller.positive, vb)
    magnitudes = [
        _implicit_sqrt!(
            m, vre[k]^2 + vim[k]^2;
            start=voltage_start, scale=voltage_start) for k in 1:3]
    u0, u1, u2 = _sequence_pair(vre, vim)

    p_raw, q_raw, vmin, vmax = _positive_smooth(
        ctx, controller.positive, magnitudes, u1, request, sb, vb)
    smax = ratings.s_max / sb
    seps = controller.limiter.power_epsilon / sb
    if _volt_watt_basis(controller.positive) === :rated
        p_raw = _smooth_min_implicit!(
            m, request.p_available/sb, p_raw, seps;
            start=request.p_available/sb)
    end
    p, q, power_scale = _limit_positive_power_smooth!(
        m, p_raw, q_raw, smax, seps, controller.limiter.pq_priority)

    floor = controller.power_voltage_floor / vb
    denominator = 3(u1[1]^2 + u1[2]^2 + floor^2)
    i1_request = ((p*u1[1] + q*u1[2]) / denominator,
                  (p*u1[2] - q*u1[1]) / denominator)
    i2v, i2r, i2_request, eta = _unbalance_smooth(
        ctx, controller.unbalance, u1, u2, i1_request, vb, ib, voltage_start)
    phase_request = _phase_pairs(i1_request, i2_request)

    ieps = controller.limiter.current_epsilon / ib
    imax = ratings.i_max / ib
    imags = [_implicit_sqrt!(
        m, pair[1]^2 + pair[2]^2;
        start=max(0.5imax, ieps), scale=imax)
        for pair in phase_request]
    current_denominator = imax
    for imag in imags
        current_denominator = _smooth_max_implicit!(
            m, current_denominator, imag, ieps; start=imax)
    end
    target_scale = imax / current_denominator
    i1_pre = (target_scale*i1_request[1], target_scale*i1_request[2])
    i2_pre = (target_scale*i2_request[1], target_scale*i2_request[2])
    phase_pre = _phase_pairs(i1_pre, i2_pre)

    # The current command is finally backed off against the physical currents at
    # both filter locations and converter-terminal apparent power. For an LCL
    # filter, I_conv = I_grid + I_shunt. Holding the solved local shunt phasor
    # fixed makes both currents affine in one scalar command. The triangle-bound
    # allocator below is conservative, continuous, and fixed-structure; it
    # prevents a protection limit from turning ordinary saturation into an
    # infeasible controller equality.
    zero_pair = (0.0, 0.0)
    shunt = ntuple(k -> handles.cri[k] === handles.gri[k] ? zero_pair : (
        handles.cri[k] - handles.gri[k],
        handles.cii[k] - handles.gii[k]), 3)
    converter_offset = controller.current_target isa ConverterCurrentTarget ?
        ntuple(_ -> zero_pair, 3) : shunt
    grid_offset = controller.current_target isa ConverterCurrentTarget ?
        ntuple(k -> (-shunt[k][1], -shunt[k][2]), 3) :
        ntuple(_ -> zero_pair, 3)

    capability_scale = 1.0
    scale_eps = max(ieps, eps(Float64))
    if inverter.i_max !== nothing
        for k in 1:3
            reserve = inverter.pwm_ac_converter_reserve[k] / ib
            limit = sqrt(max((inverter.i_max/ib)^2 - reserve^2, 0.0))
            candidate = _safe_direction_scale_implicit!(
                m, converter_offset[k], phase_pre[k], limit, scale_eps;
                magnitude_start=max(0.25limit, scale_eps), scale=max(limit, scale_eps))
            capability_scale = _combine_safe_scale_implicit!(
                m, capability_scale, candidate, scale_eps)
        end
    end
    if inverter.i_grid_max !== nothing
        for k in 1:3
            reserve = inverter.pwm_ac_grid_reserve[k] / ib
            limit = sqrt(max((inverter.i_grid_max/ib)^2 - reserve^2, 0.0))
            candidate = _safe_direction_scale_implicit!(
                m, grid_offset[k], phase_pre[k], limit, scale_eps;
                magnitude_start=max(0.25limit, scale_eps), scale=max(limit, scale_eps))
            capability_scale = _combine_safe_scale_implicit!(
                m, capability_scale, candidate, scale_eps)
        end
    end

    converter_direction = phase_pre
    converter_power_offset = (0.0, 0.0)
    converter_power_direction = (0.0, 0.0)
    for k in 1:3
        vo = (handles.vrint[k], handles.viint[k])
        so = _cmul(
            vo[1], vo[2], converter_offset[k][1], -converter_offset[k][2])
        sd = _cmul(
            vo[1], vo[2], converter_direction[k][1], -converter_direction[k][2])
        converter_power_offset = (
            converter_power_offset[1] + so[1],
            converter_power_offset[2] + so[2])
        converter_power_direction = (
            converter_power_direction[1] + sd[1],
            converter_power_direction[2] + sd[2])
    end
    power_candidate = _safe_direction_scale_implicit!(
        m, converter_power_offset, converter_power_direction, smax, seps;
        magnitude_start=max(0.25smax, seps), scale=smax)
    capability_scale = _combine_safe_scale_implicit!(
        m, capability_scale, power_candidate, max(seps, eps(Float64)))
    if inverter.dv2_max !== nothing
        ripple_offset = (0.0, 0.0)
        ripple_direction = (0.0, 0.0)
        for k in 1:3
            vo = (handles.vrint[k], handles.viint[k])
            ro = _cmul(
                vo[1], vo[2], converter_offset[k][1], converter_offset[k][2])
            rd = _cmul(
                vo[1], vo[2], converter_direction[k][1],
                converter_direction[k][2])
            ripple_offset = (
                ripple_offset[1] + ro[1], ripple_offset[2] + ro[2])
            ripple_direction = (
                ripple_direction[1] + rd[1], ripple_direction[2] + rd[2])
        end
        ripple_limit = 2(2pi*inverter.f)*inverter.c_dc*inverter.v_dc *
                       inverter.dv2_max / sb
        ripple_candidate = _safe_direction_scale_implicit!(
            m, ripple_offset, ripple_direction, ripple_limit, seps;
            magnitude_start=max(0.25ripple_limit, seps),
            scale=max(ripple_limit, seps))
        capability_scale = _combine_safe_scale_implicit!(
            m, capability_scale, ripple_candidate,
            max(seps, eps(Float64)))
    end

    current_scale = target_scale*capability_scale
    i1 = (capability_scale*i1_pre[1], capability_scale*i1_pre[2])
    i2 = (capability_scale*i2_pre[1], capability_scale*i2_pre[2])
    phase_current = _phase_pairs(i1, i2)
    target_re, target_im = controller.current_target isa ConverterCurrentTarget ?
        (handles.cri, handles.cii) : (handles.gri, handles.gii)
    for k in 1:3
        JuMP.@constraint(m, target_re[k] == phase_current[k][1])
        JuMP.@constraint(m, target_im[k] == phase_current[k][2])
    end

    u1_conj_i1 = _cmul(u1[1], u1[2], i1[1], -i1[2])
    u2_conj_i2 = _cmul(u2[1], u2[2], i2[1], -i2[2])
    sequence_power = (3u1_conj_i1[1], 3u1_conj_i1[2])
    total_power = (3(u1_conj_i1[1] + u2_conj_i2[1]),
                   3(u1_conj_i1[2] + u2_conj_i2[2]))
    u1_i2 = _cmul(u1[1], u1[2], i2[1], i2[2])
    u2_i1 = _cmul(u2[1], u2[2], i1[1], i1[2])
    ripple = (3(u1_i2[1] + u2_i1[1]), 3(u1_i2[2] + u2_i1[2]))
    return _ControlHandles(
        phase_voltage, (u0, u1, u2), vmin, vmax, eta, p, q, power_scale,
        i1_request, i2v, i2r, i2_request, current_scale, i1, i2,
        phase_current, sequence_power, total_power, ripple, sb, vb, ib,
        controller, request, ratings)
end

struct _ControlledHandles
    plant::_InvHandles
    control::_ControlHandles
end

function stamp_device!(ctx, controlled::ControlledDevice{<:AdvancedInverter};
                       request::InverterControlRequest)
    inverter = inverter_spec(controlled)
    plant = stamp_device!(ctx, inverter)
    control = stamp_smooth_control!(ctx, controlled.controller, request,
                                    inverter,
                                    inverter_handles(controlled, plant))
    _ControlledHandles(plant, control)
end

function _extract_control(handles::_ControlHandles, status::SolveStatus)
    status.publishable || return InverterControlResult(
        ntuple(_ -> ComplexF64(NaN, NaN), 3),
        ntuple(_ -> ComplexF64(NaN, NaN), 3),
        NaN, NaN, NaN, NaN, NaN, NaN,
        ComplexF64(NaN, NaN), ComplexF64(NaN, NaN),
        ComplexF64(NaN, NaN), ComplexF64(NaN, NaN), NaN,
        ComplexF64(NaN, NaN), ComplexF64(NaN, NaN),
        ntuple(_ -> ComplexF64(NaN, NaN), 3),
        ComplexF64(NaN, NaN), ComplexF64(NaN, NaN),
        ComplexF64(NaN, NaN))
    scalar(x) = Float64(JuMP.value(x))
    voltage_pair(pair) = ComplexF64(
        scalar(pair[1])*handles.vb, scalar(pair[2])*handles.vb)
    current_pair(pair) = ComplexF64(
        scalar(pair[1])*handles.ib, scalar(pair[2])*handles.ib)
    power_pair(pair) = ComplexF64(
        scalar(pair[1])*handles.sb, scalar(pair[2])*handles.sb)
    phase_voltage = ntuple(k -> voltage_pair(handles.phase_voltage[k]), 3)
    voltage_sequence = ntuple(
        k -> voltage_pair(handles.voltage_sequence[k]), 3)
    phase_current = ntuple(k -> current_pair(handles.phase_current[k]), 3)
    return InverterControlResult(
        phase_voltage, voltage_sequence,
        scalar(handles.voltage_min)*handles.vb,
        scalar(handles.voltage_max)*handles.vb,
        scalar(handles.eta),
        scalar(handles.p_request)*handles.sb,
        scalar(handles.q_request)*handles.sb,
        scalar(handles.power_scale),
        current_pair(handles.i1_request),
        current_pair(handles.i2_voltage),
        current_pair(handles.i2_ripple),
        current_pair(handles.i2_request),
        scalar(handles.current_scale),
        current_pair(handles.i1), current_pair(handles.i2), phase_current,
        power_pair(handles.sequence_power), power_pair(handles.total_power),
        power_pair(handles.ripple_power))
end

function _extract_converter_terminal(
        handles::_InvHandles, status::SolveStatus)
    nan_phasors = ntuple(_ -> ComplexF64(NaN, NaN), 3)
    status.publishable || return ConverterTerminalResult(
        nan_phasors, nan_phasors, nan_phasors, nan_phasors, nan_phasors,
        ComplexF64(NaN, NaN), ComplexF64(NaN, NaN),
        ComplexF64(NaN, NaN), ComplexF64(NaN, NaN),
        ComplexF64(NaN, NaN))
    phase_voltage = ntuple(k -> ComplexF64(
        JuMP.value(handles.vrint[k])*handles.vb,
        JuMP.value(handles.viint[k])*handles.vb), 3)
    phase_current = ntuple(k -> ComplexF64(
        JuMP.value(handles.cri[k])*handles.ib,
        JuMP.value(handles.cii[k])*handles.ib), 3)
    return _converter_terminal_result(phase_voltage, phase_current)
end

function _extract_grid_current(handles::_InvHandles, status::SolveStatus)
    status.publishable || return ntuple(_ -> ComplexF64(NaN, NaN), 3)
    return ntuple(k -> ComplexF64(
        JuMP.value(handles.gri[k])*handles.ib,
        JuMP.value(handles.gii[k])*handles.ib), 3)
end

function _safe_direction_scale_exact(offset, direction, limit)
    available = limit - abs(offset)
    available <= 0 && return 0.0
    magnitude = abs(direction)
    return iszero(magnitude) ? 1.0 : min(1.0, available/magnitude)
end

function _apply_plant_capability_exact(
        result::InverterControlResult,
        converter_voltage::NTuple{3,ComplexF64},
        converter_current::NTuple{3,ComplexF64},
        grid_current::NTuple{3,ComplexF64},
        inverter::AdvancedInverter,
        target::AbstractCurrentTarget)
    shunt = converter_current .- grid_current
    zeros3 = ntuple(_ -> 0.0im, 3)
    converter_offset = target isa ConverterCurrentTarget ? zeros3 : shunt
    grid_offset = target isa ConverterCurrentTarget ? .-shunt : zeros3
    direction = result.phase_current
    scale = 1.0
    if inverter.i_max !== nothing
        for k in 1:3
            reserve = inverter.pwm_ac_converter_reserve[k]
            limit = sqrt(max(inverter.i_max^2-reserve^2, 0.0))
            scale = min(scale, _safe_direction_scale_exact(
                converter_offset[k], direction[k], limit))
        end
    end
    if inverter.i_grid_max !== nothing
        for k in 1:3
            reserve = inverter.pwm_ac_grid_reserve[k]
            limit = sqrt(max(inverter.i_grid_max^2-reserve^2, 0.0))
            scale = min(scale, _safe_direction_scale_exact(
                grid_offset[k], direction[k], limit))
        end
    end
    power_offset = sum(
        converter_voltage[k]*conj(converter_offset[k]) for k in 1:3)
    power_direction = sum(
        converter_voltage[k]*conj(direction[k]) for k in 1:3)
    scale = min(scale, _safe_direction_scale_exact(
        power_offset, power_direction, inverter.s_max))
    if inverter.dv2_max !== nothing
        ripple_offset = sum(
            converter_voltage[k]*converter_offset[k] for k in 1:3)
        ripple_direction = sum(
            converter_voltage[k]*direction[k] for k in 1:3)
        ripple_limit = 2(2pi*inverter.f)*inverter.c_dc*inverter.v_dc *
                       inverter.dv2_max
        scale = min(scale, _safe_direction_scale_exact(
            ripple_offset, ripple_direction, ripple_limit))
    end

    return InverterControlResult(
        result.phase_voltage, result.voltage_sequence,
        result.voltage_min, result.voltage_max, result.eta,
        result.p_request, result.q_request, result.power_scale,
        result.i1_request, result.i2_voltage, result.i2_ripple,
        result.i2_request, result.current_scale*scale,
        scale*result.i1, scale*result.i2,
        ntuple(k -> scale*result.phase_current[k], 3),
        scale*result.sequence_power, scale*result.total_power,
        scale*result.ripple_power)
end

function extract_device(controlled::ControlledDevice{<:AdvancedInverter},
                        handles::_ControlledHandles, status::SolveStatus)
    return (plant=extract_device(controlled.device, handles.plant, status),
            control=_extract_control(handles.control, status),
            converter_terminal=_extract_converter_terminal(
                handles.plant, status))
end

"""Result of [`solve_controlled_inverter`](@ref)."""
struct ControlledInverterResult <: AbstractSolveResult
    termination_status::String
    plant::NamedTuple
    control::InverterControlResult
    converter_terminal::ConverterTerminalResult
    grid_phase_current::NTuple{3,ComplexF64}
    exact_control::InverterControlResult
    exact_smooth_current_residual::Float64
    bus::Dict{String,Any}
    solve::SolveStatus
end

solve_status(result::ControlledInverterResult) = result.solve
solve_diagnostics(result::ControlledInverterResult) = (
    current_scale=result.control.current_scale,
    power_scale=result.control.power_scale,
    converter_total_power=result.converter_terminal.total_power,
    converter_positive_sequence_power=
        result.converter_terminal.positive_sequence_power,
    exact_smooth_current_residual=result.exact_smooth_current_residual,
)

"""
    solve_controlled_inverter(net, controlled, request; kwargs...)
        -> ControlledInverterResult

Solve one network snapshot with a local phase-aware controller composed around
an [`AdvancedInverter`](@ref). The local law uses only the inverter's POC
voltage phasors. `per_unit=true` changes numerical coordinates only; controller
configuration and returned results remain SI. This is a controlled power flow:
the controller equalities determine the command, while `selection_objective`
only selects remaining plant allocation freedom (`:loss`, the default, or
`:zero` for an objective-invariance check).
"""
function solve_controlled_inverter(
        net::Dict{String,Any},
        controlled::ControlledDevice{<:AdvancedInverter},
        request::InverterControlRequest;
        per_unit::Bool=true,
        s_base::Real=1e6,
        optimizer=Ipopt.Optimizer,
        verbose::Bool=false,
        selection_objective::Symbol=:loss,
        solver_options=())
    selection_objective in (:loss, :zero) || throw(ArgumentError(
        "selection_objective must be :loss or :zero"))
    validate_device(controlled, (net,); periods=1)
    handles = Ref{_ControlledHandles}()
    hook! = ctx -> begin
        handles[] = stamp_device!(ctx, controlled; request=request)
        if selection_objective === :loss
            JuMP.@objective(_opf_model(ctx), Min,
                handles[].plant.p_loss + handles[].plant.p_cap_loss)
        else
            JuMP.@objective(_opf_model(ctx), Min, 0.0)
        end
    end
    ctx = build_opf_model(net; per_unit=per_unit, s_base=Float64(s_base),
                          add_objective=false, model_hook! = hook!,
                          optimizer=optimizer, verbose=verbose)
    _set_solver_options!(_opf_model(ctx), solver_options)
    enforce_kcl!(ctx)
    JuMP.optimize!(_opf_model(ctx))
    outcome = _solve_outcome(_opf_model(ctx))
    status = SolveStatus(outcome)
    extracted = extract_device(controlled, handles[], status)
    result = _extract_result(ctx, outcome)
    exact_control = status.publishable ? evaluate_exact(
        controlled.controller,
        InverterControlMeasurement(collect(extracted.control.phase_voltage)),
        request, InverterControlRatings(
            controlled.device, controlled.controller.current_target)) :
        extracted.control
    grid_current = _extract_grid_current(handles[].plant, status)
    if status.publishable
        exact_control = _apply_plant_capability_exact(
            exact_control, extracted.converter_terminal.phase_voltage,
            extracted.converter_terminal.phase_current, grid_current,
            controlled.device, controlled.controller.current_target)
    end
    residual = status.publishable ? maximum(abs,
        exact_control.phase_current .- extracted.control.phase_current) : NaN
    return ControlledInverterResult(string(outcome.termination_status),
                                    extracted.plant, extracted.control,
                                    extracted.converter_terminal,
                                    grid_current,
                                    exact_control, residual, result["bus"], status)
end
