const _CTRL_A = cis(2pi / 3)

_ctrl_phases(u1::Complex, u2::Complex=0.0im) = ComplexF64[
    u1 + u2,
    _CTRL_A^2*u1 + _CTRL_A*u2,
    _CTRL_A*u1 + _CTRL_A^2*u2,
]

function _study_controller(; ripple_blend=0.5)
    volt_watt = PiecewiseLinearLaw(
        [230.0, 240.0, 250.0], [1.0, 1.0, 0.2];
        smoothing_epsilon=0.05)
    volt_var = PiecewiseLinearLaw(
        [210.0, 220.0, 240.0, 250.0], [0.3, 0.0, 0.0, -0.3];
        smoothing_epsilon=0.05)
    negative_gain = PiecewiseLinearLaw(
        [0.0, 0.01, 0.10], [0.0, 0.0, 0.08];
        smoothing_epsilon=1e-4)
    positive = WorstPhaseVoltVarWatt(
        volt_watt=volt_watt, volt_var=volt_var,
        conflict_policy=:net,
        extrema_epsilon=1e-3, conflict_epsilon=1e-6)
    unbalance = NegativeSequenceAdmittanceDroop(
        negative_gain; impedance_angle=0.0, ripple_blend=ripple_blend,
        voltage_floor=1.0)
    return SequenceController(positive, unbalance, CommonScaleLimiter())
end

const _CTRL_REQUEST = InverterControlRequest(
    p_available=12e3, q_scale=8e3)
const _CTRL_RATINGS = InverterControlRatings(s_max=20e3, i_max=40.0)

@testset "Inverter controls: typed curves and validation" begin
    @test_throws ArgumentError PiecewiseLinearLaw(
        [1.0, 1.0], [0.0, 1.0]; smoothing_epsilon=0.1)
    @test_throws ArgumentError PiecewiseLinearLaw(
        [1.0, 2.0], [0.0]; smoothing_epsilon=0.1)
    @test_throws ArgumentError PiecewiseLinearLaw(
        [1.0, 2.0], [0.0, 1.0]; smoothing_epsilon=0.0)
    @test_throws ArgumentError WorstPhaseVoltVarWatt(
        volt_watt=PiecewiseLinearLaw(
            [220.0, 240.0], [1.0, 1.1]; smoothing_epsilon=0.1))
    @test_throws ArgumentError WorstPhaseVoltVarWatt(
        volt_var=PiecewiseLinearLaw(
            [220.0, 240.0], [-0.1, 0.1]; smoothing_epsilon=0.1))
    @test_throws ArgumentError NegativeSequenceAdmittanceDroop(
        PiecewiseLinearLaw([0.0, 1.0], [0.0, -1.0];
                           smoothing_epsilon=0.1))
    @test_throws ArgumentError InverterControlMeasurement(ComplexF64[230, 230])
    @test_throws ArgumentError InverterControlRequest(
        p_available=-1.0, q_scale=1.0)
    @test_throws ArgumentError InverterControlRequest(
        p_available=1.0, p_rated=0.0, q_scale=1.0)
    @test_throws ArgumentError InverterControlRatings(s_max=1.0, i_max=0.0)

    policy = WorstPhaseVoltVarWatt(
        volt_watt=PiecewiseLinearLaw(
            [220.0, 240.0], [1.0, 0.5]; smoothing_epsilon=0.1))
    controller = SequenceController(policy)
    @test controller.unbalance isa NoUnbalanceControl
    @test controller.limiter isa CommonScaleLimiter
    @test controller.current_target isa ConverterCurrentTarget
    @test controller.power_voltage_floor == 1.0
    @test SequenceController(policy; current_target=:grid).current_target isa
          GridCurrentTarget
    @test_throws ArgumentError SequenceController(
        policy; current_target=:bogus)
    @test_throws ArgumentError SequenceController(
        policy; power_voltage_floor=0.0)
    @test_throws ArgumentError WorstPhaseVoltVarWatt(
        volt_watt_basis=:invalid)
    @test_throws ArgumentError CommonScaleLimiter(pq_priority=:invalid)
    @test_throws ArgumentError CommonScaleLimiter(
        priority_headroom_fraction=0.0)
    @test_throws ArgumentError CommonScaleLimiter(
        priority_headroom_fraction=1.0)
    @test_throws ArgumentError CommonScaleLimiter(
        current_epsilon_fraction=0.0)
    @test CommonScaleLimiter().current_epsilon === nothing
    @test CommonScaleLimiter().power_epsilon === nothing
    @test CommonScaleLimiter(current_epsilon=2.0).current_epsilon == 2.0
    reserve_inverter = AdvancedInverter(
        id="reserve", bus="poc", s_max=10e3, i_max=10.0,
        pwm_ac_converter_reserve=(6.0, 8.0, 0.0))
    @test InverterControlRatings(reserve_inverter).i_max ≈ 6.0
end

@testset "Inverter controls: exact worst-phase and sequence laws" begin
    # Each phase occupies a different part of both curves. The active request is
    # governed by the high-voltage phase, while opposing var requests net.
    phasors = ComplexF64[
        245cis(0.05), 215cis(-2.15), 230cis(2.0)]
    result = evaluate_exact(
        _study_controller(), InverterControlMeasurement(phasors),
        _CTRL_REQUEST, _CTRL_RATINGS)
    @test result.voltage_min ≈ 215.0
    @test result.voltage_max ≈ 245.0
    @test result.p_request ≈ 0.6*_CTRL_REQUEST.p_available
    @test result.q_request ≈ 0.0 atol=1e-10
    @test abs(sum(result.phase_current)) < 1e-12
    @test maximum(abs, result.phase_current) <= _CTRL_RATINGS.i_max + 1e-12
    @test abs(result.i2) > 0.0

    low_priority = SequenceController(WorstPhaseVoltVarWatt(
        volt_var=_study_controller().positive.volt_var,
        conflict_policy=:low_voltage))
    high_priority = SequenceController(WorstPhaseVoltVarWatt(
        volt_var=_study_controller().positive.volt_var,
        conflict_policy=:high_voltage))
    lo = evaluate_exact(low_priority, InverterControlMeasurement(phasors),
                        _CTRL_REQUEST, _CTRL_RATINGS)
    hi = evaluate_exact(high_priority, InverterControlMeasurement(phasors),
                        _CTRL_REQUEST, _CTRL_RATINGS)
    @test lo.q_request > 0
    @test hi.q_request < 0

    balanced = evaluate_exact(
        SequenceController(_study_controller().positive),
        InverterControlMeasurement(_ctrl_phases(230.0 + 0im)),
        _CTRL_REQUEST, _CTRL_RATINGS)
    @test abs(balanced.voltage_sequence[1]) < 1e-12
    @test abs(balanced.voltage_sequence[3]) < 1e-12
    @test balanced.i2 == 0.0im

    # Exact curve oracle: clamping, both slopes, and both deadband endpoints.
    regimes = (
        (205.0, 1.0, 0.30),
        (215.0, 1.0, 0.15),
        (220.0, 1.0, 0.00),
        (240.0, 1.0, 0.00),
        (245.0, 0.6, -0.15),
        (250.0, 0.2, -0.30),
        (255.0, 0.2, -0.30),
    )
    for (voltage, p_fraction, q_fraction) in regimes
        curve = evaluate_exact(
            _study_controller(),
            InverterControlMeasurement(_ctrl_phases(voltage + 0im)),
            _CTRL_REQUEST, _CTRL_RATINGS)
        @test curve.p_request ≈
              p_fraction*_CTRL_REQUEST.p_available atol=1e-9
        @test curve.q_request ≈
              q_fraction*_CTRL_REQUEST.q_scale atol=1e-9
    end
end

@testset "Inverter controls: voltage references and conflict policies" begin
    measurement = InverterControlMeasurement(ComplexF64[
        245cis(0.05), 215cis(-2.15), 230cis(2.0)])
    curves = _study_controller().positive
    average = SequenceController(AverageVoltageVoltVarWatt(
        volt_watt=curves.volt_watt, volt_var=curves.volt_var))
    guarded = SequenceController(PositiveSequenceVoltVarWatt(
        volt_watt=curves.volt_watt, volt_var=curves.volt_var))
    unguarded = SequenceController(PositiveSequenceVoltVarWatt(
        volt_watt=curves.volt_watt, volt_var=curves.volt_var,
        worst_phase_watt_guard=false))
    avg = evaluate_exact(average, measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    guard = evaluate_exact(guarded, measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    no_guard = evaluate_exact(
        unguarded, measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    smooth_avg = evaluate_smooth(
        average, measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    smooth_guard = evaluate_smooth(
        guarded, measurement, _CTRL_REQUEST, _CTRL_RATINGS)

    @test avg.p_request ≈ _CTRL_REQUEST.p_available atol=1e-9
    @test avg.q_request ≈ 0.0 atol=1e-9
    @test avg.eta == 0.0
    @test guard.p_request ≈ 0.6*_CTRL_REQUEST.p_available atol=1e-9
    @test no_guard.p_request > guard.p_request

    # The numeric smooth evaluator uses the same pairwise smooth extrema as the
    # stamped graph; the exact evaluator deliberately retains hard firmware
    # comparisons. `voltage_min/max` always mean physical phase extrema, not a
    # policy's curve input.
    smooth_max(a, b, epsilon) =
        (a+b+sqrt((a-b)^2+epsilon^2))/2
    smooth_min(a, b, epsilon) =
        (a+b-sqrt((a-b)^2+epsilon^2))/2
    magnitudes = abs.(measurement.phase_voltage)
    expected_avg_min = smooth_min(
        smooth_min(magnitudes[1], magnitudes[2],
                   average.positive.extrema_epsilon),
        magnitudes[3], average.positive.extrema_epsilon)
    expected_guard_max = smooth_max(
        smooth_max(magnitudes[1], magnitudes[2],
                   guarded.positive.guard_epsilon),
        magnitudes[3], guarded.positive.guard_epsilon)
    @test smooth_avg.voltage_min ≈ expected_avg_min atol=1e-12
    @test smooth_guard.voltage_max ≈ expected_guard_max atol=1e-12
    @test smooth_guard.voltage_max != guard.voltage_max

    dominant = SequenceController(WorstPhaseVoltVarWatt(
        volt_var=curves.volt_var, conflict_policy=:dominant))
    balanced_tie = evaluate_exact(
        dominant, measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    @test abs(balanced_tie.q_request) < 1e-8

    low_dominant_measurement = InverterControlMeasurement(ComplexF64[
        245cis(0.05), 212cis(-2.15), 230cis(2.0)])
    low_dominant = evaluate_exact(
        dominant, low_dominant_measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    @test low_dominant.q_request > 0.0

    smooth_dominant = evaluate_smooth(
        dominant, low_dominant_measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    @test smooth_dominant.q_request > 0.0
    @test abs(smooth_dominant.q_request-low_dominant.q_request) < 2.0

    asymmetric = SequenceController(WorstPhaseVoltVarWatt(
        volt_var=PiecewiseLinearLaw(
            [210.0, 220.0, 240.0, 250.0], [0.44, 0.0, 0.0, -0.60];
            smoothing_epsilon=0.05),
        conflict_policy=:dominant, conflict_epsilon=1e-6))
    asymmetric_tie = evaluate_exact(
        asymmetric, measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    @test abs(asymmetric_tie.q_request) < 1e-6

    # The deployed dominant law itself is continuous. It is not a discontinuous
    # winner-take-all firmware law hidden behind a model-only smoothing.
    tie_voltages = collect(214.8:0.01:215.2)
    tie_q = [evaluate_exact(
        dominant,
        InverterControlMeasurement(ComplexF64[
            245cis(0.05), v*cis(-2.15), 230cis(2.0)]),
        _CTRL_REQUEST, _CTRL_RATINGS).q_request for v in tie_voltages]
    @test maximum(abs, diff(tie_q)) < 130.0
    @test tie_q[1] > 0 > tie_q[end]
end

@testset "Inverter controls: stamped P/Q priority under PV oversizing" begin
    inverter = AdvancedInverter(
        id="priority_saturation", bus="poc",
        phase_terminals=["a", "b", "c"], neutral="n",
        topology=:THREE_LEG, s_max=20e3, i_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, r_filter=0.05, x_filter=0.15)
    # Sweep the local capability law cheaply across DC/AC ratios 0.9--1.4 with
    # a non-zero var request, including the exact 1.0 transition.
    var_law = PiecewiseLinearLaw(
        [180.0, 280.0], [0.2, 0.2]; smoothing_epsilon=0.05)
    measurement = InverterControlMeasurement(_ctrl_phases(230.0 + 0im))
    swept = Dict{Tuple{Symbol,Float64},InverterControlResult}()
    for priority in (:proportional, :watt, :var),
        ratio in (0.9, 1.0, 1.1, 1.25, 1.4)
        controller = SequenceController(
            AverageVoltageVoltVarWatt(volt_var=var_law);
            limiter=CommonScaleLimiter(pq_priority=priority))
        request = InverterControlRequest(
            p_available=ratio*inverter.s_max,
            p_rated=ratio*inverter.s_max, q_scale=5e3)
        exact = evaluate_exact(
            controller, measurement, request, InverterControlRatings(inverter))
        smooth = evaluate_smooth(
            controller, measurement, request, InverterControlRatings(inverter))
        swept[(priority, ratio)] = smooth
        @test abs(complex(exact.p_request, exact.q_request)) <=
              inverter.s_max + 1e-8
        @test abs(complex(smooth.p_request, smooth.q_request)) <=
              inverter.s_max + 0.01
        @test smooth.p_request ≈ exact.p_request atol=0.01
        @test smooth.q_request ≈ exact.q_request atol=0.01
    end
    for ratio in (1.0, 1.1, 1.25, 1.4)
        @test swept[(:watt, ratio)].p_request ≈
              (1-CommonScaleLimiter().priority_headroom_fraction)*
              inverter.s_max atol=0.01
        @test swept[(:var, ratio)].q_request ≈ 1e3 atol=0.01
        @test swept[(:proportional, ratio)].q_request /
              swept[(:proportional, ratio)].p_request ≈
              1e3/(ratio*inverter.s_max) atol=1e-10
    end

    # Keep the stamped regression compact: all priorities at one saturated
    # point in the supported per-unit formulation.
    for priority in (:proportional, :watt, :var)
        controller = SequenceController(
            AverageVoltageVoltVarWatt();
            limiter=CommonScaleLimiter(pq_priority=priority))
        request = InverterControlRequest(
            p_available=22e3, p_rated=22e3, q_scale=0.0)
        result = solve_controlled_inverter(
            inv_grid3_bal(), ControlledDevice(inverter, controller), request;
            per_unit=true,
            solver_options=("max_iter" => 500, "tol" => 1e-8))
        @test result.solve.publishable
        @test abs(result.converter_terminal.total_power) <=
              inverter.s_max + 0.2
        @test result.control.current_scale < 1.0
    end
end

@testset "Inverter controls: Volt-watt basis and P/Q priority" begin
    law = PiecewiseLinearLaw(
        [240.0, 250.0], [1.0, 0.2]; smoothing_epsilon=0.05)
    measurement = InverterControlMeasurement(_ctrl_phases(247.5 + 0im))
    request = InverterControlRequest(
        p_available=3e3, p_rated=12e3, q_scale=15e3)
    available = SequenceController(AverageVoltageVoltVarWatt(
        volt_watt=law, volt_watt_basis=:available))
    rated = SequenceController(AverageVoltageVoltVarWatt(
        volt_watt=law, volt_watt_basis=:rated))
    @test evaluate_exact(
        available, measurement, request, _CTRL_RATINGS).p_request ≈ 1.2e3
    @test evaluate_exact(
        rated, measurement, request, _CTRL_RATINGS).p_request ≈ 3e3
    @test evaluate_smooth(
        available, measurement, request, _CTRL_RATINGS).p_request ≈ 1.2e3 atol=0.01
    @test evaluate_smooth(
        rated, measurement, request, _CTRL_RATINGS).p_request ≈ 3e3 atol=0.01

    pq_request = InverterControlRequest(
        p_available=18e3, p_rated=18e3, q_scale=30e3)
    pq_measurement = InverterControlMeasurement(_ctrl_phases(210.0 + 0im))
    proportional = evaluate_exact(
        SequenceController(AverageVoltageVoltVarWatt(
            volt_var=_study_controller().positive.volt_var);
            limiter=CommonScaleLimiter(pq_priority=:proportional)),
        pq_measurement, pq_request, _CTRL_RATINGS)
    watt = evaluate_exact(
        SequenceController(AverageVoltageVoltVarWatt(
            volt_var=_study_controller().positive.volt_var);
            limiter=CommonScaleLimiter(pq_priority=:watt)),
        pq_measurement, pq_request, _CTRL_RATINGS)
    var = evaluate_exact(
        SequenceController(AverageVoltageVoltVarWatt(
            volt_var=_study_controller().positive.volt_var);
            limiter=CommonScaleLimiter(pq_priority=:var)),
        pq_measurement, pq_request, _CTRL_RATINGS)
    @test watt.p_request > proportional.p_request > var.p_request
    @test abs(var.q_request) > abs(proportional.q_request) > abs(watt.q_request)
    @test all(r -> hypot(r.p_request, r.q_request) <=
                   _CTRL_RATINGS.s_max + 1e-8,
              (proportional, watt, var))
end

@testset "Inverter controls: PWL smoothing refinement" begin
    measurement = InverterControlMeasurement(_ctrl_phases(240.0 + 0im))
    errors = Float64[]
    for epsilon in (0.1, 0.05, 0.025)
        law = PiecewiseLinearLaw(
            [230.0, 240.0, 250.0], [1.0, 1.0, 0.2];
            smoothing_epsilon=epsilon)
        controller = SequenceController(
            AverageVoltageVoltVarWatt(volt_watt=law))
        exact = evaluate_exact(
            controller, measurement, _CTRL_REQUEST, _CTRL_RATINGS)
        smooth = evaluate_smooth(
            controller, measurement, _CTRL_REQUEST, _CTRL_RATINGS)
        push!(errors, exact.p_request-smooth.p_request)
    end
    @test all(>(0.0), errors)
    @test errors[1]/errors[2] ≈ 2.0 rtol=2e-3
    @test errors[2]/errors[3] ≈ 2.0 rtol=2e-3
    @test errors[2]/_CTRL_REQUEST.p_available < 5e-3
end

@testset "Inverter controls: smooth dominant conflict composition" begin
    curves = _study_controller().positive
    controller = SequenceController(WorstPhaseVoltVarWatt(
        volt_var=curves.volt_var, conflict_policy=:dominant))
    inverter = AdvancedInverter(
        id="dominant", bus="poc", phase_terminals=["a", "b", "c"],
        neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, r_filter=0.05, x_filter=0.15)
    request = InverterControlRequest(p_available=2e3, q_scale=4e3)
    result = solve_controlled_inverter(
        inv_grid3_src(
            mags=[245.0, 212.0, 230.0], angs=[0.05, -2.15, 2.0]),
        ControlledDevice(inverter, controller), request;
        per_unit=true,
        solver_options=("max_iter" => 400, "tol" => 1e-8))

    @test result.solve.publishable
    @test result.control.q_request > 0.0
    @test result.exact_control.q_request > 0.0
    @test result.exact_smooth_current_residual < 0.01
    @test maximum(result.plant.i_mag) <= inverter.i_max + 1e-6
end

@testset "Inverter controls: declared power-voltage floor residual" begin
    voltage = 230.0
    measurement = InverterControlMeasurement(
        _ctrl_phases(voltage + 0im))
    policy = AverageVoltageVoltVarWatt()
    residuals = Float64[]
    for floor in (0.1, 1.0, 10.0)
        controller = SequenceController(
            policy; power_voltage_floor=floor)
        result = evaluate_exact(
            controller, measurement, _CTRL_REQUEST, _CTRL_RATINGS)
        expected = result.p_request * floor^2/(voltage^2 + floor^2)
        residual = result.p_request - real(result.sequence_power)
        @test residual ≈ expected atol=1e-10
        push!(residuals, residual)
    end
    @test issorted(residuals)
    @test residuals[2]/_CTRL_REQUEST.p_available < 2e-5
end

@testset "Inverter controls: exact limiter, ripple blend, and frame invariance" begin
    measurement = InverterControlMeasurement(
        _ctrl_phases(230.0 + 0im, 10cis(0.4)))
    limited = evaluate_exact(
        _study_controller(), measurement, _CTRL_REQUEST,
        InverterControlRatings(s_max=20e3, i_max=5.0))
    @test limited.current_scale < 1.0
    @test maximum(abs, limited.phase_current) ≈ 5.0 atol=1e-10
    @test abs(sum(limited.phase_current)) < 1e-12

    ripple = evaluate_exact(
        _study_controller(ripple_blend=1.0), measurement,
        _CTRL_REQUEST, _CTRL_RATINGS)
    @test abs(ripple.ripple_power) < 0.1

    angle = 0.37
    rotated_measurement = InverterControlMeasurement(
        collect(measurement.phase_voltage .* cis(angle)))
    original = evaluate_exact(
        _study_controller(), measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    rotated = evaluate_exact(
        _study_controller(), rotated_measurement, _CTRL_REQUEST, _CTRL_RATINGS)
    @test rotated.p_request ≈ original.p_request atol=1e-10
    @test rotated.q_request ≈ original.q_request atol=1e-10
    @test rotated.current_scale ≈ original.current_scale atol=1e-12
    @test collect(rotated.phase_current) ≈
          collect(original.phase_current .* cis(angle)) atol=1e-10
end

@testset "Inverter controls: zero-sequence reference sensitivity" begin
    # A common-mode offset on all three phase phasors moves U0 only, so a
    # sequence-referred law must be invariant to it while a phase-magnitude law
    # must not be. This pins the sensing-reference contract: with
    # `neutral=nothing` the POC measurement is referred to ground, not to a
    # neutral conductor, and therefore carries the local zero-sequence
    # displacement that a three-leg bridge cannot control.
    u1 = 243cis(0.1)
    u2 = 4cis(-0.6)
    offset = 3cis(0.9)
    base = InverterControlMeasurement(_ctrl_phases(u1, u2))
    shifted = InverterControlMeasurement(
        collect(base.phase_voltage) .+ offset)

    volt_var = PiecewiseLinearLaw(
        [210.0, 220.0, 240.0, 250.0], [0.3, 0.0, 0.0, -0.3];
        smoothing_epsilon=0.05)
    negative_gain = PiecewiseLinearLaw(
        [0.0, 0.01, 0.10], [0.0, 0.0, 0.08]; smoothing_epsilon=1e-4)
    unbalance = NegativeSequenceAdmittanceDroop(
        negative_gain; impedance_angle=0.0, ripple_blend=0.0)

    sequence = SequenceController(
        PositiveSequenceVoltVarWatt(
            volt_var=volt_var, worst_phase_watt_guard=false),
        unbalance, CommonScaleLimiter())
    seq_base = evaluate_exact(
        sequence, base, _CTRL_REQUEST, _CTRL_RATINGS)
    seq_shifted = evaluate_exact(
        sequence, shifted, _CTRL_REQUEST, _CTRL_RATINGS)

    # Only the zero-sequence channel absorbs the offset.
    @test seq_shifted.voltage_sequence[1] ≈
          seq_base.voltage_sequence[1] + offset atol=1e-10
    @test seq_shifted.voltage_sequence[2] ≈
          seq_base.voltage_sequence[2] atol=1e-10
    @test seq_shifted.voltage_sequence[3] ≈
          seq_base.voltage_sequence[3] atol=1e-10

    # Sequence-referred decisions and the resulting command are invariant.
    @test seq_shifted.p_request ≈ seq_base.p_request atol=1e-9
    @test seq_shifted.q_request ≈ seq_base.q_request atol=1e-9
    @test seq_shifted.eta ≈ seq_base.eta atol=1e-12
    @test collect(seq_shifted.phase_current) ≈
          collect(seq_base.phase_current) atol=1e-9

    # The reported extrema are physical measured magnitudes, so they do move.
    @test !isapprox(seq_shifted.voltage_max, seq_base.voltage_max; atol=1e-2)

    # The recommended phase-magnitude law is reference-dependent by construction.
    worst = SequenceController(
        WorstPhaseVoltVarWatt(volt_var=volt_var, conflict_policy=:dominant),
        unbalance, CommonScaleLimiter())
    worst_base = evaluate_exact(
        worst, base, _CTRL_REQUEST, _CTRL_RATINGS)
    worst_shifted = evaluate_exact(
        worst, shifted, _CTRL_REQUEST, _CTRL_RATINGS)
    @test abs(worst_shifted.q_request - worst_base.q_request) > 100.0
end

@testset "Inverter controls: capability backoff cannot reverse the command" begin
    # An offset magnitude that already exceeds the per-leg limit leaves no
    # headroom. The exact allocator returns exactly zero; the smooth allocator
    # must agree to within its declared selector width and must never return a
    # negative factor, which would reverse the commanded current.
    epsilon = 1e-3
    limit = 5.0
    for (offset, direction) in (
            (6.0 + 0im, 4.0 + 0im),      # no headroom at all
            (5.0 + 0im, 4.0 + 0im),      # exactly no headroom
            (0.0 + 0im, 4.0 + 0im),      # unconstrained
            (1.0 + 0im, 9.0 + 0im))      # ordinary saturation
        model = Model(Ipopt.Optimizer)
        set_silent(model)
        smooth = PowerOptLab._safe_direction_scale_implicit!(
            model, (real(offset), imag(offset)),
            (real(direction), imag(direction)), limit, epsilon;
            magnitude_start=limit, scale=limit)
        @objective(model, Min, 0.0)
        optimize!(model)
        @test is_solved_and_feasible(model; allow_almost=true)
        exact = PowerOptLab._safe_direction_scale_exact(
            offset, direction, limit)
        @test 0.0 <= value(smooth) <= 1.0
        @test value(smooth) ≈ exact atol=epsilon
    end
end

@testset "Inverter controls: closed-form sequence and limiter identities" begin
    u1 = 228cis(0.2)
    u2 = 13cis(-0.4)
    measurement = InverterControlMeasurement(_ctrl_phases(u1, u2))
    ratings = InverterControlRatings(s_max=5e3, i_max=4.0)
    controller = _study_controller(ripple_blend=1.0)
    result = evaluate_exact(controller, measurement, _CTRL_REQUEST, ratings)

    @test result.voltage_sequence[1] ≈ 0.0im atol=1e-12
    @test result.voltage_sequence[2] ≈ u1 atol=1e-12
    @test result.voltage_sequence[3] ≈ u2 atol=1e-12

    # Exact apparent-power and conductor-current common scales.
    p_fraction = clamp(
        1.0 - 0.8*(result.voltage_max-240.0)/10.0, 0.2, 1.0)
    p_raw = _CTRL_REQUEST.p_available*p_fraction
    q_raw = _CTRL_REQUEST.q_scale * (
        clamp(0.3*(220.0-result.voltage_min)/10.0, 0.0, 0.3) +
        clamp(-0.3*(result.voltage_max-240.0)/10.0, -0.3, 0.0))
    expected_power_scale = min(1.0, ratings.s_max/hypot(p_raw, q_raw))
    @test result.power_scale ≈ expected_power_scale atol=1e-12
    phase_request = _ctrl_phases(result.i1_request, result.i2_request)
    expected_current_scale = min(
        1.0, ratings.i_max/maximum(abs, phase_request))
    @test result.current_scale ≈ expected_current_scale atol=1e-12
    @test maximum(abs, result.phase_current) ≈ ratings.i_max atol=1e-10

    # Fortescue power identities and the exact residual introduced only by the
    # declared voltage_floor regularization in the ripple-cancelling target.
    @test result.sequence_power ≈ 3u1*conj(result.i1) atol=1e-10
    @test result.total_power ≈
          3(u1*conj(result.i1) + u2*conj(result.i2)) atol=1e-10
    floor = controller.unbalance.voltage_floor
    expected_ripple = 3u2*result.i1 * floor^2/(abs2(u1) + floor^2)
    @test result.ripple_power ≈ expected_ripple atol=1e-10

    # The voltage-oriented endpoint has a closed-form magnitude and angle.
    voltage_only = evaluate_exact(
        _study_controller(ripple_blend=0.0), measurement,
        _CTRL_REQUEST, _CTRL_RATINGS)
    eta = abs(u2)/sqrt(abs2(u1) + floor^2)
    kappa = (eta-0.01)/(0.10-0.01) * 0.08
    expected_i2v = -kappa*u2
    @test voltage_only.i2_voltage ≈ expected_i2v atol=1e-12

    # A cyclic relabelling is a coordinate change, not a physical change.
    cyclic_measurement = InverterControlMeasurement(
        collect(measurement.phase_voltage)[[2, 3, 1]])
    cyclic = evaluate_exact(
        controller, cyclic_measurement, _CTRL_REQUEST, ratings)
    @test cyclic.p_request ≈ result.p_request atol=1e-10
    @test cyclic.q_request ≈ result.q_request atol=1e-10
    @test collect(cyclic.phase_current) ≈
          collect(result.phase_current)[[2, 3, 1]] atol=1e-10
end

@testset "Inverter controls: smooth AdvancedInverter composition" begin
    inverter = AdvancedInverter(
        id="controlled", bus="poc", phase_terminals=["a", "b", "c"],
        neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, r_filter=0.05, x_filter=0.15,
        m_max=0.96)
    controlled = ControlledDevice(inverter, _study_controller())

    @test_throws ArgumentError solve_controlled_inverter(
        inv_grid3_unbal(), controlled, _CTRL_REQUEST; per_unit=false)
    result = solve_controlled_inverter(
        inv_grid3_unbal(), controlled, _CTRL_REQUEST;
        per_unit=true, s_base=1e6,
        solver_options=("max_iter" => 500, "tol" => 1e-8))

    @test result.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test result.solve.publishable
    @test abs(sum(result.control.phase_current)) < 1e-7
    @test maximum(abs, result.control.phase_current) <= 40.0 + 1e-5
    @test result.plant.i_mag ≈
          collect(abs.(result.control.phase_current)) atol=1e-6
    @test result.plant.i_positive ≈ abs(result.control.i1) atol=1e-6
    @test abs(result.control.i2) > 0.1
    @test result.plant.i_negative ≈ abs(result.control.i2) atol=1e-6
    terminal = result.converter_terminal
    @test real(terminal.total_power) ≈ result.plant.p_conv atol=1e-6
    @test imag(terminal.total_power) ≈ result.plant.q_conv atol=1e-6
    @test terminal.total_power ≈
          terminal.zero_sequence_power +
          terminal.positive_sequence_power +
          terminal.negative_sequence_power atol=1e-6
    @test abs(terminal.ripple_power) ≈ result.plant.ripple atol=1e-6
    @test result.exact_smooth_current_residual <=
          0.4*_study_controller().positive.volt_watt.smoothing_epsilon
    @test solve_diagnostics(result).exact_smooth_current_residual ==
          result.exact_smooth_current_residual
end

@testset "Inverter controls: converter apparent-power saturation" begin
    inverter = AdvancedInverter(
        id="saturation", bus="poc", phase_terminals=["a", "b", "c"],
        neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, r_filter=0.05, x_filter=0.15)
    controller = SequenceController(AverageVoltageVoltVarWatt())
    request = InverterControlRequest(
        p_available=22e3, p_rated=22e3, q_scale=0.0)
    result = solve_controlled_inverter(
        inv_grid3_bal(), ControlledDevice(inverter, controller), request;
        per_unit=true,
        solver_options=("max_iter" => 500, "tol" => 1e-8))
    @test result.solve.publishable
    @test result.control.power_scale < 1.0
    @test result.control.current_scale < 1.0
    @test abs(result.converter_terminal.total_power) <=
          inverter.s_max + 0.1
    @test result.plant.p_poc < inverter.s_max
    @test result.exact_smooth_current_residual < 1e-3
end

@testset "Inverter controls: DC-ripple voltage saturation" begin
    # dv2_max sets the ripple budget |S_tilde| = 2*w*c_dc*v_dc*dv2_max, so
    # 0.1 V allows only 48 VA here and forces an ~87% current backoff. That
    # drives the solve into a near-zero-current corner where the capability
    # limiter's implicit square roots are worst conditioned, and the result
    # then depends on the solver tolerance and the platform's linear algebra.
    # 0.7 V keeps the backoff clearly active (current_scale ~ 0.88) while
    # solving at tol = 1e-6, 1e-8, and 1e-10.
    inverter = AdvancedInverter(
        id="ripple_limit", bus="poc", phase_terminals=["a", "b", "c"],
        neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, dv2_max=0.7,
        r_filter=0.05, x_filter=0.15)
    controller = SequenceController(AverageVoltageVoltVarWatt())
    result = solve_controlled_inverter(
        inv_grid3_unbal(), ControlledDevice(inverter, controller),
        _CTRL_REQUEST; per_unit=true,
        solver_options=("max_iter" => 500, "tol" => 1e-8))
    @test result.solve.publishable
    @test result.plant.dv2 <= inverter.dv2_max + 1e-6
    @test result.control.current_scale < 1.0
    @test result.exact_smooth_current_residual < 1e-3
end

@testset "Inverter controls: converter-side and grid-side current targets" begin
    inverter = AdvancedInverter(
        id="target", bus="poc", phase_terminals=["a", "b", "c"],
        neutral="n", topology=:THREE_LEG, s_max=20e3,
        i_max=40.0, i_grid_max=35.0,
        v_dc=750.0, c_dc=1.1e-3,
        r_filter=0.02, x_filter=0.06,
        r_filter_grid=0.03, x_filter_grid=0.09,
        c_filter_mid=30e-6, r_filter_damping=0.5)
    policy = AverageVoltageVoltVarWatt()
    converter_controller = SequenceController(
        policy; current_target=:converter)
    grid_controller = SequenceController(
        PositiveSequenceVoltVarWatt(); current_target=:grid)
    request = InverterControlRequest(p_available=9e3, q_scale=0.0)

    converter = solve_controlled_inverter(
        inv_grid3_bal(), ControlledDevice(inverter, converter_controller),
        request; per_unit=true,
        solver_options=("max_iter" => 400, "tol" => 1e-8))
    grid = solve_controlled_inverter(
        inv_grid3_bal(), ControlledDevice(inverter, grid_controller),
        request; per_unit=true,
        solver_options=("max_iter" => 400, "tol" => 1e-8))

    for result in (converter, grid)
        @test result.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
        @test real(result.converter_terminal.total_power) ≈
              result.plant.p_conv atol=1e-6
        @test imag(result.converter_terminal.total_power) ≈
              result.plant.q_conv atol=1e-6
        @test maximum(result.plant.i_mag) <= inverter.i_max + 1e-6
        @test maximum(result.plant.i_grid_mag) <=
              inverter.i_grid_max + 1e-6
    end
    @test converter.plant.i_mag ≈
          collect(abs.(converter.control.phase_current)) atol=1e-6
    @test collect(converter.converter_terminal.phase_current) ≈
          collect(converter.control.phase_current) atol=1e-6
    @test grid.plant.i_grid_mag ≈
          collect(abs.(grid.control.phase_current)) atol=5e-6
    @test collect(grid.grid_phase_current) ≈
          collect(grid.control.phase_current) atol=5e-6
    converter_numeric = evaluate_smooth(
        converter_controller,
        InverterControlMeasurement(collect(converter.control.phase_voltage)),
        request, InverterControlRatings(s_max=20e3, i_max=40.0))
    grid_numeric = evaluate_smooth(
        grid_controller,
        InverterControlMeasurement(collect(grid.control.phase_voltage)),
        request, InverterControlRatings(s_max=20e3, i_max=35.0))
    @test converter.control.voltage_min ≈
          converter_numeric.voltage_min atol=1e-6
    @test converter.control.voltage_max ≈
          converter_numeric.voltage_max atol=1e-6
    @test grid.control.voltage_min ≈ grid_numeric.voltage_min atol=1e-6
    @test grid.control.voltage_max ≈ grid_numeric.voltage_max atol=1e-6
    @test maximum(abs.(converter.plant.i_mag .-
                       converter.plant.i_grid_mag)) > 0.01
    @test maximum(abs.(grid.plant.i_mag .- grid.plant.i_grid_mag)) > 0.01
end


@testset "Inverter controls: controlled-power-flow objective invariance" begin
    inverter = AdvancedInverter(
        id="objective", bus="poc", phase_terminals=["a", "b", "c"],
        neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40.0,
        v_dc=700.0, c_dc=1.1e-3, r_filter=0.05, x_filter=0.15)
    controlled = ControlledDevice(
        inverter, SequenceController(AverageVoltageVoltVarWatt()))
    request = InverterControlRequest(p_available=8e3, q_scale=0.0)
    loss = solve_controlled_inverter(
        inv_grid3_bal(), controlled, request; selection_objective=:loss)
    zero = solve_controlled_inverter(
        inv_grid3_bal(), controlled, request; selection_objective=:zero)
    @test loss.solve.publishable && zero.solve.publishable
    @test collect(loss.control.phase_current) ≈
          collect(zero.control.phase_current) atol=1e-5
    @test loss.control.p_request ≈ zero.control.p_request atol=1e-4
end


function _ods_balanced_volt_watt(vsrc; irradiance=1.0)
    rated = 230.0
    vwx = join([0.5, 240.0/rated, 250.0/rated, 2.0], ",")
    cmds = String[
        "Clear",
        "Set DefaultBaseFreq=50",
        "New Circuit.c basekv=$(vsrc*sqrt(3)/1000) pu=1 angle=0 " *
            "phases=3 bus1=src.1.2.3 model=ideal",
        "New Reactor.gs bus1=src.4 bus2=src.0 R=0.001 X=0",
        "New Linecode.lc nphases=4 units=km " *
            "Rmatrix=[0.05|0 0.05|0 0 0.05|0 0 0 0.05] " *
            "Xmatrix=[0|0 0|0 0 0|0 0 0 0] Cmatrix=[0|0 0|0 0 0|0 0 0 0]",
        "New Line.l1 bus1=src.1.2.3.4 bus2=lb.1.2.3.4 " *
            "linecode=lc length=1 units=km",
        "New Reactor.gl bus1=lb.4 bus2=lb.0 R=0.001 X=0",
        "New PVSystem.pv phases=3 conn=wye bus1=lb.1.2.3.4 " *
            "kv=$(rated*sqrt(3)/1000) kVA=20 Pmpp=12 irradiance=$irradiance pf=1 " *
            "%cutin=0 %cutout=0 VarFollowInverter=yes WattPriority=No " *
            "Vminpu=0.1 Vmaxpu=2.0",
        "New XYcurve.vw npts=4 Xarray=[$vwx] Yarray=[1,1,0.2,0.2]",
        "New InvControl.ic mode=VOLTWATT voltwatt_curve=vw " *
            "voltage_curvex_ref=rated voltwattyaxis=PMPPPU " *
            "MonVoltageCalc=MAX ActivePChangeTolerance=0.0002 " *
            "VoltageChangeTolerance=1e-6",
        "Set VoltageBases=[$(rated*sqrt(3)/1000)]",
        "CalcVoltageBases",
        "Set Tolerance=1e-9",
        "Set MaxIterations=100",
        "Set ControlMode=Static",
        "Set MaxControlIter=100",
        "Solve",
    ]
    OpenDSSDirect.dss("Clear")
    foreach(OpenDSSDirect.dss, cmds)
    names = lowercase.(OpenDSSDirect.Circuit.AllNodeNames())
    volts = Dict(zip(names, OpenDSSDirect.Circuit.AllBusVolts()))
    OpenDSSDirect.Circuit.SetActiveElement("PVSystem.pv")
    powers = OpenDSSDirect.CktElement.Powers()
    p = -sum(real(powers[k]) for k in 1:3) * 1e3
    return volts, p
end

@testset "Inverter controls: balanced Volt-watt fixed point vs OpenDSS" begin
    if !isdefined(@__MODULE__, :_HAS_ODS) || !_HAS_ODS
        @test_skip "Requires OpenDSSDirect"
    else
        law = PiecewiseLinearLaw(
            [240.0, 250.0], [1.0, 0.2]; smoothing_epsilon=0.01)
        controller = SequenceController(
            WorstPhaseVoltVarWatt(volt_watt=law))
        inverter = AdvancedInverter(
            id="ods", bus="poc", phase_terminals=["a", "b", "c"],
            neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40.0,
            v_dc=700.0, c_dc=1.1e-3)
        controlled = ControlledDevice(inverter, controller)
        request = InverterControlRequest(p_available=12e3, q_scale=0.0)

        for vsrc in (245.0, 255.0)
            ods_voltage, ods_p = _ods_balanced_volt_watt(vsrc)
            result = solve_controlled_inverter(
                inv_grid3_bal(vsrc), controlled, request;
                per_unit=true,
                solver_options=("max_iter" => 300, "tol" => 1e-8))
            @test result.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
            @test result.plant.p_poc ≈ ods_p atol=50.0
            for (phase, node) in zip(("a", "b", "c"), ("1", "2", "3"))
                @test result.bus["poc"][phase]["vm"] ≈
                      abs(ods_voltage["lb.$node"]-ods_voltage["lb.4"]) atol=0.1
            end
        end

        # PMPPPU is a rated-power curve basis: at partial irradiance the
        # availability cap applies after the curve, rather than multiplying
        # availability by the curve ordinate a second time.
        irradiance = 0.15
        vsrc = 245.0
        ods_voltage, ods_p = _ods_balanced_volt_watt(
            vsrc; irradiance=irradiance)
        rated_controller = SequenceController(WorstPhaseVoltVarWatt(
            volt_watt=law, volt_watt_basis=:rated))
        rated_request = InverterControlRequest(
            p_available=12e3*irradiance, p_rated=12e3, q_scale=0.0)
        rated_result = solve_controlled_inverter(
            inv_grid3_bal(vsrc), ControlledDevice(inverter, rated_controller),
            rated_request; per_unit=true,
            solver_options=("max_iter" => 300, "tol" => 1e-8))
        @test rated_result.solve.publishable
        @test rated_result.plant.p_poc ≈ ods_p atol=50.0
        for (phase, node) in zip(("a", "b", "c"), ("1", "2", "3"))
            @test rated_result.bus["poc"][phase]["vm"] ≈
                  abs(ods_voltage["lb.$node"]-ods_voltage["lb.4"]) atol=0.1
        end
    end
end
