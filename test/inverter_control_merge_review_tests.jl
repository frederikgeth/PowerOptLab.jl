# Reproductions of the three merge-review findings on PR #61.
using LinearAlgebra: I
const _MERGE_REVIEW_A = cis(2pi/3)
_merge_review_phases(u1, u2=0im, u0=0im) = ComplexF64[
    u0 + u1 + u2,
    u0 + _MERGE_REVIEW_A^2*u1 + _MERGE_REVIEW_A*u2,
    u0 + _MERGE_REVIEW_A*u1 + _MERGE_REVIEW_A^2*u2]

_merge_review_inverter(; kwargs...) = AdvancedInverter(;
    id="pv", bus="poc", phase_terminals=["a", "b", "c"], neutral=nothing,
    topology=:THREE_LEG, s_max=6000.0, i_max=10.0, v_dc=750.0, c_dc=2e-3,
    r_filter=0.05, x_filter=0.15, kwargs...)

@testset "Merge review: stiff DC source with a ripple limit" begin
    request = InverterControlRequest(p_available=3000.0, p_rated=6000.0, q_scale=0.0)
    for target in (:converter, :grid), limiter in
            (CommonScaleLimiter(), SequencePriorityLimiter(priority=:power))
        controller = SequenceController(AverageVoltageVoltVarWatt();
            current_target=target, limiter=limiter)
        limited = solve_controlled_inverter(inv_grid3_unbal(),
            ControlledDevice(_merge_review_inverter(dc_source_r=0.0, dv2_max=2.0),
                controller), request; solver_options=("tol" => 1e-8,))
        reference = solve_controlled_inverter(inv_grid3_unbal(),
            ControlledDevice(_merge_review_inverter(dc_source_r=0.0), controller),
            request; solver_options=("tol" => 1e-8,))
        @test limited.solve.publishable && reference.solve.publishable
        @test limited.plant.dv2 ≈ 0.0 atol=1e-8
        @test collect(limited.control.phase_current) ≈ collect(reference.control.phase_current) atol=1e-5
        # This also exercises the exact plant-aware backoff, not just stamping.
        @test all(isfinite, limited.exact_control.phase_current)
        @test limited.exact_smooth_current_residual < 0.01
    end
end

@testset "Merge review: physical loss import is independent of voltage floor" begin
    request = InverterControlRequest(p_available=0.0, p_rated=6000.0, q_scale=0.0)
    measurement = InverterControlMeasurement(_merge_review_phases(230.0 + 0im))
    ratings = InverterControlRatings(s_max=6000.0, i_max=10.0)
    for floor in (1.0, 10.0, 20.0), basis in (:available, :rated)
        controller = SequenceController(AverageVoltageVoltVarWatt(volt_watt_basis=basis);
            power_voltage_floor=floor)
        exact = evaluate_exact(controller, measurement, request, ratings;
            dc_power_offset=100.0)
        smooth = evaluate_smooth(controller, measurement, request, ratings;
            dc_power_offset=100.0)
        @test real(exact.sequence_power) ≈ -100.0 atol=1e-9
        @test real(smooth.sequence_power) <= -100.0 + 1e-8
        @test maximum(abs, smooth.phase_current .- exact.phase_current) < 0.001
    end
    # A required import cannot be delivered at zero positive-sequence voltage.
    @test_throws ArgumentError evaluate_exact(
        SequenceController(AverageVoltageVoltVarWatt()),
        InverterControlMeasurement(zeros(ComplexF64, 3)), request, ratings;
        dc_power_offset=100.0)

    for floor in (1.0, 20.0), target in (:converter, :grid), limiter in
            (CommonScaleLimiter(), SequencePriorityLimiter(priority=:power))
        @testset "floor=$floor target=$target limiter=$(nameof(typeof(limiter)))" begin
            controller = SequenceController(AverageVoltageVoltVarWatt();
                power_voltage_floor=floor, current_target=target, limiter=limiter)
            # Match the existing inverter tests' solver tolerance; keep the
            # physical DC-power assertion independent of solver coordinates.
            result = solve_controlled_inverter(inv_grid3_bal(),
                ControlledDevice(_merge_review_inverter(p_loss_fixed=100.0), controller),
                request; solver_options=("tol" => 1e-8,))
            @test result.solve.publishable
            @test result.plant.p_dc <= 1e-5
            @test result.plant.p_dc >= -0.01
            @test result.exact_smooth_current_residual < 0.001
        end
    end
end

@testset "Merge review: reduced loop gain retains sensed zero-sequence voltage" begin
    vv = PiecewiseLinearLaw([207.0, 220.0, 240.0, 258.0],
        [0.44, 0.0, 0.0, -0.60]; smoothing_epsilon=0.05)
    request = InverterControlRequest(p_available=3000.0, p_rated=6000.0, q_scale=6000.0)
    inverter = _merge_review_inverter(i_max=20.0)
    z = 0.05Matrix{Float64}(I, 4, 4)
    point = [230.0, 0.0, 0.0, 0.0]
    for target in (:converter, :grid), u0 in (20.0 + 0im, -15.0 + 4im)
        controller = SequenceController(WorstPhaseVoltVarWatt(volt_var=vv);
            current_target=target, sensing_reference="n")
        measurement = InverterControlMeasurement(_merge_review_phases(230.0 + 0im, 0im, u0))
        ratings = InverterControlRatings(inverter, controller.current_target)
        # Independent inverse Fortescue reconstruction, retaining measured U0.
        function current_at(state)
            phase_voltage = _merge_review_phases(
                complex(state[1], state[2]), complex(state[3], state[4]), u0)
            result = evaluate_exact(controller, InverterControlMeasurement(phase_voltage),
                request, ratings)
            [real(result.i1), imag(result.i1), real(result.i2), imag(result.i2)]
        end
        jacobian = zeros(4, 4)
        for k in 1:4
            delta = zeros(4); delta[k] = 1e-3
            jacobian[:, k] = (current_at(point + delta) - current_at(point - delta))/0.002
        end
        @test abs(jacobian[2, 1]) > 0.01
        single = inverter_control_loop_gain(controller, measurement, request, ratings, z;
            coordinate_system=:positive_negative_sequence)
        fleet = ControlledInverterFleetSpec(
            Dict("pv" => ControlledDevice(inverter, controller)), Dict("pv" => request))
        simultaneous = controlled_inverter_fleet_loop_gain(
            fleet, Dict("pv" => measurement), z)
        @test single.jacobian ≈ z*jacobian atol=1e-7 rtol=1e-5
        @test simultaneous.jacobian ≈ z*jacobian atol=1e-7 rtol=1e-5
    end
end
