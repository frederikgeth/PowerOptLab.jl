@testset "Closed-loop evidence diagnostics" begin
    @testset "fixed-point oracle and local gain screen" begin
        A = [0.2 0.0; 0.0 0.4]
        b = [1.0, -1.0]
        affine(x) = A*x + b
        result = fixed_point_oracle(
            affine, zeros(2); max_iterations=100, atol=1e-10, rtol=1e-10)
        @test result.converged
        @test !result.cycled
        @test result.cycle_period == 0
        @test result.iterations < 50
        @test result.final_state ≈ [1.25, -1.6666666666666667] atol=1e-8
        @test size(result.trajectory) == (2, result.iterations + 1)

        screen = screen_fixed_point_gain(
            affine, result.final_state; step=1e-6)
        @test screen.jacobian ≈ A atol=1e-8
        @test screen.spectral_radius ≈ 0.4 atol=1e-8
        @test screen.induced_inf_norm ≈ 0.4 atol=1e-8
        @test screen.local_contractive
        @test screen.margin ≈ 0.6 atol=1e-8
    end

    @testset "cycle detection" begin
        result = fixed_point_oracle(
            x -> -x, [1.0]; max_iterations=20, atol=0.0, rtol=0.0,
            cycle_window=4, cycle_tolerance=1e-12)
        @test !result.converged
        @test result.cycled
        @test result.cycle_period == 2
        @test result.iterations == 2

        no_trajectory = fixed_point_oracle(
            x -> -x, [1.0]; max_iterations=20, atol=0.0, rtol=0.0,
            cycle_window=4, cycle_tolerance=1e-12,
            store_trajectory=false)
        @test no_trajectory.cycled
        @test size(no_trajectory.trajectory) == (1, 0)
    end

    @testset "sequence-controller loop adapter" begin
        controller = SequenceController(AverageVoltageVoltVarWatt())
        measurement = InverterControlMeasurement(ComplexF64[
            230.0 + 0im, 230.0cis(-2pi/3), 230.0cis(2pi/3)])
        request = InverterControlRequest(
            p_available=5e3, p_rated=5e3, q_scale=2e3)
        ratings = InverterControlRatings(s_max=10e3, i_max=20.0)
        jacobian = inverter_control_current_jacobian(
            controller, measurement, request, ratings; step=1e-5)
        @test size(jacobian) == (6, 6)
        @test all(isfinite, jacobian)

        zero_sensitivity = zeros(6, 6)
        screen = inverter_control_loop_gain(
            controller, measurement, request, ratings, zero_sensitivity;
            step=1e-5)
        @test screen.jacobian ≈ zeros(6, 6) atol=1e-10
        @test screen.spectral_radius == 0.0
        @test screen.local_contractive
    end

    @testset "scaling audit and exact equilibrium adapter" begin
        controller = SequenceController(AverageVoltageVoltVarWatt())
        ratings = InverterControlRatings(s_max=20e3, i_max=40.0)
        audit = inverter_control_scaling_audit(controller, ratings)
        @test audit.voltage_start_V == 230.0
        @test audit.voltage_scale_V == 230.0
        @test audit.sequence_voltage_start_V == 11.5
        @test audit.power_start_VA == ratings.s_max
        @test audit.power_scale_VA == ratings.s_max
        @test audit.current_scale_A == ratings.i_max
        @test audit.capacity_aux_start_VA ≈ ratings.s_max *
              sqrt(2e-3 - 1e-6) atol=1e-10
        @test audit.power_epsilon_VA ≈
              CommonScaleLimiter().power_epsilon_fraction*ratings.s_max

        reference = InverterControlMeasurement(ComplexF64[
            230.0 + 0im, 230.0cis(-2pi/3), 230.0cis(2pi/3)])
        initial = InverterControlMeasurement(ComplexF64[
            240.0 + 0im, 240.0cis(-2pi/3), 240.0cis(2pi/3)])
        request = InverterControlRequest(
            p_available=5e3, p_rated=5e3, q_scale=2e3)
        result = inverter_control_fixed_point_oracle(
            controller, initial, request, ratings, reference, zeros(6, 6);
            max_iterations=10, atol=1e-10, rtol=1e-10)
        @test result.oracle.converged
        @test result.oracle.iterations == 2
        @test result.final_measurement.phase_voltage == reference.phase_voltage
        @test collect(result.final_control.phase_current) ≈ collect(
              evaluate_exact(controller, reference, request, ratings).phase_current)
        @test result.voltage_sensitivity == zeros(6, 6)
    end

    @testset "plant-backed exact equilibrium" begin
        inverter = AdvancedInverter(
            id="fixed_point", bus="poc", phase_terminals=["a", "b", "c"],
            neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40.0,
            v_dc=700.0, c_dc=1.1e-3, r_filter=0.05, x_filter=0.15,
            m_max=0.96)
        controlled = ControlledDevice(
            inverter, SequenceController(AverageVoltageVoltVarWatt()))
        request = InverterControlRequest(
            p_available=5e3, p_rated=5e3, q_scale=0.0)
        result = solve_inverter_control_network_fixed_point(
            inv_grid3_bal(), controlled, request;
            max_iterations=6, solver_options=("max_iter" => 500, "tol" => 1e-8))
        @test result.smooth.solve.publishable
        @test result.solve.publishable
        @test result.converged
        @test !result.cycled
        @test result.plant !== nothing
        @test result.residual_voltage_inf_V < 1e-3
        @test result.residual_current_inf_A < 1e-3
        @test size(result.voltage_trajectory, 1) == 6
        @test size(result.voltage_trajectory, 2) >= 2
        @test result.bus["poc"]["a"]["vm"] > 0

        failed = solve_inverter_control_network_fixed_point(
            inv_grid3_bal(), controlled, request;
            max_iterations=2, solver_options=("max_iter" => 0,))
        @test !failed.converged
        @test !failed.solve.publishable
        @test size(failed.voltage_trajectory) == (6, 0)
    end
end
