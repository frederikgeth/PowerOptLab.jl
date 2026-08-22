using LinearAlgebra

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
        @test screen.maximum_real_eigenvalue ≈ 0.4 atol=1e-8
        @test screen.continuous_time_margin ≈ 0.6 atol=1e-8
        @test screen.alpha_max ≈ 2 / 0.8 atol=1e-8

        negative_gain = screen_fixed_point_gain(
            x -> -2.0 .* x, [1.0]; step=1e-6)
        @test !negative_gain.local_contractive
        @test negative_gain.maximum_real_eigenvalue ≈ -2.0 atol=1e-8
        @test negative_gain.continuous_time_margin ≈ 3.0 atol=1e-8
        @test negative_gain.alpha_max ≈ 2 / 3 atol=1e-8

        rectangular(x) = [x[1] + 2x[2], x[1] - x[2], 3x[2]]
        rectangular_jacobian = finite_difference_jacobian(
            rectangular, [2.0, -1.0]; step=1e-6)
        @test size(rectangular_jacobian) == (3, 2)
        @test rectangular_jacobian ≈ [1.0 2.0; 1.0 -1.0; 0.0 3.0] atol=1e-8
        @test_throws DimensionMismatch screen_fixed_point_gain(
            rectangular, [2.0, -1.0]; step=1e-6)
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

        volt_watt = PiecewiseLinearLaw(
            [200.0, 220.0, 240.0, 260.0], [1.0, 1.0, 0.5, 0.0];
            smoothing_epsilon=0.05)
        volt_var = PiecewiseLinearLaw(
            [200.0, 220.0, 240.0, 260.0], [0.4, 0.2, -0.2, -0.4];
            smoothing_epsilon=0.05)
        curved = ControlledDevice(
            inverter,
            SequenceController(AverageVoltageVoltVarWatt(
                volt_watt=volt_watt, volt_var=volt_var)))
        curved_result = solve_inverter_control_network_fixed_point(
            inv_grid3_bal(), curved, request;
            max_iterations=4,
            solver_options=("max_iter" => 2000, "tol" => 1e-7))
        @test curved_result.smooth.solve.publishable
        @test curved_result.solve.publishable
        @test curved_result.converged || curved_result.cycled
        @test isfinite(curved_result.exact_control.p_request)

        floating_result = solve_inverter_control_network_fixed_point(
            inv_grid3_floating_neutral(), curved, request;
            max_iterations=4,
            solver_options=("max_iter" => 2000, "tol" => 1e-7))
        @test floating_result.smooth.solve.publishable
        @test floating_result.solve.publishable
        @test floating_result.converged || floating_result.cycled
        neutral = floating_result.bus["poc"]["n"]
        @test hypot(neutral["vr"], neutral["vi"]) > 1e-4
        sensed = ntuple(3) do k
            phase = ("a", "b", "c")[k]
            floating_result.bus["poc"][phase]["vr"] +
                im*floating_result.bus["poc"][phase]["vi"] -
                (neutral["vr"] + im*neutral["vi"])
        end
        @test maximum(abs, sensed .-
            floating_result.exact_control.phase_voltage) < 1e-6

        failed = solve_inverter_control_network_fixed_point(
            inv_grid3_bal(), controlled, request;
            max_iterations=2, solver_options=("max_iter" => 0,))
        @test !failed.converged
        @test !failed.solve.publishable
        @test size(failed.voltage_trajectory) == (6, 0)
    end

    @testset "simultaneous fleet exact equilibrium" begin
        inverter = AdvancedInverter(
            id="pv", bus="poc", phase_terminals=["a", "b", "c"],
            neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40.0,
            v_dc=700.0, c_dc=1.1e-3, r_filter=0.05, x_filter=0.15,
            m_max=0.96)
        controlled = ControlledDevice(
            inverter, SequenceController(AverageVoltageVoltVarWatt()))
        spec = ControlledInverterFleetSpec(
            Dict("pv" => controlled),
            Dict("pv" => InverterControlRequest(
                p_available=5e3, p_rated=5e3, q_scale=0.0)))
        net = inv_grid3_bal()
        net["ibr"] = Dict("pv" => Dict{String,Any}(
            "bus" => "poc", "terminal_map" => ["a", "b", "c", "n"],
            "topology" => "FOUR_LEG", "prime_mover" => "PV",
            "s_max" => fill(20e3, 3), "p_min" => zeros(3),
            "p_max" => fill(20e3, 3), "q_min" => zeros(3),
            "q_max" => zeros(3)))
        result = solve_controlled_inverter_fleet_network_fixed_point(
            net, spec; max_iterations=6,
            solver_options=("max_iter" => 500, "tol" => 1e-8))
        @test result.smooth.solve.publishable
        @test result.solve.publishable
        @test result.converged
        @test !result.cycled
        @test result.cycle_period == 0
        @test result.residual_voltage_inf_V < 1e-3
        @test result.residual_current_inf_A < 1e-3
        rows = controlled_inverter_network_fixed_point_rows(result)
        @test length(rows) == 1
        @test rows[1].device_id == "pv"
        @test rows[1].converged
        @test rows[1].exact_p_poc_W ≈ result.plants["pv"].p_poc atol=1e-8

        starts = Dict(
            "high" => Dict("pv" => ComplexF64[
                232.0 + 0im, 232.0cis(-2pi/3), 232.0cis(2pi/3)]),
            "low" => Dict("pv" => ComplexF64[
                228.0 + 0im, 228.0cis(-2pi/3), 228.0cis(2pi/3)]))
        multi = solve_controlled_inverter_fleet_multistart(
            net, spec, starts; smooth_result=result.smooth, max_iterations=6,
            solver_options=("max_iter" => 500, "tol" => 1e-8))
        @test multi.start_ids == ["high", "low"]
        @test length(multi.results) == 2
        multi_rows = controlled_inverter_fleet_multistart_rows(multi)
        @test [row.start_id for row in multi_rows] == ["high", "low"]
        @test all(row.converged for row in multi_rows)
        @test all(isfinite(row.residual_voltage_inf_V) for row in multi_rows)
        @test all(isfinite(
            row.maximum_pairwise_equilibrium_voltage_spread_inf_V)
            for row in multi_rows)

        @test_throws ArgumentError inverter_control_network_voltage_sensitivity(
            net, controlled, zeros(3); step=0.0)
        sensitivity_failure = inverter_control_network_voltage_sensitivity(
            net, controlled, zeros(3); step=1e-3,
            solver_options=("max_iter" => 0,))
        @test !sensitivity_failure.base_status.publishable
        @test all(isnan, sensitivity_failure.sensitivity)
        @test isempty(sensitivity_failure.plus_status)
        @test isempty(sensitivity_failure.minus_status)

        sensitivity_success = inverter_control_network_voltage_sensitivity(
            inv_grid3_bal(), controlled,
            result.smooth.devices["pv"].control.phase_current;
            step=1e-3,
            solver_options=("max_iter" => 500, "tol" => 1e-8))
        @test sensitivity_success.coordinate_system ==
              :positive_negative_sequence
        @test size(sensitivity_success.sensitivity) == (4, 4)
        @test sensitivity_success.base_status.publishable
        @test all(s -> s.publishable, sensitivity_success.plus_status)
        @test all(s -> s.publishable, sensitivity_success.minus_status)
        @test all(isfinite, sensitivity_success.sensitivity)
        @test sensitivity_success.sensitivity ≈
              0.05 .* Matrix{Float64}(I, 4, 4) atol=1e-5
        sequence_screen = inverter_control_loop_gain(
            controlled.controller,
            InverterControlMeasurement(
                result.smooth.devices["pv"].control.phase_voltage),
            spec.requests["pv"],
            InverterControlRatings(inverter, controlled.controller.current_target),
            sensitivity_success;
            step=1e-5)
        @test size(sequence_screen.jacobian) == (4, 4)
        @test all(isfinite, sequence_screen.jacobian)

        fleet_sensitivity_failure =
            controlled_inverter_fleet_network_voltage_sensitivity(
                net, spec, Dict("pv" => zeros(3)); step=1e-3,
                solver_options=("max_iter" => 0,))
        @test fleet_sensitivity_failure.device_ids == ["pv"]
        @test fleet_sensitivity_failure.coordinate_system == :positive_negative_sequence
        @test size(fleet_sensitivity_failure.sensitivity) == (4, 4)
        @test all(isnan, fleet_sensitivity_failure.sensitivity)
        @test isempty(fleet_sensitivity_failure.plus_status)
        @test isempty(fleet_sensitivity_failure.minus_status)
        sensitivity_rows = controlled_inverter_fleet_network_sensitivity_rows(
            fleet_sensitivity_failure)
        @test length(sensitivity_rows) == 16
        @test all(isnan, getfield.(sensitivity_rows,
            :voltage_sensitivity_V_per_A))
        @test all(!row.plus_publishable && !row.minus_publishable
            for row in sensitivity_rows)

        diagonal_sensitivity = zeros(4, 4)
        for k in 1:4
            diagonal_sensitivity[k, k] = 0.01
        end
        fleet_screen = controlled_inverter_fleet_loop_gain(
            spec,
            Dict("pv" => result.smooth.devices["pv"].control.phase_voltage),
            diagonal_sensitivity; step=1e-5)
        @test all(isfinite, fleet_screen.jacobian)
        @test fleet_screen.spectral_radius >= 0.0
        @test fleet_screen.local_contractive
    end
end
