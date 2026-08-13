@testset "Closed-loop evidence diagnostics" begin
    @testset "fixed-point oracle and local gain screen" begin
        A = [0.2 0.0; 0.0 0.4]
        b = [1.0, -1.0]
        affine(x) = A*x + b
        result = fixed_point_oracle(
            affine, zeros(2); max_iterations=100, atol=1e-10, rtol=1e-10)
        @test result.converged
        @test !result.cycled
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
end
