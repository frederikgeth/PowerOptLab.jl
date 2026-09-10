# Regressions for the September 2026 numerical/scientific review.
using SparseArrays: SparseMatrixCSC, nnz
using Random: MersenneTwister, randn

@testset "Compiled SE review: loaded four-wire derivatives" begin
    net = compiled_se_net()
    # Include complex series impedance and phase/neutral coupling.
    net["linecode"]["lc"]["X_series_1_1"] = 0.04
    net["linecode"]["lc"]["R_series_1_4"] = 0.01
    net["linecode"]["lc"]["X_series_1_4"] = 0.005
    ms = Any[Measurement(kind=kind, bus="b1", terminal=phase, value=0., sigma=3.)
             for kind in (:vr, :vi, :vmag, :pinj, :qinj), phase in ("1", "2", "3")][:]
    append!(ms, [BranchMeasurement(kind=kind, line="l1", side=side,
        terminal="1", value=0., sigma=2.) for kind in (:ire, :iim, :imag, :pflow, :qflow)
        for side in (:from, :to)])
    device = ExactDeviceEquation(ZIPDevice(
        [TerminalConnection(("b1", "1"), ("b1", "n"))],
        ComplexF64[4000 + 1000im], ComplexF64[2 - .5im], ComplexF64[.01 - .002im]))
    s = compile_state_estimator(net, ms; exact_devices=[device], zero_injection=[("b2", "1")])
    p = SEParameters(s, ms; exact_devices=[device])
    rng = MersenneTwister(741)
    for _ in 1:4
        x = flat_compiled_state(s) + 4randn(rng, 2length(s.free_state_map))
        H = residual_jacobian(s, p, x)
        @test H isa SparseMatrixCSC
        @test H ≈ central_jacobian(y -> evaluate_state_estimator(s, p, y).residual, x) rtol=1e-5 atol=1e-5
        C = constraint_jacobian(s, p, x)
        @test C isa SparseMatrixCSC
        @test C ≈ central_jacobian(y -> evaluate_state_estimator(s, p, y).constraints, x) rtol=1e-5 atol=1e-5
    end
end

@testset "Compiled SE review: regularised Hachtel and singular LS" begin
    for damping in (0., 0.5, 2.)
        step, _ = PowerOptLab._se_hachtel_step(ones(1, 1), zeros(0, 1), [1.], Float64[], [1.], damping)
        @test step[1] ≈ -1 / (1 + damping)
    end
    ms = [Measurement(kind=:vr, bus="b", reference=nothing, value=v, sigma=1.) for v in (990., 992.)]
    s = compile_state_estimator(constant_power_test_net(), ms); p = SEParameters(s, ms)
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator)
        r = solver(s, p, [1000., 0.])
        @test r.status == :converged_underobserved
        @test r.state[1] ≈ 991.
    end
end

@testset "Compiled SE review: exact constraints dominate conflicting meters" begin
    ms = [Measurement(kind=:vr, bus="b", reference=nothing, value=900., sigma=.1),
          Measurement(kind=:vi, bus="b", reference=nothing, value=0., sigma=1.)]
    s = compile_state_estimator(constant_power_test_net(), ms; zero_injection=["b"])
    p = SEParameters(s, ms)
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator)
        r = solver(s, p, [900., 0.])
        @test r.status == :converged_unique
        @test r.state ≈ [1000., 0.]
        @test norm(r.evaluation.constraints) < 1e-7
        @test last(r.history).penalty > 10
    end
end

@testset "Compiled SE review: reject maxima and invalid initial domains" begin
    ms = [Measurement(kind=:vmag, bus="b", reference=nothing, value=230., sigma=1.),
          Measurement(kind=:vr, bus="b", reference=nothing, value=0., sigma=100.),
          Measurement(kind=:vi, bus="b", reference=nothing, value=0., sigma=100.)]
    s = compile_state_estimator(constant_power_test_net(), ms)
    p = SEParameters(s, ms; magnitude_epsilon=.1)
    current = [BranchMeasurement(kind=:imag, line="l", terminal="1", value=5., sigma=1.)]
    sc = compile_state_estimator(constant_power_test_net(), current)
    for epsilon in (.001, 1.), solver in (solve_compiled_state_estimator, solve_sparse_state_estimator)
        pc = SEParameters(sc, current; current_epsilon=epsilon)
        @test solver(sc, pc, [1000., 0.]).status == :stationary_not_minimum
    end
    device = ExactDeviceEquation(ConstantPowerDevice(
        [TerminalConnection(("b", "1"), nothing)], ComplexF64[1000]))
    sd = compile_state_estimator(constant_power_test_net(); exact_devices=[device])
    pd = SEParameters(sd; exact_devices=[device])
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator)
        r = solver(s, p, [0., 0.])
        @test r.status == :stationary_not_minimum
        @test !solve_status(r).publishable
        bad = solver(sd, pd, [0., 0.])
        @test bad.status == :invalid_initial_domain
        @test !solve_status(bad).publishable
        @test all(isnan, bad.evaluation.constraints)
        @test all(z -> isnan(real(z)), bad.evaluation.voltage)
    end
end

@testset "Compiled SE review: nonzero residual optimum with reactive telemetry" begin
    ms = [Measurement(kind=:vr, bus="b", reference=nothing, value=990., sigma=1.),
          Measurement(kind=:vi, bus="b", reference=nothing, value=5., sigma=1.),
          Measurement(kind=:qinj, bus="b", reference=nothing, value=-45000., sigma=1000.)]
    s = compile_state_estimator(constant_power_test_net(), ms); p = SEParameters(s, ms)
    # Q = -10000*Vi, so the optimum is an independent scalar weighted mean.
    expected = [990., 455 / 101]
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator)
        r = solver(s, p, [1000., 0.])
        @test r.status == :converged_unique
        @test r.state ≈ expected atol=1e-7
        @test norm(central_jacobian(y -> evaluate_state_estimator(s,p,y).residual, r.state)' * r.evaluation.residual) < 1e-6
    end
end

@testset "Compiled SE review: identifiable covariance and prior provenance" begin
    ms = [Measurement(kind=:vr, bus="b", reference=nothing, value=990., sigma=2.)]
    s = compile_state_estimator(constant_power_test_net(), ms); p = SEParameters(s, ms)
    @test selected_state_covariance(s, p, [990., 0.], [1]) ≈ fill(4., 1, 1)
    @test_throws ArgumentError selected_state_covariance(s, p, [990., 0.], [2])
    @test derived_covariance(s, p, [990., 0.], reshape([3., 0.], 1, 2)) ≈ fill(36., 1, 1)
    pp = SEParameters(s, ms; prior=StatePrior([2], [0.], [10.]))
    @test observability_diagnostics(s, pp, [990., 0.]).unobservable_dimension == 0
    @test observability_diagnostics(s, pp, [990., 0.]; include_prior=false).unobservable_dimension == 1
    # Linear independent measurements: Monte Carlo variance validates the scale
    # without relying on the estimator's own derivative as an oracle.
    rng = MersenneTwister(882)
    draws = Float64[]
    for _ in 1:1000
        p.measurement_values[1] = 990 + 2randn(rng)
        result = solve_compiled_state_estimator(s, p, [990., 0.])
        @test result.status == :converged_underobserved
        push!(draws, result.state[1])
    end
    empirical = sum((draws .- sum(draws)/length(draws)).^2) / (length(draws)-1)
    @test empirical ≈ selected_state_covariance(s,p,[990.,0.],[1])[1,1] rtol=.15
end

@testset "Compiled SE review: sparse allocation growth" begin
    include(joinpath(@__DIR__, "..", "scripts", "benchmark_compiled_state_estimation.jl"))
    allocations = Tuple{Int,Int}[]
    for n in (50, 200)
        s, p, x = benchmark_se_chain(n)
        H = residual_jacobian(s,p,x)
        evaluate_state_estimator(s,p,x)
        @test nnz(H) == 2n
        push!(allocations, (@allocated(residual_jacobian(s,p,x)),
                            @allocated(evaluate_state_estimator(s,p,x))))
    end
    # Four times as many states should not restore the former quadratic arrays.
    @test allocations[2][1] < 6allocations[1][1]
    @test allocations[2][2] < 6allocations[1][2]
end
