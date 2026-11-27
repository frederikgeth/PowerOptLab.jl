# Regressions for the November 2026 review of the two state-estimation
# approaches: convergence at the residual noise floor, the size-independent
# rank threshold, an effective continuation homotopy, the bounded curvature
# check, statistical consistency of both estimators, and derivative coverage
# for connection types that previously had none.
using Test, PowerOptLab, LinearAlgebra, Logging
using StableRNGs
using BMOPFTools: parse_bmopf

# Sample moments without a Statistics test dependency.
sr_mean(v) = sum(v) / length(v)
sr_var(v) = sum(abs2, v .- sr_mean(v)) / (length(v) - 1)
sr_column_means(X) = vec(sum(X; dims=1)) ./ size(X, 1)
sr_column_vars(X) = [sr_var(view(X, :, j)) for j in axes(X, 2)]

# Four-wire feeder with a bonded source neutral and floating load neutrals.
function solver_regression_feeder()
    parse_bmopf("""
    {"bus":{
        "src":{"terminal_names":["1","2","3","n"],"perfectly_grounded_terminals":["n"]},
        "b1": {"terminal_names":["1","2","3","n"]},
        "b2": {"terminal_names":["1","2","3","n"]}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1","2","3"],
        "v_magnitude":[230.0,230.0,230.0],"v_angle":[0.0,-2.0943951023931953,2.0943951023931953]}},
     "linecode":{"lc":{"R_series_1_1":0.1,"R_series_2_2":0.1,"R_series_3_3":0.1,"R_series_4_4":0.1,
        "X_series_1_1":0.05,"X_series_2_2":0.05,"X_series_3_3":0.05,"X_series_4_4":0.05}},
     "line":{
        "l1":{"bus_from":"src","bus_to":"b1","terminal_map_from":["1","2","3","n"],"terminal_map_to":["1","2","3","n"],"linecode":"lc","length":1.0},
        "l2":{"bus_from":"b1","bus_to":"b2","terminal_map_from":["1","2","3","n"],"terminal_map_to":["1","2","3","n"],"linecode":"lc","length":1.0}}}
    """; from_string=true)
end

# A radial four-wire chain, long enough for cond(H*Z) to exceed the old
# dimension-scaled rank threshold.
function solver_regression_chain(n)
    tm = "[\"1\",\"2\",\"3\",\"n\"]"
    buses = ["\"src\":{\"terminal_names\":$tm,\"perfectly_grounded_terminals\":[\"n\"]}"]
    lines = String[]
    previous = "src"
    for k in 1:n
        push!(buses, "\"b$k\":{\"terminal_names\":$tm}")
        push!(lines, "\"l$k\":{\"bus_from\":\"$previous\",\"bus_to\":\"b$k\",\"terminal_map_from\":$tm,\"terminal_map_to\":$tm,\"linecode\":\"lc\",\"length\":0.05}")
        previous = "b$k"
    end
    parse_bmopf("{\"bus\":{" * join(buses, ",") * "},\"line\":{" * join(lines, ",") *
        "},\"voltage_source\":{\"s\":{\"bus\":\"src\",\"terminal_map\":[\"1\",\"2\",\"3\"]," *
        "\"v_magnitude\":[230.0,230.0,230.0],\"v_angle\":[0.0,-2.0943951023931953,2.0943951023931953]}}," *
        "\"linecode\":{\"lc\":{\"R_series_1_1\":0.1,\"R_series_2_2\":0.1,\"R_series_3_3\":0.1,\"R_series_4_4\":0.1}}}";
        from_string=true)
end

# Exact-device power flow through the compiled model, then phase-to-neutral
# |V|, P and Q telemetry at every load. Returns the truth and the meter list.
function solver_regression_problem(net, loads; sigma_v=1.0, sigma_s=50.0)
    nodes = sort!(collect(keys(loads)))
    device = ExactDeviceEquation(ConstantPowerDevice(
        [TerminalConnection(node, (node[1], "n")) for node in nodes],
        ComplexF64[loads[node] for node in nodes]))
    s0 = compile_state_estimator(net, Measurement[]; exact_devices=[device])
    p0 = SEParameters(s0, Measurement[]; exact_devices=[device])
    flow = solve_sparse_state_estimator(s0, p0, initial_state_estimator(s0, p0))
    voltage = Dict(s0.nodes[i] => flow.evaluation.voltage[i] for i in eachindex(s0.nodes))
    rows = NamedTuple[]
    for node in nodes
        drop = voltage[node] - voltage[(node[1], "n")]
        push!(rows, (kind=:vmag, node=node, value=abs(drop), sigma=sigma_v))
        push!(rows, (kind=:pinj, node=node, value=-real(loads[node]), sigma=sigma_s))
        push!(rows, (kind=:qinj, node=node, value=-imag(loads[node]), sigma=sigma_s))
    end
    meters(values) = [Measurement(kind=row.kind, bus=row.node[1], terminal=row.node[2],
                                  value=value, sigma=row.sigma) for (row, value) in zip(rows, values)]
    (; flow, voltage, meters, truth=[row.value for row in rows], sigmas=[row.sigma for row in rows])
end

function solver_regression_state(s, voltage)
    nf = length(s.free_state_map)
    x = zeros(length(s.nominal_state))
    for (node, k) in s.free_state_map
        x[k] = real(voltage[node]); x[nf + k] = imag(voltage[node])
    end
    x
end

function solver_regression_fd(f, x)
    y = f(x)
    J = zeros(length(y), length(x))
    for j in eachindex(x)
        h = 1e-6 * max(1.0, abs(x[j]))
        up = copy(x); down = copy(x)
        up[j] += h; down[j] -= h
        J[:, j] .= (f(up) .- f(down)) ./ (2h)
    end
    J
end

const SOLVER_REGRESSION_LOADS = Dict(
    ("b1", "1") => 4000 + 800im, ("b1", "2") => 1500 + 200im, ("b1", "3") => 500 + 0im,
    ("b2", "1") => 2500 + 500im, ("b2", "2") => 3500 - 300im, ("b2", "3") => 1000 + 400im)

@testset "SE solvers: convergence at the residual noise floor" begin
    problem = solver_regression_problem(solver_regression_feeder(), SOLVER_REGRESSION_LOADS)
    @test problem.flow.status == :converged_unique
    ms = problem.meters(problem.truth)
    s = compile_state_estimator(solver_regression_feeder(), ms)
    p = SEParameters(s, ms)
    xtrue = solver_regression_state(s, problem.voltage)
    @test norm(evaluate_state_estimator(s, p, xtrue).residual) < 1e-6

    reference = solve_sparse_state_estimator(s, p, initial_state_estimator(s, p))
    @test reference.status == :converged_unique
    dof = solve_diagnostics(reference).dof
    @test dof == 6
    covariance = selected_state_covariance(s, p, reference.state, collect(eachindex(xtrue)))

    # With noisy data the residual is nonzero at the optimum. Defaults used to
    # report `:trust_region_stalled` for about half of such draws.
    rng = StableRNG(20261127)
    draws = 240
    states = zeros(draws, length(xtrue)); objectives = zeros(draws)
    for draw in 1:draws
        p.measurement_values .= problem.truth .+ problem.sigmas .* randn(rng, length(problem.truth))
        x0 = initial_state_estimator(s, p)
        sparse_fit = solve_sparse_state_estimator(s, p, x0)
        @test sparse_fit.status == :converged_unique
        @test sparse_fit.iterations <= 8
        states[draw, :] = sparse_fit.state
        objectives[draw] = solve_diagnostics(sparse_fit).objective
        if draw <= 60
            dense_fit = solve_compiled_state_estimator(s, p, x0)
            @test dense_fit.status == :converged_unique
            @test dense_fit.state ≈ sparse_fit.state atol=1e-5
            # A request far below the achievable gradient must not turn a
            # correct estimate into a failure.
            for solve in (solve_sparse_state_estimator, solve_compiled_state_estimator)
                strict = solve(s, p, x0; optimality_tolerance=1e-13)
                @test strict.status == :converged_unique
                @test strict.state ≈ sparse_fit.state atol=1e-5
            end
        end
        if draw <= 8
            ipopt_fit = solve_ipopt_state_estimator(s, p, x0)
            @test ipopt_fit.status == :converged_unique
            @test ipopt_fit.state ≈ sparse_fit.state atol=1e-5
        end
    end

    # Statistical consistency: J ~ chi-square(dof), unbiased, and the reported
    # covariance matches the empirical one (asserted, not merely logged).
    @test abs(sr_mean(objectives) - dof) <= 4 * sqrt(2dof / draws)
    @test 0.6 * 2dof <= sr_var(objectives) <= 1.5 * 2dof
    empirical = sr_column_vars(states)
    ratio = sqrt.(empirical) ./ sqrt.(diag(covariance))
    @test all(0.82 .<= ratio .<= 1.18)
    bias = sr_column_means(states) .- xtrue
    @test all(abs.(bias) .<= 4.5 .* sqrt.(empirical ./ draws))

    # The time-series driver stops at the first non-converged snapshot, so the
    # same defect ended a default run at its second snapshot.
    snapshots = [deepcopy(p) for _ in 1:6]
    for snapshot in snapshots
        snapshot.measurement_values .= problem.truth .+ problem.sigmas .* randn(rng, length(problem.truth))
    end
    for solver in (:sparse, :dense)
        series = solve_time_series_state_estimator(s, snapshots, initial_state_estimator(s, p); solver)
        @test series.status == :converged
        @test series.completed_snapshots == 6
    end
end

@testset "SE solvers: genuine failures keep failure statuses" begin
    problem = solver_regression_problem(solver_regression_feeder(), SOLVER_REGRESSION_LOADS)
    ms = problem.meters(problem.truth .+ problem.sigmas .* randn(StableRNG(5), length(problem.truth)))
    s = compile_state_estimator(solver_regression_feeder(), ms)
    p = SEParameters(s, ms)
    far = initial_state_estimator(s, p) .* 0.6
    for solve in (solve_sparse_state_estimator, solve_compiled_state_estimator)
        limited = solve(s, p, far; max_iterations=1)
        @test limited.status in (:max_iterations, :constraint_not_satisfied, :constraint_restoration_failed)
        @test !solve_status(limited).publishable
        @test ismissing(solve_diagnostics(limited).dof)
        # A trust region too small to move is a stall, not convergence: the
        # model still predicts a large reduction.
        pinned = solve(s, p, far; initial_radius=1e-9, max_radius=1e-9, min_radius=5e-10,
                       constraint_tolerance=1e6)
        @test pinned.status in (:trust_region_stalled, :max_iterations)
        @test !solve_status(pinned).publishable
    end
end

@testset "SE legacy WLS: statistical consistency and agreement" begin
    problem = solver_regression_problem(solver_regression_feeder(), SOLVER_REGRESSION_LOADS)
    ms = problem.meters(problem.truth)
    s = compile_state_estimator(solver_regression_feeder(), ms)
    p = SEParameters(s, ms)
    nf = length(s.free_state_map)
    rng = StableRNG(418)
    draws = 40
    objectives = zeros(draws)
    for draw in 1:draws
        values = problem.truth .+ problem.sigmas .* randn(rng, length(problem.truth))
        p.measurement_values .= values
        compiled = solve_sparse_state_estimator(s, p, initial_state_estimator(s, p))
        legacy = with_logger(NullLogger()) do
            solve_state_estimation(solver_regression_feeder(), problem.meters(values))
        end
        @test solve_status(legacy).publishable
        @test legacy.observability.observable
        @test solve_diagnostics(legacy).dof == solve_diagnostics(compiled).dof == 6
        @test legacy.objective ≈ solve_diagnostics(compiled).objective rtol=1e-6 atol=1e-8
        for (node, k) in s.free_state_map
            estimate = complex(legacy.bus[node[1]][node[2]]["vr"], legacy.bus[node[1]][node[2]]["vi"])
            @test estimate ≈ complex(compiled.state[k], compiled.state[nf + k]) atol=1e-5
        end
        # The two APIs report opposite residual signs; pin that convention.
        @test [r.standardized for r in legacy.residuals] ≈
              -compiled.evaluation.residual atol=1e-5
        objectives[draw] = legacy.objective
    end
    @test abs(sr_mean(objectives) - 6) <= 4 * sqrt(12 / draws)

    disabled = with_logger(NullLogger()) do
        solve_state_estimation(solver_regression_feeder(), ms; check_observability=false)
    end
    enabled = with_logger(NullLogger()) do
        solve_state_estimation(solver_regression_feeder(), ms)
    end
    @test Set(keys(disabled.observability)) == Set(keys(enabled.observability))
    @test ismissing(disabled.observability.observable)
    @test ismissing(solve_diagnostics(disabled).dof)
end

@testset "SE rank threshold does not grow with network size" begin
    # A conditioning threshold is a property of the singular-value ratio.
    n = 1500
    wide = Matrix(Diagonal(exp10.(range(0, -6; length=n))))
    rank, basis = PowerOptLab._se_rank_nullspace(wide)
    @test rank == n
    @test size(basis, 2) == 0
    # Exact deficiency is still found after an orthogonal change of basis.
    rng = StableRNG(77)
    Q = Matrix(qr(randn(rng, 60, 60)).Q)
    deficient = Q * Diagonal(vcat(exp10.(range(0, -5; length=57)), zeros(3))) * Q'
    rank, basis = PowerOptLab._se_rank_nullspace(deficient)
    @test rank == 57
    @test norm(deficient * basis) < 1e-10
    # The least-squares helper truncates at the same ratio.
    rhs = randn(rng, n)
    @test PowerOptLab._se_least_squares(wide, rhs) ≈ rhs ./ diag(wide) rtol=1e-8

    buses = 200
    net = solver_regression_chain(buses)
    rng = StableRNG(buses)
    loads = Dict((("b$k", t) => complex(200 + 300rand(rng), 50rand(rng)) * 25 / buses)
                 for k in 1:buses for t in ("1", "2", "3"))
    problem = solver_regression_problem(net, loads; sigma_s=20.0)
    # 1200 exact device equations leave no tangent direction.
    @test problem.flow.status == :converged_unique
    @test problem.flow.tangent_dimension == 0
    ms = problem.meters(problem.truth .+ problem.sigmas .* randn(rng, length(problem.truth)))
    s = compile_state_estimator(net, ms)
    p = SEParameters(s, ms)
    fit = solve_sparse_state_estimator(s, p, initial_state_estimator(s, p))
    @test fit.status == :converged_unique
    @test fit.tangent_dimension == 6buses
    @test fit.observable_dimension == 6buses
    diagnostic = observability_diagnostics(s, p, fit.state)
    @test diagnostic.unobservable_dimension == 0
    @test 1e4 < diagnostic.condition_number < 1e7
    @test maximum(abs(fit.evaluation.voltage[i] - problem.voltage[s.nodes[i]])
                  for i in eachindex(s.nodes)) < 5.0
end

@testset "SE curvature check: Lanczos estimate for large tangent spaces" begin
    n = 300
    @test n > PowerOptLab._SE_CURVATURE_DENSE_LIMIT
    rng = StableRNG(9)
    Q = Matrix(qr(randn(rng, n, n)).Q)
    positive = Q * Diagonal(range(1.0, 2.0; length=n)) * Q'
    indefinite = Q * Diagonal(vcat(-0.5, range(1.0, 2.0; length=n - 1))) * Q'
    calls = Ref(0)
    counted(A) = v -> (calls[] += 1; A * v)
    @test PowerOptLab._se_lanczos_curvature_status(counted(positive), n, :converged_unique) ==
          :converged_unique
    @test calls[] <= PowerOptLab._SE_CURVATURE_LANCZOS_STEPS
    @test PowerOptLab._se_lanczos_curvature_status(counted(indefinite), n, :converged_unique) ==
          :stationary_not_minimum
    @test PowerOptLab._se_lanczos_curvature_status(v -> nothing, n, :converged_unique) ==
          :curvature_check_failed
    # Negligible negative round-off is not negative curvature.
    nearly = Q * Diagonal(vcat(-1e-9, range(1.0, 2.0; length=n - 1))) * Q'
    @test PowerOptLab._se_lanczos_curvature_status(counted(nearly), n, :converged_underobserved) ==
          :converged_underobserved
end

@testset "SE continuation: the default homotopy changes the model" begin
    net = solver_regression_feeder()
    nodes = sort!(collect(keys(SOLVER_REGRESSION_LOADS)))
    device = ExactDeviceEquation(ConstantPowerDevice(
        [TerminalConnection(node, (node[1], "n")) for node in nodes],
        ComplexF64[SOLVER_REGRESSION_LOADS[node] for node in nodes]))
    s = compile_state_estimator(net, Measurement[]; exact_devices=[device])
    p = SEParameters(s, Measurement[]; exact_devices=[device])
    # The regularisation voltage defaults to the network's nominal level; a
    # millivolt value made every stage identical to the physical law.
    @test p.regularization_voltage ≈ 230.0
    @test SEParameters(s, Measurement[]; exact_devices=[device],
                       regularization_voltage=12.5).regularization_voltage == 12.5
    x0 = initial_state_estimator(s, p)
    continued = solve_with_continuation(s, p, x0)
    @test continued.status == :converged_unique
    @test p.continuation_alpha == 1.0
    @test all(stage.status == :converged_regularized for stage in continued.stages[1:end-1])
    # Each stage has to move: the internal law differs from the physical one.
    @test all(stage.iterations >= 1 for stage in continued.stages)
    first_stage, last_stage = continued.stages[1].state, continued.stages[end].state
    @test norm(first_stage - last_stage, Inf) > 0.5
    # The lightly loaded first stage sits between no load and full load.
    @test norm(first_stage - x0, Inf) < norm(last_stage - x0, Inf)
    direct = solve_compiled_state_estimator(s, SEParameters(s, Measurement[]; exact_devices=[device]), x0)
    @test direct.status == :converged_unique
    @test last_stage ≈ direct.state atol=1e-6
end

@testset "SE current-magnitude page: zero-current statuses" begin
    net = parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1"]},"b":{"terminal_names":["1"]}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.5}},
     "line":{"l":{"bus_from":"src","bus_to":"b","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
    """; from_string=true)
    meter = BranchMeasurement(kind=:imag, line="l", side=:from, terminal="1", value=1.0, sigma=0.01)
    s = compile_state_estimator(net, [meter])
    flat = [230.0, 0.0]
    expected = ((0.0, :undefined_derivative), (1e-3, :stationary_not_minimum),
                (1.0, :converged_underobserved))
    for (epsilon, status) in expected
        p = SEParameters(s, [meter]; current_epsilon=epsilon)
        @test evaluate_state_estimator(s, p, flat).predicted[1] ≈ epsilon
        for solve in (solve_compiled_state_estimator, solve_sparse_state_estimator,
                      solve_ipopt_state_estimator)
            @test solve(s, p, flat).status == status
        end
    end
end

@testset "SE derivatives: phase-to-phase, delta and voltage-following devices" begin
    net = solver_regression_feeder()
    net["linecode"]["lc"]["R_series_1_4"] = 0.01
    net["linecode"]["lc"]["X_series_1_4"] = 0.005
    # Phase-to-phase node telemetry of every kind, next to phase-to-neutral rows.
    ms = Any[Measurement(kind=kind, bus="b1", terminal="1", reference="2", value=0.0, sigma=2.0)
             for kind in (:vr, :vi, :vmag, :pinj, :qinj)]
    append!(ms, [Measurement(kind=kind, bus="b2", terminal="3", value=0.0, sigma=3.0)
                 for kind in (:vmag, :pinj, :qinj)])
    append!(ms, [BranchMeasurement(kind=kind, line="l2", side=:to, terminal="2", reference="n",
                                   value=0.0, sigma=2.0) for kind in (:pflow, :qflow, :imag)])
    delta = [TerminalConnection(("b2", "1"), ("b2", "2")), TerminalConnection(("b2", "2"), ("b2", "3"))]
    devices = [
        ExactDeviceEquation(ConstantPowerDevice(delta, ComplexF64[3000 + 600im, 1200 - 300im])),
        ExactDeviceEquation(ZIPDevice([TerminalConnection(("b1", "3"), ("b1", "1"))],
            ComplexF64[900 + 200im], ComplexF64[1.5 - 0.4im], ComplexF64[0.004 - 0.001im])),
        ExactDeviceEquation(ConstantCurrentDevice([TerminalConnection(("b1", "2"), ("b1", "n"))],
            ComplexF64[6 - 2im]; angle_reference=:voltage)),
        ExactDeviceEquation(ConstantCurrentDevice([TerminalConnection(("b1", "2"), ("b1", "3"))],
            ComplexF64[2 + 1im])),
        ExactDeviceEquation(ConstantPowerDevice([TerminalConnection(("b2", "n"), nothing)],
            ComplexF64[40 + 10im])),
    ]
    rng = StableRNG(31)
    for formulation in (:auto, :current)
        s = compile_state_estimator(net, ms; exact_devices=devices, formulation)
        for (alpha, regularization) in ((1.0, nothing), (0.5, 100.0), (0.0, 100.0))
            p = SEParameters(s, ms; exact_devices=devices, continuation_alpha=alpha,
                             regularization_voltage=regularization, voltage_min_model=1e-3)
            base = initial_state_estimator(s, p)
            for _ in 1:3
                x = base + 3 .* randn(rng, length(base))
                H = residual_jacobian(s, p, x)
                @test H ≈ solver_regression_fd(y -> evaluate_state_estimator(s, p, y).residual, x) rtol=2e-5 atol=2e-5
                C = constraint_jacobian(s, p, x)
                @test C ≈ solver_regression_fd(y -> evaluate_state_estimator(s, p, y).constraints, x) rtol=2e-5 atol=2e-5
            end
        end
    end
end
