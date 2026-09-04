using BMOPFTools: parse_bmopf, solve_pf
using LinearAlgebra: norm, I

# A four-conductor feeder with an explicitly modelled (not perfectly grounded)
# neutral.  The source phase phasors pin the electrical angle; the neutral stays
# in the state, which is exactly the behaviour needed for four-wire estimation.
function compiled_se_net()
    parse_bmopf("""
    {"bus":{
        "src":{"terminal_names":["1","2","3","n"]},
        "b1": {"terminal_names":["1","2","3","n"]},
        "b2": {"terminal_names":["1","2","3","n"]}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1","2","3"],
        "v_magnitude":[230.0,230.0,230.0],"v_angle":[0.0,-2.0943951023931953,2.0943951023931953]}},
     "linecode":{"lc":{"R_series_1_1":0.1,"R_series_2_2":0.1,
        "R_series_3_3":0.1,"R_series_4_4":0.1}},
     "line":{
        "l1":{"bus_from":"src","bus_to":"b1","terminal_map_from":["1","2","3","n"],"terminal_map_to":["1","2","3","n"],"linecode":"lc","length":1.0},
        "l2":{"bus_from":"b1","bus_to":"b2","terminal_map_from":["1","2","3","n"],"terminal_map_to":["1","2","3","n"],"linecode":"lc","length":1.0}}}
    """; from_string=true)
end

function constant_power_test_net()
    parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1"]},"b":{"terminal_names":["1"]}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1"],"v_magnitude":[1000.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.1}},
     "line":{"l":{"bus_from":"src","bus_to":"b","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
    """; from_string=true)
end

function flat_compiled_state(s)
    nf = length(s.free_state_map)
    x = zeros(2nf)
    phase = Dict("1" => 230.0 * cis(0.0),
                 "2" => 230.0 * cis(-2.0943951023931953),
                 "3" => 230.0 * cis(2.0943951023931953),
                 "n" => 0.0 + 0.0im)
    for ((_, terminal), k) in s.free_state_map
        x[k] = real(phase[terminal])
        x[nf + k] = imag(phase[terminal])
    end
    x
end

function central_jacobian(f, x)
    y = f(x)
    J = zeros(length(y), length(x))
    for j in eachindex(x)
        h = 1e-6 * max(1.0, abs(x[j]))
        xp = copy(x); xm = copy(x)
        xp[j] += h; xm[j] -= h
        J[:, j] .= (f(xp) .- f(xm)) ./ (2h)
    end
    J
end

@testset "Compiled constrained state estimator: four-wire evaluator" begin
    measurements = [
        Measurement(kind=:vr,   bus="b1", terminal="1", value=230.0, sigma=0.5),
        Measurement(kind=:vi,   bus="b1", terminal="2", value=-230.0 * sqrt(3) / 2, sigma=0.5),
        Measurement(kind=:vmag, bus="b1", terminal="3", value=230.0, sigma=1.0),
        Measurement(kind=:pinj, bus="b1", terminal="1", value=0.0, sigma=10.0),
        Measurement(kind=:qinj, bus="b1", terminal="1", value=0.0, sigma=10.0),
    ]
    s = compile_state_estimator(compiled_se_net(), measurements;
                                zero_injection=[("b2", "1")])
    p = SEParameters(s, measurements)
    x = flat_compiled_state(s)
    e = evaluate_state_estimator(s, p, x)

    @test length(s.nodes) == 12                 # all four conductors remain explicit
    @test length(x) == 2length(s.free_state_map)
    @test s.free_state_map[("b1", "n")] > 0    # neutral was not silently grounded
    @test e.residual ≈ zeros(length(measurements)) atol=1e-9
    @test e.constraints ≈ zeros(2) atol=1e-9

    Hr = residual_jacobian(s, p, x)
    Hfd = central_jacobian(y -> evaluate_state_estimator(s, p, y).residual, x)
    C = constraint_jacobian(s, p, x)
    Cfd = central_jacobian(y -> evaluate_state_estimator(s, p, y).constraints, x)
    @test Hr ≈ Hfd rtol=1e-5 atol=1e-6
    @test C ≈ Cfd rtol=1e-8 atol=1e-8
end

@testset "Compiled constrained state estimator: parameter updates retain structure" begin
    measurements = [Measurement(kind=:vr, bus="b1", terminal="1", value=230.0, sigma=1.0)]
    s = compile_state_estimator(compiled_se_net(), measurements)
    p = SEParameters(s, measurements)
    x = flat_compiled_state(s)
    e1 = evaluate_state_estimator(s, p, x)
    p.measurement_values[1] = 229.0
    p.covariance_values[1] = 0.5
    e2 = evaluate_state_estimator(s, p, x)
    @test e1.predicted == e2.predicted == [230.0]
    @test e2.residual == [2.0]
    @test size(s.passive_pattern) == (12, 12)
end

@testset "Compiled constrained state estimator: source updates preserve I = YV" begin
    measurements = [Measurement(kind=:vr, bus="b1", terminal="1", value=230.0, sigma=1.0)]
    s = compile_state_estimator(compiled_se_net(), measurements)
    p = SEParameters(s, measurements)
    x = flat_compiled_state(s)
    source = s.node_index[("src", "1")]
    p.fixed_voltages[source] = 235.0 + 0.0im
    e = evaluate_state_estimator(s, p, x)
    @test e.current ≈ s.passive_pattern * e.voltage atol=1e-10
    @test abs(e.current[source]) > 1.0
end

@testset "Compiled constrained state estimator: branch telemetry" begin
    measurements = [
        BranchMeasurement(kind=:ire,   line="l1", side=:from, terminal="1", value=0.0, sigma=0.1),
        BranchMeasurement(kind=:iim,   line="l1", side=:from, terminal="1", value=0.0, sigma=0.1),
        BranchMeasurement(kind=:imag,  line="l1", side=:from, terminal="1", value=0.0, sigma=0.1),
        BranchMeasurement(kind=:pflow, line="l1", side=:from, terminal="1", value=0.0, sigma=10.0),
        BranchMeasurement(kind=:qflow, line="l1", side=:from, terminal="1", value=0.0, sigma=10.0),
    ]
    s = compile_state_estimator(compiled_se_net(), measurements)
    p = SEParameters(s, measurements)
    x = flat_compiled_state(s)
    x[s.free_state_map[("b1", "1")]] -= 1.0  # non-zero l1 current for |I| derivative
    p.measurement_values .= evaluate_state_estimator(s, p, x).predicted
    @test evaluate_state_estimator(s, p, x).residual ≈ zeros(length(measurements)) atol=1e-10
    @test residual_jacobian(s, p, x) ≈
          central_jacobian(y -> evaluate_state_estimator(s, p, y).residual, x) rtol=1e-5 atol=1e-6
end

@testset "Compiled constrained state estimator: dense composite-step reference solver" begin
    net = compiled_se_net()
    seed = compile_state_estimator(net)
    xtrue = flat_compiled_state(seed)
    nf = length(seed.free_state_map)
    # Rectangular phasors to ground fully observe every free conductor, including
    # the source and feeder neutrals.  This isolates the solver test from an
    # observability ambiguity rather than hiding one behind a prior.
    measurements = Measurement[]
    for ((bus, terminal), k) in seed.free_state_map
        push!(measurements, Measurement(kind=:vr, bus=bus, terminal=terminal,
                                        reference=nothing, value=xtrue[k], sigma=1.0))
        push!(measurements, Measurement(kind=:vi, bus=bus, terminal=terminal,
                                        reference=nothing, value=xtrue[nf + k], sigma=1.0))
    end
    s = compile_state_estimator(net, measurements)
    result = solve_compiled_state_estimator(s, SEParameters(s, measurements), zeros(length(xtrue));
                                            initial_radius=2.0)
    @test result.status == :converged_unique
    @test result.state ≈ xtrue atol=1e-9
    @test norm(result.evaluation.residual) ≤ 1e-9
    @test result.constraint_rank == 0

    sparse_result = solve_sparse_state_estimator(s, SEParameters(s, measurements), zeros(length(xtrue));
                                                initial_radius=2.0)
    @test sparse_result.status == :converged_unique
    @test sparse_result.state ≈ xtrue atol=1e-9
    @test isempty(sparse_result.constraint_multipliers)

    diagnostic = observability_diagnostics(s, SEParameters(s, measurements), xtrue)
    @test diagnostic.observable_dimension == length(xtrue)
    @test diagnostic.unobservable_dimension == 0
    covariance = selected_state_covariance(s, SEParameters(s, measurements), xtrue, [1, nf + 1])
    @test covariance ≈ Matrix{Float64}(I, 2, 2) atol=1e-10
    @test derived_covariance(s, SEParameters(s, measurements), xtrue,
                             reshape([1.0; zeros(length(xtrue) - 1)], 1, :)) ≈ ones(1, 1)

    # One magnitude row cannot identify the full four-wire state.  The
    # unobservable basis is generated only when explicitly requested, and a
    # finite covariance is refused rather than invented.
    weak = [Measurement(kind=:vmag, bus="b1", terminal="1", value=230.0, sigma=1.0)]
    sw = compile_state_estimator(net, weak)
    pw = SEParameters(sw, weak)
    ow = observability_diagnostics(sw, pw, xtrue)
    @test ow.unobservable_dimension > 0
    direction = unobservable_directions(sw, pw, xtrue; count=1)
    @test norm(residual_jacobian(sw, pw, xtrue) * direction) ≤ 1e-8
    @test_throws ArgumentError selected_state_covariance(sw, pw, xtrue, [1])

    # Priors are ordinary whitened residuals.  They add rows to H without
    # changing the compiled topology or measurement incidence.
    prior = StatePrior([1], [xtrue[1] - 1.0], [2.0])
    pp = SEParameters(s, measurements; prior=prior)
    ep = evaluate_state_estimator(s, pp, xtrue)
    @test length(ep.prior_residual) == 1
    @test ep.prior_residual == [0.5]
    @test size(residual_jacobian(s, pp, xtrue), 1) == length(measurements) + 1

    # A time-series solve reuses `s`, warm-starts snapshot 2, and makes the
    # previous state a stochastic (not exact) prior.
    p1 = SEParameters(s, measurements)
    p2 = SEParameters(s, measurements)
    p2.measurement_values[1] += 2.0
    series = solve_time_series_state_estimator(s, [p1, p2], zeros(length(xtrue));
                                                previous_state_sigma=10.0, initial_radius=2.0)
    @test series.status == :converged
    @test series.completed_snapshots == 2
    @test length(series.snapshots) == 2
    @test length(p1.prior_indices) == 0
    @test p2.prior_values ≈ series.snapshots[1].state
    measured_state = s.free_state_map[(measurements[1].bus, measurements[1].terminal)]
    @test series.snapshots[2].state[measured_state] > series.snapshots[1].state[measured_state]
    @test series.snapshots[2].state[measured_state] < p2.measurement_values[1]
end

@testset "Compiled constrained state estimator: exact nonlinear devices" begin
    connection = TerminalConnection(("b1", "1"), ("b1", "n"))
    devices = ExactDeviceEquation[
        ExactDeviceEquation(ConstantPowerDevice([connection], ComplexF64[4_000 + 1_000im])),
        ExactDeviceEquation(ConstantCurrentDevice([connection], ComplexF64[2.0 - 0.5im])),
        ExactDeviceEquation(ZIPDevice([connection], ComplexF64[1_000 + 200im],
                                      ComplexF64[0.5 + 0.1im], ComplexF64[0.01 - 0.002im])),
    ]
    for device in devices
        s = compile_state_estimator(compiled_se_net(); exact_devices=[device])
        p = SEParameters(s; exact_devices=[device])
        x = flat_compiled_state(s)
        C = constraint_jacobian(s, p, x)
        Cfd = central_jacobian(y -> evaluate_state_estimator(s, p, y).constraints, x)
        @test C ≈ Cfd rtol=1e-5 atol=1e-6
        @test any(x -> !iszero(x), evaluate_state_estimator(s, p, x).device_current)
    end

    # A phase-to-phase connection injects equal and opposite conductor current:
    # this is the same primitive used for a delta winding or a signed generator.
    delta = ExactDeviceEquation(ConstantPowerDevice(
        [TerminalConnection(("b1", "1"), ("b1", "2"))], ComplexF64[2_000 + 0im]))
    sδ = compile_state_estimator(compiled_se_net(); exact_devices=[delta])
    pδ = SEParameters(sδ; exact_devices=[delta])
    eδ = evaluate_state_estimator(sδ, pδ, flat_compiled_state(sδ))
    i1, i2 = sδ.node_index[("b1", "1")], sδ.node_index[("b1", "2")]
    @test eδ.device_current[i1] ≈ -eδ.device_current[i2]

    cp = devices[1]
    sguard = compile_state_estimator(compiled_se_net(); exact_devices=[cp])
    pguard = SEParameters(sguard; exact_devices=[cp], voltage_min_model=1.0)
    @test_throws DomainError evaluate_state_estimator(sguard, pguard, zeros(2length(sguard.free_state_map)))
    pguard.continuation_alpha = 0.0
    @test isfinite(norm(evaluate_state_estimator(sguard, pguard, zeros(2length(sguard.free_state_map))).constraints))

    load = ExactDeviceEquation(ConstantPowerDevice(
        [TerminalConnection(("b", "1"), nothing)], ComplexF64[1_000 + 0im]))
    ssolve = compile_state_estimator(constant_power_test_net(); exact_devices=[load])
    psolve = SEParameters(ssolve; exact_devices=[load])
    direct = solve_compiled_state_estimator(ssolve, psolve, [1000.0, 0.0];
                                             initial_radius=0.1, constraint_tolerance=1e-7)
    @test direct.status == :converged_unique
    @test norm(direct.evaluation.constraints) ≤ 1e-7
    continuation = solve_with_continuation(ssolve, psolve, [1000.0, 0.0];
                                           alphas=[0.0, 0.5, 1.0], initial_radius=0.1,
                                           constraint_tolerance=1e-7)
    @test continuation.status == :converged_unique
    @test psolve.continuation_alpha == 1.0

    # The sparse Hachtel path exposes the multiplier directly.  Here the
    # voltage reading conflicts with the exact load equation, so the multiplier
    # is nonzero rather than hiding the conflict in a soft compromise.
    conflicting = [Measurement(kind=:vr, bus="b", terminal="1", reference=nothing,
                               value=999.0, sigma=1.0)]
    smult = compile_state_estimator(constant_power_test_net(), conflicting; exact_devices=[load])
    pmult = SEParameters(smult, conflicting; exact_devices=[load])
    sparse = solve_sparse_state_estimator(smult, pmult, [1000.0, 0.0];
                                          initial_radius=0.1, constraint_tolerance=1e-7)
    @test sparse.status == :converged_unique
    @test norm(sparse.evaluation.constraints) ≤ 1e-7
    @test abs(sparse.constraint_multipliers[1]) > 1e-3
end

# ── regression coverage for the review fixes ────────────────────────────────

# Single-phase feeder src—b—c, `r` ohm per segment, scaled to `vsrc` volts.
function scaled_feeder(vsrc, r)
    parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1"]},"b":{"terminal_names":["1"]},"c":{"terminal_names":["1"]}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1"],"v_magnitude":[$vsrc],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":$r}},
     "line":{"l1":{"bus_from":"src","bus_to":"b","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
             "l2":{"bus_from":"b","bus_to":"c","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
    """; from_string=true)
end

@testset "Compiled SE: magnitude_epsilon reaches branch current magnitudes" begin
    net = compiled_se_net()
    bm = [BranchMeasurement(kind=:imag, line="l1", side=:from, terminal="1",
                            value=0.0, sigma=0.1)]
    s = compile_state_estimator(net, bm)
    x = flat_compiled_state(s)          # no load anywhere => l1 carries no current

    # Without smoothing the prediction is a hard zero and the derivative does not
    # exist; with it, both the value and the analytic Jacobian are the smoothed
    # magnitude. The option used to be accepted and then ignored for :imag.
    p0 = SEParameters(s, bm; magnitude_epsilon=0.0)
    @test evaluate_state_estimator(s, p0, x).predicted == [0.0]
    @test_throws DomainError residual_jacobian(s, p0, x)

    p5 = SEParameters(s, bm; magnitude_epsilon=5.0)
    @test evaluate_state_estimator(s, p5, x).predicted ≈ [5.0]
    @test residual_jacobian(s, p5, x) ≈
          central_jacobian(y -> evaluate_state_estimator(s, p5, y).residual, x) rtol=1e-5 atol=1e-6

    # The same smoothing must reach a node voltage magnitude, which it always did.
    vm = [Measurement(kind=:vmag, bus="b1", terminal="n", reference=nothing,
                      value=0.0, sigma=1.0)]
    sv = compile_state_estimator(net, vm)
    @test evaluate_state_estimator(sv, SEParameters(sv, vm; magnitude_epsilon=5.0),
                                   flat_compiled_state(sv)).predicted ≈ [5.0]
end

@testset "Compiled SE: an undefined magnitude derivative is a status, not a throw" begin
    # A line carrying exactly zero current is the textbook flat-start failure of
    # ampere measurements. Both solvers must diagnose it and return.
    net = compiled_se_net()
    seed = compile_state_estimator(net)
    xtrue = flat_compiled_state(seed)
    nf = length(seed.free_state_map)
    meas = Any[BranchMeasurement(kind=:imag, line="l1", side=:from, terminal="1",
                                 value=1.0, sigma=0.1)]
    for ((bus, terminal), k) in seed.free_state_map
        push!(meas, Measurement(kind=:vr, bus=bus, terminal=terminal, reference=nothing,
                                value=xtrue[k], sigma=1.0))
        push!(meas, Measurement(kind=:vi, bus=bus, terminal=terminal, reference=nothing,
                                value=xtrue[nf + k], sigma=1.0))
    end
    s = compile_state_estimator(net, meas)
    p = SEParameters(s, meas)
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator)
        result = solver(s, p, xtrue)
        @test result.status == :undefined_derivative
        @test !solve_status(result).publishable
    end
    # Smoothing the magnitude makes the same start solvable.
    ps = SEParameters(s, meas; magnitude_epsilon=1e-3)
    @test solve_compiled_state_estimator(s, ps, xtrue).status != :undefined_derivative
end

@testset "Compiled SE: a bare zero-injection bus covers its neutral conductor" begin
    net = compiled_se_net()
    bare = compile_state_estimator(net; zero_injection=["b2"])
    tuples = compile_state_estimator(net;
        zero_injection=[("b2","1"), ("b2","2"), ("b2","3"), ("b2","n")])
    covered(s) = Set(s.nodes[i] for i in s.constraint_pattern)
    # Four-wire KCL is per conductor against earth, so "nothing is attached to
    # b2" constrains its neutral too. Omitting it left the neutral voltage free.
    @test covered(bare) == covered(tuples)
    @test ("b2", "n") in covered(bare)
    # The WLS estimator states zero injection per phase against a return
    # terminal, so its expansion deliberately excludes the neutral.
    @test !(("b2","n") in PowerOptLab._zero_injection_set(net, ["b2"], "n"))
end

@testset "Compiled SE: time series preserves a caller-supplied prior" begin
    net = compiled_se_net()
    seed = compile_state_estimator(net)
    xtrue = flat_compiled_state(seed)
    nf = length(seed.free_state_map)
    meas = Measurement[]
    for ((bus, terminal), k) in seed.free_state_map
        push!(meas, Measurement(kind=:vr, bus=bus, terminal=terminal, reference=nothing,
                                value=xtrue[k], sigma=1.0))
        push!(meas, Measurement(kind=:vi, bus=bus, terminal=terminal, reference=nothing,
                                value=xtrue[nf + k], sigma=1.0))
    end
    s = compile_state_estimator(net, meas)
    prior = StatePrior([1], [xtrue[1] - 4.0], [2.0])

    # Without a previous-state prior the caller's prior must survive untouched.
    p1 = SEParameters(s, meas; prior=prior)
    p2 = SEParameters(s, meas; prior=prior)
    series = solve_time_series_state_estimator(s, [p1, p2], zeros(length(xtrue));
                                               initial_radius=2.0)
    @test series.status == :converged
    @test p1.prior_indices == [1] && p1.prior_values == [xtrue[1] - 4.0]
    @test p2.prior_indices == [1] && p2.prior_values == [xtrue[1] - 4.0]

    # With one, it is LAYERED ON TOP of the caller's prior, not substituted for
    # it: two independent observations of the same state, so both rows survive.
    q1 = SEParameters(s, meas; prior=prior)
    q2 = SEParameters(s, meas; prior=prior)
    solve_time_series_state_estimator(s, [q1, q2], zeros(length(xtrue));
                                      previous_state_sigma=10.0, initial_radius=2.0)
    @test q1.prior_indices == [1]
    @test length(q2.prior_indices) == 1 + length(xtrue)
    @test q2.prior_values[1] == xtrue[1] - 4.0
end

@testset "Compiled SE: recovers a BMOPFTools power flow from exact telemetry" begin
    # Ground truth from the engine's own power flow. The compiled estimator had
    # no test tying it to a solved network state at all.
    loaded = parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1"]},"b":{"terminal_names":["1"]},"c":{"terminal_names":["1"]}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1"],"v_magnitude":[1000.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.5}},
     "line":{"l1":{"bus_from":"src","bus_to":"b","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
             "l2":{"bus_from":"b","bus_to":"c","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}},
     "load":{"d":{"bus":"c","terminal_map":["1"],"configuration":"SINGLE_PHASE","p_nom":[20000.0],"q_nom":[5000.0]}}}
    """; from_string=true)
    pf = solve_pf(loaded; per_unit=false)
    vtrue = Dict(b => pf["bus"][b]["1"]["vr"] + im * pf["bus"][b]["1"]["vi"] for b in ("b","c"))

    est_net = scaled_feeder(1000.0, 0.5)   # same topology, load removed
    meas = [Measurement(kind=:vmag, bus="b", reference=nothing, value=abs(vtrue["b"]), sigma=0.1),
            Measurement(kind=:pinj, bus="c", reference=nothing, value=-20_000.0, sigma=1.0),
            Measurement(kind=:qinj, bus="c", reference=nothing, value=-5_000.0, sigma=1.0)]
    s = compile_state_estimator(est_net, meas; zero_injection=["b"])
    p = SEParameters(s, meas)
    x0 = zeros(2length(s.free_state_map))
    for ((_, _), k) in s.free_state_map; x0[k] = 1000.0; end

    dense = solve_compiled_state_estimator(s, p, x0; initial_radius=1.0)
    sparse_r = solve_sparse_state_estimator(s, p, x0; initial_radius=1.0)
    @test dense.status == :converged_unique
    @test sparse_r.status == :converged_unique
    for r in (dense, sparse_r), b in ("b", "c")
        k = s.free_state_map[(b, "1")]
        v = r.state[k] + im * r.state[k + length(s.free_state_map)]
        @test v ≈ vtrue[b] rtol=1e-6
    end
    # The two solvers implement different steps; on the same problem they must
    # still land on the same estimate.
    @test dense.state ≈ sparse_r.state rtol=1e-6
end

@testset "Compiled SE: feasibility tolerance is scale invariant" begin
    # The same dimensionless problem at two voltage levels. `norm(c)` is a
    # difference of Y*V products, so its double-precision floor scales with the
    # network's SI magnitude; a purely absolute tolerance made the high-voltage
    # case report failure at a perfectly converged point.
    for (vsrc, r) in ((230.0, 0.1), (132_000.0, 1e-4))
        load = ExactDeviceEquation(ConstantPowerDevice(
            [TerminalConnection(("c", "1"), nothing)], ComplexF64[vsrc * 10 + 0im]))
        net = scaled_feeder(vsrc, r)
        s = compile_state_estimator(net; zero_injection=["b"], exact_devices=[load])
        p = SEParameters(s; exact_devices=[load])
        x0 = zeros(2length(s.free_state_map))
        for ((_, _), k) in s.free_state_map; x0[k] = vsrc; end
        for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator)
            result = solver(s, p, x0; initial_radius=0.5)
            @test result.status == :converged_unique
            # Feasible to round-off RELATIVE to the network's own current scale.
            @test norm(result.evaluation.constraints) <=
                  1e-11 * PowerOptLab._se_current_scale(s, p)
        end
    end
end

@testset "Compiled SE: |I| admits mirror solutions that |I_re|,|I_im| do not" begin
    # The classic difficulty with ampere measurements: a magnitude carries no
    # direction, so |V| plus |I| is satisfied exactly by a conjugate pair of
    # states. Both are reported :converged_unique — that status is a LOCAL rank
    # statement about one point, never a global uniqueness claim.
    net = scaled_feeder(230.0, 0.5)
    net["line"] = Dict("l" => net["line"]["l1"])       # keep a single segment
    delete!(net["bus"], "c")
    vtrue = 220.0 * cis(-0.05)
    itrue = (230.0 - vtrue) / 0.5

    ambiguous = Any[Measurement(kind=:vmag, bus="b", reference=nothing,
                                value=abs(vtrue), sigma=0.01),
                    BranchMeasurement(kind=:imag, line="l", side=:from, terminal="1",
                                      value=abs(itrue), sigma=0.01)]
    s = compile_state_estimator(net, ambiguous)
    p = SEParameters(s, ambiguous)
    angles = Float64[]
    for x0 in ([219.0, -11.0], [219.0, +11.0])
        r = solve_compiled_state_estimator(s, p, x0; initial_radius=1.0)
        @test r.status == :converged_unique
        @test norm(r.evaluation.residual) < 1e-6      # both fit the data exactly
        @test abs(r.state[1] + im * r.state[2]) ≈ abs(vtrue) rtol=1e-6
        push!(angles, angle(r.state[1] + im * r.state[2]))
    end
    @test angles[1] ≈ -angles[2] rtol=1e-4            # a mirror pair, not one answer
    @test abs(angles[1] - angles[2]) > 1e-3

    # The same current as a rectangular pair fixes the direction, and both starts
    # converge to the one true state.
    resolved = Any[Measurement(kind=:vmag, bus="b", reference=nothing,
                               value=abs(vtrue), sigma=0.01),
                   BranchMeasurement(kind=:ire, line="l", side=:from, terminal="1",
                                     value=real(itrue), sigma=0.01),
                   BranchMeasurement(kind=:iim, line="l", side=:from, terminal="1",
                                     value=imag(itrue), sigma=0.01)]
    s2 = compile_state_estimator(net, resolved)
    p2 = SEParameters(s2, resolved)
    for x0 in ([219.0, -11.0], [219.0, +11.0])
        r = solve_compiled_state_estimator(s2, p2, x0; initial_radius=1.0)
        @test r.status == :converged_unique
        @test r.state[1] + im * r.state[2] ≈ vtrue rtol=1e-6
    end
end
