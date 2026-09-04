using BMOPFTools: parse_bmopf, solve_pf
using LinearAlgebra: norm, svdvals, I

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
    @test length(p2.prior_indices) == 0     # driver scratch is removed on exit
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

@testset "Compiled SE: voltage and current magnitudes smooth independently" begin
    net = compiled_se_net()
    bm = [BranchMeasurement(kind=:imag, line="l1", side=:from, terminal="1",
                            value=0.0, sigma=0.1)]
    s = compile_state_estimator(net, bm)
    x = flat_compiled_state(s)          # no load anywhere => l1 carries no current

    # Without smoothing the prediction is a hard zero and the derivative does not
    # exist; with it, both the value and the analytic Jacobian are the smoothed
    # magnitude. The option used to be accepted and then ignored for :imag.
    p0 = SEParameters(s, bm)
    @test evaluate_state_estimator(s, p0, x).predicted == [0.0]
    @test_throws DomainError residual_jacobian(s, p0, x)

    # `current_epsilon` is in AMPERES and smooths |I|. `magnitude_epsilon` is in
    # VOLTS and must not touch a branch current: one scalar cannot carry two
    # physical units, so a mixed :vmag/:imag dataset needs both knobs.
    pv = SEParameters(s, bm; magnitude_epsilon=5.0)
    @test evaluate_state_estimator(s, pv, x).predicted == [0.0]
    @test_throws DomainError residual_jacobian(s, pv, x)

    pc = SEParameters(s, bm; current_epsilon=5.0)
    @test evaluate_state_estimator(s, pc, x).predicted ≈ [5.0]
    @test residual_jacobian(s, pc, x) ≈
          central_jacobian(y -> evaluate_state_estimator(s, pc, y).residual, x) rtol=1e-5 atol=1e-6

    # Symmetrically, a node voltage magnitude answers only to `magnitude_epsilon`.
    vm = [Measurement(kind=:vmag, bus="b1", terminal="n", reference=nothing,
                      value=0.0, sigma=1.0)]
    sv = compile_state_estimator(net, vm)
    xv = flat_compiled_state(sv)
    @test evaluate_state_estimator(sv, SEParameters(sv, vm; magnitude_epsilon=5.0),
                                   xv).predicted ≈ [5.0]
    @test evaluate_state_estimator(sv, SEParameters(sv, vm; current_epsilon=5.0),
                                   xv).predicted == [0.0]
    @test_throws ArgumentError SEParameters(sv, vm; current_epsilon=-1.0)
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
    # Smoothing the CURRENT magnitude makes the same start solvable.
    ps = SEParameters(s, meas; current_epsilon=1e-3)
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

    # With one, the previous-state prior is LAYERED ON TOP of the caller's while
    # the snapshot solves, then removed again: it is driver scratch, not caller
    # input. Repeated calls on the same objects must therefore be idempotent --
    # otherwise each invocation would stack another previous-state block and
    # count the same historical estimate again, shrinking reported uncertainty.
    q1 = SEParameters(s, meas; prior=prior)
    q2 = SEParameters(s, meas; prior=prior)
    for call in 1:3
        series2 = solve_time_series_state_estimator(s, [q1, q2], zeros(length(xtrue));
                                                    previous_state_sigma=10.0, initial_radius=2.0)
        @test series2.status == :converged
        @test q1.prior_indices == [1]
        @test q2.prior_indices == [1]
        @test q2.prior_values == [xtrue[1] - 4.0]
    end

    # A stalled run must restore the caller's priors too.
    bad = SEParameters(s, meas; prior=prior)
    bad.covariance_values[1] = 1e-14        # force a hopeless snapshot
    solve_time_series_state_estimator(s, [q1, bad], zeros(length(xtrue));
                                      previous_state_sigma=10.0, initial_radius=2.0,
                                      max_iterations=3)
    @test q1.prior_indices == [1]
    @test bad.prior_indices == [1]
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

@testset "Compiled SE: convergence tests are scale invariant" begin
    # DIMENSIONALLY EQUIVALENT networks: holding Z fixed and scaling V by k
    # scales I by k and S by k^2, so the per-unit problem is IDENTICAL at every
    # k and only the SI magnitudes move. Anything that changes across this sweep
    # is an artefact of units, not of the estimation problem.
    #
    # Both stopping tests are exercised: the exact device drives the FEASIBILITY
    # test, and the |V| measurement supplies a stochastic residual row so the
    # STATIONARITY test has a non-empty H to normalise (the earlier version of
    # this test had no measurements at all, so it never touched that path).
    reference = nothing
    for k in (1.0, 1e2, 1e4, 5.74e5)          # 230 V, 23 kV, 2.3 MV, 132 MV
        vsrc = 230.0 * k
        load = ExactDeviceEquation(ConstantPowerDevice(
            [TerminalConnection(("c", "1"), nothing)], ComplexF64[(4_000 + 1_000im) * k^2]))
        net = scaled_feeder(vsrc, 0.1)        # Z fixed across the sweep
        meas = [Measurement(kind=:vmag, bus="b", reference=nothing,
                            value=0.999 * vsrc, sigma=1e-3 * vsrc)]
        s = compile_state_estimator(net, meas; zero_injection=["b"], exact_devices=[load])
        p = SEParameters(s, meas; exact_devices=[load], voltage_min_model=1e-3 * k)
        x0 = zeros(2length(s.free_state_map))
        for ((_, _), kk) in s.free_state_map; x0[kk] = vsrc; end

        scale = PowerOptLab._se_current_scale(s, p)
        for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator)
            result = solver(s, p, x0; initial_radius=0.5)
            @test result.status == :converged_unique
            # Feasible to round-off RELATIVE to the network's own current scale.
            @test norm(result.evaluation.constraints) <= 1e-11 * scale
            # The per-unit answer must not depend on k.
            pu = result.state ./ vsrc
            reference === nothing && (reference = pu)
            @test pu ≈ reference rtol=1e-6
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

# ── observability: the numbers quoted in docs/src/estimation/observability.md ──

# Two-bus radial feeder; state order is [Vre_b, Vre_c, Vim_b, Vim_c].
obs_net() = scaled_feeder(230.0, 0.5)
const OBS_VB = 211.58 + 7.79im
const OBS_VC = 204.90 + 10.90im
obs_state() = [real(OBS_VB), real(OBS_VC), imag(OBS_VB), imag(OBS_VC)]
obs_m(kind, bus, value; sigma=0.1) =
    Measurement(kind=kind, bus=bus, reference=nothing, value=value, sigma=sigma)

function obs_diag(ms; zi=String[])
    s = compile_state_estimator(obs_net(), ms; zero_injection=zi)
    p = SEParameters(s, ms)
    (s, p, observability_diagnostics(s, p, obs_state()))
end

@testset "Observability: redundancy is not observability" begin
    vb = obs_m(:vmag, "b", abs(OBS_VB))
    vc = obs_m(:vmag, "c", abs(OBS_VC))
    # Three copies of one reading: redundant against noise, rank one.
    _, _, d3 = obs_diag([vb, vb, vb])
    @test d3.tangent_dimension == 4
    @test d3.observable_dimension == 1
    @test d3.unobservable_dimension == 3

    # Same count as the determined set below, opposite verdict: placement decides.
    _, _, dm = obs_diag([vb, vc])
    @test dm.observable_dimension == 2
    @test dm.unobservable_dimension == 2

    rect = [obs_m(:vr, "b", real(OBS_VB)), obs_m(:vi, "b", imag(OBS_VB)),
            obs_m(:vr, "c", real(OBS_VC)), obs_m(:vi, "c", imag(OBS_VC))]
    _, _, dr = obs_diag(rect)
    @test dr.observable_dimension == 4
    @test dr.unobservable_dimension == 0
    @test dr.min_singular ≈ 10.0 rtol=1e-6      # 1/sigma with sigma = 0.1
end

@testset "Observability: the unobservable direction of a magnitude-only set" begin
    # |V| at both buses leaves both phase angles free. The null space should be
    # exactly the tangent to each constant-magnitude circle -- the analytic
    # direction (-Im V, Re V)/|V| -- which is what makes this diagnostic worth
    # reading rather than just counting.
    ms = [obs_m(:vmag, "b", abs(OBS_VB)), obs_m(:vmag, "c", abs(OBS_VC))]
    s, p, d = obs_diag(ms)
    @test d.unobservable_dimension == 2
    x = obs_state()
    U = unobservable_directions(s, p, x)
    @test size(U) == (4, 2)
    H = residual_jacobian(s, p, x)
    for k in axes(U, 2)
        @test norm(H * U[:, k]) < 1e-10
    end
    # `unobservable_directions` returns an orthonormal basis, so projecting an
    # analytic phase direction onto it must reproduce it exactly.
    for (V, ire, iim) in ((OBS_VB, 1, 3), (OBS_VC, 2, 4))
        u = zeros(4)
        u[ire] = -imag(V) / abs(V)
        u[iim] =  real(V) / abs(V)
        @test norm(U * (U' * u) - u) < 1e-8     # u is already in the span
    end
end

@testset "Observability: an exact constraint can make a set sufficient" begin
    # The measurements do not change; the constraint removes two directions
    # they would otherwise have had to explain. rank(H) alone is the wrong test.
    ms = [obs_m(:vr, "b", real(OBS_VB)), obs_m(:vi, "b", imag(OBS_VB))]
    _, _, free = obs_diag(ms)
    @test free.tangent_dimension == 4
    @test free.observable_dimension == 2
    @test free.unobservable_dimension == 2

    _, _, held = obs_diag(ms; zi=["c"])
    @test held.tangent_dimension == 2
    @test held.observable_dimension == 2
    @test held.unobservable_dimension == 0
end

@testset "Observability: a critical measurement hides bad data completely" begin
    # Exactly determined => every measurement is critical => every residual is
    # identically zero, so a gross error leaves no signature at all.
    rect = [obs_m(:vr, "b", real(OBS_VB)), obs_m(:vi, "b", imag(OBS_VB)),
            obs_m(:vr, "c", real(OBS_VC)), obs_m(:vi, "c", imag(OBS_VC))]
    s = compile_state_estimator(obs_net(), rect)
    p = SEParameters(s, rect)
    x0 = zeros(4); x0[1] = x0[2] = 230.0
    clean = solve_compiled_state_estimator(s, p, x0; initial_radius=2.0)
    @test clean.status == :converged_unique
    @test norm(clean.evaluation.residual) < 1e-10

    p.measurement_values[1] += 5.0                     # 50 sigma
    bad = solve_compiled_state_estimator(s, p, x0; initial_radius=2.0)
    @test bad.status == :converged_unique
    @test norm(bad.evaluation.residual) < 1e-10        # NO residual signature
    # The estimate absorbed the whole error.
    @test bad.state[1] - clean.state[1] ≈ 5.0 rtol=1e-6

    # One redundant row makes the same error unmistakable.
    red = vcat(rect, [obs_m(:vmag, "b", abs(OBS_VB))])
    s2 = compile_state_estimator(obs_net(), red)
    p2 = SEParameters(s2, red)
    @test norm(solve_compiled_state_estimator(s2, p2, x0; initial_radius=2.0
              ).evaluation.residual) < 1e-8
    p2.measurement_values[1] += 5.0
    exposed = solve_compiled_state_estimator(s2, p2, x0; initial_radius=2.0)
    @test norm(exposed.evaluation.residual) > 10.0
    @test maximum(abs, exposed.evaluation.residual) > 20.0
end

@testset "Observability: covariance refuses to invent information" begin
    ms = [obs_m(:vr, "b", real(OBS_VB)), obs_m(:vi, "b", imag(OBS_VB))]
    s = compile_state_estimator(obs_net(), ms)
    x = obs_state()

    # Rank deficient: an unobservable direction has infinite variance, and a
    # pseudo-inverse would silently report confidence that does not exist.
    p = SEParameters(s, ms)
    @test observability_diagnostics(s, p, x).unobservable_dimension == 2
    @test_throws ArgumentError selected_state_covariance(s, p, x, [1])

    # A prior restores rank -- and the covariance says where the answer came
    # from: the measured bus keeps its meter's sigma, the unmeasured bus
    # reports exactly the prior's.
    for psig in (1.0, 10.0, 100.0)
        pp = SEParameters(s, ms;
            prior=StatePrior([2, 4], [real(OBS_VC), imag(OBS_VC)], [psig, psig]))
        @test observability_diagnostics(s, pp, x).unobservable_dimension == 0
        C = selected_state_covariance(s, pp, x, [1, 2])
        @test sqrt(C[1, 1]) ≈ 0.1  rtol=1e-8      # the meter's sigma
        @test sqrt(C[2, 2]) ≈ psig rtol=1e-8      # entirely the prior
    end
end

# ── current-magnitude conditioning: BOTH the benign and the degenerate path ──

# Two-bus feeder, 0.5 ohm, 230 V source. State is [Vre_b, Vim_b].
cm_net() = parse_bmopf("""
{"bus":{"src":{"terminal_names":["1"]},"b":{"terminal_names":["1"]}},
 "voltage_source":{"s":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.5}},
 "line":{"l":{"bus_from":"src","bus_to":"b","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

# |V_b| and |I_l| at an operating point with current `Ic`; returns (H, x, s, p).
function cm_setup(Ic; sigma=0.01)
    V = 230.0 - 0.5 * Ic
    ms = Any[Measurement(kind=:vmag, bus="b", reference=nothing, value=abs(V), sigma=sigma),
             BranchMeasurement(kind=:imag, line="l", side=:from, terminal="1",
                               value=abs(Ic), sigma=sigma)]
    s = compile_state_estimator(cm_net(), ms)
    p = SEParameters(s, ms)
    ([real(V), imag(V)], s, p)
end

@testset "Current magnitude: conditioning is bounded as |I| falls at fixed angle" begin
    # Shrinking the current at a FIXED power-factor angle does not degrade the
    # conditioning: the |I| row is a unit vector along the current direction, so
    # its magnitude does not vanish with |I|. This is the path that makes the
    # blanket claim "ampere measurements are ill-conditioned" wrong.
    conds = Float64[]
    for imag_a in (60.0, 20.0, 5.0, 1.0, 1e-2, 1e-4)
        x, s, p = cm_setup(imag_a * cis(-0.451))
        sv = svdvals(residual_jacobian(s, p, x))
        push!(conds, maximum(sv) / minimum(sv))
    end
    @test all(c -> 4.0 < c < 6.0, conds)          # converges to about 5.6
    @test conds[end] ≈ conds[end - 1] rtol=1e-3   # it settles rather than diverging
end

@testset "Current magnitude: |V| and |I| go collinear at unity power factor" begin
    # ...but the conditioning IS operating-point dependent, and this is the
    # degeneracy that matters in practice. |V_b| senses the direction V-hat and
    # |I| senses I-hat. When the current is in phase with the bus voltage the two
    # rows are parallel and the pair loses rank -- at ANY current magnitude.
    # Unity power factor is common (resistive load, PV at unity pf), so this is
    # not a corner case.
    conds = Float64[]
    for phi in (-0.6, -0.3, -0.1, -0.03, -0.01, -0.003)
        x, s, p = cm_setup(40.0 * cis(phi))
        sv = svdvals(residual_jacobian(s, p, x))
        push!(conds, maximum(sv) / minimum(sv))
        @test observability_diagnostics(s, p, x).unobservable_dimension == 0
    end
    # Monotone blow-up, roughly inversely with the angle between V and I.
    @test issorted(conds)
    @test conds[1] < 5.0
    @test conds[end] > 500.0

    # Exactly at unity power factor the rows are parallel: rank 1 of 2, at 40 A.
    x, s, p = cm_setup(40.0 + 0.0im)
    sv = svdvals(residual_jacobian(s, p, x))
    @test minimum(sv) < 1e-10
    d = observability_diagnostics(s, p, x)
    @test d.unobservable_dimension == 1
    @test_throws ArgumentError selected_state_covariance(s, p, x, [1])

    # A rectangular voltage pair is completely insensitive to the power factor.
    for phi in (-0.6, -0.01, 0.0)
        V = 230.0 - 0.5 * (40.0 * cis(phi))
        ms = [Measurement(kind=:vr, bus="b", reference=nothing, value=real(V), sigma=0.01),
              Measurement(kind=:vi, bus="b", reference=nothing, value=imag(V), sigma=0.01)]
        sr = compile_state_estimator(cm_net(), ms)
        sv = svdvals(residual_jacobian(sr, SEParameters(sr, ms), [real(V), imag(V)]))
        @test maximum(sv) / minimum(sv) ≈ 1.0 rtol=1e-9
    end
end

@testset "Current magnitude: uncertainty tracks the conditioning" begin
    # sigma_min is not decoration: 1/sigma_min bounds the worst-case standard
    # deviation, so the estimated angle degrades in lockstep with cond(H).
    sds = Float64[]
    for phi in (-0.6, -0.1, -0.01, -0.003)
        x, s, p = cm_setup(40.0 * cis(phi))
        V = x[1] + im * x[2]
        J = zeros(1, 2)
        J[1, 1] = -imag(V) / abs(V)^2       # d(angle V)/d(Vre)
        J[1, 2] =  real(V) / abs(V)^2       # d(angle V)/d(Vim)
        push!(sds, sqrt(derived_covariance(s, p, x, J)[1, 1]))
    end
    @test issorted(sds)
    @test sds[1] < 1e-4                     # 0.6 rad power factor: tight
    @test sds[end] > 1e-2                   # 0.003 rad: two orders worse
end
