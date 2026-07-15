using BMOPFTools: parse_bmopf
using LinearAlgebra: norm

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
end
