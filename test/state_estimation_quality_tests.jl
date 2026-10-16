using Test, PowerOptLab, LinearAlgebra
using BMOPFTools: parse_bmopf

function quality_se_net(; grounded=false)
    bond = grounded ? ",\"perfectly_grounded_terminals\":[\"n\"]" : ""
    parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
            "b":{"terminal_names":["1","n"]$bond}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.5,"R_series_2_2":0.5}},
     "line":{"l":{"bus_from":"src","bus_to":"b","terminal_map_from":["1","n"],
                    "terminal_map_to":["1","n"],"linecode":"lc","length":1.0}}}
    """; from_string=true)
end

@testset "DSSE review: exact devices respect ideal-source boundaries" begin
    net = quality_se_net(; grounded=true)
    downstream = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("b", "1"), nothing)], ComplexF64[2]))
    at_source = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("src", "1"), nothing)], ComplexF64[3]))
    devices = [downstream, at_source]
    s = compile_state_estimator(net; exact_devices=devices)
    p = SEParameters(s; exact_devices=devices)
    # The ideal source supplies both devices. Its current is a free boundary
    # reaction, not zero: only downstream KCL constrains the voltage state.
    x = [229., 0.]
    @test length(evaluate_state_estimator(s, p, x).constraints) == 2
    @test norm(evaluate_state_estimator(s, p, x).constraints) < 1e-10
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)
        r = solver(s, p, [230., 0.])
        @test solve_status(r).publishable
        @test r.state ≈ x atol=1e-6
    end
    floating = quality_se_net()
    through_source = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("src", "1"), ("b", "n"))], ComplexF64[3]))
    load = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("b", "1"), ("b", "n"))], ComplexF64[2]))
    sf = compile_state_estimator(floating; exact_devices=[load, through_source])
    pf = SEParameters(sf; exact_devices=[load, through_source])
    expected = zeros(4)
    expected[sf.free_state_map[("b", "1")]] = 229.
    expected[sf.free_state_map[("b", "n")]] = 2.5
    @test norm(evaluate_state_estimator(sf, pf, expected).constraints) < 1e-10
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)
        r = solver(sf, pf, initial_state_estimator(sf, pf))
        @test solve_status(r).publishable
        @test r.state ≈ expected atol=1e-6
    end
end

@testset "DSSE review: support partially exact shared injection nodes" begin
    net = quality_se_net()
    known = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("b", "1"), ("b", "n"))], ComplexF64[2]))
    # A second customer's unknown current at an aliased terminal would also
    # enter the exact-device KCL rows. The voltage-only model cannot pretend
    # that the first customer's exact law fixes the whole node injection.
    net["bus"]["z"] = deepcopy(net["bus"]["b"])
    net["line"]["alias"] = Dict{String,Any}("bus_from"=>"b", "bus_to"=>"z",
        "terminal_map_from"=>["1", "n"], "terminal_map_to"=>["1", "n"],
        "linecode"=>"lc", "length"=>0.)
    other = [Measurement(kind=:pinj, bus="z", value=-1760., sigma=10.),
             Measurement(kind=:qinj, bus="z", value=0., sigma=10.)]
    @test state_estimator_preflight(net, other; exact_devices=[known]).supported
        se=compile_state_estimator(net,other;exact_devices=[known])
    pe=SEParameters(se,other;exact_devices=[known])
    for solver in (solve_compiled_state_estimator,solve_sparse_state_estimator,solve_ipopt_state_estimator)
        result=solver(se,pe,initial_state_estimator(se,pe))
        @test solve_status(result).publishable
        @test result.evaluation.voltage[se.node_index[("b","1")]]≈225 atol=1e-7
        @test result.evaluation.voltage[se.node_index[("b","n")]]≈5 atol=1e-7
    end
    # The same protection applies to connections declared without telemetry.
    @test state_estimator_preflight(net; exact_devices=[known],
        injection_connections=[TerminalConnection(("z", "1"), ("z", "n"))]).supported
    # Metering the exact device itself is supported.
    own = [Measurement(kind=:pinj, bus="b", value=-456., sigma=10.),
           Measurement(kind=:qinj, bus="b", value=0., sigma=10.)]
    @test state_estimator_preflight(net, own; exact_devices=[known]).supported
    known_alias = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("z", "1"), ("z", "n"))], ComplexF64[8]))
    @test state_estimator_preflight(net, vcat(own, other);
        exact_devices=[known, known_alias]).supported

    # A shared ungrounded return also couples distinct phase conductors.
    for bus in ("src", "b")
        push!(net["bus"][bus]["terminal_names"], "2")
    end
    push!(net["voltage_source"]["s"]["terminal_map"], "2")
    push!(net["voltage_source"]["s"]["v_magnitude"], 230.)
    push!(net["voltage_source"]["s"]["v_angle"], -2pi/3)
    push!(net["line"]["l"]["terminal_map_from"], "2")
    push!(net["line"]["l"]["terminal_map_to"], "2")
    net["linecode"]["lc"]["R_series_3_3"] = .5
    unknown = [TerminalConnection(("b", "2"), ("b", "n"))]
    @test state_estimator_preflight(net; exact_devices=[known], injection_connections=unknown).supported
    # Sharing ideal earth is safe: it has no KCL/voltage state to overconstrain.
    for bus in ("b", "z")
        net["bus"][bus]["perfectly_grounded_terminals"] = ["n"]
    end
    @test state_estimator_preflight(net; exact_devices=[known], injection_connections=unknown).supported
end

@testset "DSSE quality: parameter updates preserve measurement identity" begin
    net = quality_se_net()
    ms = [Measurement(kind=:vr, bus="b", reference=nothing, value=225., sigma=1.),
          Measurement(kind=:vi, bus="b", reference=nothing, value=1., sigma=2.)]
    s = compile_state_estimator(net, ms)
    @test_throws ArgumentError SEParameters(s, reverse(ms))
    different_ref = [Measurement(kind=:vr, bus="b", value=225., sigma=1.), ms[2]]
    @test_throws ArgumentError SEParameters(s, different_ref)
    moved = [Measurement(kind=:vr, bus="src", reference=nothing, value=230., sigma=1.), ms[2]]
    @test_throws ArgumentError SEParameters(s, moved)
    changed = [Measurement(kind=m.kind, bus=m.bus, terminal=m.terminal, reference=m.reference,
                           value=m.value+3, sigma=2m.sigma) for m in ms]
    p = SEParameters(s, changed)
    @test p.measurement_values == [228., 4.]
    @test p.covariance_values == [2., 4.]
    branches = [BranchMeasurement(kind=:ire, line="l", side=side, terminal="1", value=0., sigma=1.) for side in (:from, :to)]
    sb = compile_state_estimator(net, branches)
    @test_throws ArgumentError SEParameters(sb, reverse(branches))
end

@testset "DSSE quality: exact device topology and finite parameters" begin
    net = quality_se_net()
    connection = TerminalConnection(("b", "1"), ("b", "n"))
    device = ExactDeviceEquation(ConstantCurrentDevice([connection], ComplexF64[2+im]))
    s = compile_state_estimator(net; exact_devices=[device])
    updated = ExactDeviceEquation(ConstantCurrentDevice([connection], ComplexF64[3+2im]))
    @test SEParameters(s; exact_devices=[updated]).device_currents == [3+2im]
    wrong_connection = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("b", "1"), nothing)], ComplexF64[2+im]))
    @test_throws ArgumentError SEParameters(s; exact_devices=[wrong_connection])
    wrong_kind = ExactDeviceEquation(ConstantPowerDevice([connection], ComplexF64[1000]))
    @test_throws ArgumentError SEParameters(s; exact_devices=[wrong_kind])
    @test_throws ArgumentError ConstantPowerDevice([connection], ComplexF64[NaN])
    @test_throws ArgumentError ConstantCurrentDevice([connection], ComplexF64[Inf])
    @test_throws ArgumentError ZIPDevice([connection], ComplexF64[0], ComplexF64[0], ComplexF64[NaN])
    p = SEParameters(s; exact_devices=[device])
    x = initial_state_estimator(s, p)
    p.device_currents[1] = Inf
    @test_throws ArgumentError evaluate_state_estimator(s, p, x)
    @test_throws ArgumentError constraint_jacobian(s, p, x)
end

@testset "DSSE quality: named earth returns and zero source neutrals" begin
    net = quality_se_net(; grounded=true)
    explicit = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("b", "1"), ("b", "n"))], ComplexF64[2]))
    implicit = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("b", "1"), nothing)], ComplexF64[2]))
    report = state_estimator_preflight(net; exact_devices=[explicit])
    @test report.supported
    if report.supported
        se = compile_state_estimator(net; exact_devices=[explicit])
        si = compile_state_estimator(net; exact_devices=[implicit])
        pe = SEParameters(se; exact_devices=[explicit]); pi_ = SEParameters(si; exact_devices=[implicit])
        x = [229., 0.]
        @test norm(evaluate_state_estimator(se, pe, x).constraints) < 1e-10
        @test evaluate_state_estimator(se, pe, x).constraints ≈ evaluate_state_estimator(si, pi_, x).constraints
        @test constraint_jacobian(se, pe, x) ≈ constraint_jacobian(si, pi_, x)
    end
    source = net["voltage_source"]["s"]
    push!(source["terminal_map"], "n"); push!(source["v_magnitude"], 0.); push!(source["v_angle"], 0.)
    @test state_estimator_preflight(net).supported
    source["v_magnitude"][2] = 1.
    @test !state_estimator_preflight(net).supported
    absent = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("b", "1"), ("b", "absent"))], ComplexF64[2]))
    @test !state_estimator_preflight(quality_se_net(); exact_devices=[absent]).supported
end

@testset "DSSE quality: continuation stages are not physical estimates" begin
    net = quality_se_net(; grounded=true)
    device = ExactDeviceEquation(ConstantPowerDevice(
        [TerminalConnection(("b", "1"), nothing)], ComplexF64[1000]))
    s = compile_state_estimator(net; exact_devices=[device])
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)
        p = SEParameters(s; exact_devices=[device], regularization_voltage=100., continuation_alpha=0.)
        r = solver(s, p, [230., 0.])
        @test r.status == :converged_regularized
        @test !solve_status(r).publishable
        physical = SEParameters(s; exact_devices=[device])
        @test norm(evaluate_state_estimator(s, physical, r.state).constraints) > .01
    end
    p = SEParameters(s; exact_devices=[device], regularization_voltage=100.)
    failed = solve_with_continuation(s, p, [230., 0.]; alphas=[0., 1.], max_iterations=1, initial_radius=1e-5)
    @test failed.status == :power_flow_initialization_failed
    @test p.continuation_alpha == 1. # failed continuation must restore caller data
    @test !solve_status(failed).publishable
    p.continuation_alpha = .4
    @test_throws ArgumentError solve_with_continuation(s, p, [230., 0.]; max_iterations=0)
    @test p.continuation_alpha == .4
    result = solve_with_continuation(s, p, [230., 0.]; alphas=[0., .5, 1.])
    @test result.status == :converged_unique
    @test all(r -> r.status == :converged_regularized && !solve_status(r).publishable, result.stages[1:end-1])
    @test solve_status(result).publishable
    @test p.continuation_alpha == 1.
    regularized = SEParameters(s; exact_devices=[device], continuation_alpha=0.)
    series = solve_time_series_state_estimator(s, [regularized], [230., 0.])
    @test series.status == :time_series_stalled
    @test series.completed_snapshots == 0
    # Alpha has no effect on a pure constant-current law, which remains physical.
    current = ExactDeviceEquation(ConstantCurrentDevice(
        [TerminalConnection(("b", "1"), nothing)], ComplexF64[2]))
    sc = compile_state_estimator(net; exact_devices=[current])
    pc = SEParameters(sc; exact_devices=[current], continuation_alpha=0.)
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)
        @test solve_status(solver(sc, pc, [230., 0.])).publishable
    end
end

@testset "DSSE quality: solver options cannot disable feasibility silently" begin
    net = quality_se_net(; grounded=true)
    device = ExactDeviceEquation(ConstantPowerDevice(
        [TerminalConnection(("b", "1"), nothing)], ComplexF64[1000]))
    s = compile_state_estimator(net; exact_devices=[device])
    p = SEParameters(s; exact_devices=[device])
    # At 230 V the exact load equation is violated by over four amperes.
    # Infinite feasibility tolerance used to publish that initial point.
    @test norm(evaluate_state_estimator(s, p, [230., 0.]).constraints) > 4.
    invalid = [(; constraint_tolerance=Inf), (; constraint_tolerance=-1.),
               (; constraint_rtol=Inf), (; constraint_rtol=NaN),
               (; optimality_tolerance=Inf), (; optimality_tolerance=-1.),
               (; rank_rtol=NaN), (; rank_rtol=-1.), (; min_radius=0.),
               (; initial_radius=Inf, max_radius=Inf), (; penalty=Inf)]
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator), options in invalid
        @test_throws ArgumentError solver(s, p, [230., 0.]; options...)
    end
    @test_throws ArgumentError solve_sparse_state_estimator(s, p, [230., 0.]; damping=Inf)
    ms = [Measurement(kind=:vr, bus="b", reference=nothing, value=230., sigma=1.),
          Measurement(kind=:vi, bus="b", reference=nothing, value=0., sigma=1.)]
    sl = compile_state_estimator(net, ms); pl = SEParameters(sl, ms)
    for solver in (solve_compiled_state_estimator, solve_sparse_state_estimator, solve_ipopt_state_estimator)
        @test solve_status(solver(sl, pl, [230., 0.]; constraint_tolerance=0., constraint_rtol=0.,
                                  optimality_tolerance=0., rank_rtol=0.)).publishable
    end
end
