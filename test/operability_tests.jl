using BMOPFTools

@testset "post-OPF operability checker" begin
    net = single_bus_net(pload=100.0)
    pf = solve_pf(net; per_unit=false)

    @test_throws ArgumentError check_opf_operability(net, pf)

    spec = OperabilitySpec(
        scaling_policy=SIUnitsScaling(),
        voltage_min=800.0,
        voltage_max=1100.0,
    )
    report = check_opf_operability(net, pf; spec)
    @test report.status == :pass
    @test report.checks["endpoint"].status == :pass
    @test report.endpoint_residual < 1e-6
    @test report.checks["jacobian_regular"].status == :pass
    @test report.checks["load_scale_sensitivity"].status == :pass
    connection = report.load_connections["ld1/1"]
    @test connection["requested_power"] ≈ connection["realized_power"] atol=1e-8
    load_scale_connection = report.sensitivities["load_scale"]["load_connections"]["ld1/1"]
    @test load_scale_connection["magnitude_derivative"] < 0.0
    @test isfinite(load_scale_connection["path_dP_dV"])
    @test load_scale_connection["path_dP_dV"] < 0.0
    @test haskey(report.sensitivities["directions"], "P")
    @test haskey(report.sensitivities["directions"], "Q")
    p_direction = report.sensitivities["directions"]["P"]["ld1/1"]
    @test p_direction["units"] == "W"
    @test p_direction["load_connections"]["ld1/1"]["magnitude_derivative"] < 0.0
    @test report.provenance["operability"]["scope"] == "static_ybus_linearized"
    critical = report.branch_evidence["critical_mode"]
    @test critical["status"] == :pass
    @test length(critical["left_vector"]) == 2 * length(report.state_nodes)
    @test length(critical["right_node_participation"]) == length(report.state_nodes)

    bad = deepcopy(pf)
    bad["bus"]["bus1"]["1"]["vr"] += 10.0
    bad_report = check_opf_operability(net, bad; spec)
    @test bad_report.status == :fail
    @test bad_report.checks["endpoint"].status == :fail

    # The balanced fixture uses the phase ordering/convention shared with the
    # existing inverter symmetrical-component helper.
    net3 = inv_grid3_bal()
    pf3 = solve_pf(net3; per_unit=false)
    report3 = check_opf_operability(net3, pf3;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(), vuf_max=0.01))
    @test report3.status == :pass
    @test report3.checks["sequence_unbalance"].status == :pass
    @test report3.sequences["poc"]["vuf"] < 1e-4
    @test abs(report3.sequences["poc"]["v1"]) > 200.0

    delta_net = inv_grid3_bal()
    delta_net["load"] = Dict("dΔ" => Dict{String,Any}(
        "bus" => "poc", "terminal_map" => ["a", "b", "c"],
        "configuration" => "DELTA", "p_nom" => [100.0, 100.0, 100.0],
        "q_nom" => [0.0, 0.0, 0.0]))
    delta_pf = solve_pf(delta_net; per_unit=false)
    delta_report = check_opf_operability(delta_net, delta_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()))
    @test delta_report.status == :pass
    @test length(delta_report.load_connections) == 3
    @test all(r["positive"] !== r["negative"] for r in values(delta_report.load_connections))

    helm_report = check_opf_operability(net, pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(), compute_helm=true))
    @test helm_report.checks["helm_reachability"].status == :pass
    @test helm_report.branch_evidence["reachability"]["base"] == "energized_no_load_germ"
    @test helm_report.branch_evidence["reachability"]["endpoint_mismatch"] < 1e-4

    current_net = single_bus_net(pload=100.0)
    current_net["load"]["ld1"]["model"] = "constant_current"
    current_net["load"]["ld1"]["v_nom"] = [1000.0]
    current_pf = solve_pf(current_net; per_unit=false)
    current_report = check_opf_operability(current_net, current_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling(), compute_helm=true))
    @test current_report.checks["helm_reachability"].status == :not_applicable
    @test any(occursin("constant_current", reason) for reason in
        current_report.branch_evidence["reachability"]["reasons"])

    trace = continue_opf_operability(net, pf;
        spec=spec,
        continuation=OperabilityContinuationSpec(initial_step=0.2, max_step=0.25))
    @test trace.status == :pass
    @test first(trace.lambdas) == 0.0
    @test last(trace.lambdas) == 1.0
    @test trace.endpoint_match === true
    @test trace.endpoint_distance < 1e-4
    @test maximum(trace.residuals) < 1e-5
    @test trace.provenance["continuation"]["pseudo_arclength"] === false
    @test_throws ArgumentError OperabilityContinuationSpec(initial_step=0.01, min_step=0.1)

    # Scope exclusions are explicit evidence rather than a silent partial pass.
    ibr_net = doe_ibr_feeder()
    ibr_pf = solve_pf(ibr_net; per_unit=false)
    ibr_report = check_opf_operability(ibr_net, ibr_pf; spec)
    @test ibr_report.status == :not_applicable
    @test !isempty(ibr_report.unsupported)
    @test ibr_report.checks["scope"].status == :not_applicable
    ibr_trace = continue_opf_operability(ibr_net, ibr_pf; spec=spec)
    @test ibr_trace.status == :not_applicable
    @test ibr_trace.events[1]["kind"] == "unsupported_physics"
end
