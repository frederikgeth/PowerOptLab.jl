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

    pseudo_trace = continue_opf_operability_pseudo_arclength(net, pf;
        spec=spec,
        continuation=OperabilityPseudoArclengthSpec(
            initial_step=0.05, max_step=0.1, target_lambda_tol=0.01))
    @test pseudo_trace.status == :pass
    @test first(pseudo_trace.lambdas) == 0.0
    @test abs(last(pseudo_trace.lambdas) - 1.0) <= 0.01
    @test pseudo_trace.endpoint_match === true
    @test pseudo_trace.endpoint_distance < 1e-4
    @test maximum(pseudo_trace.residuals) < 1e-5
    @test pseudo_trace.provenance["continuation"]["pseudo_arclength"] === true
    @test all(>(0.0), pseudo_trace.provenance["continuation"]["arclength_state_scale"])
    @test_throws ArgumentError OperabilityPseudoArclengthSpec(initial_step=0.01,
                                                               min_step=0.1)

    # Analytic two-bus resistive feeder: P = V(E - V) / R has a high branch
    # at 600 V and a low branch at 400 V for P = 2.4 MW.  Both are valid
    # equilibria, but only the high branch is connected to the energized germ
    # before the λ=1 target crossing.
    nose_net = single_bus_net(pload=2.4e6)
    nose_pf = solve_pf(nose_net; per_unit=false)
    @test nose_pf["bus"]["bus1"]["1"]["vr"] ≈ 600.0 atol=1e-6
    nose_trace = continue_opf_operability_pseudo_arclength(nose_net, nose_pf;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
        continuation=OperabilityPseudoArclengthSpec(
            initial_step=0.05, max_step=0.1, max_steps=40, target_lambda_tol=0.01))
    @test nose_trace.status == :pass
    @test nose_trace.endpoint_match === true
    @test any(get(event, "kind", "") == "target_crossing" for event in nose_trace.events)
    @test any(get(event, "kind", "") == "target_refinement" &&
              get(event, "status", nothing) == :pass for event in nose_trace.events)

    low_solution = deepcopy(nose_pf)
    low_solution["bus"]["bus1"]["1"]["vr"] = 400.0
    low_report = check_opf_operability(nose_net, low_solution;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()))
    @test low_report.status == :pass
    low_trace = continue_opf_operability_pseudo_arclength(nose_net, low_solution;
        spec=OperabilitySpec(scaling_policy=SIUnitsScaling()),
        continuation=OperabilityPseudoArclengthSpec(
            initial_step=0.05, max_step=0.1, max_steps=40, target_lambda_tol=0.01))
    @test low_trace.status == :fail
    @test low_trace.endpoint_match === false
    @test low_trace.endpoint_distance ≈ 200.0 atol=1e-6
    @test any(get(event, "kind", "") == "target_refinement" &&
              get(event, "status", nothing) == :pass for event in low_trace.events)

    fold_guess = deepcopy(nose_pf)
    fold_guess["bus"]["bus1"]["1"]["vr"] = 500.0
    fold = locate_opf_operability_fold(nose_net, fold_guess;
        lambda=1.04, spec=OperabilitySpec(scaling_policy=SIUnitsScaling()))
    @test fold.status == :pass
    @test fold.lambda ≈ 25 / 24 atol=1e-7
    @test fold.node_voltages[("bus1", "1")] ≈ 500.0 + 0.0im atol=1e-5
    @test fold.residual_norm < 1e-8
    @test fold.sigma_min < 1e-8
    @test fold.critical_mode["coordinate_order"] == "[real(state_nodes); imag(state_nodes)]"
    @test fold.provenance["fold_localization"]["equations"] == "F=0,Jv=0,norm(v)=1"

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
    ibr_pseudo_trace = continue_opf_operability_pseudo_arclength(
        ibr_net, ibr_pf; spec=spec)
    @test ibr_pseudo_trace.status == :not_applicable
    @test ibr_pseudo_trace.events[1]["kind"] == "unsupported_physics"
    ibr_fold = locate_opf_operability_fold(ibr_net, ibr_pf;
        spec=spec, lambda=1.0)
    @test ibr_fold.status == :not_applicable
    @test_throws ArgumentError locate_opf_operability_fold(
        net, pf; spec=spec, lambda=0.0)
end
