function _batch_case(scenario_id, variant_id, net, fleet;
                     duration_h=0.5, weight=2.0, penetration=50)
    InverterControlStudyCase(
        scenario_id=scenario_id, variant_id=variant_id,
        network=net, fleet=fleet, duration_h=duration_h, weight=weight,
        metadata=Dict("penetration" => penetration, "feeder" => "tiny"))
end

@testset "Inverter-control study runner: contracts and retained failures" begin
    request = InverterControlRequest(p_available=8e3, q_scale=0.0)
    baseline_device = _fleet_controlled(
        "pv", AverageVoltageVoltVarWatt())
    comparison_device = _fleet_controlled(
        "pv", PositiveSequenceVoltVarWatt())
    baseline_fleet = ControlledInverterFleetSpec(
        Dict("pv" => baseline_device), Dict("pv" => request))
    comparison_fleet = ControlledInverterFleetSpec(
        Dict("pv" => comparison_device), Dict("pv" => request))
    net = inv_grid3_bal()
    net["ibr"] = Dict("pv" => _fleet_native_ibr())

    @test_throws ArgumentError InverterControlStudyCase(
        scenario_id="", variant_id="baseline", network=net,
        fleet=baseline_fleet)
    @test_throws ArgumentError InverterControlStudyCase(
        scenario_id="s", variant_id="baseline", network=net,
        fleet=baseline_fleet, duration_h=0.0)
    @test_throws ArgumentError InverterControlStudyCase(
        scenario_id="s", variant_id="baseline", network=net,
        fleet=baseline_fleet, weight=0.0)
    @test_throws ArgumentError run_inverter_control_study(
        InverterControlStudyCase[])
    @test_throws ArgumentError run_inverter_control_study(Any["not a case"])
    source_metadata = Dict("nested" => [1])
    reference_case = InverterControlStudyCase(
        scenario_id="reference", variant_id="baseline", network=net,
        fleet=baseline_fleet, metadata=source_metadata)
    source_metadata["nested"][1] = 2
    @test reference_case.network === net
    @test reference_case.metadata !== source_metadata
    @test reference_case.metadata["nested"] == [1]
    @test_throws ArgumentError InverterControlStudyCase(
        scenario_id="s", variant_id="baseline", network=net,
        fleet=baseline_fleet,
        metadata=Dict{Any,Any}("duplicate" => 1, :duplicate => 2))

    good_baseline = _batch_case(
        "good", "baseline", net, baseline_fleet)
    good_comparison = _batch_case(
        "good", "sequence", net, comparison_fleet)
    bad_net = inv_grid3_bal() # no native IBR for the selected fleet id
    bad_baseline = _batch_case(
        "bad", "baseline", bad_net, baseline_fleet)
    @test_throws ArgumentError run_inverter_control_study(
        [good_baseline, good_baseline])
    @test_throws ArgumentError run_inverter_control_study(
        [bad_baseline]; continue_on_error=false)

    result = run_inverter_control_study(
        [good_comparison, bad_baseline, good_baseline];
        solver_options=("max_iter" => 500, "tol" => 1e-8))
    @test !solve_status(result).publishable
    @test getproperty.(inverter_control_study_case_rows(result), :scenario_id) ==
          ["bad", "good", "good"]
    @test getproperty.(inverter_control_study_case_rows(result), :variant_id) ==
          ["baseline", "baseline", "sequence"]
    diagnostics = solve_diagnostics(result)
    @test diagnostics.case_count == 3
    @test diagnostics.publishable_case_count == 2
    @test diagnostics.error_case_count == 1
    @test diagnostics.unpublished_case_count == 0
    @test result.settings.selection_objective == :loss
    @test result.settings.solver_options ==
          Dict("max_iter" => 500, "tol" => 1e-8)
    @test occursin("Ipopt", result.settings.optimizer)
    @test result.cases[1].error_type !== nothing
    @test occursin("not present in network", result.cases[1].error_message)

    device_rows = inverter_control_study_device_rows(result)
    phase_rows = inverter_control_study_phase_rows(result)
    @test length(device_rows) == 2
    @test length(phase_rows) == 6
    @test getproperty.(device_rows, :variant_id) == ["baseline", "sequence"]
    @test all(row.metadata["penetration"] == 50 for row in device_rows)
    @test all(row.p_available_W == 8e3 for row in device_rows)

    summaries = inverter_control_study_summary_rows(
        result; group_by=["penetration"])
    @test inverter_control_study_summary_rows(
        result; group_by="penetration") == summaries
    @test length(summaries) == 2
    baseline_summary = only(filter(
        row -> row["variant_id"] == "baseline", summaries))
    sequence_summary = only(filter(
        row -> row["variant_id"] == "sequence", summaries))
    @test baseline_summary["penetration"] == 50
    @test baseline_summary["case_count"] == 2
    @test baseline_summary["publishable_case_count"] == 1
    @test baseline_summary["error_case_count"] == 1
    @test baseline_summary["weighted_available_energy_kWh"] ≈ 8.0
    @test sequence_summary["case_count"] == 1
    @test sequence_summary["weighted_available_energy_kWh"] ≈ 8.0
    @test isfinite(sequence_summary["converter_current_A_p95"])
    @test_throws ArgumentError inverter_control_study_summary_rows(
        result; group_by=["missing_metadata"])

    failed_only = run_inverter_control_study([bad_baseline])
    failed_summary = only(inverter_control_study_summary_rows(failed_only))
    @test failed_summary["publishable_device_points"] == 0
    @test failed_summary["weighted_available_energy_kWh"] == 0.0
    @test isnan(failed_summary["converter_current_A_p95"])

    unpublished = run_inverter_control_study(
        [good_baseline]; solver_options=("max_iter" => 0,))
    @test unpublished.cases[1].result !== nothing
    @test !inverter_control_study_case_rows(unpublished)[1].publishable
    @test solve_diagnostics(unpublished).unpublished_case_count == 1
    @test solve_diagnostics(unpublished).error_case_count == 0
    @test length(inverter_control_study_device_rows(unpublished)) == 1
    @test isnan(inverter_control_study_device_rows(unpublished)[1].p_poc_W)

    paired = inverter_control_paired_rows(result, "baseline")
    @test length(paired) == 2
    @test paired[1].scenario_id == "bad"
    @test !paired[1].publishable_pair
    @test isnan(paired[1].delta_converter_current_A)
    @test paired[2].scenario_id == "good"
    @test paired[2].matched_case_definition
    @test paired[2].publishable_pair
    @test paired[2].duration_h == 0.5
    @test paired[2].weight == 2.0
    @test paired[2].metadata === result.cases[2].case.metadata
    @test abs(paired[2].delta_converter_current_A) < 0.1
    @test_throws ArgumentError inverter_control_paired_rows(
        result, "not_a_variant")

    @test PowerOptLab._weighted_quantile(
        [1.0, 2.0, 3.0], [1.0, 8.0, 1.0], 0.50) == 2.0

    extra_device = _fleet_controlled(
        "extra", AverageVoltageVoltVarWatt())
    extra_fleet = ControlledInverterFleetSpec(
        Dict("extra" => extra_device), Dict("extra" => request))
    extra_net = inv_grid3_bal()
    extra_net["ibr"] = Dict("extra" => _fleet_native_ibr())
    mismatched = run_inverter_control_study([
        _batch_case("mismatch", "baseline", net, baseline_fleet),
        _batch_case("mismatch", "sequence", extra_net, extra_fleet),
    ]; solver_options=("max_iter" => 500, "tol" => 1e-8))
    mismatched_pairs = inverter_control_paired_rows(mismatched, "baseline")
    @test getproperty.(mismatched_pairs, :device_id) == ["extra", "pv"]
    @test all(!row.publishable_pair for row in mismatched_pairs)

    exposure_case = _batch_case(
        "good", "sequence", net, comparison_fleet; duration_h=0.25)
    exposure_outcome = InverterControlStudyCaseResult(
        exposure_case, result.cases[3].result, nothing, nothing, 0.0)
    exposure_result = InverterControlStudyResult(
        [result.cases[2], exposure_outcome], result.solve, result.settings)
    exposure_pair = only(inverter_control_paired_rows(
        exposure_result, "baseline"))
    @test !exposure_pair.matched_case_definition
    @test !exposure_pair.publishable_pair
    @test inverter_control_study_device_rows(result)[1].metadata ===
          result.cases[2].case.metadata
end
