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
    @test diagnostics.validation_error_case_count == 1
    @test diagnostics.unexpected_error_case_count == 0
    @test diagnostics.unpublished_case_count == 0
    @test result.settings.selection_objective == :loss
    @test result.settings.continue_on_error == :validation
    @test result.settings.solver_options ==
          Dict("max_iter" => 500, "tol" => 1e-8)
    @test occursin("Ipopt", result.settings.optimizer)
    @test result.cases[1].error_type !== nothing
    @test result.cases[1].error_class == :validation
    @test inverter_control_study_case_rows(result)[1].error_class == :validation
    @test occursin("not present in network", result.cases[1].error_message)
    @test_throws ErrorException run_inverter_control_study(
        [good_baseline]; optimizer=Int)
    unexpected = run_inverter_control_study(
        [good_baseline]; optimizer=Int, continue_on_error=:all)
    @test unexpected.cases[1].error_class == :unexpected
    @test solve_diagnostics(unexpected).unexpected_error_case_count == 1
    @test only(inverter_control_study_summary_rows(unexpected))[
        "unexpected_error_case_count"] == 1
    @test_throws ArgumentError run_inverter_control_study(
        [good_baseline]; continue_on_error=:unsupported)
    @test_throws ArgumentError run_inverter_control_study(
        [good_baseline]; selection_objective=:invalid)

    device_rows = inverter_control_study_device_rows(result)
    phase_rows = inverter_control_study_phase_rows(result)
    @test length(device_rows) == 2
    @test length(phase_rows) == 6
    @test isconcretetype(eltype(device_rows))
    @test isconcretetype(eltype(phase_rows))
    @test getproperty.(device_rows, :variant_id) == ["baseline", "sequence"]
    @test all(row.metadata["penetration"] == 50 for row in device_rows)
    @test all(row.p_available_W == 8e3 for row in device_rows)

    summaries = inverter_control_study_summary_rows(
        result; group_by=["penetration"])
    @test isequal(inverter_control_study_summary_rows(
        result; group_by="penetration"), summaries)
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
    @test sequence_summary["converter_current_A_finite_points"] == 1
    @test sequence_summary[
        "capacitor_current_utilization_finite_points"] == 0
    @test isnan(sequence_summary["capacitor_current_utilization_p50"])
    # A failure fraction is uninterpretable without knowing whether the
    # non-publishable cases were infeasible or merely unconverged.
    @test baseline_summary["termination_status_counts"]["ERROR"] == 1
    @test sum(values(baseline_summary["termination_status_counts"])) ==
          baseline_summary["case_count"]
    @test baseline_summary["iteration_limit_case_count"] == 0
    @test baseline_summary["locally_infeasible_case_count"] == 0
    @test baseline_summary["numerical_error_case_count"] == 0
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
    unpublished_summary = only(
        inverter_control_study_summary_rows(unpublished))
    @test unpublished_summary["iteration_limit_case_count"] == 1
    @test unpublished_summary["error_case_count"] == 0
    @test unpublished_summary["termination_status_counts"] ==
          Dict("ITERATION_LIMIT" => 1)

    paired = inverter_control_paired_rows(result, "baseline")
    @test length(paired) == 2
    @test isconcretetype(eltype(paired))
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
    paired_summary = only(inverter_control_paired_summary_rows(
        result, "baseline"))
    @test paired_summary["candidate_pair_count"] == 2
    @test paired_summary["matched_definition_pair_count"] == 1
    @test paired_summary["publishable_pair_count"] == 1
    @test paired_summary["dropped_pair_count"] == 1
    @test paired_summary["publishable_pair_fraction"] == 0.5
    @test isfinite(paired_summary["mean_delta_converter_current_A"])
    @test paired_summary["converter_current_A_finite_pair_count"] == 1
    @test paired_summary["negative_delta_fraction_converter_current_A"] +
          paired_summary["positive_delta_fraction_converter_current_A"] +
          paired_summary[
              "indistinguishable_delta_fraction_converter_current_A"] ≈ 1.0
    @test isfinite(paired_summary["delta_converter_current_A_p50"])
    @test_throws ArgumentError inverter_control_paired_summary_rows(
        result, "baseline"; delta_relative_tolerance=-1)
    @test length(inverter_control_paired_summary_rows(
        result, "baseline"; group_by="scenario_id")) == 2
    @test length(inverter_control_paired_summary_rows(
        result, "baseline"; group_by="duration_h")) == 1
    @test length(inverter_control_paired_summary_rows(
        result, "baseline"; group_by="weight")) == 1
    @test_throws ArgumentError inverter_control_paired_rows(
        result, "not_a_variant")

    identical_variant_case = _batch_case(
        "tie", "identical", net, baseline_fleet)
    identical_baseline_case = _batch_case(
        "tie", "baseline", net, baseline_fleet)
    identical_result = InverterControlStudyResult([
        InverterControlStudyCaseResult(
            identical_baseline_case, result.cases[2].result,
            nothing, nothing, nothing, 0.0),
        InverterControlStudyCaseResult(
            identical_variant_case, result.cases[2].result,
            nothing, nothing, nothing, 0.0),
    ], result.solve, result.settings)
    tie_summary = only(inverter_control_paired_summary_rows(
        identical_result, "baseline"))
    @test tie_summary["mean_delta_converter_current_A"] == 0.0
    @test tie_summary["negative_delta_fraction_converter_current_A"] == 0.0
    @test tie_summary["positive_delta_fraction_converter_current_A"] == 0.0
    @test tie_summary[
        "indistinguishable_delta_fraction_converter_current_A"] == 1.0
    @test tie_summary["delta_converter_current_A_p05"] == 0.0
    @test tie_summary["delta_converter_current_A_p50"] == 0.0
    @test tie_summary["delta_converter_current_A_p95"] == 0.0

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
    # Both arms solved; only the pairing is invalid. Each arm therefore keeps
    # its own achieved value on the device it actually owns, so a dropped pair
    # can be inspected rather than merely counted. The difference stays NaN.
    extra_row = only(filter(
        row -> row.device_id == "extra", mismatched_pairs))
    pv_row = only(filter(row -> row.device_id == "pv", mismatched_pairs))
    @test extra_row.variant_published && !extra_row.baseline_published
    @test pv_row.baseline_published && !pv_row.variant_published
    @test isfinite(extra_row.variant_converter_current_A)
    @test isnan(extra_row.baseline_converter_current_A)
    @test isfinite(pv_row.baseline_converter_current_A)
    @test isnan(pv_row.variant_converter_current_A)
    @test all(isnan(row.delta_converter_current_A)
              for row in mismatched_pairs)
    mismatched_summary = only(inverter_control_paired_summary_rows(
        mismatched, "baseline"))
    @test mismatched_summary["dropped_pair_count"] == 2
    @test mismatched_summary["dropped_pair_known_baseline_count"] == 1

    scaled_comparison = resize_controlled_inverter_fleet(
        comparison_fleet,
        InverterHardwareSweepPoint(
            id="confounded", converter_current_scale=2.0))
    scaled_case = InverterControlStudyCase(
        scenario_id="good", variant_id="sequence", network=net,
        fleet=scaled_comparison, duration_h=good_baseline.duration_h,
        weight=good_baseline.weight, metadata=good_baseline.metadata)
    scaled_outcome = InverterControlStudyCaseResult(
        scaled_case, result.cases[3].result, nothing, nothing, nothing, 0.0)
    scaled_result = InverterControlStudyResult(
        [result.cases[2], scaled_outcome], result.solve, result.settings)
    scaled_pair = only(inverter_control_paired_rows(
        scaled_result, "baseline"))
    # A confounded pair: both arms published the same device, so both values
    # are retained for inspection, but the difference is refused.
    @test scaled_pair.baseline_published && scaled_pair.variant_published
    @test isfinite(scaled_pair.baseline_converter_current_A)
    @test isfinite(scaled_pair.variant_converter_current_A)
    @test isnan(scaled_pair.delta_converter_current_A)
    @test !scaled_pair.matched_case_definition
    @test !scaled_pair.publishable_pair

    network_case = InverterControlStudyCase(
        scenario_id="good", variant_id="sequence", network=deepcopy(net),
        fleet=comparison_fleet, duration_h=good_baseline.duration_h,
        weight=good_baseline.weight, metadata=good_baseline.metadata)
    network_outcome = InverterControlStudyCaseResult(
        network_case, result.cases[3].result, nothing, nothing, nothing, 0.0)
    network_result = InverterControlStudyResult(
        [result.cases[2], network_outcome], result.solve, result.settings)
    @test !only(inverter_control_paired_rows(
        network_result, "baseline")).matched_case_definition
    @test only(inverter_control_paired_rows(
        network_result, "baseline";
        require_shared_network=false)).matched_case_definition

    numeric_outcomes = InverterControlStudyCaseResult[]
    for penetration in (100, 0, 20, 10)
        numeric_case = InverterControlStudyCase(
            scenario_id="numeric_$penetration", variant_id="baseline",
            network=net, fleet=baseline_fleet,
            metadata=Dict("penetration" => penetration))
        push!(numeric_outcomes, InverterControlStudyCaseResult(
            numeric_case, result.cases[2].result,
            nothing, nothing, nothing, 0.0))
    end
    numeric_result = InverterControlStudyResult(
        numeric_outcomes, result.solve, result.settings)
    @test [row["penetration"] for row in
           inverter_control_study_summary_rows(
               numeric_result; group_by="penetration")] == [0, 10, 20, 100]

    long_error_case = InverterControlStudyCase(
        scenario_id="long_error", variant_id="baseline", network=bad_net,
        fleet=baseline_fleet, duration_h=9.0, weight=1.0)
    long_error = InverterControlStudyCaseResult(
        long_error_case, nothing, :validation, "ArgumentError", "bad", 0.0)
    short_good_case = InverterControlStudyCase(
        scenario_id="short_good", variant_id="baseline", network=net,
        fleet=baseline_fleet, duration_h=1.0, weight=1.0)
    short_good = InverterControlStudyCaseResult(
        short_good_case, result.cases[2].result,
        nothing, nothing, nothing, 0.0)
    exposure_result = InverterControlStudyResult(
        [long_error, short_good], result.solve, result.settings)
    @test only(inverter_control_study_summary_rows(exposure_result))[
        "weighted_publishable_fraction"] == 0.1

    exposure_case = _batch_case(
        "good", "sequence", net, comparison_fleet; duration_h=0.25)
    exposure_outcome = InverterControlStudyCaseResult(
        exposure_case, result.cases[3].result,
        nothing, nothing, nothing, 0.0)
    exposure_result = InverterControlStudyResult(
        [result.cases[2], exposure_outcome], result.solve, result.settings)
    exposure_pair = only(inverter_control_paired_rows(
        exposure_result, "baseline"))
    @test !exposure_pair.matched_case_definition
    @test !exposure_pair.publishable_pair
    @test inverter_control_study_device_rows(result)[1].metadata ===
          result.cases[2].case.metadata
end
