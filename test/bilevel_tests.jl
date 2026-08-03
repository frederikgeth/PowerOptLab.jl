@testset "Bilevel PV/tap DiffOpt proof of concept" begin
    net = bilevel_demo_network()

    declared_neutral_net = Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "lv1" => Dict{String,Any}(
                "terminal_names" => ["1", "ret"])),
        "terminal_conventions" => Dict{String,Any}(
            "phase" => ["1"], "neutral" => ["ret"], "earth" => String[]))
    @test PowerOptLab._bilevel_phase_terminals(
        declared_neutral_net, "lv1"; neutral_labels=Set(["ret"])) ==
        (["1"], "ret")

    three_wire_net = Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "mv" => Dict{String,Any}(
                "terminal_names" => ["a", "b", "c"])))
    @test PowerOptLab._bilevel_phase_terminals(
        three_wire_net, "mv"; voltage_measurement=:pg) ==
        (["a", "b", "c"], nothing)
    @test_throws ArgumentError PowerOptLab._bilevel_phase_terminals(
        three_wire_net, "mv"; voltage_measurement=:pn)

    no_rating = deepcopy(net)
    delete!(no_rating["ibr"]["pv1"], "s_max")
    no_rating_response = solve_bilevel_pv_response(no_rating;
        transformer_id="reg", pv_ids=["pv1", "pv2"],
        monitored_buses=["lv1", "lv2"], tap=1.0)
    @test solve_status(no_rating_response).publishable

    result = solve_bilevel_pv_tap(net;
        transformer_id="reg", pv_ids=["pv1", "pv2"],
        monitored_buses=["lv1", "lv2"], max_iterations=3)

    @test solve_status(result).publishable
    @test result.termination_status == "LOCALLY_SOLVED"
    @test result.lower_level == :aggregate
    @test 0.95 <= result.tap <= 1.05
    @test result.exported_power_W > 0.0
    @test all(216.0 - 1e-4 <= v <= 244.0 + 1e-4
              for v in values(result.voltages_V))
    @test !isempty(result.history)
    @test result.differentiability_report.ready
    @test result.differentiability_report.kkt_diagnostic.status == :accepted
    @test result.converged
    @test result.termination_reason == :at_bound

    response = solve_bilevel_pv_response(net;
        transformer_id="reg", pv_ids=["pv1", "pv2"],
        monitored_buses=["lv1", "lv2"], tap=1.0,
        lower_level=:local_controller)
    @test solve_status(response).publishable
    @test response.termination_status == "LOCALLY_SOLVED"
    @test all(isfinite, values(response.voltage_sensitivity_V_per_tap))
    @test response.differentiability_report.kkt_diagnostic.status == :accepted

    pg_response = solve_bilevel_pv_response(net;
        transformer_id="reg", pv_ids=["pv1", "pv2"],
        monitored_buses=["lv1", "lv2"], tap=1.0,
        voltage_measurement=:pg)
    @test solve_status(pg_response).publishable
    @test Set(keys(pg_response.voltages_V)) ==
          Set([("lv1", "1"), ("lv2", "1")])

    for lower_level in (:aggregate, :local_controller)
        nominal = lower_level == :local_controller ? response :
            solve_bilevel_pv_response(net;
                transformer_id="reg", pv_ids=["pv1", "pv2"],
                monitored_buses=["lv1", "lv2"], tap=1.0,
                lower_level=lower_level)
        h = 1e-5
        plus = solve_bilevel_pv_response(net;
            transformer_id="reg", pv_ids=["pv1", "pv2"],
            monitored_buses=["lv1", "lv2"], tap=1.0 + h,
            lower_level=lower_level)
        minus = solve_bilevel_pv_response(net;
            transformer_id="reg", pv_ids=["pv1", "pv2"],
            monitored_buses=["lv1", "lv2"], tap=1.0 - h,
            lower_level=lower_level)
        for bus in (("lv1", "1"), ("lv2", "1"))
            fd = (plus.voltages_V[bus] - minus.voltages_V[bus]) / (2h)
            @test isapprox(nominal.voltage_sensitivity_V_per_tap[bus], fd;
                           atol=2e-3, rtol=2e-5)
        end
    end

    base, triples = PowerOptLab._bilevel_breakpoints_to_triples(
        [216.0, 225.0, 235.0, 244.0], [0.44, 0.0, 0.0, -0.44])
    curve_value(u) = base + sum(a * max(0.0, u - x) for (a, x) in triples)
    @test curve_value(200.0) ≈ 0.44
    @test curve_value(230.0) ≈ 0.0 atol=1e-12
    @test curve_value(250.0) ≈ -0.44

    low_start = solve_bilevel_pv_tap(net;
        transformer_id="reg", pv_ids=["pv1", "pv2"],
        monitored_buses=["lv1", "lv2"], tap_initial=0.97, max_iterations=8)
    @test low_start.converged
    @test low_start.termination_reason == :at_bound
    @test solve_status(low_start).publishable

    interior = solve_bilevel_pv_tap(net;
        transformer_id="reg", pv_ids=["pv1", "pv2"],
        monitored_buses=["lv1", "lv2"], voltage_reference=234.0,
        tap_penalty=0.0, max_iterations=20)
    @test interior.converged
    @test interior.termination_reason == :tap_tolerance
    @test 0.95 < interior.tap < 1.05

    local_controller = solve_bilevel_pv_tap(net;
        transformer_id="reg", pv_ids=["pv1", "pv2"],
        monitored_buses=["lv1", "lv2"], lower_level=:local_controller,
        max_iterations=2)
    @test solve_status(local_controller).publishable
    @test local_controller.termination_status == "LOCALLY_SOLVED"
    @test local_controller.lower_level == :local_controller
    @test local_controller.exported_power_W > 0.0
    @test all(216.0 - 1e-4 <= v <= 244.0 + 1e-4
              for v in values(local_controller.voltages_V))
    @test local_controller.differentiability_report.kkt_diagnostic.status == :accepted
    @test isfinite(local_controller.exported_power_W)
    @test local_controller.exported_power_W <=
          result.exported_power_W + 2e-5 * 12_000.0
    @test all(local_controller.lower_result["ibr"][id]["1"]["pg"] <=
              net["ibr"][id]["p_max"][1] + 1e-3
              for id in ("pv1", "pv2"))
    @test net["control_profile"]["vvw"]["volt_watt"]["p_limits"] == [0.0, 1.0]
    @test_throws ArgumentError solve_bilevel_pv_tap(net;
        transformer_id="reg", pv_ids=["pv1", "pv2"],
        monitored_buses=["lv1", "lv2"], lower_level=:not_a_mode,
        max_iterations=1)

    centralized = solve_single_level_pv_tap(net;
        transformer_id="reg", pv_ids=["pv1", "pv2"],
        monitored_buses=["lv1", "lv2"])
    @test solve_status(centralized).publishable
    @test centralized.termination_status == "LOCALLY_SOLVED"
    @test 0.95 <= centralized.tap <= 1.05
    @test centralized.exported_power_W >= 0.0
    @test all(216.0 - 1e-4 <= v <= 244.0 + 1e-4
              for v in values(centralized.voltages_V))
    @test isfinite(centralized.objective)
end
