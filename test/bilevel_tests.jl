@testset "Bilevel PV/tap DiffOpt proof of concept" begin
    net = bilevel_demo_network()
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
