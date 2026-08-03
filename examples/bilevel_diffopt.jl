# Bilevel distribution-network POC.
#
# Run from the package root:
#   julia --project=. examples/bilevel_diffopt.jl

using PowerOptLab

net = bilevel_demo_network()
aggregate = solve_bilevel_pv_tap(net;
    transformer_id="reg", pv_ids=["pv1", "pv2"],
    monitored_buses=["lv1", "lv2"], max_iterations=12)

local_controller = solve_bilevel_pv_tap(net;
    transformer_id="reg", pv_ids=["pv1", "pv2"],
    monitored_buses=["lv1", "lv2"], lower_level=:local_controller,
    max_iterations=12)

centralized = solve_single_level_pv_tap(net;
    transformer_id="reg", pv_ids=["pv1", "pv2"],
    monitored_buses=["lv1", "lv2"])

println("Aggregate lower-level bilevel response")
println("  tap:       ", round(aggregate.tap; digits=6))
println("  PV export: ", round(aggregate.exported_power_W; digits=2), " W")
println("  voltages:  ", aggregate.voltages_V)
println("  report:    ", aggregate.differentiability_report.ready,
        " / ", aggregate.differentiability_report.kkt_diagnostic.status)

println("Local-controller lower-level bilevel response")
println("  tap:       ", round(local_controller.tap; digits=6))
println("  PV export: ", round(local_controller.exported_power_W; digits=2), " W")
println("  voltages:  ", local_controller.voltages_V)
println("  report:    ", local_controller.differentiability_report.ready,
        " / ", local_controller.differentiability_report.kkt_diagnostic.status)

println("Centralized single-level benchmark")
println("  tap:       ", round(centralized.tap; digits=6))
println("  PV export: ", round(centralized.exported_power_W; digits=2), " W")
println("  voltages:  ", centralized.voltages_V)
