using BMOPFTools: add_statcom!
using Dates

@testset "Operating envelope: equal allocation is uniform, voltage-limited, load-dependent" begin
    cps = [ConnectionPoint(id="der1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="der2", bus="bus2", export_max=10e3)]
    # Two intervals: low baseline load, then high. Higher local load absorbs
    # export and frees voltage headroom, so the envelope should grow.
    nets = [doe_feeder(p1=200.0, p2=200.0), doe_feeder(p1=5000.0, p2=5000.0)]

    res = solve_operating_envelope(nets, cps; fairness=:equal)
    @test all(s in ("LOCALLY_SOLVED", "OPTIMAL") for s in res.termination_status)

    for t in 1:2
        e1 = res.envelope["der1"][t]; e2 = res.envelope["der2"][t]
        @test e1 ≈ e2  rtol=1e-3                       # equal allocation
        @test 0.0 <= e1 <= 10e3                         # within the inverter cap
        @test e1 < 10e3 - 1.0                           # bound by voltage, not the cap
        # The far bus sits at its v_max — the binding operational constraint.
        vmax_bus = max(res.snapshots[t]["bus"]["bus1"]["1"]["vm"],
                       res.snapshots[t]["bus"]["bus2"]["1"]["vm"])
        @test vmax_bus ≈ 245.0  atol=0.05
    end

    # Dynamic response: the high-load interval admits a larger envelope.
    @test res.envelope["der1"][2] > res.envelope["der1"][1] + 100.0
    @test res.total_export[2] > res.total_export[1]
end

@testset "Operating envelope: :sum is more efficient but less equitable than :equal" begin
    cps = [ConnectionPoint(id="der1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="der2", bus="bus2", export_max=10e3)]
    nets = [doe_feeder(p1=200.0, p2=200.0), doe_feeder(p1=5000.0, p2=5000.0)]

    eq  = solve_operating_envelope(nets, cps; fairness=:equal)
    sm  = solve_operating_envelope(nets, cps; fairness=:sum)
    @test all(s in ("LOCALLY_SOLVED", "OPTIMAL") for s in sm.termination_status)

    for t in 1:2
        # Maximising the total allocates at least as much in aggregate as the
        # equitable rule.
        @test sm.total_export[t] >= eq.total_export[t] - 1.0
        # …and every allocation still respects the inverter cap.
        @test sm.envelope["der1"][t] <= 10e3 + 1.0
        @test sm.envelope["der2"][t] <= 10e3 + 1.0
    end
    # At low load the efficient rule is visibly uneven (the electrically stronger
    # point near the source gets most of the allocation).
    @test sm.envelope["der1"][1] > sm.envelope["der2"][1] + 1000.0
end

@testset "Operating envelope: proportional fairness is the middle ground" begin
    cps = [ConnectionPoint(id="der1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="der2", bus="bus2", export_max=10e3)]
    nets = [doe_feeder(p1=200.0, p2=200.0), doe_feeder(p1=5000.0, p2=5000.0)]

    eq = solve_operating_envelope(nets, cps; fairness=:equal)
    pr = solve_operating_envelope(nets, cps; fairness=:proportional)
    sm = solve_operating_envelope(nets, cps; fairness=:sum)
    @test all(s in ("LOCALLY_SOLVED", "OPTIMAL") for s in pr.termination_status)

    for t in 1:2
        e1 = pr.envelope["der1"][t]; e2 = pr.envelope["der2"][t]
        @test e1 > 100.0 && e2 > 100.0                 # no point is starved
        @test e1 >= e2 - 1.0                           # stronger point gets ≥ weaker
        @test e1 <= 10e3 + 1.0 && e2 <= 10e3 + 1.0     # within the inverter cap
        # Total sits between the equitable and efficient extremes.
        @test pr.total_export[t] >= eq.total_export[t] - 1.0
        @test pr.total_export[t] <= sm.total_export[t] + 1.0
    end

    # Where :sum starves the weak point (low load), proportional keeps it well
    # above zero — the defining property of proportional fairness.
    @test pr.envelope["der2"][1] > sm.envelope["der2"][1] + 1000.0
end

@testset "Operating envelope: single-net convenience + input validation" begin
    cps = [ConnectionPoint(id="der1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="der2", bus="bus2", export_max=10e3)]
    res = solve_operating_envelope(doe_feeder(p1=1000.0, p2=1000.0), cps; fairness=:equal)
    @test length(res.termination_status) == 1
    @test length(res.envelope["der1"]) == 1
    @test res.total_export[1] ≈ res.envelope["der1"][1] + res.envelope["der2"][1]  rtol=1e-9

    @test_throws ArgumentError solve_operating_envelope(
        [doe_feeder(p1=1000.0, p2=1000.0)], cps; fairness=:bogus)
    dup = [ConnectionPoint(id="x", bus="bus1", export_max=1e3),
           ConnectionPoint(id="x", bus="bus2", export_max=1e3)]
    @test_throws ArgumentError solve_operating_envelope(
        [doe_feeder(p1=1000.0, p2=1000.0)], dup)
end

@testset "Operating envelope: failure-safe extraction and strict validation" begin
    net = doe_feeder(p1=200.0, p2=200.0)
    @test_throws ArgumentError solve_operating_envelope(net, ConnectionPoint[])
    @test_throws ArgumentError solve_operating_envelope(net,
        [ConnectionPoint(id="bad", bus="bus1", export_max=-1.0)])
    @test_throws ArgumentError solve_operating_envelope(net,
        [ConnectionPoint(id="bad", bus="missing", export_max=1.0)])
    @test_throws ArgumentError solve_operating_envelope(net,
        [ConnectionPoint(id="bad", bus="bus1", phase_terminals=String[], export_max=1.0)])
    @test_throws ArgumentError solve_operating_envelope(net,
        [ConnectionPoint(id="bad", bus="bus1", export_max=1.0)];
        fairness=FairnessPolicy(kind=:equal, normalization=:request))
    @test_throws ArgumentError solve_operating_envelope(net,
        [ConnectionPoint(id="bad", bus="bus1", export_max=1.0)];
        fairness=FairnessPolicy(weights=Dict("unknown"=>1.0)))

    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    infeasible = solve_operating_envelope(
        doe_feeder(p1=200.0, p2=200.0, vmax=200.0), cps)
    @test !infeasible.diagnostics[1]["feasible"]
    @test all(isnan(infeasible.envelope[id][1]) for id in ("d1", "d2"))
    @test isnan(infeasible.total_capacity[1])
    @test infeasible.snapshots[1]["primal_status"] == "INFEASIBLE_POINT"
end

@testset "Operating envelope: import direction and result semantics" begin
    cp = ConnectionPoint(id="battery", bus="bus2", import_max=20e3)
    r = solve_operating_envelope(doe_feeder(p1=200.0, p2=200.0), [cp];
                                 direction=:import)
    @test r.direction == :import
    @test r.diagnostics[1]["direction"] == :import
    @test 0.0 < r.envelope["battery"][1] < 20e3
    @test r.total_capacity == r.total_export  # compatibility alias
end

@testset "Operating envelope: scenario sharing and all-corner security" begin
    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    scenarios = [[doe_feeder(p1=200.0, p2=200.0),
                  doe_feeder(p1=1000.0, p2=1000.0)]]
    r = solve_operating_envelope(scenarios, cps; security=:corners)
    @test all(isfinite(r.envelope[id][1]) for id in ("d1", "d2"))
    @test r.diagnostics[1]["scenario_count"] == 2
    @test r.diagnostics[1]["dispatch_points_per_scenario"] == 4
    @test r.diagnostics[1]["security_scope"] == :all_box_corners
    @test r.diagnostics[1]["guarantee"] == :local_ac_feasibility_at_tested_dispatches
    @test_throws ArgumentError solve_operating_envelope(scenarios, cps;
        security=:corners, max_exact_corners=1)

    # A forecast scenario whose zero-DER corner is infeasible must invalidate the
    # whole box instead of publishing the feasible all-export endpoint.
    bad_scenarios = [[doe_feeder(p1=200.0, p2=200.0),
                      doe_feeder(p1=5000.0, p2=5000.0)]]
    bad = solve_operating_envelope(bad_scenarios, cps; security=:corners)
    @test !bad.diagnostics[1]["feasible"]
    @test all(isnan(bad.envelope[id][1]) for id in ("d1", "d2"))
end

@testset "Operating envelope: parameterized and normalized fairness" begin
    cps = [ConnectionPoint(id="large", bus="bus1", export_max=10e3,
                           requested=8e3, normalization=5e3),
           ConnectionPoint(id="small", bus="bus2", export_max=5e3,
                           requested=4e3, normalization=2.5e3)]
    net = doe_feeder(p1=200.0, p2=200.0)

    flat = solve_operating_envelope(net, cps;
        fairness=FairnessPolicy(kind=:equal))
    proportional_capacity = solve_operating_envelope(net, cps;
        fairness=FairnessPolicy(kind=:equal, normalization=:capacity))
    maxmin = solve_operating_envelope(net, cps;
        fairness=FairnessPolicy(kind=:max_min, normalization=:capacity))
    requested = solve_operating_envelope(net, cps;
        fairness=FairnessPolicy(kind=:equal, normalization=:request))
    custom = solve_operating_envelope(net, cps;
        fairness=FairnessPolicy(kind=:equal, normalization=:custom))

    @test flat.envelope["large"][1] ≈ flat.envelope["small"][1] rtol=1e-4
    for r in (proportional_capacity, maxmin, requested, custom)
        @test r.envelope["large"][1] / r.envelope["small"][1] ≈ 2.0 rtol=2e-3
    end
    @test proportional_capacity.total_capacity[1] > flat.total_capacity[1] + 100.0

    weighted = solve_operating_envelope(net, cps; fairness=FairnessPolicy(
        kind=:proportional, weights=Dict("small"=>2.0)))
    @test weighted.envelope["small"][1] > weighted.envelope["large"][1]

    alpha0 = solve_operating_envelope(net, cps;
        fairness=FairnessPolicy(kind=:alpha, alpha=0.0))
    efficient = solve_operating_envelope(net, cps; fairness=:sum)
    @test alpha0.total_capacity[1] ≈ efficient.total_capacity[1] rtol=1e-4
end

@testset "Operating envelope: prescribed IBR Q-V law is retained" begin
    cp = ConnectionPoint(id="pv", bus="b1", ibr_id="pv1", export_max=10e3)
    unity = solve_operating_envelope(doe_ibr_feeder(volt_var=false), [cp])
    qv = solve_operating_envelope(doe_ibr_feeder(volt_var=true), [cp];
                                  security=:corners)
    @test all(s in ("LOCALLY_SOLVED", "OPTIMAL") for s in
              vcat(unity.termination_status, qv.termination_status))
    @test abs(unity.snapshots[1]["ibr"]["pv1"]["1"]["qg"]) < 1.0
    @test qv.snapshots[1]["ibr"]["pv1"]["1"]["qg"] < -1000.0
    @test qv.envelope["pv"][1] > unity.envelope["pv"][1] + 100.0
    @test qv.envelope["pv"][1] <= 10e3
    @test qv.diagnostics[1]["dispatch_points_per_scenario"] == 2

    @test_throws ArgumentError solve_operating_envelope(doe_ibr_feeder(),
        [ConnectionPoint(id="pv", bus="b1", ibr_id="missing", export_max=10e3)])
end

@testset "Operating envelope: optional STATCOM expands active-power DOE" begin
    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    base = doe_feeder_rx()
    with_statcom = deepcopy(base)
    add_statcom!(with_statcom, "bus2"; s_max=5000.0)
    stat = with_statcom["ibr"]["statcom_bus2"]
    stat["p_min"] = [0.0]; stat["p_max"] = [0.0]
    stat["q_min"] = [-5000.0]; stat["q_max"] = [5000.0]

    r0 = solve_operating_envelope(base, cps)
    rs = solve_operating_envelope(with_statcom, cps)
    @test rs.total_capacity[1] > r0.total_capacity[1] + 1000.0
    @test rs.snapshots[1]["ibr"]["statcom_bus2"]["1"]["qg"] < -1000.0
end

@testset "Operating envelope: inherited thermal and unbalance constraints" begin
    thermal = solve_operating_envelope(doe_thermal_feeder(),
        [ConnectionPoint(id="d", bus="b1", export_max=20e3)])
    @test 4000.0 < thermal.envelope["d"][1] < 5000.0
    @test thermal.snapshots[1]["line"]["l1"]["1"]["cm_fr"] ≈ 20.0 atol=1e-3
    @test "line:l1:1:i_max" in thermal.diagnostics[1]["binding_constraints"]

    cp = ConnectionPoint(id="single_phase", bus="b1", phase_terminals=["1"],
                         neutral="n", export_max=20e3)
    loose = solve_operating_envelope(doe_unbalanced_feeder(vneg_max=20.0), [cp])
    tight = solve_operating_envelope(doe_unbalanced_feeder(vneg_max=1.0), [cp])
    @test tight.envelope["single_phase"][1] < loose.envelope["single_phase"][1] - 1000.0
    @test tight.envelope["single_phase"][1] > 0.0
    @test "bus:b1:vneg_max" in tight.diagnostics[1]["binding_constraints"]
end

@testset "Operating envelope: operational publication and verification" begin
    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    nets = [doe_feeder(p1=200.0, p2=200.0), doe_feeder(p1=300.0, p2=300.0)]
    issued = DateTime(2026, 7, 16, 9, 0)
    r = solve_operating_envelope(nets, cps;
        temporal_fairness=:cumulative_max_min, temporal_dt_h=0.25,
        issued_at=issued, interval_seconds=300.0, validity_seconds=600.0)
    @test r.fairness_metrics[1]["available"]
    @test 0.0 <= r.fairness_metrics[1]["jain_index"] <= 1.0
    @test haskey(r.fairness_metrics[2], "cumulative_normalized")
    @test r.schedule[2]["valid_from"] == issued + Minute(5)
    @test r.schedule[2]["valid_until"] == issued + Minute(15)

    compared = compare_operating_envelope_policies(nets[1], cps,
        ["equal"=>:equal, "sum"=>:sum])
    @test Set(keys(compared)) == Set(["equal", "sum"])
    verified = verify_operating_envelope(nets[1], cps, compared["equal"];
        utilizations=:bound_point)
    @test verified.feasible == [true]
    @test verified.diagnostics[1]["verification"]

    fallback = solve_operating_envelope([nets[1], doe_feeder(p1=5e3, p2=5e3)], cps;
        security=:corners, fallback=:last_feasible)
    @test fallback.schedule[2]["publication_source"] in (:optimized, :last_feasible_fallback)
end
