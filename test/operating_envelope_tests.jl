using BMOPFTools: OpfModelKey, add_statcom!, opf_model, register_opf_object!
using Dates
using JuMP

@testset "Operating envelope: equal allocation is uniform, voltage-limited, load-dependent" begin
    cps = [ConnectionPoint(id="der1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="der2", bus="bus2", export_max=10e3)]
    # Two intervals: low baseline load, then high. Higher local load absorbs
    # export and frees voltage headroom, so the envelope should grow.
    nets = [doe_feeder(p1=200.0, p2=200.0), doe_feeder(p1=5000.0, p2=5000.0)]

    res = solve_operating_envelope(nets, cps; fairness=:equal)
    @test all(s in ("LOCALLY_SOLVED", "OPTIMAL") for s in res.termination_status)
    @test solve_status(res).publishable
    @test solve_status(res).optimal
    @test solve_diagnostics(res).published_count == 2
    @test all(diag["primal_residual_available"] for diag in res.diagnostics)
    @test all(diag["maximum_primal_constraint_violation"] < 1e-5
              for diag in res.diagnostics)

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
    @test !solve_status(infeasible).publishable
    @test !solve_status(infeasible).feasible
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
    @test !r.diagnostics[1]["global_certificate"]
    @test r.diagnostics[1]["solver_class"] == :local_nonlinear
    @test r.diagnostics[1]["uncertainty_semantics"] == :finite_scenario_set
    @test r.diagnostics[1]["control_recourse"] == :perfect_recourse
    @test r.diagnostics[1]["control_policy"] == :perfect_recourse
    @test r.diagnostics[1]["control_policy_source"] == :legacy_default
    @test r.diagnostics[1]["prescribed_ibr_controls"] == :retained
    @test !r.diagnostics[1]["control_nonanticipativity"]
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

@testset "Operating envelope: explicit control-recourse policies" begin
    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    net = doe_feeder_rx()
    add_statcom!(net, "bus2"; s_max=5000.0)
    stat = net["ibr"]["statcom_bus2"]
    stat["p_min"] = [0.0]; stat["p_max"] = [0.0]
    stat["q_min"] = [-5000.0]; stat["q_max"] = [5000.0]

    perfect = solve_operating_envelope(net, cps; security=:corners,
        control_policy=PerfectRecourse())
    issued = solve_operating_envelope(net, cps; security=:corners,
        control_policy=IssuePlusLocalLaws())
    @test all(s in ("LOCALLY_SOLVED", "OPTIMAL") for s in
              vcat(perfect.termination_status, issued.termination_status))
    @test perfect.total_capacity[1] > issued.total_capacity[1] + 100.0
    @test perfect.diagnostics[1]["control_policy"] == :perfect_recourse
    @test perfect.diagnostics[1]["control_policy_source"] == :explicit
    @test perfect.diagnostics[1]["perfect_recourse_controls_present"]
    @test !perfect.diagnostics[1]["control_nonanticipativity"]
    @test issued.diagnostics[1]["control_policy"] == :issue_plus_local_laws
    @test length(issued.diagnostics[1]["control_policy_signature"]) == 64
    @test issued.diagnostics[1]["control_nonanticipativity"]
    @test issued.diagnostics[1]["control_link_constraints"] == 3
    @test length(issued.diagnostics[1]["issued_control_values"]) == 1
    q_audit = only(filter(item -> item["id"] == "statcom_bus2" &&
                                  item["quantity"] == :reactive_power,
                          issued.diagnostics[1]["control_audit"]))
    @test q_audit["stage"] == :issue
    @test q_audit["contexts_present"] == 4
    @test q_audit["link_constraints"] == 3

    operational_replay = verify_operating_envelope(net, cps, perfect;
        utilizations=:corners, control_policy=IssuePlusLocalLaws())
    @test operational_replay.feasible == [false]
    @test operational_replay.diagnostics[1]["control_policy"] ==
          :issue_plus_local_laws
    @test length(operational_replay.context_results[1]) == 4
    @test all(context.feasible === true for context in
              operational_replay.context_results[1])
    @test operational_replay.diagnostics[1]["infeasibility_interpretation"] ==
          :shared_control_incompatibility_or_joint_nlp_failure

    issued_replay = verify_operating_envelope(net, cps, issued;
        utilizations=:corners, control_policy=IssuePlusLocalLaws())
    @test issued_replay.feasible == [true]
    @test issued_replay.diagnostics[1]["issued_control_replay_source"] ==
          :operating_envelope_result
    @test issued_replay.diagnostics[1]["issued_control_replay_count"] == 1

    scenario_policy = PerfectRecourse(rules=[DOEControlRule(
        component=:ibr, id="statcom_bus2", quantity=:reactive_power,
        stage=:scenario)])
    scenario = solve_operating_envelope([[net, deepcopy(net)]], cps;
        security=:corners, control_policy=scenario_policy)
    @test scenario.diagnostics[1]["control_link_constraints"] == 6
    @test length(scenario.diagnostics[1]["issued_control_values"]) == 2
    @test_throws ArgumentError DOEControlRule(
        component=:ibr, id="x", quantity=:reactive_power, stage=:unknown)
    @test_throws ArgumentError DOEControlPolicy(
        rules=[DOEControlRule(component=:ibr, id="x",
                              quantity=:reactive_power, stage=:issue),
               DOEControlRule(component=:ibr, id="x",
                              quantity=:reactive_power, stage=:context)])
    @test_throws ArgumentError solve_operating_envelope(net, cps;
        control_policy=PerfectRecourse(rules=[DOEControlRule(
            component=:ibr, id="missing", quantity=:reactive_power,
            stage=:context)]))

    qv = solve_operating_envelope(doe_ibr_feeder(volt_var=true),
        [ConnectionPoint(id="pv", bus="b1", ibr_id="pv1", export_max=10e3)];
        security=:corners, control_policy=IssuePlusLocalLaws())
    qv_q = only(filter(item -> item["id"] == "pv1" &&
                               item["quantity"] == :reactive_power,
                       qv.diagnostics[1]["control_audit"]))
    @test qv_q["stage"] == :local_law
    @test qv_q["automatic_laws"] == [:volt_var]

    tap_net = bilevel_demo_network()
    tap = tap_net["transformer"]["single_phase"]["reg"]
    tap["tap_min"] = 0.95; tap["tap_max"] = 1.05
    tap_cps = [ConnectionPoint(id="p1", bus="lv1", ibr_id="pv1",
                               export_max=6000.0),
               ConnectionPoint(id="p2", bus="lv2", ibr_id="pv2",
                               export_max=6000.0)]
    tap_result = solve_operating_envelope(tap_net, tap_cps;
        security=:corners, control_policy=IssuePlusLocalLaws())
    tap_audit = only(filter(item -> item["component"] == :transformer,
                            tap_result.diagnostics[1]["control_audit"]))
    @test tap_audit["id"] == "reg"
    @test tap_audit["stage"] == :issue
    @test tap_audit["link_constraints"] == 3

    generator_net = doe_feeder(p1=200.0, p2=200.0)
    generator_net["generator"] = Dict("g1" => Dict(
        "bus" => "bus1", "terminal_map" => ["1", "n"],
        "configuration" => "SINGLE_PHASE",
        "p_min" => [0.0], "p_max" => [100.0],
        "q_min" => [0.0], "q_max" => [0.0]))
    generator_cp = [ConnectionPoint(id="d", bus="bus2", export_max=1000.0)]
    generator_ideal = solve_operating_envelope(generator_net, generator_cp;
        control_policy=PerfectRecourse())
    generator_audit = only(filter(item -> item["component"] == :generator,
                                  generator_ideal.diagnostics[1]["control_audit"]))
    @test generator_audit["stage"] == :context
    @test !generator_audit["linkage_supported"]
    @test_throws ArgumentError solve_operating_envelope(
        generator_net, generator_cp; control_policy=IssuePlusLocalLaws())

    extension_key = OpfModelKey(:variable, :research_support, "setting")
    extension_hook! = ctx -> begin
        setting = @variable(opf_model(ctx), lower_bound=-1.0,
                            upper_bound=1.0,
                            base_name="research_support_setting")
        register_opf_object!(ctx, extension_key, setting)
    end
    registration = DOEControlRegistration(
        component=:research_device, id="support1", quantity=:setting,
        handle=extension_key, canonical_unit=:per_unit_setting,
        metadata=Dict("controller" => "laboratory fixture"))
    extension_policy = IssuePlusLocalLaws(registrations=[registration])
    extension_result = solve_operating_envelope(
        doe_feeder(p1=200.0, p2=200.0), cps;
        security=:corners, control_policy=extension_policy,
        context_hook! = extension_hook!)
    extension_audit = only(filter(
        item -> item["component"] == :research_device,
        extension_result.diagnostics[1]["control_audit"]))
    @test extension_audit["stage"] == :issue
    @test extension_audit["link_constraints"] == 3
    @test extension_audit["canonical_unit"] == :per_unit_setting
    @test extension_audit["metadata"]["controller"] == "laboratory fixture"
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
    @test solve_status(verified).publishable
    @test solve_diagnostics(verified).verification
    @test length(verified.context_results[1]) == 1
    replay = verified.context_results[1][1].diagnostics["independent_replay"]
    @test replay["feasible"]
    @test replay["control_replay_complete"]
    @test replay["maximum_voltage_difference_V"] < 1e-3

    custom = verify_operating_envelope(nets[1], cps, compared["equal"];
        utilizations=[[0.25, 0.75], [0.75, 0.25]])
    @test custom.feasible == [true]
    @test custom.diagnostics[1]["security_scope"] == :explicit_utilization_points
    @test custom.diagnostics[1]["dispatch_points_per_scenario"] == 2
    @test [context.utilization for context in custom.context_results[1]] ==
          [[0.25, 0.75], [0.75, 0.25]]

    fallback = solve_operating_envelope([nets[1], doe_feeder(p1=5e3, p2=5e3)], cps;
        security=:corners, fallback=:last_feasible)
    @test fallback.schedule[2]["publication_source"] in (:optimized, :last_feasible_fallback)
    @test solve_status(fallback).publishable
end

@testset "Operating envelope: explicit utilization and search-stable screening" begin
    net = doe_feeder(p1=200.0, p2=200.0)
    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    explicit = solve_operating_envelope(net, cps;
        utilizations=[[0.0, 0.0], [0.5, 0.25], [1.0, 1.0]])
    @test explicit.diagnostics[1]["security_scope"] ==
          :explicit_utilization_points
    @test explicit.diagnostics[1]["dispatch_points_per_scenario"] == 3
    @test haskey(explicit.diagnostics[1], "minimum_normalized_margins")
    @test_throws ArgumentError solve_operating_envelope(net, cps;
        security=:corners, utilizations=[[1.0, 1.0]])

    conservative = Dict(id => 0.5 * explicit.envelope[id][1]
                        for id in keys(explicit.envelope))
    stable = search_operating_envelope_utilizations(
        net, cps, conservative; samples=3)
    @test stable.outcome == :search_stable
    @test stable.diagnostics["global_certificate"] == false
    @test stable.utilization_points[1] == [0.0, 0.0]
    @test stable.utilization_points[2] == [1.0, 1.0]
    @test stable.utilization_points[3] == [0.5, 1 / 3]
    @test solve_diagnostics(stable).tested_point_count == 5

    unsafe = search_operating_envelope_utilizations(
        net, cps, Dict("d1" => 10e3, "d2" => 10e3); samples=0)
    @test unsafe.outcome == :candidate_counterexample
    @test !isempty(unsafe.candidate_contexts)
    @test unsafe.diagnostics["claim"] ==
          :candidate_violation_requires_confirmation

    adversarial = search_operating_envelope_adversarial(
        net, cps, conservative; seed_samples=1,
        refinement_rounds=1, restarts=1, initial_step=0.2)
    @test adversarial.outcome == :search_stable
    @test length(adversarial.verifications) == 2
    @test length(adversarial.utilization_points) > 3
    @test all(isfinite, adversarial.point_scores)
    @test adversarial.worst_interval == 1
    @test adversarial.worst_context !== nothing
    @test adversarial.diagnostics["score_definition"] ==
          :negative_minimum_declared_normalized_constraint_margin
    @test !solve_diagnostics(adversarial).global_certificate

    adversarial_unsafe = search_operating_envelope_adversarial(
        net, cps, Dict("d1" => 10e3, "d2" => 10e3);
        seed_samples=0, refinement_rounds=2)
    @test adversarial_unsafe.outcome == :candidate_counterexample
    @test length(adversarial_unsafe.verifications) == 1
    @test !isempty(adversarial_unsafe.candidate_contexts)
    @test adversarial_unsafe.diagnostics["global_certificate"] == false
    @test_throws ArgumentError search_operating_envelope_adversarial(
        net, cps, conservative; initial_step=0.0)

    # A balanced three-phase upper point satisfies a tight negative-sequence
    # limit, while partial utilization breaks the balance. This is the key
    # range-safety failure that a bound-point allocation cannot reveal.
    unbalanced_net = doe_unbalanced_feeder(vneg_max=1.0)
    phase_cps = [ConnectionPoint(
        id=string(phase), bus="b1", phase_terminals=[string(phase)],
        neutral="n", export_max=20e3) for phase in 1:3]
    balanced_bound = solve_operating_envelope(
        unbalanced_net, phase_cps; security=:bound_point,
        control_policy=PerfectRecourse())
    @test all(isfinite, balanced_bound.total_capacity)
    @test verify_operating_envelope(
        unbalanced_net, phase_cps, balanced_bound;
        utilizations=:bound_point,
        control_policy=PerfectRecourse()).feasible == [true]
    range_failure = search_operating_envelope_adversarial(
        unbalanced_net, phase_cps, balanced_bound;
        seed_samples=0, refinement_rounds=1, restarts=2,
        initial_step=0.5, control_policy=PerfectRecourse())
    @test range_failure.outcome == :candidate_counterexample
    @test length(range_failure.verifications) == 2
    @test any(context -> any(x -> 0 < x < 1, context.utilization),
              range_failure.candidate_contexts)
    candidate = first(range_failure.candidate_contexts)
    confirmation = confirm_operating_envelope_counterexample(
        unbalanced_net, phase_cps, balanced_bound, candidate.utilization;
        start_scales=(1.0, 0.95), control_policy=PerfectRecourse())
    @test confirmation.outcome == :repeated_candidate
    @test confirmation.diagnostics["candidate_reproduced_runs"] ==
          [true, true]
    @test !solve_diagnostics(confirmation).global_certificate

    adaptive_range = solve_adversarial_search_stable_operating_envelope(
        unbalanced_net, phase_cps;
        max_rounds=2,
        control_policy=PerfectRecourse(),
        search_keywords=(seed_samples=0, refinement_rounds=1,
                         restarts=2, initial_step=0.5))
    @test adaptive_range.outcome == :search_stable
    @test adaptive_range.rounds == 2
    @test length(adaptive_range.allocations) == 2
    @test adaptive_range.allocations[2].total_capacity[1] <
          adaptive_range.allocations[1].total_capacity[1]
    @test adaptive_range.diagnostics["final_allocation_screened"]

    adaptive = solve_search_stable_operating_envelope(
        net, cps; samples_per_round=2, max_rounds=2)
    @test adaptive.outcome == :search_stable
    @test adaptive.rounds == 1
    @test length(adaptive.searches) == 1
    @test !adaptive.diagnostics["global_certificate"]

    multistart = solve_operating_envelope_multistart(
        net, cps; start_scales=(1.0, 0.98))
    @test length(multistart.runs) == 2
    @test multistart.selected_index in (1, 2)
    @test all(count > 0 for count in
              multistart.diagnostics["start_changed_variable_counts"])
    @test isfinite(multistart.diagnostics["maximum_capacity_spread_W"])
    @test !multistart.diagnostics["global_certificate"]
end

@testset "Operating envelope: reproducible study manifest" begin
    net = doe_feeder(p1=200.0, p2=300.0)
    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3)]
    first_spec = DOEStudySpec(net, cps;
        security=:corners,
        control_policy=IssuePlusLocalLaws(),
        fairness=FairnessPolicy(kind=:max_min, normalization=:capacity),
        seeds=Dict("scenario_generation" => 42),
        metadata=Dict("case" => "manifest fixture"))
    second_spec = DOEStudySpec(deepcopy(net), deepcopy(cps);
        security=:corners,
        control_policy=IssuePlusLocalLaws(),
        fairness=FairnessPolicy(kind=:max_min, normalization=:capacity),
        seeds=Dict("scenario_generation" => 42),
        metadata=Dict("case" => "manifest fixture"))
    @test first_spec.study_id == second_spec.study_id
    @test first_spec.network_hashes == second_spec.network_hashes
    @test doe_study_manifest(first_spec)["study_id"] == first_spec.study_id
    @test first_spec.software_versions["PowerOptLab"] == "0.1.0"

    changed = deepcopy(net)
    changed["load"]["d1"]["p_nom"][1] += 1.0
    changed_spec = DOEStudySpec(changed, cps;
        security=:corners, control_policy=IssuePlusLocalLaws(),
        fairness=FairnessPolicy(kind=:max_min, normalization=:capacity),
        seeds=Dict("scenario_generation" => 42),
        metadata=Dict("case" => "manifest fixture"))
    @test changed_spec.study_id != first_spec.study_id
    @test changed_spec.network_hashes != first_spec.network_hashes

    benchmark_spec = DOEStudySpec(net, cps)
    result = solve_operating_envelope(net, cps)
    rows = doe_benchmark_rows(benchmark_spec, result; method_label="baseline")
    @test length(rows) == 1
    @test rows[1].study_id == benchmark_spec.study_id
    @test rows[1].method == "baseline"
    @test length(rows[1].control_policy_signature) == 64
    @test rows[1].issued_control_count == 0
    @test !rows[1].global_certificate
    verification = verify_operating_envelope(net, cps, result)
    context_rows = doe_context_benchmark_rows(
        benchmark_spec, verification)
    @test length(context_rows) == 1
    @test context_rows[1].replay_feasible === true
    @test !isempty(context_rows[1].minimum_normalized_margins)
    @test context_rows[1].issued_control_replay_source == :no_issued_controls
end

@testset "Operating envelope: typed scenarios and held-out coverage" begin
    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    calibration = DOEScenario(
        id="calibration-high-load",
        network=doe_feeder(p1=5000.0, p2=5000.0),
        role=:calibration,
        weight=1.0,
        source="synthetic fixture",
        generation_method=:deterministic_fixture,
        seed=11,
        timestamp=DateTime(2026, 1, 1, 12))
    restrictive_test = DOEScenario(
        id="test-low-load",
        network=doe_feeder(p1=200.0, p2=200.0),
        role=:test,
        weight=0.7,
        source="synthetic fixture",
        generation_method=:held_out_fixture,
        seed=12)
    permissive_test = DOEScenario(
        id="test-high-load",
        network=doe_feeder(p1=5000.0, p2=5000.0),
        role=:test,
        weight=0.3,
        source="synthetic fixture",
        generation_method=:held_out_fixture,
        seed=13)
    scenarios = DOEScenarioSet(
        [calibration, restrictive_test, permissive_test];
        dataset_id="doe-held-out-fixture-v1",
        metadata=Dict("license" => "synthetic"))
    calibration_set = select_doe_scenarios(
        scenarios; roles=:calibration)
    test_set = select_doe_scenarios(scenarios; roles=:test)
    @test length(calibration_set.intervals[1]) == 1
    @test length(test_set.intervals[1]) == 2
    @test test_set.metadata["selected_roles"] == [:test]

    allocation = solve_operating_envelope(
        calibration_set, cps; control_policy=PerfectRecourse())
    @test allocation.diagnostics[1]["uncertainty_semantics"] ==
          :typed_finite_scenario_set
    @test allocation.diagnostics[1]["scenario_provenance"][1]["id"] ==
          "calibration-high-load"
    typed_spec = DOEStudySpec(
        calibration_set, cps; control_policy=PerfectRecourse())
    manifest = doe_study_manifest(typed_spec)
    @test manifest["scenario_set_metadata"]["dataset_id"] ==
          "doe-held-out-fixture-v1"
    @test manifest["scenario_provenance"][1][1]["role"] == :calibration

    curve = evaluate_operating_envelope_coverage_curve(
        scenarios, cps, allocation;
        scales=(1.0, 0.5, 1.0),
        roles=:test, iid_assumption=true, confidence=0.95,
        control_policy=PerfectRecourse())
    @test curve.scales == [0.5, 1.0]
    @test length(curve.coverages) == 2
    @test curve.diagnostics["issue_control_treatment"] ==
          :retained_from_issued_result
    @test curve.diagnostics["control_reoptimization"] == :none
    @test curve.diagnostics["continuous_threshold_estimated"] == false
    @test haskey(curve.diagnostics, "candidate_count_reversals")
    @test last(curve.rows).total_capacity_W ≈ allocation.total_capacity
    coverage = last(curve.coverages)
    @test coverage.outcome == :candidate_violations_observed
    @test coverage.metrics["scenario_count"] == 2
    @test coverage.metrics["candidate_scenario_count"] == 1
    @test coverage.metrics["passed_scenario_count"] == 1
    @test coverage.metrics["weighted_candidate_scenario_frequency"] ≈ 0.7
    @test 0.5 <= coverage.metrics["one_sided_hoeffding_upper_bound"] <= 1.0
    @test coverage.diagnostics["iid_assumption"]
    @test coverage.diagnostics["distribution_shift_assessed"] == false
    @test Set(row.scenario_id for row in coverage.scenario_rows) ==
          Set(["test-low-load", "test-high-load"])

    shift = compare_doe_coverage_shift(
        first(curve.coverages), last(curve.coverages);
        reference_label="half capacity", shifted_label="issued capacity")
    @test shift.outcome in (:higher_candidate_frequency,
                            :no_observed_change)
    @test shift.metric_deltas["conservative_scenario_frequency"] ==
          last(curve.rows).conservative_scenario_frequency -
          first(curve.rows).conservative_scenario_frequency
    @test shift.diagnostics["performance_shift_observed"] ==
          (shift.outcome == :higher_candidate_frequency)
    @test shift.diagnostics["distribution_shift_test"] == :not_performed
    @test shift.diagnostics["distribution_shift_detected"] == false

    historical = DOEScenarioSet([
        DOEScenario(id="history-1", network=doe_feeder(p1=1000.0, p2=1000.0),
            role=:unspecified, weight=1.0, timestamp=DateTime(2026, 1, 1),
            source="chronological fixture"),
        DOEScenario(id="history-2", network=doe_feeder(p1=1100.0, p2=1100.0),
            role=:train, weight=1.0, timestamp=DateTime(2026, 1, 2),
            source="chronological fixture"),
        DOEScenario(id="gap", network=doe_feeder(p1=1200.0, p2=1200.0),
            role=:validation, weight=1.0, timestamp=DateTime(2026, 1, 3),
            source="chronological fixture"),
        DOEScenario(id="future", network=doe_feeder(p1=1300.0, p2=1300.0),
            role=:stress, weight=1.0, timestamp=DateTime(2026, 1, 4),
            source="chronological fixture"),
    ]; dataset_id="chronological-fixture")
    time_split = split_doe_scenarios_by_time(
        historical;
        calibration_end=DateTime(2026, 1, 3),
        test_start=DateTime(2026, 1, 4),
        split_name="blocked-holdout")
    @test [scenario.id for scenario in time_split.calibration.intervals[1]] ==
          ["history-1", "history-2"]
    @test [scenario.role for scenario in time_split.calibration.intervals[1]] ==
          [:calibration, :calibration]
    @test only(time_split.test.intervals[1]).id == "future"
    @test only(time_split.test.intervals[1]).role == :test
    @test only(time_split.test.intervals[1]).metadata["original_role"] == :stress
    @test time_split.excluded_scenario_ids == ["gap"]
    @test time_split.diagnostics["temporal_overlap"] == false
    @test time_split.diagnostics["group_or_site_leakage_assessed"] == false

    grouped_history = DOEScenarioSet([
        DOEScenario(id="cal-site-a",
            network=doe_feeder(p1=1400.0, p2=1400.0),
            role=:unspecified, weight=1.0,
            timestamp=DateTime(2026, 2, 1), seed=21,
            source="grouped chronological fixture",
            generation_method=:deterministic_fixture,
            metadata=Dict("site_id" => "site-a", "load_kw" => 1.4)),
        DOEScenario(id="cal-site-b",
            network=doe_feeder(p1=1500.0, p2=1500.0),
            role=:unspecified, weight=1.0,
            timestamp=DateTime(2026, 2, 2), seed=22,
            source="grouped chronological fixture",
            generation_method=:deterministic_fixture,
            metadata=Dict("site_id" => "site-b", "load_kw" => 1.5)),
        DOEScenario(id="test-site-a",
            network=doe_feeder(p1=1600.0, p2=1600.0),
            role=:unspecified, weight=1.0,
            timestamp=DateTime(2026, 2, 4), seed=23,
            source="grouped chronological fixture",
            generation_method=:deterministic_fixture,
            metadata=Dict("site_id" => "site-a", "load_kw" => 1.6)),
        DOEScenario(id="test-site-c",
            network=doe_feeder(p1=1700.0, p2=1700.0),
            role=:unspecified, weight=1.0,
            timestamp=DateTime(2026, 2, 5), seed=24,
            source="grouped chronological fixture",
            generation_method=:deterministic_fixture,
            metadata=Dict("site_id" => "site-c", "load_kw" => 1.7)),
    ]; dataset_id="grouped-chronological-fixture")
    @test_throws ArgumentError split_doe_scenarios_by_time(
        grouped_history;
        calibration_end=DateTime(2026, 2, 3),
        group_key="site_id")
    allowed_group_split = split_doe_scenarios_by_time(
        grouped_history;
        calibration_end=DateTime(2026, 2, 3),
        group_key="site_id", group_overlap_policy=:allow)
    @test allowed_group_split.diagnostics["group_overlap_present"]
    @test allowed_group_split.diagnostics["group_overlap_retained"]
    @test allowed_group_split.diagnostics["overlapping_groups"] == ["site-a"]
    excluded_group_split = split_doe_scenarios_by_time(
        grouped_history;
        calibration_end=DateTime(2026, 2, 3),
        group_key="site_id", group_overlap_policy=:exclude_test)
    @test [scenario.id for scenario in
           excluded_group_split.test.intervals[1]] == ["test-site-c"]
    @test "test-site-a" in excluded_group_split.excluded_scenario_ids
    @test excluded_group_split.diagnostics[
        "group_overlap_excluded_count"] == 1
    @test excluded_group_split.diagnostics["group_overlap_detected"]
    @test !excluded_group_split.diagnostics["group_overlap_present"]

    clean_audit = audit_doe_scenario_calibration(
        excluded_group_split.calibration, excluded_group_split.test;
        group_key="site_id", require_chronological_order=true)
    @test clean_audit.outcome == :no_declared_leakage_detected
    @test clean_audit.leakage_checks["chronological_order"] === true
    @test isempty(clean_audit.leakage_checks["group_overlap"])
    @test clean_audit.calibration_summary[
        "effective_sample_size_by_interval"] == [2.0]
    @test clean_audit.evaluation_summary[
        "effective_sample_size_by_interval"] == [1.0]
    @test clean_audit.diagnostics["probabilistic_calibration_assessed"] == false

    overlap_audit = audit_doe_scenario_calibration(
        allowed_group_split.calibration, allowed_group_split.test;
        group_key="site_id", require_chronological_order=true)
    @test overlap_audit.outcome == :leakage_candidates_detected
    @test "metadata_group_overlap" in
          overlap_audit.leakage_checks["violated_requirements"]
    longitudinal_audit = audit_doe_scenario_calibration(
        allowed_group_split.calibration, allowed_group_split.test;
        group_key="site_id", require_group_disjoint=false,
        require_chronological_order=true)
    @test longitudinal_audit.outcome == :no_declared_leakage_detected
    incomplete_audit = audit_doe_scenario_calibration(
        time_split.calibration, time_split.test;
        group_key="site_id", require_chronological_order=true)
    @test incomplete_audit.outcome == :required_metadata_incomplete
    @test "group_disjointness_missing_metadata" in
          incomplete_audit.leakage_checks["unassessed_requirements"]
    reuse_audit = audit_doe_scenario_calibration(
        excluded_group_split.calibration,
        excluded_group_split.calibration)
    @test reuse_audit.outcome == :leakage_candidates_detected
    @test !isempty(reuse_audit.leakage_checks["scenario_id_overlap"])
    @test !isempty(reuse_audit.leakage_checks["exact_network_overlap"])

    descriptive_shift = test_doe_covariate_shift(
        allowed_group_split.calibration, allowed_group_split.test;
        features="load_kw")
    @test descriptive_shift.outcome == :descriptive_difference_only
    @test descriptive_shift.p_value === missing
    @test descriptive_shift.energy_distance >= 0
    @test only(descriptive_shift.feature_rows).mean_difference ≈ 0.2
    @test descriptive_shift.diagnostics["permutations"] == 0
    @test descriptive_shift.diagnostics[
        "general_distribution_shift_assessed"] == false
    @test_throws ArgumentError test_doe_covariate_shift(
        allowed_group_split.calibration, allowed_group_split.test;
        features="load_kw", exchangeability_assumption=true,
        permutation_unit=:group, group_key="site_id", permutations=19)

    feature_net = doe_feeder(p1=1000.0, p2=1000.0)
    reference_features = DOEScenarioSet([
        DOEScenario(
            id="feature-reference-$index", network=feature_net,
            role=:calibration, source="feature-shift fixture",
            metadata=Dict(
                "load_kw" => Float64(index),
                "temperature_c" => 20.0 + 0.1 * index,
                "site_id" => "reference-site-$index"))
        for index in 1:8]; dataset_id="feature-reference")
    shifted_features = DOEScenarioSet([
        DOEScenario(
            id="feature-shifted-$index", network=feature_net,
            role=:test, source="feature-shift fixture",
            metadata=Dict(
                "load_kw" => 20.0 + index,
                "temperature_c" => 30.0 + 0.1 * index,
                "site_id" => "shifted-site-$index"))
        for index in 1:8]; dataset_id="feature-shifted")
    tested_shift = test_doe_covariate_shift(
        reference_features, shifted_features;
        features=("load_kw", "temperature_c"),
        exchangeability_assumption=true,
        permutations=199, seed=31, alpha=0.05)
    @test tested_shift.outcome == :declared_covariate_shift_detected
    @test 0 < tested_shift.p_value <= 0.05
    @test tested_shift.features == ["load_kw", "temperature_c"]
    @test tested_shift.diagnostics["seed"] == 31
    repeated_shift = test_doe_covariate_shift(
        reference_features, shifted_features;
        features=("load_kw", "temperature_c"),
        exchangeability_assumption=true,
        permutations=199, seed=31, alpha=0.05)
    @test repeated_shift.p_value == tested_shift.p_value
    grouped_shift = test_doe_covariate_shift(
        reference_features, shifted_features;
        features="load_kw", exchangeability_assumption=true,
        permutation_unit=:group, group_key="site_id",
        permutations=99, seed=32, alpha=0.1)
    @test grouped_shift.outcome == :declared_covariate_shift_detected
    @test grouped_shift.diagnostics["reference_group_count"] == 8
    @test grouped_shift.diagnostics["shifted_group_count"] == 8
    @test_throws ArgumentError test_doe_covariate_shift(
        reference_features, shifted_features; features="missing")
    @test_throws ArgumentError test_doe_covariate_shift(
        reference_features, shifted_features;
        features="load_kw", permutations=0)

    @test_throws ArgumentError DOEScenario(
        id="bad", network=doe_feeder(p1=1.0, p2=1.0), role=:unknown)
    @test_throws ArgumentError DOEScenarioSet(
        [restrictive_test, DOEScenario(
            id="unweighted", network=doe_feeder(p1=1.0, p2=1.0), role=:test)];
        dataset_id="incomplete-weights")
    @test_throws ArgumentError select_doe_scenarios(
        scenarios; roles=:validation)
    @test_throws ArgumentError split_doe_scenarios_by_time(
        DOEScenarioSet([historical.intervals[1][1:2],
                        historical.intervals[1][3:4]];
                       dataset_id="multi-interval");
        calibration_end=DateTime(2026, 1, 3))
    @test_throws ArgumentError split_doe_scenarios_by_time(
        DOEScenarioSet([DOEScenario(
            id="missing-time", network=doe_feeder(p1=1.0, p2=1.0))];
            dataset_id="missing-time");
        calibration_end=DateTime(2026, 1, 3))
    @test_throws ArgumentError split_doe_scenarios_by_time(
        historical;
        calibration_end=DateTime(2026, 1, 4),
        test_start=DateTime(2026, 1, 3))
    @test_throws ArgumentError split_doe_scenarios_by_time(
        historical;
        calibration_end=DateTime(2026, 1, 3),
        group_overlap_policy=:unknown)
    @test_throws ArgumentError evaluate_operating_envelope_coverage_curve(
        scenarios, cps, allocation; scales=(-0.1,))
end
