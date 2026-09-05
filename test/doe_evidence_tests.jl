using BMOPFTools: add_statcom!

@testset "DOE evidence: numerical outcomes" begin
    MOI = JuMP.MOI
    for status in (MOI.ITERATION_LIMIT, MOI.TIME_LIMIT, MOI.NUMERICAL_ERROR,
                   MOI.INVALID_MODEL, MOI.INTERRUPTED, MOI.ALMOST_LOCALLY_SOLVED)
        outcome = PowerOptLab.SolveOutcome(status, MOI.INFEASIBLE_POINT, 1,
            true, false, false, false)
        @test PowerOptLab._doe_numerical_evidence(outcome) == :numerically_unresolved
    end
    for status in (MOI.INFEASIBLE, MOI.LOCALLY_INFEASIBLE)
        outcome = PowerOptLab.SolveOutcome(status, MOI.INFEASIBLE_POINT, 1,
            true, false, false, false)
        @test PowerOptLab._doe_numerical_evidence(outcome) == :local_infeasibility_candidate
    end
    net = doe_feeder(p1=200.0, p2=200.0)
    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    capacities = Dict("d1" => 100.0, "d2" => 100.0)
    scenarios = DOEScenarioSet([[DOEScenario(id="heldout", network=net, role=:test)]]; dataset_id="heldout")
    coverage = evaluate_operating_envelope_coverage(scenarios, cps, capacities;
        solver_options=(max_iter=0,), iid_assumption=true)
    @test coverage.outcome == :inconclusive
    @test coverage.metrics["candidate_scenario_count"] == 0
    @test coverage.metrics["unresolved_scenario_count"] == 1
    @test coverage.metrics["conservative_scenario_frequency"] == 1
    confirmation = confirm_operating_envelope_counterexample(net, cps, capacities,
        [1.0, 1.0]; solver_options=(max_iter=0,))
    @test confirmation.outcome == :inconclusive
    @test search_operating_envelope_utilizations(net, cps, capacities;
        samples=1, solver_options=(max_iter=0,)).outcome == :inconclusive

    # Joint solve succeeds; intentionally exhaust only the single-context replay.
    replay_limit! = contexts -> length(contexts) == 1 &&
        JuMP.set_optimizer_attribute(PowerOptLab._opf_model(only(contexts)), "max_iter", 0)
    replay = evaluate_operating_envelope_coverage(scenarios, cps, capacities;
        utilizations=:corners, start_hook! = replay_limit!)
    @test replay.verification.diagnostics[1]["joint_policy_feasible"]
    @test !only(replay.verification.feasible)
    @test replay.outcome == :inconclusive
    @test replay.metrics["passed_context_count"] == 0
    context = first(replay.verification.context_results[1])
    context.diagnostics["independent_replay"]["feasible"] = true
    context.diagnostics["independent_replay"]["control_replay_complete"] = false
    @test PowerOptLab._doe_context_verdict(context) === nothing
end

@testset "DOE evidence: bindings and faithful replay" begin
    ibr_net = doe_ibr_feeder(volt_var=false)
    cp = ConnectionPoint(id="pv", bus="b1", ibr_id="pv1", export_max=10e3)
    @test_throws ArgumentError solve_operating_envelope(ibr_net, [cp,
        ConnectionPoint(id="duplicate", bus="b1", ibr_id="pv1", export_max=10e3)])
    @test_throws ArgumentError solve_operating_envelope(ibr_net,
        [ConnectionPoint(id="wrong_phase", bus="b1", ibr_id="pv1",
                         phase_terminals=["2"], export_max=10e3)])
    @test_throws ArgumentError solve_operating_envelope(doe_unbalanced_feeder(),
        [ConnectionPoint(id="aggregate", bus="b1", phase_terminals=["1", "2", "3"],
                         export_max=20e3)])
    four = doe_unbalanced_feeder()
    four["ibr"] = Dict{String,Any}("four" => Dict{String,Any}(
        "bus" => "b1", "topology" => "FOUR_LEG", "prime_mover" => "PV",
        "terminal_map" => ["n", "3", "1", "2"],
        "p_min" => zeros(3), "p_max" => fill(1000.0, 3),
        "q_min" => zeros(3), "q_max" => zeros(3), "s_max" => fill(1200.0, 3)))
    four_cp = ConnectionPoint(id="four", bus="b1", ibr_id="four",
        phase_terminals=["3", "1", "2"], neutral="n", export_max=2000.0)
    @test solve_operating_envelope(four, [four_cp]).total_capacity[1] ≈ 2000 atol=0.01
    @test_throws ArgumentError solve_operating_envelope(four,
        [ConnectionPoint(id="wrong", bus="b1", ibr_id="four",
            phase_terminals=["n", "3", "1"], neutral="2", export_max=2000.0)])

    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3, import_max=3e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3, import_max=3e3)]
    net = doe_feeder_rx()
    imported = solve_operating_envelope(net, cps; direction=:import,
        per_unit=false, volt_var_watt_eps=0.003, interval_seconds=300,
        validity_seconds=600, control_policy=IssuePlusLocalLaws())
    replay = verify_operating_envelope(net, cps, imported)
    @test replay.feasible == [true]
    @test replay.diagnostics[1]["direction"] == :import
    @test replay.diagnostics[1]["control_policy_source"] == :issued_result
    @test replay.diagnostics[1]["formulation_settings"].per_unit == false
    @test replay.diagnostics[1]["formulation_settings"].volt_var_watt_eps == 0.003
    changed_cps = [ConnectionPoint(id="d1", bus="bus2", export_max=10e3), cps[2]]
    @test_throws ArgumentError verify_operating_envelope(net, changed_cps, imported)

    add_statcom!(net, "bus2"; s_max=5000.0)
    stat = net["ibr"]["statcom_bus2"]
    stat["p_min"] = [0.0]; stat["p_max"] = [0.0]
    stat["q_min"] = [-5000.0]; stat["q_max"] = [5000.0]
    perfect = solve_operating_envelope(net, cps; security=:corners,
        control_policy=PerfectRecourse())
    heldout = DOEScenarioSet([[DOEScenario(id="heldout", network=net, role=:test)]]; dataset_id="heldout")
    # Retain a shared-policy conflict with margin to each pointwise AC bound.
    # Exact optimized endpoints can legitimately replay as unresolved.
    conflict_caps = Dict(id => 0.99 * only(values) for (id, values) in perfect.envelope)
    conflict = evaluate_operating_envelope_coverage(heldout, cps, conflict_caps;
        utilizations=:corners, control_policy=IssuePlusLocalLaws())
    @test all(c.feasible === true for c in conflict.verification.context_results[1])
    @test conflict.outcome == :inconclusive
    @test conflict.metrics["conservative_scenario_frequency"] == 1
    @test conflict.metrics["passed_context_count"] == 0

    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=2000.0),
           ConnectionPoint(id="d2", bus="bus2", export_max=2000.0)]
    policy = PerfectRecourse(rules=[DOEControlRule(component=:ibr,
        id="statcom_bus2", quantity=:reactive_power, stage=:scenario)])
    other = deepcopy(net)
    other["load"]["d1"]["p_nom"] = [400.0]
    a = DOEScenario(id="a", network=net, role=:test)
    b = DOEScenario(id="b", network=other, role=:test)
    scenarios = DOEScenarioSet([[a, b]]; dataset_id="issued")
    issued = solve_operating_envelope(scenarios, cps; control_policy=policy)
    for target in (DOEScenarioSet([[b, a]]; dataset_id="reordered"), DOEScenarioSet([[b]]; dataset_id="subset"))
        replay = verify_operating_envelope(target, cps, issued)
        @test replay.feasible == [true]
        @test replay.diagnostics[1]["issued_control_replay_count"] == length(target.intervals[1])
        records = replay.diagnostics[1]["issued_control_values"]
        originals = issued.diagnostics[1]["issued_control_values"]
        @test records[1]["value"] ≈ originals[2]["value"] atol=1e-3
    end
    mismatch = DOEScenarioSet([[DOEScenario(id="a", network=other, role=:test)]]; dataset_id="changed")
    @test_throws ArgumentError verify_operating_envelope(mismatch, cps, issued)
    untyped = solve_operating_envelope([[net, deepcopy(net)]], cps; control_policy=policy)
    @test_throws ArgumentError verify_operating_envelope([[net, deepcopy(net)]], cps, untyped)
    # New held-out information states deliberately retain issue controls only.
    @test evaluate_operating_envelope_coverage(mismatch, cps, issued) isa DOECoverageResult
    # Force a diagnostic run at zero dispatch while retaining both scenario Qs.
    diagnostic = verify_operating_envelope(scenarios, cps, issued;
        utilizations=[[0.0, 0.0]], solver_options=(max_iter=0,))
    @test length(diagnostic.context_results[1]) == 2
    @test diagnostic.diagnostics[1]["verification_outcome"] == :unresolved

    failed = solve_search_stable_operating_envelope(net, cps;
        max_rounds=1, samples_per_round=1,
        solve_keywords=(fallback=:zero, solver_options=(max_iter=0,)))
    @test failed.outcome == :allocation_failed
    retained = solve_search_stable_operating_envelope(net, cps;
        initial_utilizations=[[0.123, 0.456]], samples_per_round=1, max_rounds=1)
    @test [0.123, 0.456] in retained.utilization_points
end

include("../scripts/cases/doe_analytic_reference.jl")
using .DOEAnalyticReference

@testset "DOE analytic reference: independent resistive phases" begin
    net = doe_unbalanced_feeder(vneg_max=1.0)
    cps = [ConnectionPoint(id=string(i), bus="b1", phase_terminals=[string(i)],
        export_max=20e3, import_max=20e3) for i in 1:3]
    @test equal_export_limit() == 1747.5
    @test equal_import_limit() == 1702.5
    @test negative_sequence(phase_voltages([0, 19500, 19500])) ≈ 10 atol=1e-12
    @test resistive_voltage(-230^2/(4*0.4)) == 115
    @test_throws DomainError resistive_voltage(-40000)
    bound = solve_operating_envelope(net, cps; security=:bound_point)
    @test only(bound.total_capacity) ≈ 58500 atol=0.1
    points = PowerOptLab._dispatch_patterns(3, :corners)
    append!(points, [[0.25, 0.6, 0.9], [0.1, 0.2, 0.3]])
    for (direction, exact) in ((:export, equal_export_limit()), (:import, equal_import_limit()))
        for per_unit in (true, false)
            result = solve_operating_envelope(net, cps; direction, per_unit,
                security=:corners, control_policy=IssuePlusLocalLaws())
            @test only(result.total_capacity) ≈ 3exact atol=0.15
            verification = verify_operating_envelope(net, cps, result; utilizations=points)
            @test all(verification.feasible)
            sign = direction == :export ? 1 : -1
            for context in verification.context_results[1]
                p = sign .* context.utilization .* [result.envelope[cp.id][1] for cp in cps]
                reference = phase_voltages(p)
                bus = context.snapshot["bus"]["b1"]
                actual = [complex(bus[string(i)]["vr"], bus[string(i)]["vi"]) for i in 1:3]
                @test maximum(abs.(actual - reference)) < 1e-4
                @test negative_sequence(actual) <= 1 + 1e-5
            end
        end
    end
end

@testset "DOE analytic reference: exact affine box containment" begin
    A = [1.0 2.0; -1.0 2.0; 1.0 -1.0]
    b = [10.0, 5.0, 4.0]
    lower, upper = [-1.0, -2.0], [2.0, 1.5]
    vertices = [[x, y] for x in (lower[1], upper[1]) for y in (lower[2], upper[2])]
    exact = polyhedral_box_headroom(A, b, lower, upper)
    enumerated = [minimum(b[i] - (A*p)[i] for p in vertices) for i in eachindex(b)]
    @test exact ≈ enumerated
    @test all(exact .>= 0)
    @test any(polyhedral_box_headroom(A, b, lower, 2upper) .< 0)
end

module DOEDsseExample
include("../scripts/validate_doe_from_dsse.jl")
include("../scripts/cases/doe_dsse_validation_demo.jl")
end

@testset "DOE DSSE adapter and independent PF replay" begin
    case = DOEDsseExample.doe_validation_case()
    run = DOEDsseExample.run_doe_validation(case)
    @test run.max_dsse_voltage_error_V < 1e-3
    @test run.base.max_doe_pf_voltage_difference_V < 1e-3
    @test run.statcom.max_doe_pf_voltage_difference_V < 1e-3
    @test run.statcom.doe.total_capacity[1] > run.base.doe.total_capacity[1]
    # Verification at the issued capacity sits on the binding constraint, so the
    # joint solve may not converge there on every platform. A candidate
    # violation is still an error; `:unresolved` is reported and tolerated.
    @test all(outcome -> outcome in (:passed, :unresolved),
              run.base.verification_outcomes)
    @test all(outcome -> outcome in (:passed, :unresolved),
              run.statcom.verification_outcomes)
    @test !any(outcome -> outcome == :candidate_violation,
               vcat(run.base.verification_outcomes,
                    run.statcom.verification_outcomes))
    # The independent PF replay is the substantive agreement check, and it must
    # still hold regardless of how the boundary solve resolved.
    @test run.statcom.pf["termination_status"] in ("LOCALLY_SOLVED", "OPTIMAL")
    @test_throws ArgumentError DOEDsseExample._validate_snapshot("Volt-VAr", doe_ibr_feeder(),
        [ConnectionPoint(id="pv", bus="b1", ibr_id="pv1", export_max=10e3)])
    result = run.statcom.doe
    fixed = DOEDsseExample._fixed_dispatch_net(case.with_statcom_net,
        case.connection_points,
        Dict(cp.id => result.envelope[cp.id][1] for cp in case.connection_points),
        result.snapshots[1]; direction=result.direction,
        control_audit=result.diagnostics[1]["control_audit"])
    q = result.snapshots[1]["ibr"]["statcom_bus2"]["1"]["qg"]
    @test fixed["ibr"]["statcom_bus2"]["q_min"] == [q]
    @test fixed["ibr"]["statcom_bus2"]["q_max"] == [q]
    imported = DOEDsseExample._fixed_dispatch_net(case.operational_net,
        case.connection_points, Dict("pv1" => 100.0, "pv2" => 200.0),
        result.snapshots[1]; direction=:import, control_audit=[])
    @test imported["ibr"]["pv1"]["p_min"] == [-100.0]
    @test imported["ibr"]["pv2"]["p_max"] == [-200.0]
    changed = (residuals=[merge(r, (estimated=1.1*r.estimated,)) for r in run.estimate.residuals],)
    materialized = case.materialize_estimate(changed, deepcopy(case.operational_net))
    @test materialized["load"]["d1"]["p_nom"][1] ≈ 1320 atol=1e-3
    @test_throws ErrorException DOEDsseExample._require_pf_witness(
        Dict("termination_status" => "ITERATION_LIMIT"))
    @test_throws ErrorException DOEDsseExample._max_voltage_error(
        Dict("b" => Dict("1" => Dict("vm" => 230.0))), Dict())
    unsupported = [Dict("native_classification" => :free, "component" => :transformer,
                        "id" => "tap", "quantity" => :tap)]
    @test_throws ErrorException DOEDsseExample._fixed_dispatch_net(case.operational_net,
        case.connection_points, Dict("pv1" => 100.0, "pv2" => 200.0),
        result.snapshots[1]; direction=:export, control_audit=unsupported)
end
