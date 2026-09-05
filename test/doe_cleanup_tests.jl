@testset "DOE scrambled Halton and bounded dropout generation" begin
    inverse = PowerOptLab._doe_radical_inverse
    # Complement every binary digit, including the infinite trailing zeros:
    # phi_complement(n) = 1 - phi(n). The truncated implementation collided.
    expected = [1//2, 3//4, 1//4, 7//8, 3//8, 5//8, 1//8, 15//16]
    @test [inverse(n, 2, [1, 0]) for n in 1:8] ≈ expected
    @test length(unique(inverse(n, 2, [1, 0]) for n in 1:1024)) == 1024
    @test all(inverse(n, 3, [2, 1, 0]) ≈ 1 - inverse(n, 3) for n in 1:100)
    @test all(inverse(n, 5, collect(0:4)) == inverse(n, 5) for n in 1:100)
    diagnostics = PowerOptLab._doe_search_point_diagnostics
    @test !diagnostics(20, 2; dropout_depth=0,
        halton_scramble_seed=7)["halton_seed_stratified"]
    @test !diagnostics(2, 0; dropout_depth=0,
        halton_scramble_seed=7)["halton_seed_stratified"]
    @test diagnostics(1, 2; dropout_depth=0, sequence_offset=2,
        halton_scramble_seed=nothing)["halton_occupied_first_digit_strata"] == [2]
    # Seed 8 permutes the base-2 digits, which maps both samples into one
    # first-digit stratum: scrambling redistributes coverage, it does not
    # create it. (Seed 7 leaves base 2 unpermuted under the package stream.)
    @test PowerOptLab._doe_stable_permutation!(
        PowerOptLab._doe_stable_stream(8), 0:1) == [1, 0]
    @test diagnostics(1, 2; dropout_depth=0, sequence_offset=0,
        halton_scramble_seed=8)["halton_occupied_first_digit_strata"] == [1]
    @test PowerOptLab._doe_dropout_count(64, 64) == big(2)^64 - 1
    @test_throws ArgumentError PowerOptLab._doe_search_points(64, 0, 0;
        include_zero=false, include_bound=false, include_corners=false,
        max_exact_corners=10, dropout_depth=64, max_dropout_points=10)
end

@testset "DOE weighted rolling fairness and objective selection" begin
    cps = [ConnectionPoint(id="a", bus="a", export_max=10.0),
           ConnectionPoint(id="b", bus="b", export_max=10.0)]
    weighted = FairnessPolicy(kind=:max_min, weights=Dict("a"=>2.0, "b"=>1.0))
    function solve_service(history, dt)
        model = JuMP.Model(Ipopt.Optimizer)
        JuMP.set_silent(model)
        cap = Dict(id => JuMP.@variable(model, lower_bound=0.0, upper_bound=10.0)
                   for id in ("a", "b"))
        JuMP.@constraint(model, cap["a"] + cap["b"] <= 3.0)
        stage = PowerOptLab._set_fairness_objective!(
            model, cap, cps, weighted, :export, 1.0;
            temporal_history=history, temporal_dt_h=dt)
        PowerOptLab._optimize_fairness!(model, stage; max_min_tolerance=1e-8)
        allocation = Dict(id => JuMP.value(x) for (id, x) in cap)
        key = PowerOptLab._doe_fairness_selection_key(model, stage;
            temporal_history=history, temporal_dt_h=dt)
        return allocation, key
    end
    # A shared budget of 3 and entitlements 2:1 have the analytic allocation
    # (2,1). A balanced accumulated history must preserve this at every step.
    history = Dict("a"=>2.0, "b"=>1.0)
    for step in 1:3
        allocation, key = solve_service(history, 0.5)
        @test allocation["a"] ≈ 2.0 atol=1e-6
        @test allocation["b"] ≈ 1.0 atol=1e-6
        @test key[1] ≈ 1 + step / 2 atol=1e-6
        for id in keys(history)
            history[id] += 0.5 * allocation[id]
        end
    end
    whole, _ = solve_service(Dict("a"=>2.0, "b"=>1.0), 1.5)
    @test history["a"] ≈ 2 + 1.5 * whole["a"] atol=1e-6
    @test history["b"] ≈ 1 + 1.5 * whole["b"] atol=1e-6

    function fixed_key(values_, policy)
        model = JuMP.Model(Ipopt.Optimizer)
        JuMP.set_silent(model)
        cap = Dict(cp.id => JuMP.@variable(model) for cp in cps)
        for (cp, value) in zip(cps, values_)
            JuMP.fix(cap[cp.id], value)
        end
        stage = PowerOptLab._set_fairness_objective!(
            model, cap, cps, policy, :export, 1.0)
        PowerOptLab._optimize_fairness!(model, stage; max_min_tolerance=1e-8)
        return PowerOptLab._doe_fairness_selection_key(model, stage)
    end
    better(candidate, incumbent, tolerances) =
        PowerOptLab._doe_select_fairness_run([incumbent, candidate], [1, 2], tolerances) == 2
    @test PowerOptLab._doe_select_fairness_run(
        [[[1.0, 0.0]], [[0.94, 1.0]], [[0.88, 2.0]]], [1, 2, 3], [[0.1, 0.0]]) == 2
    # Lower total capacity can have higher log utility (product 16 versus 9).
    proportional = FairnessPolicy(kind=:proportional)
    balanced = fixed_key([4.0, 4.0], proportional)
    unequal = fixed_key([9.0, 1.0], proportional)
    @test balanced[1] ≈ 2log(4 + proportional.epsilon)
    @test better([balanced], [unequal], [[0.0]])
    @test !better([unequal], [balanced], [[0.0]])
    @test fixed_key([4.0, 4.0], FairnessPolicy(kind=:equal))[1] ≈ 4 atol=1e-6
    for (kind, alpha) in ((:max_total, 0.0), (:alpha, 0.0), (:alpha, 2.0))
        policy_ = FairnessPolicy(kind=kind, alpha=alpha, weights=Dict("a"=>2.0))
        expected_ = alpha == 2.0 ? -2 / (9 + policy_.epsilon) - 1 / (1 + policy_.epsilon) : 19.0
        @test fixed_key([9.0, 1.0], policy_)[1] ≈ expected_ atol=1e-6
    end
    maxmin = FairnessPolicy(kind=:max_min)
    @test better([fixed_key([4.0, 4.0], maxmin)],
        [fixed_key([9.0, 1.0], maxmin)], [[1e-7, 0.0]])
    @test better([[4.0, 9.0]], [[4.0, 8.0]], [[1e-7, 0.0]])
    @test !better([[3.0, 20.0]], [[4.0, 8.0]], [[1e-7, 0.0]])
    curtailment = FairnessPolicy(kind=:equal_curtailment)
    @test fixed_key([8.0, 8.0], curtailment)[1] ≈ -2.0 atol=1e-6
    @test better([fixed_key([8.0, 8.0], curtailment)],
        [fixed_key([5.0, 5.0], curtailment)], [[0.0]])
    # Earlier interval objectives take priority over later service.
    @test better([[2.0], [1.0]], [[1.0], [100.0]], [[0.0], [0.0]])
    @test !better([[2.0]], [[2.0]], [[0.0]])

    net = doe_feeder(p1=200.0, p2=200.0)
    physical_cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
                    ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    policy = FairnessPolicy(kind=:proportional, normalization=:capacity)
    result = solve_operating_envelope_multistart(net, physical_cps;
        fairness=policy, start_scales=(1.0, 0.9))
    scores = [sum(log(run.envelope[cp.id][1] / cp.export_max + policy.epsilon)
                  for cp in physical_cps) for run in result.runs]
    @test scores[result.selected_index] ≈ maximum(scores) atol=1e-7
    @test all(only(only(key)) ≈ score for (key, score) in
        zip(result.diagnostics["fairness_selection_keys"], scores))
end

@testset "DOE pairing requires experimental identity" begin
    net = doe_feeder(p1=300.0, p2=300.0)
    cps = [ConnectionPoint(id="d1", bus="bus1", export_max=10e3),
           ConnectionPoint(id="d2", bus="bus2", export_max=10e3)]
    draws = [sample_doe_gaussian_uncertainty([300.0], fill(1.0, 1, 1);
        parameter_names=["demand"], count=3, seed=seed) for seed in (71, 72)]
    @test [s.id for s in draws[1].samples] == [s.id for s in draws[2].samples]
    @test draws[1].samples[1].parameters != draws[2].samples[1].parameters
    materialize = (_, s) -> doe_feeder(p1=s.parameters["demand"], p2=300.0)
    models = [materialize_doe_scenarios(net, samples, materialize;
        dataset_id="ensemble-$i", materializer_id="demand-fixture", role=:test)
        for (i, samples) in enumerate(draws)]
    independent = compare_doe_uncertainty_models(
        ["a"=>models[1], "b"=>models[2]], cps, Dict("d1"=>100.0, "d2"=>100.0);
        scales=(0.5, 1.0))
    @test all(row.evidence_design == :unpaired && row.matched_scenario_count == 0
              for row in independent.pairwise_rows)
    @test independent.diagnostics["pairing_basis"] == :declared_shared_namespace

    left, right = [first(independent.curves[id].coverages) for id in ("a", "b")]
    function scoped(coverage, namespaces)
        rows = NamedTuple[merge(row, (metadata=merge(row.metadata,
            Dict("pairing_namespace"=>namespace)),))
            for (row, namespace) in zip(coverage.scenario_rows, namespaces)]
        return DOECoverageResult(coverage.outcome, coverage.selected_roles,
            coverage.verification, rows, coverage.metrics, coverage.diagnostics)
    end
    pair(a, b; method=:auto, key=nothing) = PowerOptLab._doe_coverage_pairing(
        a, b; pairing=method, pair_key=key)
    @test pair(left, right).method == :none
    @test pair(left, right; method=:uncertainty_sample_id).evidence_design == :paired
    @test pair(left, right; method=:scenario_id).evidence_design == :paired
    @test pair(left, right; method=:none).evidence_design == :unpaired
    @test pair(scoped(left, fill("experiment-a", 3)),
        scoped(right, fill("experiment-b", 3))).evidence_design == :unpaired
    # Identical evaluations of the same draws are valid matching units.
    scoped_left = scoped(left, fill("shared-draws", 3))
    @test pair(scoped_left, scoped(left, fill("shared-draws", 3))).evidence_design == :paired
    partial = pair(scoped_left, scoped(left, ["shared-draws", "shared-draws", "another-study"]))
    @test partial.evidence_design == :partially_paired
    @test length(partial.keys) == 2
    @test pair(scoped_left, scoped(left, ["shared-draws", "shared-draws", ""])).method == :none
    duplicate = deepcopy(left)
    duplicate.scenario_rows[2].metadata["uncertainty_sample_id"] =
        duplicate.scenario_rows[1].metadata["uncertainty_sample_id"]
    @test_throws ArgumentError pair(left, duplicate; method=:uncertainty_sample_id)

    reused = [materialize_doe_scenarios(net, draws[1], materialize;
        dataset_id="reused-$i", materializer_id="demand-fixture", role=:test,
        scenario_metadata=Dict("pairing_namespace"=>"shared-draws")) for i in 1:2]
    paired = compare_doe_uncertainty_models(["a"=>reused[1], "b"=>reused[2]], cps,
        Dict("d1"=>100.0, "d2"=>100.0); scales=(0.5, 1.0))
    @test all(row.evidence_design == :paired && row.matched_scenario_count == 3
              for row in paired.pairwise_rows)
    @test all(row.paired_conservative_agreement_count == 3 for row in paired.pairwise_rows)
end

@testset "DOE original and binned Brier identities" begin
    observation(id, p, y; weight=1.0) = DOEProbabilityObservation(
        id=id, predicted_violation_probability=p, observed_violation=y, weight=weight)
    result = evaluate_doe_probability_calibration(
        [observation("a", 0.1, false), observation("b", 0.9, true)]; bins=1)
    m = result.metrics
    @test m["brier_score"] ≈ 0.01
    @test m["binned_forecast_brier_score"] ≈ 0.25
    @test m["within_bin_brier_remainder"] ≈ -0.24
    @test m["within_bin_forecast_variance"] ≈ 0.16
    @test m["within_bin_forecast_outcome_covariance"] ≈ 0.2
    @test result.diagnostics["brier_decomposition_target"] == :bin_mean_forecasts
    @test !result.diagnostics["conditional_outcome_independence_asserted"]

    # Unequal weights and forecasts within each bin exercise the full identity.
    probabilities = [0.1, 0.4, 0.6, 0.9]
    outcomes = [0.0, 1.0, 0.0, 1.0]
    weights = [1.0, 3.0, 2.0, 4.0]
    rows = [observation(string(i), probabilities[i], Bool(outcomes[i]); weight=weights[i])
            for i in 1:4]
    weighted = evaluate_doe_probability_calibration(rows; bins=2, independence_assumption=true)
    m = weighted.metrics
    # Bin means (0.325, 0.8) computed independently from the four observations.
    bin_forecasts = [0.325, 0.325, 0.8, 0.8]
    @test m["brier_score"] ≈ sum(weights .* (probabilities .- outcomes).^2) / 10
    @test m["binned_forecast_brier_score"] ≈ sum(weights .* (bin_forecasts .- outcomes).^2) / 10
    @test m["binned_brier_reconstruction"] ≈ m["binned_forecast_brier_score"]
    @test m["brier_score"] ≈ m["binned_brier_reconstruction"] + m["within_bin_brier_remainder"]
    @test m["within_bin_brier_remainder"] ≈ m["within_bin_forecast_variance"] -
        2m["within_bin_forecast_outcome_covariance"]
    @test weighted.diagnostics["outcome_independent_design_asserted"]
    @test !weighted.diagnostics["overall_and_bin_intervals_joint"]
    equal_forecasts = evaluate_doe_probability_calibration(
        [observation("a", 0.2, false), observation("b", 0.2, true)]; bins=3)
    @test equal_forecasts.metrics["within_bin_brier_remainder"] ≈ 0 atol=1e-14
    @test equal_forecasts.metrics["within_bin_forecast_variance"] ≈ 0 atol=1e-14
    for treatment in (:exclude, :pass, :violation)
        unresolved = evaluate_doe_probability_calibration(
            [rows; observation("unresolved", 0.3, nothing)]; bins=2, unresolved=treatment)
        m = unresolved.metrics
        @test m["brier_score"] ≈ m["binned_brier_reconstruction"] + m["within_bin_brier_remainder"]
    end
end

@testset "DOE seeded draws are Julia-version independent" begin
    stream = PowerOptLab._doe_stable_stream
    index! = PowerOptLab._doe_stable_index!
    permutation! = PowerOptLab._doe_stable_permutation!

    # `rand(rng, ::UnitRange)`, `randperm` and `shuffle!` map random bits to
    # indices through Julia implementation details that have changed between
    # releases, so a recorded seed alone did not reproduce a study on a
    # different Julia. The stream below is defined inside the package. These
    # pins are a compatibility contract: changing them changes every recorded
    # study identity that depends on a seeded draw.
    @test [index!(stream(47), 3) for _ in 1:1] == [1]
    @test [index!(stream(47), 3), index!(stream(47), 3)] == [1, 1]
    let s = stream(47)
        @test [index!(s, 3) for _ in 1:4] == [1, 3, 2, 3]
    end
    let s = stream(47)
        @test [index!(s, 5) for _ in 1:4] == [1, 2, 5, 1]
    end
    @test permutation!(stream(7), 0:1) == [0, 1]
    @test permutation!(stream(3), 4) == [3, 4, 1, 2]

    # A single-outcome draw must not consume a word, so that streams stay
    # aligned when a degenerate dimension appears.
    let s = stream(11)
        @test index!(s, 1) == 1
        @test index!(s, 1) == 1
        @test index!(s, 9) == index!(stream(11), 9)
    end
    @test_throws ArgumentError index!(stream(1), 0)

    # Rejection sampling must not bias the low indices of a non-power-of-two
    # range; a modulo-only draw fails this at these counts.
    let s = stream(1234), counts = zeros(Int, 3)
        for _ in 1:60_000
            counts[index!(s, 3)] += 1
        end
        @test all(count -> abs(count - 20_000) < 700, counts)
    end
end
