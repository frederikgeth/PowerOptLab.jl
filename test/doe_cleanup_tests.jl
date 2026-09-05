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
    @test diagnostics(1, 2; dropout_depth=0, sequence_offset=0,
        halton_scramble_seed=7)["halton_occupied_first_digit_strata"] == [1]
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
