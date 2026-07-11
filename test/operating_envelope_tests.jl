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
