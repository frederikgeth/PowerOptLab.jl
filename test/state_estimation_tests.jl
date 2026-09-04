using BMOPFTools: solve_pf, parse_bmopf

# A physics-only radial feeder: src—bus1—bus2, grounded neutrals, no loads/limits.
se_net() = parse_bmopf("""
{"bus":{
    "src": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus1":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
    "bus2":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[1000.0],"v_angle":[0.0]}},
 "linecode":{"lc":{"R_series_1_1":0.5}},
 "line":{
    "l1":{"bus_from":"src","bus_to":"bus1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
    "l2":{"bus_from":"bus1","bus_to":"bus2","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
"""; from_string=true)

# The same feeder with two nominal loads (an OPERATIONAL net — rejected by SE).
se_net_loaded() = two_bus_net(; load=true)

# Full, well-posed measurement set from a determined power flow of the loaded net.
function se_full_meas()
    pf = solve_pf(se_net_loaded(); per_unit=false)
    tv(b) = hypot(pf["bus"][b]["1"]["vr"], pf["bus"][b]["1"]["vi"])
    m = Measurement[]
    for b in ("bus1","bus2")
        push!(m, Measurement(kind=:vmag, bus=b, value=tv(b), sigma=2.0))
        push!(m, Measurement(kind=:pinj, bus=b, value=-20000.0, sigma=400.0))
        push!(m, Measurement(kind=:qinj, bus=b, value=0.0, sigma=400.0))
    end
    m, Dict(b => tv(b) for b in ("bus1","bus2"))
end

@testset "Shared measurement-model equations" begin
    @test measurement_prediction(:vr, 3.0, 4.0) == 3.0
    @test measurement_prediction(:vi, 3.0, 4.0) == 4.0
    @test measurement_prediction(:vmag, 3.0, 4.0) == 5.0
    @test measurement_prediction(:vmag, 3.0, 4.0; magnitude=7.0) == 7.0
    @test measurement_prediction(:pinj, 3.0, 4.0; ir=2.0, ii=-1.0) == 2.0
    @test measurement_prediction(:qinj, 3.0, 4.0; ir=2.0, ii=-1.0) == 11.0
    @test measurement_prediction(:pflow, 3.0, 4.0; ir=2.0, ii=-1.0) == 2.0
    @test measurement_prediction(:qflow, 3.0, 4.0; ir=2.0, ii=-1.0) == 11.0
    @test measurement_prediction(:ire, 0.0, 0.0; ir=3.0, ii=4.0) == 3.0
    @test measurement_prediction(:iim, 0.0, 0.0; ir=3.0, ii=4.0) == 4.0
    @test measurement_prediction(:imag, 0.0, 0.0; ir=3.0, ii=4.0) == 5.0
    @test_throws ArgumentError measurement_prediction(:pinj, 3.0, 4.0)
    @test_throws ArgumentError measurement_prediction(:bogus, 3.0, 4.0)
    @test_throws ArgumentError measurement_prediction(:vmag, 3.0, 4.0;
                                                       magnitude_epsilon=-1.0)
end

@testset "State estimation: exact recovery from noiseless measurements" begin
    meas, true_vm = se_full_meas()
    @test all(m -> m isa AbstractMeasurement, meas)
    @test measurement_kind(meas[1]) == meas[1].kind
    @test measurement_value(meas[1]) == meas[1].value
    @test measurement_sigma(meas[1]) == meas[1].sigma
    se = solve_state_estimation(se_net(), meas)
    @test se.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test se.primal_status == "FEASIBLE_POINT"
    for b in ("bus1","bus2")
        @test se.bus[b]["1"]["vm"] ≈ true_vm[b]  atol=1e-2
    end
    @test all(abs(r.residual) < 1e-1 for r in se.residuals if r.kind == :vmag)
    @test se.objective < 1e-3
    @test se.observability.observable === true
    @test se.observability.rank == se.observability.n_states == 4
    @test solve_status(se).publishable
    @test solve_diagnostics(se).residual_count == length(meas)
end

@testset "State estimation: noise reduction vs raw measurements" begin
    _, true_vm = se_full_meas()
    dv = Dict("bus1" => 1.5, "bus2" => -2.5)
    dp = Dict("bus1" => -300.0, "bus2" => 350.0)
    meas = Measurement[]
    for b in ("bus1","bus2")
        push!(meas, Measurement(kind=:vmag, bus=b, value=true_vm[b] + dv[b], sigma=2.0))
        push!(meas, Measurement(kind=:pinj, bus=b, value=-20000.0 + dp[b], sigma=400.0))
        push!(meas, Measurement(kind=:qinj, bus=b, value=0.0, sigma=400.0))
    end
    se = solve_state_estimation(se_net(), meas)
    @test se.primal_status == "FEASIBLE_POINT"
    raw_rms = sqrt(sum(dv[b]^2 for b in ("bus1","bus2")) / 2)
    est_rms = sqrt(sum((se.bus[b]["1"]["vm"] - true_vm[b])^2 for b in ("bus1","bus2")) / 2)
    @test est_rms < raw_rms
    for b in ("bus1","bus2")
        @test abs(se.bus[b]["1"]["vm"] - true_vm[b]) <= 3*2.0
    end
    for r in se.residuals
        @test r.residual ≈ r.measured - r.estimated  atol=1e-6
        @test r.standardized ≈ r.residual / (r.kind == :vmag ? 2.0 : 400.0)  rtol=1e-6
    end
end

@testset "State estimation: SI and per-unit agree" begin
    meas, true_vm = se_full_meas()
    se_pu = solve_state_estimation(se_net(), meas; per_unit=true)
    se_si = solve_state_estimation(se_net(), meas; per_unit=false)
    @test se_pu.primal_status == "FEASIBLE_POINT"
    @test se_si.primal_status == "FEASIBLE_POINT"
    for b in ("bus1","bus2")
        @test se_pu.bus[b]["1"]["vm"] ≈ true_vm[b]  atol=1e-2
        @test se_si.bus[b]["1"]["vm"] ≈ true_vm[b]  atol=1e-2   # was NaN before init fix
        @test se_si.bus[b]["1"]["vm"] ≈ se_pu.bus[b]["1"]["vm"] atol=1e-2
    end
    # s_base should not move the SI-referenced estimate.
    se_sb = solve_state_estimation(se_net(), meas; per_unit=true, s_base=1e3)
    for b in ("bus1","bus2")
        @test se_sb.bus[b]["1"]["vm"] ≈ se_pu.bus[b]["1"]["vm"] atol=1e-2
    end
end

@testset "State estimation: network contract is enforced" begin
    meas, _ = se_full_meas()
    # Loads/generators make it an operational model — rejected.
    @test_throws ArgumentError solve_state_estimation(se_net_loaded(), meas)
    # Operational voltage bound on a bus — rejected.
    bounded = se_net()
    bounded["bus"]["bus1"]["v_max"] = [950.0]
    @test_throws ArgumentError solve_state_estimation(bounded, meas)
    # Thermal limit on a line — rejected.
    thermal = se_net()
    thermal["line"]["l1"]["i_max"] = [100.0]
    @test_throws ArgumentError solve_state_estimation(thermal, meas)
    # allow_operational downgrades to a warning and still solves.
    se = solve_state_estimation(se_net_loaded(), meas; allow_operational=true)
    @test se.primal_status == "FEASIBLE_POINT"
end

@testset "State estimation: injection coverage contract" begin
    meas, true_vm = se_full_meas()
    # Drop bus2's injection pair and do NOT declare it — un-covered ⇒ error.
    m_gap = filter(m -> !(m.bus == "bus2" && m.kind in (:pinj, :qinj)), meas)
    @test_throws ArgumentError solve_state_estimation(se_net(), m_gap)

    # P without Q (unpaired) is ill-posed ⇒ error.
    m_ponly = [Measurement(kind=:vmag, bus="bus1", value=979.0, sigma=2.0),
               Measurement(kind=:pinj, bus="bus1", value=-20000.0, sigma=400.0),
               Measurement(kind=:vmag, bus="bus2", value=969.0, sigma=2.0)]
    @test_throws ArgumentError solve_state_estimation(se_net(), m_ponly; zero_injection=["bus2"])

    # A declared zero-injection bus that is also measured ⇒ error.
    @test_throws ArgumentError solve_state_estimation(se_net(), meas; zero_injection=["bus1"])
end

@testset "State estimation: zero-injection declaration recovers the state" begin
    # bus1 loaded, bus2 genuinely zero-injection; declare bus2, measure |V| there.
    net = parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
        "bus1":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
        "bus2":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
     "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[1000.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.5}},
     "line":{"l1":{"bus_from":"src","bus_to":"bus1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
             "l2":{"bus_from":"bus1","bus_to":"bus2","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}},
     "load":{"d1":{"bus":"bus1","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[20000.0],"q_nom":[0.0]}}}
    """; from_string=true)
    pf = solve_pf(net; per_unit=false)
    tv(b) = hypot(pf["bus"][b]["1"]["vr"], pf["bus"][b]["1"]["vi"])
    est_net = se_net()   # same topology, no load
    meas = [Measurement(kind=:vmag, bus="bus1", value=tv("bus1"), sigma=2.0),
            Measurement(kind=:pinj, bus="bus1", value=-20000.0, sigma=400.0),
            Measurement(kind=:qinj, bus="bus1", value=0.0, sigma=400.0),
            Measurement(kind=:vmag, bus="bus2", value=tv("bus2"), sigma=2.0)]
    se = solve_state_estimation(est_net, meas; zero_injection=["bus2"])
    @test se.primal_status == "FEASIBLE_POINT"
    @test se.observability.observable === true
    @test se.bus["bus2"]["1"]["vm"] ≈ tv("bus2") atol=1e-2
end

# A four-wire feeder whose neutral conductor is explicitly modelled. `ground`
# selects whether the SOURCE neutral is bonded to earth.
function se_net_fourwire(; ground::Bool)
    bond = ground ? ",\"perfectly_grounded_terminals\":[\"n\"]" : ""
    parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1","n"]$bond},
        "bus1":{"terminal_names":["1","n"]},
        "bus2":{"terminal_names":["1","n"]}},
     "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[1000.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.5,"R_series_2_2":0.5}},
     "line":{"l1":{"bus_from":"src","bus_to":"bus1","terminal_map_from":["1","n"],"terminal_map_to":["1","n"],"linecode":"lc","length":1.0},
             "l2":{"bus_from":"bus1","bus_to":"bus2","terminal_map_from":["1","n"],"terminal_map_to":["1","n"],"linecode":"lc","length":1.0}}}
    """; from_string=true)
end

@testset "State estimation: observability diagnostic" begin
    net = se_net()
    pf = solve_pf(se_net_loaded(); per_unit=false)
    estbus = pf["bus"]
    zi = Set{Tuple{String,String}}()

    # A well-posed set is observable, with redundancy from the surplus |V| rows.
    full, _ = se_full_meas()
    o2 = PowerOptLab._observability(net, full, "n", zi, estbus)
    @test o2.observable === true
    @test o2.redundancy == length(full) - o2.n_states

    # The diagnostic reproduces the equation set the SOLVER imposes: KCL holds at
    # every non-source node, and a measured (bus, terminal) absorbs its own KCL
    # row into a free injection. A node with no free injection is therefore a
    # zero-injection node to the diagnostic, so a |V|-only set describes a
    # determined passive network. (The public API separately refuses that set,
    # because absence of telemetry is not evidence of zero injection.)
    m_vonly = [Measurement(kind=:vmag, bus="bus1", value=1000.0, sigma=1.0),
               Measurement(kind=:vmag, bus="bus2", value=1000.0, sigma=1.0)]
    o_v = PowerOptLab._observability(net, m_vonly, "n", zi, estbus)
    @test o_v.observable === true
    @test_throws ArgumentError solve_state_estimation(net, m_vonly)
end

@testset "State estimation: exact-row scaling cannot change numerical rank" begin
    # Exact constraints and whitened measurements have different physical units.
    # Ranking their stacked matrix lets the larger block set the tolerance and
    # can erase an otherwise independent direction.
    Hm = [1.0 0.0]
    C = [0.0 1.0]
    for scale in (1e-30, 1.0, 1e30)
        d = PowerOptLab._se_observability_metrics(Hm, scale .* C)
        @test d.observable
        @test d.rank == d.n_states == 2
        @test d.tangent_dimension == 1
        @test d.redundancy == 0
        @test d.min_singular ≈ 1.0
        @test d.cond ≈ 1.0
    end
end

@testset "State estimation: observability covers non-grounded return terminals" begin
    # Regression: the diagnostic used to emit rows only for measurements and for
    # DECLARED zero-injection nodes. An explicitly modelled neutral is neither,
    # so its KCL rows went missing, the state was reported UNOBSERVABLE, and
    # `redundancy` went negative — a structural giveaway that equations, not
    # measurements, were absent.
    meas = Measurement[]
    for b in ("bus1", "bus2")
        push!(meas, Measurement(kind=:vmag, bus=b, value=980.0, sigma=2.0))
        push!(meas, Measurement(kind=:pinj, bus=b, value=-20_000.0, sigma=400.0))
        push!(meas, Measurement(kind=:qinj, bus=b, value=0.0, sigma=400.0))
    end
    zi = Set{Tuple{String,String}}()

    bonded = se_net_fourwire(; ground=true)
    ob = PowerOptLab._observability(bonded, meas, "n", zi,
                                    solve_pf(bonded; per_unit=false)["bus"])
    @test ob.observable === true
    @test ob.rank == ob.n_states == 8      # both conductors of both buses
    @test ob.tangent_dimension == 4        # four real KCL rows are exact
    # `redundancy` counts measurement rows beyond the tangent directions they
    # identify: 6 readings, 4 identified. It was -2 while the KCL rows were
    # missing from the diagnostic altogether.
    @test ob.redundancy == 2
    @test ob.min_singular > 1.0

    # With the source neutral unbonded the whole neutral system floats: a genuine
    # gauge freedom, and the one case here that IS rank deficient.
    # The rank is a LOCAL property, so the exact deficiency depends on the point
    # the diagnostic is evaluated at; that the state is not identifiable does not.
    floating = se_net_fourwire(; ground=false)
    of = PowerOptLab._observability(floating, meas, "n", zi, Dict{String,Any}())
    @test of.observable === false
    @test of.rank < of.n_states
    @test of.min_singular < 1e-10
    # More redundancy AND less observability than the bonded network: the
    # missing direction is a gauge freedom, not a measurement gap.
    @test of.redundancy > ob.redundancy
end

@testset "State estimation: Measurement validation" begin
    @test_throws ArgumentError Measurement(kind=:bogus, bus="b", value=1.0, sigma=1.0)
    @test_throws ArgumentError Measurement(kind=:vmag, bus="b", value=NaN, sigma=1.0)
    @test_throws ArgumentError Measurement(kind=:vmag, bus="b", value=Inf, sigma=1.0)
    @test_throws ArgumentError Measurement(kind=:vmag, bus="b", value=1.0, sigma=0.0)
    @test_throws ArgumentError Measurement(kind=:vmag, bus="b", value=1.0, sigma=-1.0)
    @test_throws ArgumentError Measurement(kind=:vmag, bus="b", value=1.0, sigma=NaN)
    @test_throws ArgumentError Measurement(kind=:vmag, bus="",  value=1.0, sigma=1.0)
    @test_throws ArgumentError Measurement(kind=:vmag, bus="b", value=1.0, sigma=1.0, terminal="")
    m = Measurement(kind=:vmag, bus="b", value=230.0, sigma=1.0)
    @test m.reference === missing            # inherits the solve neutral
    m2 = Measurement(kind=:vmag, bus="b", value=230.0, sigma=1.0, reference="n2")
    @test m2.reference == "n2"
end

@testset "State estimation: unconverged solve is not published as an estimate" begin
    meas, _ = se_full_meas()
    # One iteration cannot converge this nonlinear fit ⇒ no feasible point.
    se = solve_state_estimation(se_net(), meas; solver_options=["max_iter" => 0])
    @test !(se.termination_status in ("LOCALLY_SOLVED", "OPTIMAL"))
    @test !solve_status(se).publishable
    @test isnan(se.objective)
    @test all(isnan(se.bus[b]["1"]["vm"]) for b in ("bus1","bus2"))
    @test all(isnan(r.estimated) for r in se.residuals)
end

@testset "State estimation: empty and bad inputs" begin
    @test_throws ArgumentError solve_state_estimation(se_net(), Measurement[])
    @test_throws ArgumentError solve_state_estimation(se_net(), [1, 2, 3])
end

@testset "State estimation: unresolvable reference terminal is diagnosed" begin
    # A two-terminal-free net whose buses have no "n": the default neutral cannot
    # resolve. The error must name the solve's `neutral` default as the source of
    # the terminal, not surface an unknown-OPF-terminal error from the engine.
    net = parse_bmopf("""
    {"bus":{"src":{"terminal_names":["1"]},"b":{"terminal_names":["1"]}},
     "voltage_source":{"s":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
     "linecode":{"lc":{"R_series_1_1":0.5}},
     "line":{"l":{"bus_from":"src","bus_to":"b","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}}}
    """; from_string=true)
    meas = [Measurement(kind=:vmag, bus="b", value=229.0, sigma=1.0),
            Measurement(kind=:pinj, bus="b", value=-1_000.0, sigma=10.0),
            Measurement(kind=:qinj, bus="b", value=0.0, sigma=10.0)]
    err = try
        solve_state_estimation(net, meas); nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("neutral", sprint(showerror, err))
    # Referencing to ground instead resolves, and the estimate is recovered.
    se = solve_state_estimation(net, meas; neutral=nothing)
    @test se.primal_status == "FEASIBLE_POINT"
    @test se.observability.observable === true
end
