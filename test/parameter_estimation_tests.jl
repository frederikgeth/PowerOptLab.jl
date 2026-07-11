using BMOPFTools: parse_bmopf, solve_pf
using Random

# ── Fixtures ────────────────────────────────────────────────────────────────
const PE_R0, PE_X0 = 0.4, 0.25

# True feeder: line src─b1, tapped transformer b1─b2, line b2─b3, loads at b1,b2,b3.
pe_truenet(l; L1, L2, tau) = parse_bmopf("""
{"bus":{
  "src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
  "b1": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
  "b2": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
  "b3": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "linecode":{"lc1":{"R_series_1_1":$(PE_R0*L1),"X_series_1_1":$(PE_X0*L1)},
             "lc2":{"R_series_1_1":$(PE_R0*L2),"X_series_1_1":$(PE_X0*L2)}},
 "line":{"l1":{"bus_from":"src","bus_to":"b1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc1","length":1.0},
         "l2":{"bus_from":"b2","bus_to":"b3","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc2","length":1.0}},
 "transformer":{"single_phase":{"t1":{
     "bus_from":"b1","bus_to":"b2","terminal_map_from":["1","n"],"terminal_map_to":["1","n"],
     "v_nom_from":230.0,"v_nom_to":230.0,"tap":$tau,"s_rating":1e6,"x_series_from":0.001}}},
 "load":{"d1":{"bus":"b1","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$(l[1])],"q_nom":[$(l[2])]},
         "d2":{"bus":"b2","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$(l[3])],"q_nom":[$(l[4])]},
         "d3":{"bus":"b3","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$(l[5])],"q_nom":[$(l[6])]}}}
"""; from_string=true)

# Calibration physics net: source + KNOWN transformer (kept; its tap is freed by
# CalibTap); the uncertain lines are omitted and NO loads are present — injections
# come from the measurements.
pe_calnet() = parse_bmopf("""
{"bus":{
  "src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
  "b1": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
  "b2": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
  "b3": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "transformer":{"single_phase":{"t1":{
     "bus_from":"b1","bus_to":"b2","terminal_map_from":["1","n"],"terminal_map_to":["1","n"],
     "v_nom_from":230.0,"v_nom_to":230.0,"s_rating":1e6,"x_series_from":0.001}}}}
"""; from_string=true)

const PE_LOADS = [(2000.0,400.0,3000.0,600.0,1500.0,300.0),
                  (4000.0,800.0,1500.0,300.0,3500.0,700.0),
                  (3500.0,700.0,4000.0,800.0,2500.0,500.0),
                  (1000.0,200.0,3800.0,760.0,4500.0,900.0),
                  (4500.0,900.0,2500.0,500.0,1500.0,300.0),
                  (3000.0,600.0,2000.0,400.0,4000.0,800.0)]

# Snapshots + noisy (P, Q, |V|_pn) triples at each metered bus. `noise` scales the
# Gaussian perturbation; a fixed seed keeps the test deterministic.
function pe_dataset(; L1, L2, tau, noise=1.0, seed=1)
    rng = Random.MersenneTwister(seed)
    nets = Any[]; meas = Vector{Vector{Measurement}}()
    for l in PE_LOADS
        pf = solve_pf(pe_truenet(l; L1=L1, L2=L2, tau=tau); per_unit=false)
        ms = Measurement[]
        for (bi, b) in enumerate(("b1","b2","b3"))
            v = hypot(pf["bus"][b]["1"]["vr"] - pf["bus"][b]["n"]["vr"],
                      pf["bus"][b]["1"]["vi"] - pf["bus"][b]["n"]["vi"])
            p = l[2bi-1]; q = l[2bi]
            push!(ms, Measurement(kind=:vmag, bus=b, value=v + noise*0.3*randn(rng), sigma=0.3))
            push!(ms, Measurement(kind=:pinj, bus=b, value=-p*(1 + noise*0.01*randn(rng)), sigma=0.01*p+50))
            push!(ms, Measurement(kind=:qinj, bus=b, value=-q*(1 + noise*0.02*randn(rng)), sigma=0.02*q+50))
        end
        push!(nets, pe_calnet()); push!(meas, ms)
    end
    return nets, meas
end

pe_lines() = [CalibLine(id="l1", bus_from="src", bus_to="b1", r_per_length=PE_R0, x_per_length=PE_X0),
              CalibLine(id="l2", bus_from="b2", bus_to="b3", r_per_length=PE_R0, x_per_length=PE_X0)]
pe_taps()  = [CalibTap(id="t1", tap_min=0.9, tap_max=1.15)]

# ── Tests ───────────────────────────────────────────────────────────────────

@testset "Parameter estimation: joint lengths + tap from noisy (P,Q,|V|)" begin
    nets, meas = pe_dataset(; L1=1.7, L2=1.2, tau=1.05)
    r = solve_parameter_estimation(nets, meas; lines=pe_lines(), taps=pe_taps())
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.line_length["l1"] ≈ 1.7   rtol=3e-2
    @test r.line_length["l2"] ≈ 1.2   rtol=3e-2
    @test r.tap["t1"]         ≈ 1.05  rtol=5e-3
    @test r.residual_rms < 0.5                       # ≈ the 0.3 V meter noise floor
    @test length(r.snapshots) == length(PE_LOADS)
end

@testset "Parameter estimation: per-unit equals SI" begin
    nets, meas = pe_dataset(; L1=1.7, L2=1.2, tau=1.05)
    si = solve_parameter_estimation(nets, meas; lines=pe_lines(), taps=pe_taps(), per_unit=false)
    pu = solve_parameter_estimation(nets, meas; lines=pe_lines(), taps=pe_taps(), per_unit=true)
    @test si.line_length["l1"] ≈ pu.line_length["l1"]  rtol=1e-5
    @test si.line_length["l2"] ≈ pu.line_length["l2"]  rtol=1e-5
    @test si.tap["t1"]         ≈ pu.tap["t1"]           rtol=1e-5
    @test si.residual_rms      ≈ pu.residual_rms        rtol=1e-4
end

@testset "Parameter estimation: robust WLAV objective" begin
    nets, meas = pe_dataset(; L1=1.7, L2=1.2, tau=1.05)
    r = solve_parameter_estimation(nets, meas; lines=pe_lines(), taps=pe_taps(), objective=:wlav)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.line_length["l1"] ≈ 1.7   rtol=4e-2
    @test r.tap["t1"]         ≈ 1.05  rtol=1e-2
end

@testset "Parameter estimation: native free tap actually moved" begin
    # Estimating ONLY the tap (lengths known/absent): with the true tap off-nominal,
    # the engine's free-tap variable must move away from 1.0 to fit the data.
    nets, meas = pe_dataset(; L1=1.7, L2=1.2, tau=1.10)
    # Known lines are absent from the physics net here too, so also estimate them,
    # but the point is the tap: it should land near 1.10, not stay at nominal.
    r = solve_parameter_estimation(nets, meas; lines=pe_lines(), taps=pe_taps())
    @test r.tap["t1"] ≈ 1.10  rtol=1e-2
    @test abs(r.tap["t1"] - 1.0) > 0.05           # genuinely off-nominal
end

@testset "Parameter estimation: argument validation" begin
    nets, meas = pe_dataset(; L1=1.7, L2=1.2, tau=1.05)
    aline = CalibLine(id="l1", bus_from="src", bus_to="b1", r_per_length=PE_R0, x_per_length=PE_X0)
    @test_throws ArgumentError solve_parameter_estimation(nets, meas)                       # nothing to estimate
    @test_throws ArgumentError solve_parameter_estimation(nets, meas[1:end-1]; lines=[aline])  # not parallel
    @test_throws ArgumentError solve_parameter_estimation(nets, meas;                       # duplicate id
        lines=[aline], taps=[CalibTap(id="l1")])
    @test_throws ArgumentError solve_parameter_estimation(nets, meas;                       # bad objective
        lines=[aline], objective=:huber)
    @test_throws ArgumentError solve_parameter_estimation(nets, meas;                       # tap id not in net
        taps=[CalibTap(id="nope")])
end
