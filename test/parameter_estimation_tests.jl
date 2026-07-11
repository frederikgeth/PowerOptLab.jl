using BMOPFTools: parse_bmopf, solve_pf

# ── Fixtures ────────────────────────────────────────────────────────────────
# Per-unit-length line impedance shared by the calibration lines.
const PE_R0, PE_X0 = 0.4, 0.25

# True feeder with a line src─b1, a line b2─b3, and a tapped transformer b1─b2.
# `l` = (p1,q1, p2,q2, p3,q3) loads (W/var) at b1,b2,b3; `L1,L2` line lengths;
# `tau` the transformer tap multiplier.
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

# Calibration physics net: same source and metered loads, but the uncertain
# elements (both lines and the transformer) are omitted — they are the unknowns.
pe_calnet(l) = parse_bmopf("""
{"bus":{
  "src":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
  "b1": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
  "b2": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
  "b3": {"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
 "voltage_source":{"vs":{"bus":"src","terminal_map":["1"],"v_magnitude":[230.0],"v_angle":[0.0]}},
 "load":{"d1":{"bus":"b1","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$(l[1])],"q_nom":[$(l[2])]},
         "d2":{"bus":"b2","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$(l[3])],"q_nom":[$(l[4])]},
         "d3":{"bus":"b3","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[$(l[5])],"q_nom":[$(l[6])]}}}
"""; from_string=true)

# Six diverse-load snapshots (kept light so the forward PF stays well-conditioned).
const PE_LOADS = [(2000.0,400.0,3000.0,600.0,1500.0,300.0),
                  (4000.0,800.0,1500.0,300.0,3500.0,700.0),
                  (3500.0,700.0,4000.0,800.0,2500.0,500.0),
                  (1000.0,200.0,3800.0,760.0,4500.0,900.0),
                  (4500.0,900.0,2500.0,500.0,1500.0,300.0),
                  (3000.0,600.0,2000.0,400.0,4000.0,800.0)]

# Build snapshots + noiseless voltage-magnitude measurements at the metered buses.
function pe_dataset(buses; L1, L2, tau)
    nets = Any[]; meas = Vector{Vector{Measurement}}()
    for l in PE_LOADS
        pf = solve_pf(pe_truenet(l; L1=L1, L2=L2, tau=tau); per_unit=false)
        ms = Measurement[]
        for b in buses
            v = hypot(pf["bus"][b]["1"]["vr"], pf["bus"][b]["1"]["vi"])
            push!(ms, Measurement(kind=:vmag, bus=b, value=v, sigma=0.5))
        end
        push!(nets, pe_calnet(l)); push!(meas, ms)
    end
    return nets, meas
end

# ── Tests ───────────────────────────────────────────────────────────────────

@testset "Parameter estimation: line lengths only" begin
    # Two uncertain lines, tap known/absent (nominal transformer): recover both
    # lengths from noiseless multi-snapshot voltage data.
    nets, meas = pe_dataset(("b1", "b2", "b3"); L1=1.7, L2=1.2, tau=1.0)
    r = solve_parameter_estimation(nets, meas;
        lines=[CalibLine(id="l1", bus_from="src", bus_to="b1", r_per_length=PE_R0, x_per_length=PE_X0),
               CalibLine(id="l2", bus_from="b2", bus_to="b3", r_per_length=PE_R0, x_per_length=PE_X0)],
        taps=[CalibTap(id="t1", bus_from="b1", bus_to="b2", ratio_nom=1.0, tap_min=0.99, tap_max=1.01)])
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.line_length["l1"] ≈ 1.7  rtol=1e-2
    @test r.line_length["l2"] ≈ 1.2  rtol=1e-2
    @test r.tap["t1"] ≈ 1.0  atol=2e-3
    @test r.residual_rms < 1e-2                    # noiseless data ⇒ near-perfect fit
    @test length(r.snapshots) == length(PE_LOADS)  # per-snapshot fitted state returned
end

@testset "Parameter estimation: lengths and tap jointly" begin
    # All three unknowns at once (both lengths + an off-nominal tap).
    nets, meas = pe_dataset(("b1", "b2", "b3"); L1=1.7, L2=1.2, tau=1.05)
    r = solve_parameter_estimation(nets, meas;
        lines=[CalibLine(id="l1", bus_from="src", bus_to="b1", r_per_length=PE_R0, x_per_length=PE_X0),
               CalibLine(id="l2", bus_from="b2", bus_to="b3", r_per_length=PE_R0, x_per_length=PE_X0)],
        taps=[CalibTap(id="t1", bus_from="b1", bus_to="b2", ratio_nom=1.0, tap_min=0.9, tap_max=1.15)])
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.line_length["l1"] ≈ 1.7   rtol=1.5e-2
    @test r.line_length["l2"] ≈ 1.2   rtol=1.5e-2
    @test r.tap["t1"]         ≈ 1.05  rtol=5e-3
    @test r.residual_rms < 1e-2
end

@testset "Parameter estimation: a wrong tap shows up as misfit" begin
    # If the tap is FIXED wrong (bounds pinned away from truth) while only the
    # lengths are free, the fit cannot reach the data — residual stays large.
    nets, meas = pe_dataset(("b1", "b2", "b3"); L1=1.7, L2=1.2, tau=1.05)
    r = solve_parameter_estimation(nets, meas;
        lines=[CalibLine(id="l1", bus_from="src", bus_to="b1", r_per_length=PE_R0, x_per_length=PE_X0),
               CalibLine(id="l2", bus_from="b2", bus_to="b3", r_per_length=PE_R0, x_per_length=PE_X0)],
        taps=[CalibTap(id="t1", bus_from="b1", bus_to="b2", ratio_nom=1.0,
                       tap_init=1.0, tap_min=1.0, tap_max=1.0)])   # pinned at nominal (wrong)
    @test r.tap["t1"] ≈ 1.0  atol=1e-6
    @test r.residual_rms > 1.0    # ~5 % ratio error at ~220 V ⇒ volts of misfit
end

@testset "Parameter estimation: argument validation" begin
    nets, meas = pe_dataset(("b1",); L1=1.7, L2=1.2, tau=1.0)
    aline = CalibLine(id="l1", bus_from="src", bus_to="b1", r_per_length=PE_R0, x_per_length=PE_X0)
    # Nothing to estimate.
    @test_throws ArgumentError solve_parameter_estimation(nets, meas)
    # Measurements not parallel to nets.
    @test_throws ArgumentError solve_parameter_estimation(nets, meas[1:end-1]; lines=[aline])
    # Duplicate ids across lines/taps.
    @test_throws ArgumentError solve_parameter_estimation(nets, meas;
        lines=[aline], taps=[CalibTap(id="l1", bus_from="b1", bus_to="b2")])
    # Non-:vmag measurement.
    bad = [[Measurement(kind=:pinj, bus="b1", value=0.0, sigma=1.0)] for _ in nets]
    @test_throws ArgumentError solve_parameter_estimation(nets, bad; lines=[aline])
end
