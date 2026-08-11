@testset "Advanced inverter: plain converter exports up to s_max at unity PF" begin
    inv = AdvancedInverter(id="inv", bus="poc", s_max=5000.0)
    r = solve_advanced_inverter(inv_grid(), inv)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.p_poc ≈ 5000.0  rtol=1e-3          # max export = nameplate
    @test abs(r.q_poc) < 5.0                    # unity power factor
    @test r.p_loss ≈ 0.0  atol=1e-6             # no loss model
    @test r.p_dc ≈ r.p_conv  atol=1e-6          # P_dc = P_conv (+0 loss)
    # With no filter the internal node collapses onto the POC voltage.
    @test r.v_int_mag[1] ≈ r.bus["poc"]["1"]["vm"]  rtol=1e-3
    @test r.v_filter_mag[1] ≈ r.bus["poc"]["1"]["vm"] rtol=1e-3
    @test r.i_grid_mag ≈ r.i_mag rtol=1e-12
    @test r.i_filter_shunt_mag == [0.0]
    @test isnan(r.filter_resonance_hz)
end

@testset "Advanced inverter: output filter reduces POC power below the converter rating" begin
    inv = AdvancedInverter(id="inv", bus="poc", s_max=5000.0,
                           r_filter=0.2, x_filter=0.5, v_int_max=245.0)
    r = solve_advanced_inverter(inv_grid(), inv)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    # The converter operates on its apparent-power circle …
    @test hypot(r.p_conv, r.q_conv) ≈ 5000.0  rtol=1e-2
    # … but the series filter means less active power reaches the grid.
    @test r.p_poc < r.p_conv - 10.0
    @test r.p_poc <= 5000.0 + 1.0
    @test r.v_int_mag[1] <= 245.0 + 1e-3        # EMF within its cap
end

@testset "Advanced inverter: converter losses are non-branching (P_dc = P_conv + P_loss)" begin
    inv = AdvancedInverter(id="inv", bus="poc", s_max=5000.0, r_filter=0.2, x_filter=0.5,
                           p_loss_fixed=20.0, a_loss=0.3, c_loss=0.02)
    r = solve_advanced_inverter(inv_grid(), inv; objective=:min_loss, p_set=3000.0)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.p_poc ≈ 3000.0  rtol=1e-3
    # Three-term loss curve evaluated at the solved current magnitude.
    imag = r.i_mag[1]
    @test r.p_loss ≈ 20.0 + 0.3*imag + 0.02*imag^2  rtol=1e-4
    @test r.p_loss > 0.0
    @test r.p_dc ≈ r.p_conv + r.p_loss  rtol=1e-9   # the single non-branching equation
    @test r.p_dc > r.p_poc                           # DC side supplies more than delivered
end

@testset "Advanced inverter: single-phase double-frequency ripple bound limits export" begin
    common = (bus="poc", s_max=5000.0, r_filter=0.1, x_filter=0.2)
    ref = solve_advanced_inverter(inv_grid(), AdvancedInverter(; id="inv", common...))
    lim = solve_advanced_inverter(inv_grid(),
              AdvancedInverter(; id="inv", p_ripple_max=2500.0, common...))
    @test lim.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test lim.ripple <= 2500.0 + 1.0                 # ripple held at its bound
    @test lim.p_poc < ref.p_poc - 100.0              # export curtailed by the ripple cap
end

@testset "Advanced inverter: three-phase grid-forming holds a balanced internal EMF" begin
    inv = AdvancedInverter(id="inv", bus="poc", phase_terminals=["a","b","c"], neutral="n",
                           s_max=15000.0, r_filter=0.1, x_filter=0.3,
                           grid_forming=true, v_int_min=225.0, v_int_max=245.0)
    r = solve_advanced_inverter(inv_grid3(), inv)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    # Balanced positive-sequence EMF: equal magnitude on every phase …
    @test r.v_int_mag[1] ≈ r.v_int_mag[2]  rtol=1e-3
    @test r.v_int_mag[2] ≈ r.v_int_mag[3]  rtol=1e-3
    @test 225.0 - 1e-3 <= r.v_int_mag[1] <= 245.0 + 1e-3   # within the magnitude box
    # … and balanced three-phase produces essentially no 2ω ripple.
    @test r.ripple < 10.0
    @test hypot(r.p_conv, r.q_conv) <= 15000.0 + 1.0       # respects the rating
end

@testset "Advanced inverter: input validation" begin
    inv = AdvancedInverter(id="inv", bus="poc", s_max=5000.0)
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(), inv; objective=:min_loss)
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(), inv; objective=:bogus)
    # Three-phase topology requires 3 phases and v_dc/c_dc (and In_max for 4-wire).
    bad1 = AdvancedInverter(id="i", bus="poc", s_max=5e3, topology=:FOUR_LEG)  # 1 phase, no dc
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(), bad1)
    bad2 = AdvancedInverter(id="i", bus="poc", phase_terminals=["a","b","c"], neutral="n",
                            s_max=5e3, topology=:FOUR_LEG, v_dc=700.0, c_dc=1e-3)  # no In_max
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(), bad2)
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5e3, topology=:BOGUS))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=-5e3))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5e3, grid_forming=true))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5e3, modulation_max=1.0))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5e3, p_loss_fixed=-1.0))
    bad_dc = (bus="poc", phase_terminals=["a","b","c"], neutral="n",
              s_max=5e3, topology=:THREE_LEG, v_dc=700.0)
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", c_dc=0.0, bad_dc...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", c_dc=1e-3, f=0.0, bad_dc...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", c_dc=1e-3, n_samples=0, bad_dc...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", c_dc=1e-3, m_max=1.01, bad_dc...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5e3, i_negative_max=5.0))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=1e-3,
                         In_max=40.0, dv_mid_max=2.0, bus="poc",
                         phase_terminals=["a","b","c"], neutral="n", s_max=20e3))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", neutral=nothing, s_max=5e3,
                         r_filter_neutral=0.1))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(id="i", bus="poc", phase_terminals=["a","b","c"], neutral="n",
                         s_max=5e3, r_filter_matrix=zeros(3, 3)))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(id="i", bus="poc", phase_terminals=["a","b","c"], neutral="n",
                         s_max=5e3, r_filter=0.1, r_filter_matrix=zeros(4, 4)))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(id="i", bus="poc", phase_terminals=["a","b","c"], neutral="n",
                         s_max=5e3, r_filter_matrix=[1.0 0 0 0; 0 -0.1 0 0;
                                                   0 0 1.0 0; 0 0 0 1.0]))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(id="i", bus="poc", phase_terminals=["a","b","c"], neutral="n",
                         s_max=5e3, x_filter_matrix=[0.1 0.01 0 0; 0 0.1 0 0;
                                                   0 0 0.1 0; 0 0 0 0.1]))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5e3, c_filter_mid=-1e-6))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5e3, r_filter_damping=1.0))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5e3, r_filter_grid=-0.1))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5e3, i_grid_max=0.0))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(id="i", bus="poc", phase_terminals=["a","b","c"], neutral="n",
                         s_max=5e3, r_filter_grid_matrix=zeros(3, 3)))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(id="i", bus="poc", phase_terminals=["a","b","c"], neutral="n",
                         s_max=5e3, r_filter_grid=0.1,
                         r_filter_grid_matrix=zeros(4, 4)))
    splitbase = (bus="poc", phase_terminals=["a","b","c"], neutral="n",
                 topology=:SPLIT_DC, s_max=20e3, v_dc=800.0, c_dc=2.8e-3,
                 In_max=30.0)
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", c_dc_upper=0.0, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", esr_dc=0.01, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(id="i", bus="poc", phase_terminals=["a","b","c"], neutral="n",
            topology=:FOUR_LEG, s_max=20e3, v_dc=800.0, c_dc=2.8e-3,
            In_max=30.0, c_dc_upper=2.5e-3))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", cap_thermal_weights=(1.0, -1.0, 1.0), splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:BOGUS, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:SPWM, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:CENTERED, f_sw=10e3,
                           i_cap_max=20.0, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:SPWM, f_sw=10e3,
                           pwm_carrier_samples=15, i_cap_max=20.0, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_dc_source_r=0.1, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:SPWM, f_sw=10e3,
            pwm_dc_source_r=-0.1, i_cap_max=20.0, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:SPWM, f_sw=10e3,
            pwm_dc_source_l=-1e-3, i_cap_max=20.0, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:SPWM, f_sw=10e3,
            pwm_dc_source_r=0.0, i_cap_max=20.0, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:SPWM, f_sw=10e3,
            pwm_dc_harmonics=129, pwm_carrier_samples=256,
            pwm_dc_source_r=0.1, i_cap_max=20.0, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:SPWM, f_sw=10e3,
            i_cap_max=20.0, pwm_ac_ripple=true, splitbase...)) # no inductance
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_strategy=:SPWM, f_sw=10e3,
            i_cap_max=20.0, pwm_ac_ripple=true, pwm_ac_harmonics=129,
            pwm_carrier_samples=256, x_filter=0.1, splitbase...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", pwm_ac_converter_reserve=(1.0, 0.0, 0.0),
            splitbase...)) # reserve requires i_max
    # Natural offset = 800(2.5-3.5)/(2·6) = -66.7 V; without an actuator a
    # 1 V mean-offset limit is a contradictory parameterisation.
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", c_dc_upper=2.5e-3, c_dc_lower=3.5e-3,
            v_mid_mean_max=1.0, splitbase...))
end

# Shared knobs for the three-phase topology tests.
const _TOPO_COMMON = (bus="poc", phase_terminals=["a","b","c"], neutral="n",
                      s_max=20e3, i_max=40.0, r_filter=0.05, x_filter=0.15, m_max=0.96)

@testset "Advanced inverter: three-phase topologies on a balanced grid" begin
    net = inv_grid3_bal()
    r3 = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:THREE_LEG, v_dc=700.0, c_dc=1.1e-3, _TOPO_COMMON...))
    r4 = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3, In_max=40.0, _TOPO_COMMON...))
    # In_max = 21 A: the 12 A/half bank rating net of the 2ω allocation (the same
    # caps carry both currents). See the split-link note in the component docs.
    rs = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3, In_max=21.0, _TOPO_COMMON...))
    for r in (r3, r4, rs)
        @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
        @test r.i_neutral < 0.1                       # balanced ⇒ no neutral current
        @test r.dv2 < 0.05                            # balanced ⇒ no 2ω ripple
        @test hypot(r.p_conv, r.q_conv) <= 20e3 + 1.0 # converter rating respected
    end
    # All three deliver the same balanced-grid export (topology only differs in
    # how it purchases neutral/zero-sequence capability, absent here).
    @test r3.p_poc ≈ r4.p_poc rtol=1e-3
    @test r4.p_poc ≈ rs.p_poc rtol=1e-3
end

@testset "Advanced inverter: 4-leg draws bounded neutral current on an unbalanced grid" begin
    net = inv_grid3_unbal()
    r = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3, In_max=40.0, _TOPO_COMMON...))
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.i_neutral > 1.0                # unbalance drives real neutral current …
    @test r.i_neutral <= 40.0 + 1e-3       # … within the 4th-leg rating
    @test r.dv2 > 0.5                       # and a non-zero 2ω bus ripple
end

@testset "Advanced inverter: 3-leg carries no neutral current even when unbalanced" begin
    net = inv_grid3_unbal()
    r = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:THREE_LEG, v_dc=700.0, c_dc=1.1e-3, _TOPO_COMMON...))
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.i_neutral < 0.1                # 3-wire: Σ I_abc = 0, no zero-sequence path
    # The internal EMF is genuinely unbalanced to serve the unbalanced grid.
    @test maximum(r.v_int_mag) - minimum(r.v_int_mag) > 10.0
end

@testset "Advanced inverter: split-DC utilization penalty (needs a higher Vdc than 4-leg)" begin
    net = inv_grid3_bal(245.0)             # 245 V demands more DC utilisation
    four = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=650.0, c_dc=1.1e-3, In_max=40.0, _TOPO_COMMON...))
    split_lo = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:SPLIT_DC, v_dc=650.0, c_dc=2.8e-3, In_max=21.0, _TOPO_COMMON...))
    split_hi = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3, In_max=21.0, _TOPO_COMMON...))
    @test four.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")   # 4-leg fine at 650 V
    @test isnan(split_lo.p_poc)                                       # split infeasible at 650 V
    @test all(isnan(split_lo.bus["poc"][ph]["vm"]) for ph in ("a", "b", "c"))
    @test split_hi.termination_status in ("LOCALLY_SOLVED", "OPTIMAL") # feasible at 800 V
    @test !isnan(split_hi.p_poc)
end

@testset "Advanced inverter: 2ω ripple phasor responds to C_dc and honours dv2_max" begin
    net = inv_grid3_unbal()
    solve_kw = (per_unit=true, s_base=100e3)
    big = solve_advanced_inverter(net, AdvancedInverter(; id="i",
        topology=:FOUR_LEG, v_dc=700.0, c_dc=2.0e-3, In_max=40.0,
        _TOPO_COMMON...); solve_kw...)
    small = solve_advanced_inverter(net, AdvancedInverter(; id="i",
        topology=:FOUR_LEG, v_dc=700.0, c_dc=0.3e-3, In_max=40.0,
        _TOPO_COMMON...); solve_kw...)
    @test small.dv2 > big.dv2 + 1.0        # smaller capacitor ⇒ larger bus ripple
    capped = solve_advanced_inverter(net, AdvancedInverter(; id="i",
        topology=:FOUR_LEG, v_dc=700.0, c_dc=0.3e-3, In_max=40.0,
        dv2_max=3.0, _TOPO_COMMON...); solve_kw...)
    @test capped.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test capped.dv2 <= 3.0 + 1e-2         # the amplitude cap binds
end

@testset "Advanced inverter: split-link capacitor allocation is endogenous" begin
    # The split link's half-banks carry the 2ω bus current AND the fundamental
    # neutral current. Different frequencies ⇒ they combine in RMS and share one
    # thermal budget, so tightening the bank rating must SQUEEZE the neutral
    # current rather than leave it at a hand-computed nameplate.
    net = inv_grid3_unbal()
    split = (topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3, In_max=21.0)
    loose = solve_advanced_inverter(net, AdvancedInverter(; id="i", split..., _TOPO_COMMON...))
    tight = solve_advanced_inverter(net, AdvancedInverter(; id="i", split..., i_cap_max=2.0, _TOPO_COMMON...))
    @test tight.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test loose.i_cap > 2.0                       # unconstrained bank current …
    @test tight.i_cap <= 2.0 + 1e-3               # … the rating binds …
    @test tight.i_neutral < loose.i_neutral - 1.0 # … and neutral capability pays for it

    # The reported bank current is exactly the RMS sum of its two components:
    #   i_cap² = (k·|D|)² + (|I_n|/2)²,  k = ωC/√2  (split link: C_eq = C/2)
    k = 2pi * 50.0 * 2.8e-3 / sqrt(2)
    @test tight.i_cap^2 ≈ (k*tight.dv2)^2 + (tight.i_neutral/2)^2  rtol=1e-4

    # i_sw reserves part of the budget, leaving √(3² − 2²) = √5 for the
    # low-frequency terms — but the REPORTED i_cap is the total bank loading, so
    # it reads 3.0 and is directly comparable against i_cap_max.
    sw = solve_advanced_inverter(net,
        AdvancedInverter(; id="i", split..., i_cap_max=3.0, i_sw=2.0, _TOPO_COMMON...))
    @test sw.i_cap ≈ 3.0  rtol=1e-3                        # total, i_sw included
    @test sqrt(sw.i_cap^2 - 2.0^2) ≈ sqrt(5.0)  rtol=1e-2  # low-frequency part
    # Reserving switching headroom costs neutral capability at the SAME rating:
    # the low-frequency budget drops from 3.0 to √(3²−2²) = √5 ≈ 2.24.
    nosw = solve_advanced_inverter(net,
        AdvancedInverter(; id="i", split..., i_cap_max=3.0, _TOPO_COMMON...))
    @test sw.i_neutral < nosw.i_neutral - 0.5
end

@testset "Advanced inverter: split link can use i_cap_max instead of In_max" begin
    # The split link's neutral current flows through the capacitors, so a bank
    # rating bounds it on its own — In_max is optional (this is the documented
    # "endogenous instead of pre-computed" usage).
    net = inv_grid3_unbal()
    r = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:SPLIT_DC,
            v_dc=800.0, c_dc=2.8e-3, i_cap_max=2.0, _TOPO_COMMON...))
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.i_cap <= 2.0 + 1e-3            # the bank rating alone governs …
    @test r.i_neutral > 0.1                # … and it is what bounds |I_n|
    # The 4-leg still needs In_max: its neutral goes through the 4th leg, whose
    # device rating has nothing to do with the capacitor bank.
    @test_throws ArgumentError solve_advanced_inverter(net,
        AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3,
                         i_cap_max=5.0, _TOPO_COMMON...))
    # And a split link with neither knob is still rejected.
    @test_throws ArgumentError solve_advanced_inverter(net,
        AdvancedInverter(; id="i", topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3,
                         _TOPO_COMMON...))
end

@testset "Advanced inverter: 2ω capacitor current pins the RMS convention" begin
    # I_2ω,rms = |S̃|/(√2·V_dc) — an RMS quantity, NOT the peak harmonic
    # amplitude |S̃|/V_dc. The identity is checked at several (C, V_dc) settings;
    # since C does not appear in it, independence from C follows algebraically
    # *at equal |S̃|*. Note these are separate optima — changing C or V_dc moves
    # the operating point — so this pins the FORMULA, not a same-operating-point
    # invariance (capacitance still sets the ripple VOLTAGE, checked via dv2).
    net = inv_grid3_unbal()
    four(c, vdc) = solve_advanced_inverter(net, AdvancedInverter(; id="i",
        topology=:FOUR_LEG, v_dc=vdc, c_dc=c, In_max=40.0, _TOPO_COMMON...))
    a = four(1.1e-3, 700.0)
    @test a.i_cap ≈ a.ripple/(sqrt(2)*700.0)  rtol=1e-4    # RMS, not peak
    @test !isapprox(a.i_cap, a.ripple/700.0; rtol=1e-2)    # peak would be √2× larger

    # Independent of C at fixed operating point; scales as 1/V_dc.
    b = four(4.0e-3, 700.0)
    @test b.i_cap ≈ b.ripple/(sqrt(2)*700.0)  rtol=1e-4
    @test b.dv2 < a.dv2 - 1.0                              # …but dv2 ∝ 1/C
    c = four(1.1e-3, 1400.0)
    @test c.i_cap ≈ c.ripple/(sqrt(2)*1400.0) rtol=1e-4
end

@testset "Advanced inverter: literature ripple-phasor benchmarks" begin
    # Deakin, Heidari & Deng (arXiv:2512.18293), Table V, prescribes six
    # fundamental-current patterns for a 416 V line-line, 20 A peak test system.
    # These equation-level cases are independent of the OPF optimum and pin the
    # unconjugated sum S̃ = Σ V_x I_x used by the stamped model.
    a = cis(2pi/3)
    vph = 416.0 / sqrt(3)                 # RMS phase voltage
    iref = 20.0 / sqrt(2)                 # 20 A peak -> RMS phasor
    V = vph .* ComplexF64[1, a^-1, a^-2]
    currents = (
        iref .* ComplexF64[1, a^-1, a^-2],       # 3a: balanced P
        iref .* ComplexF64[1, 0, 0],              # 3b: phase-a P
        iref .* ComplexF64[1, -1/2, -1/2],        # 3c: a to b,c
        iref .* ComplexF64[1, a, a^2],            # 3d: negative sequence
        iref .* ComplexF64[im, im*a^-1, im*a^-2], # 3e: balanced Q
        iref .* ComplexF64[im, 0, 0],             # 3f: phase-a Q
    )
    ripple(I) = begin
        re, im_ = PowerOptLab._ripple_components(real.(V), imag.(V),
                                                  real.(I), imag.(I))
        hypot(re, im_)
    end
    predicted = collect(ripple.(currents))
    unit = vph * iref
    @test predicted ≈ [0.0, unit, 1.5unit, 3unit, 0.0, unit] atol=1e-8

    # The ideal-terminal predictions reproduce the paper's PLECS values for the
    # strongest non-cancelling cases (small differences are filter/controller
    # effects in the time-domain model): 5.122 kW and 10.200 kW in Table VII.
    @test predicted[3] / 1e3 ≈ 5.122 rtol=0.01
    @test predicted[4] / 1e3 ≈ 10.200 rtol=0.01

    # S̃ is a sinusoidal AMPLITUDE. The capacitor-current amplitude is |S̃|/Vdc;
    # PowerOptLab reports RMS loading, hence the additional 1/sqrt(2).
    vdc = 700.0
    @test predicted[4] / (sqrt(2)*vdc) ≈ (predicted[4] / vdc) / sqrt(2)

    # Table VII is also retained as a literal regression fixture. The paper's
    # reported current uses its Fourier-magnitude convention, I_2ω=P_2ω/Vdc;
    # it should not be confused with this package's time-domain RMS bank loading.
    proposed_kw = [0.0, 3.458, 5.122, 10.200, 0.0, 2.772]
    plecs_kw    = [6.68e-4, 3.466, 5.127, 10.207, 5.57e-3, 2.777]
    plecs_irms  = [1.9e-3, 4.979, 7.386, 14.582, 8.1e-3, 3.967]
    proposed_i  = [0.0, 4.940, 7.317, 14.571, 0.0, 3.960]
    active = [2, 3, 4, 6]
    @test proposed_i[active] ≈ proposed_kw[active] .* 1e3 ./ vdc rtol=2e-4
    @test maximum(abs.((plecs_kw[active] .- proposed_kw[active]) ./
                       proposed_kw[active])) < 3e-3
    @test maximum(abs.((plecs_irms[active] .- proposed_i[active]) ./
                       proposed_i[active])) < 1e-2
    @test plecs_kw[1] < 1e-3 && plecs_kw[5] < 1e-2 # numerical residue in cancelling cases
end

@testset "Advanced inverter: sequence components split neutral vs 2ω ripple" begin
    # S̃ = Σ_x U_x I_x is the UNCONJUGATED sum, so with positive-sequence voltage
    # a zero-sequence current contributes Σ_x U_x · I₀ = I₀ Σ_x U_x = 0. The two
    # unbalance mechanisms therefore load different hardware, and superpose:
    #   negative sequence ⇒ 2ω ripple (capacitor heating), NO neutral current
    #   zero sequence     ⇒ neutral current (4th leg / caps), NO 2ω ripple
    seqnet(; vp=230.0, vn=0.0, v0=0.0) = begin
        a = cis(2pi/3)
        Va = vp + vn + v0; Vb = a^2*vp + a*vn + v0; Vc = a*vp + a^2*vn + v0
        inv_grid3_src(mags=[abs(Va), abs(Vb), abs(Vc)],
                      angs=[angle(Va), angle(Vb), angle(Vc)], vmax=300.0)
    end
    run4(net) = solve_advanced_inverter(net, AdvancedInverter(; id="i",
        topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3, In_max=40.0, _TOPO_COMMON...))
    bal  = run4(seqnet())
    negs = run4(seqnet(vn=0.08*230))
    zers = run4(seqnet(v0=0.08*230))
    both = run4(seqnet(vn=0.08*230, v0=0.08*230))

    @test bal.i_neutral < 1e-2 && bal.ripple < 10.0          # balanced ⇒ neither
    # Negative sequence: large oscillating power, no zero-sequence current.
    @test negs.ripple > 3000.0
    @test negs.i_neutral < 1e-2
    @test negs.i_negative > 1.0
    @test negs.i_zero < 1e-2
    # Zero sequence: the mirror image — neutral current, negligible ripple.
    @test zers.i_neutral > 5.0
    @test zers.i_neutral ≈ 3 * zers.i_zero rtol=1e-4
    @test zers.i_negative < 1e-2
    @test zers.ripple < negs.ripple/10
    # And the two superpose rather than trading against each other.
    @test both.i_neutral > 5.0
    @test both.ripple > 3000.0
    # The bank current tracks the ripple (4-leg: neutral bypasses the caps).
    @test negs.i_cap ≈ negs.ripple/(sqrt(2)*700.0)  rtol=1e-3
end

@testset "Advanced inverter: Fortescue transform convention" begin
    a = cis(2pi/3)
    phasors = (ComplexF64[1, 1, 1],
               ComplexF64[1, a^2, a],
               ComplexF64[1, a, a^2])
    for (expected, I) in enumerate(phasors)
        seq = PowerOptLab._sequence_components(real.(I), imag.(I))
        mags = [hypot(s...) for s in seq]
        @test mags[expected] ≈ 1.0 atol=1e-12
        @test sum(mags) ≈ 1.0 atol=1e-12
    end
end

@testset "Advanced inverter: symmetrical-component limits bind independently" begin
    a = cis(2pi/3); vp = 230.0
    seqnet(; vn=0.0, v0=0.0) = begin
        Va = vp + vn + v0; Vb = a^2*vp + a*vn + v0; Vc = a*vp + a^2*vn + v0
        inv_grid3_src(mags=[abs(Va), abs(Vb), abs(Vc)],
                      angs=[angle(Va), angle(Vb), angle(Vc)], vmax=300.0)
    end
    common = (topology=:FOUR_LEG, v_dc=750.0, c_dc=1.1e-3, In_max=60.0)
    neg = solve_advanced_inverter(seqnet(vn=0.12vp),
        AdvancedInverter(; id="i", common..., i_negative_max=1.0, _TOPO_COMMON...))
    zer = solve_advanced_inverter(seqnet(v0=0.12vp),
        AdvancedInverter(; id="i", common..., i_zero_max=1.0, _TOPO_COMMON...))
    @test neg.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test zer.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test neg.i_negative <= 1.0 + 1e-3
    @test zer.i_zero <= 1.0 + 1e-3
    @test zer.i_neutral <= 3.0 + 3e-3
end

@testset "Advanced inverter: neutral filter KVL includes the shared return drop" begin
    rp, rn = 0.05, 0.20
    r = solve_advanced_inverter(inv_grid3_unbal(), AdvancedInverter(; id="i",
        topology=:FOUR_LEG, v_dc=750.0, c_dc=1.1e-3, In_max=60.0,
        r_filter_neutral=rn, merge(_TOPO_COMMON, (r_filter=rp,))...))
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.i_neutral > 0.1
    # Summing phase-to-neutral converter power counts both phase-conductor and
    # common neutral-conductor I²R losses exactly once.
    filter_loss = rp*sum(abs2, r.i_mag) + rn*r.i_neutral^2
    @test r.p_conv - r.p_poc ≈ filter_loss rtol=2e-4
end

@testset "Advanced inverter: primitive conductor-domain filter matrix" begin
    # Equation-level check with unequal self terms and mutual coupling. The
    # primitive is applied to [Ia,Ib,Ic,In], then the neutral drop is subtracted
    # to recover phase-to-neutral KVL.
    R = [0.11 0.01 0.00 0.02;
         0.01 0.13 0.01 0.01;
         0.00 0.01 0.12 0.02;
         0.02 0.01 0.02 0.20]
    X = [0.30 0.04 0.02 0.01;
         0.04 0.32 0.03 0.02;
         0.02 0.03 0.31 0.01;
         0.01 0.02 0.01 0.25]
    Iph = ComplexF64[12+3im, -7+5im, 2-4im]
    dr, di = PowerOptLab._filter_voltage_drops(
        R, X, real.(Iph), imag.(Iph), 3, true)
    J = [Iph; -sum(Iph)]
    primitive_drop = complex.(dr, di)
    expected = (R + im*X) * J
    @test primitive_drop ≈ expected[1:3] .- expected[4] atol=1e-12

    # A diagonal primitive is exactly the legacy scalar phase/neutral model.
    net = inv_grid3_unbal()
    base = (id="i", bus="poc", phase_terminals=["a","b","c"], neutral="n",
            s_max=20e3, i_max=40.0, topology=:FOUR_LEG, v_dc=750.0,
            c_dc=1.1e-3, In_max=60.0, m_max=0.96)
    scalar = solve_advanced_inverter(net, AdvancedInverter(; base...,
        r_filter=0.05, x_filter=0.15, r_filter_neutral=0.20,
        x_filter_neutral=0.25))
    Rdiag = [0.05 0 0 0; 0 0.05 0 0; 0 0 0.05 0; 0 0 0 0.20]
    Xdiag = [0.15 0 0 0; 0 0.15 0 0; 0 0 0.15 0; 0 0 0 0.25]
    matrix = solve_advanced_inverter(net, AdvancedInverter(; base...,
        r_filter_matrix=Rdiag, x_filter_matrix=Xdiag))
    @test matrix.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test matrix.p_poc ≈ scalar.p_poc rtol=2e-5
    @test matrix.q_poc ≈ scalar.q_poc rtol=2e-5
    @test matrix.i_neutral ≈ scalar.i_neutral rtol=2e-5
    @test matrix.ripple ≈ scalar.ripple rtol=2e-5
end

@testset "Advanced inverter: explicit LCL filter limiting cases and losses" begin
    # With C=0, splitting one series impedance across the two arms must reduce
    # exactly to the original single-arm circuit.
    reduced = solve_advanced_inverter(inv_grid(), AdvancedInverter(
        id="i", bus="poc", s_max=5e3, r_filter=0.10, x_filter=0.20);
        objective=:min_loss, p_set=3000.0)
    split = solve_advanced_inverter(inv_grid(), AdvancedInverter(
        id="i", bus="poc", s_max=5e3,
        r_filter=0.04, x_filter=0.08,
        r_filter_grid=0.06, x_filter_grid=0.12);
        objective=:min_loss, p_set=3000.0)
    @test split.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test split.p_conv ≈ reduced.p_conv rtol=2e-6
    @test split.q_conv ≈ reduced.q_conv atol=0.01
    @test split.i_mag ≈ reduced.i_mag rtol=2e-6
    @test split.i_grid_mag ≈ split.i_mag rtol=2e-8
    @test split.i_filter_shunt_mag[1] < 1e-10

    # A damped midpoint capacitor separates converter and grid currents. Real
    # filter loss must equal both terminal power difference and the sum of the
    # two arm I²R terms plus the physical damping-resistor loss.
    c, rd = 20e-6, 1.0
    lcl = solve_advanced_inverter(inv_grid(), AdvancedInverter(
        id="i", bus="poc", s_max=5e3,
        r_filter=0.05, x_filter=0.10,
        r_filter_grid=0.05, x_filter_grid=0.10,
        c_filter_mid=c, r_filter_damping=rd, i_grid_max=30.0);
        objective=:min_loss, p_set=3000.0)
    @test lcl.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test abs(lcl.i_grid_mag[1] - lcl.i_mag[1]) > 0.01
    xc = 1 / (2pi*50.0*c)
    ymag = 1 / hypot(rd, xc)
    @test lcl.i_filter_shunt_mag[1] ≈ ymag*lcl.v_filter_mag[1] rtol=2e-5
    expected_loss = 0.05*lcl.i_mag[1]^2 + 0.05*lcl.i_grid_mag[1]^2 +
                    rd*lcl.i_filter_shunt_mag[1]^2
    @test lcl.p_filter_loss ≈ lcl.p_conv - lcl.p_poc atol=1e-8
    @test lcl.p_filter_loss ≈ expected_loss rtol=2e-5
    l1 = 0.10/(2pi*50.0); l2 = 0.10/(2pi*50.0)
    fres = sqrt((l1+l2)/(l1*l2*c))/(2pi)
    @test lcl.filter_resonance_hz ≈ fres rtol=1e-12

    # Diagonal primitive matrices on both arms reproduce the scalar LCL case.
    R1 = [0.05 0.0; 0.0 0.0]; X1 = [0.10 0.0; 0.0 0.0]
    R2 = [0.05 0.0; 0.0 0.0]; X2 = [0.10 0.0; 0.0 0.0]
    mlcl = solve_advanced_inverter(inv_grid(), AdvancedInverter(
        id="i", bus="poc", s_max=5e3,
        r_filter_matrix=R1, x_filter_matrix=X1,
        r_filter_grid_matrix=R2, x_filter_grid_matrix=X2,
        c_filter_mid=c, r_filter_damping=rd, i_grid_max=30.0);
        objective=:min_loss, p_set=3000.0)
    @test mlcl.p_conv ≈ lcl.p_conv rtol=2e-6
    @test mlcl.i_grid_mag ≈ lcl.i_grid_mag rtol=2e-6
    @test isnan(mlcl.filter_resonance_hz) # no unique scalar estimate in matrix mode
end

@testset "Advanced inverter: LCL filter composes with a four-leg topology" begin
    common = (bus="poc", phase_terminals=["a","b","c"], neutral="n",
              topology=:FOUR_LEG, s_max=20e3, i_max=50.0, i_grid_max=50.0,
              v_dc=750.0, c_dc=1.1e-3, In_max=60.0, m_max=0.96,
              r_filter=0.02, x_filter=0.06,
              r_filter_grid=0.03, x_filter_grid=0.09,
              c_filter_mid=30e-6, r_filter_damping=0.5)
    r = solve_advanced_inverter(inv_grid3_unbal(), AdvancedInverter(; id="i", common...))
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.i_neutral > 0.1
    @test maximum(r.i_filter_shunt_mag) > 1.0
    @test maximum(abs.(r.i_mag .- r.i_grid_mag)) > 0.01
    @test r.p_filter_loss > 0.0
    @test r.switching_margin > 0.0
end

@testset "Advanced inverter: fourth-leg current participates in converter loss" begin
    af, cf = 0.4, 0.015
    r = solve_advanced_inverter(inv_grid3_unbal(), AdvancedInverter(; id="i",
        topology=:FOUR_LEG, v_dc=750.0, c_dc=1.1e-3, In_max=60.0,
        a_loss=af, c_loss=cf, _TOPO_COMMON...);
        objective=:max_export, per_unit=true)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.i_neutral > 1.0
    expected = af*(sum(r.i_mag) + r.i_neutral) +
               cf*(sum(abs2, r.i_mag) + r.i_neutral^2)
    @test r.p_loss ≈ expected rtol=2e-4 atol=0.01
end

@testset "Advanced inverter: split-link midpoint ripple is reported and bounded" begin
    w = 2pi*50.0; c = 2.8e-3
    common = (topology=:SPLIT_DC, v_dc=800.0, c_dc=c, In_max=40.0)
    loose = solve_advanced_inverter(inv_grid3_unbal(),
        AdvancedInverter(; id="i", common..., _TOPO_COMMON...))
    tight = solve_advanced_inverter(inv_grid3_unbal(),
        AdvancedInverter(; id="i", common..., dv_mid_max=0.5, _TOPO_COMMON...))
    @test loose.dv_mid ≈ loose.i_neutral/(2w*c) rtol=1e-4
    @test tight.dv_mid <= 0.5 + 1e-3
    @test tight.i_neutral <= 2w*c*0.5 + 1e-3
end

@testset "Advanced inverter: asymmetric split-link charge and current sharing" begin
    w = 2pi*50.0; vdc = 900.0
    cu, cl = 2.5e-3, 3.5e-3
    csum = cu + cl; ceq = cu*cl/csum
    common = (bus="poc", phase_terminals=["a","b","c"], neutral="n",
              topology=:SPLIT_DC, s_max=20e3, i_max=50.0, In_max=40.0,
              v_dc=vdc, c_dc=3e-3, m_max=0.96, r_filter=0.05, x_filter=0.15)

    symmetric = solve_advanced_inverter(inv_grid3_unbal(),
        AdvancedInverter(; id="i", common...))
    explicit_symmetric = solve_advanced_inverter(inv_grid3_unbal(),
        AdvancedInverter(; id="i", c_dc_upper=3e-3, c_dc_lower=3e-3, common...))
    @test explicit_symmetric.p_poc ≈ symmetric.p_poc rtol=2e-6
    @test explicit_symmetric.dv2 ≈ symmetric.dv2 rtol=2e-6
    @test explicit_symmetric.dv_mid ≈ symmetric.dv_mid rtol=2e-6
    @test explicit_symmetric.v_mid_mean ≈ 0.0 atol=1e-12
    @test explicit_symmetric.i_cap_upper ≈ explicit_symmetric.i_cap_lower rtol=1e-8
    @test explicit_symmetric.i_cap ≈ explicit_symmetric.i_cap_upper rtol=1e-12

    asym = solve_advanced_inverter(inv_grid3_unbal(),
        AdvancedInverter(; id="i", c_dc_upper=cu, c_dc_lower=cl, common...))
    @test asym.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test asym.v_mid_mean ≈ vdc*(cu-cl)/(2csum) atol=1e-10
    @test asym.q_mid_balance == 0.0
    @test asym.dv2 ≈ asym.ripple/(2w*ceq*vdc) rtol=1e-4
    @test asym.dv_mid ≈ asym.i_neutral/(w*csum) rtol=1e-4
    i2 = asym.ripple/(sqrt(2)*vdc)
    @test asym.i_cap_upper^2 ≈ i2^2 + (cu/csum*asym.i_neutral)^2 rtol=1e-4
    @test asym.i_cap_lower^2 ≈ i2^2 + (cl/csum*asym.i_neutral)^2 rtol=1e-4
    @test asym.i_cap_lower > asym.i_cap_upper + 0.5
    @test asym.i_cap ≈ asym.i_cap_lower rtol=1e-12
    @test asym.switching_margin > 0.0
end

@testset "Advanced inverter: split-link balancing and frequency-weighted heating" begin
    w = 2pi*50.0; vdc = 900.0
    cu, cl = 2.5e-3, 3.5e-3; csum = cu + cl
    weights = (2.0, 1.0, 0.2); esru, esrl = 0.03, 0.05
    common = (bus="poc", phase_terminals=["a","b","c"], neutral="n",
              topology=:SPLIT_DC, s_max=20e3, i_max=50.0, In_max=40.0,
              v_dc=vdc, c_dc=3e-3, c_dc_upper=cu, c_dc_lower=cl,
              m_max=0.96, r_filter=0.05, x_filter=0.15)
    r = solve_advanced_inverter(inv_grid3_unbal(), AdvancedInverter(; id="i",
        q_mid_balance_max=0.25, v_mid_mean_max=1.0,
        i_cap_max=15.0, cap_thermal_weights=weights,
        esr_dc_upper=esru, esr_dc_lower=esrl, common...))
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test abs(r.v_mid_mean) <= 1.0 + 1e-5
    natural = vdc*(cu-cl)/(2csum)
    @test r.v_mid_mean ≈ natural + 2r.q_mid_balance/csum rtol=1e-8 atol=1e-8
    @test abs(r.q_mid_balance) <= 0.25 + 1e-8
    i2 = r.ripple/(sqrt(2)*vdc)
    expected_tu2 = weights[2]*i2^2 +
        weights[1]*(cu/csum*r.i_neutral)^2 + weights[3]*0.0^2
    expected_tl2 = weights[2]*i2^2 +
        weights[1]*(cl/csum*r.i_neutral)^2 + weights[3]*0.0^2
    @test r.i_cap_thermal_upper^2 ≈ expected_tu2 rtol=1e-4
    @test r.i_cap_thermal_lower^2 ≈ expected_tl2 rtol=1e-4
    @test r.i_cap_thermal ≈ max(r.i_cap_thermal_upper,
                                r.i_cap_thermal_lower) rtol=1e-12
    expected_loss = esru*expected_tu2 + esrl*expected_tl2
    @test r.p_cap_loss ≈ expected_loss rtol=1e-4
    @test r.p_dc ≈ r.p_conv + r.p_loss + r.p_cap_loss rtol=1e-9

    # Insufficient charge authority cannot meet the same mean-offset bound.
    bad = solve_advanced_inverter(inv_grid3_unbal(), AdvancedInverter(; id="i",
        q_mid_balance_max=0.10, v_mid_mean_max=1.0, common...))
    @test !(bad.termination_status in ("LOCALLY_SOLVED", "OPTIMAL"))
    @test isnan(bad.v_mid_mean)
end

@testset "Advanced inverter: carrier PWM audit reproduces published RMS formulas" begin
    a = cis(2pi/3)
    Ipk = 1.0
    for m in (0.1, 0.3, 0.5), strategy in (:SPWM, :CENTERED)
        U = (m/sqrt(2)) .* ComplexF64[1, a^-1, a]
        I = (Ipk/sqrt(2)) .* ComplexF64[1, a^-1, a]
        inv = AdvancedInverter(id="pwm", bus="poc", phase_terminals=["a","b","c"],
            neutral="n", topology=:FOUR_LEG, s_max=10.0, In_max=10.0,
            v_dc=1.0, c_dc=1.0, i_cap_max=10.0, f=0.01, f_sw=1.0,
            pwm_strategy=strategy, pwm_fundamental_samples=720,
            pwm_carrier_samples=512)
        audit = PowerOptLab._pwm_ripple_audit(inv, real.(U), imag.(U),
                                               real.(I), imag.(I))
        expected = if strategy == :SPWM
            m*Ipk*sqrt(15pi - 88sqrt(3)*m + 45pi*m^2)/(8sqrt(5pi))
        else
            m*Ipk*sqrt(120pi - 704sqrt(3)*m +
                (540pi - 405sqrt(3))*m^2)/(16sqrt(10pi))
        end
        # Mandrioli et al. (2021), Eqs. (40) and (42), with Cdc=fsw=1.
        @test audit.dv_rms ≈ expected rtol=2e-3 atol=2e-5
        @test audit.i_rms > 0.0
        @test audit.dv_pp >= 2*audit.dv_rms
        @test audit.modulation_margin >= -1e-10
    end

    # Switching-current RMS is independent of fsw under the frozen-fundamental
    # assumption; voltage ripple scales as 1/(Cdc*fsw) and both scale with I.
    m = 0.3
    U = (m/sqrt(2)) .* ComplexF64[1, a^-1, a]
    I = (1/sqrt(2)) .* ComplexF64[1, a^-1, a]
    common = (id="pwm", bus="poc", phase_terminals=["a","b","c"], neutral="n",
              topology=:FOUR_LEG, s_max=10.0, In_max=10.0, v_dc=1.0,
              i_cap_max=10.0, f=0.01, pwm_strategy=:CENTERED,
              pwm_fundamental_samples=360, pwm_carrier_samples=256)
    r1 = PowerOptLab._pwm_ripple_audit(
        AdvancedInverter(; c_dc=1.0, f_sw=1.0, common...),
        real.(U), imag.(U), real.(I), imag.(I))
    r2 = PowerOptLab._pwm_ripple_audit(
        AdvancedInverter(; c_dc=2.0, f_sw=2.0, common...),
        real.(U), imag.(U), 2 .* real.(I), 2 .* imag.(I))
    @test r2.i_rms ≈ 2*r1.i_rms rtol=2e-3
    @test r2.dv_rms ≈ r1.dv_rms/2 rtol=3e-3
    @test r2.dv_pp ≈ r1.dv_pp/2 rtol=3e-3
end

@testset "Advanced inverter: finite DC-source harmonic network" begin
    omega, ceq = 2pi*2e3, 1.5e-3
    ib = 3.0 + 4.0im
    open = PowerOptLab._dc_link_harmonic_response(ib, omega, ceq, nothing, 0.0)
    @test open.cap_current == -ib
    @test open.source_current == 0.0im
    @test ib + open.cap_current + open.source_current == 0.0im
    @test open.admittance_margin == 1.0

    shared = PowerOptLab._dc_link_harmonic_response(ib, omega, ceq, 0.2, 1e-3)
    @test ib + shared.cap_current + shared.source_current ≈ 0.0im atol=1e-12
    @test shared.admittance_margin > 0.0

    # A lossless source inductor in parallel with C has a carrier-harmonic
    # antiresonance at omega^2*L*C=1. The helper reports, rather than hides, the
    # singularity so callers can reject that parameterisation.
    resonant_l = inv(omega^2 * ceq)
    resonant = PowerOptLab._dc_link_harmonic_response(
        ib, omega, ceq, nothing, resonant_l)
    @test resonant.admittance_margin <= 100eps(Float64)
    @test !isfinite(resonant.voltage)

    a = cis(2pi/3)
    U = (0.3/sqrt(2)) .* ComplexF64[1, a^-1, a]
    I = (1/sqrt(2)) .* ComplexF64[1, a^-1, a]
    common = (id="pwm", bus="poc", phase_terminals=["a","b","c"],
              neutral="n", topology=:FOUR_LEG, s_max=10.0, In_max=10.0,
              v_dc=1.0, c_dc=1.0, i_cap_max=10.0, f=0.01, f_sw=1.0,
              pwm_strategy=:CENTERED, pwm_fundamental_samples=72,
              pwm_carrier_samples=128, pwm_dc_harmonics=64)
    audit(inv) = PowerOptLab._pwm_ripple_audit(
        inv, real.(U), imag.(U), real.(I), imag.(I))
    open_audit = audit(AdvancedInverter(; common...))
    high_z = audit(AdvancedInverter(; pwm_dc_source_r=1e6, common...))
    low_z = audit(AdvancedInverter(; pwm_dc_source_r=0.01, common...))
    @test high_z.i_rms ≈ open_audit.i_rms rtol=2e-4
    @test high_z.dv_rms ≈ open_audit.dv_rms rtol=3e-2
    @test high_z.i_source_rms < 1e-5*high_z.i_bridge_rms
    @test low_z.i_source_rms > low_z.i_rms
    @test low_z.i_rms < open_audit.i_rms
    @test low_z.dv_rms < open_audit.dv_rms
    @test low_z.source_loss ≈ 0.01*low_z.i_source_rms^2 rtol=1e-10
    @test low_z.dc_network_margin > 0.0

    cu, cl = 2.5e-3, 3.5e-3
    split_parameters = merge(common, (topology=:SPLIT_DC,
        c_dc=3e-3, c_dc_upper=cu, c_dc_lower=cl,
        pwm_strategy=:SPWM, pwm_dc_source_r=0.1))
    split = AdvancedInverter(; split_parameters...)
    split_audit = audit(split)
    ceq_split = cu*cl/(cu+cl)
    @test split_audit.dv_upper_rms ≈
          ceq_split/cu*split_audit.dv_rms rtol=1e-12
    @test split_audit.dv_lower_rms ≈
          ceq_split/cl*split_audit.dv_rms rtol=1e-12
    @test split_audit.dv_upper_rms + split_audit.dv_lower_rms ≈
          split_audit.dv_rms rtol=1e-12
    @test split_audit.dv_upper_pp + split_audit.dv_lower_pp ≈
          split_audit.dv_pp rtol=1e-12
end

@testset "Advanced inverter: automatic PWM reserve closes the capacitor constraint" begin
    net = inv_grid3_bal()
    solve_kw = (per_unit=true, s_base=100e3)
    common = (topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3,
              In_max=40.0, i_cap_max=10.0, esr_dc=0.02, _TOPO_COMMON...)
    fixed = solve_advanced_inverter(net, AdvancedInverter(; id="i", common...);
        solve_kw...)
    pwm = solve_advanced_inverter(net, AdvancedInverter(; id="i", i_sw=2.0,
        pwm_strategy=:CENTERED, f_sw=10e3, pwm_fundamental_samples=72,
        pwm_carrier_samples=128, common...); solve_kw...)
    @test pwm.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test pwm.pwm_iterations >= 2
    @test pwm.i_cap_switching > 1.0
    @test pwm.i_cap_switching_reserved + 2e-2 >= pwm.i_cap_switching
    @test pwm.pwm_reserve_margin >= -2e-2
    @test pwm.i_cap_thermal <= 10.0 + 1e-3
    @test pwm.p_poc < fixed.p_poc - 100.0
    @test pwm.dv_switching_rms > 0.0
    @test pwm.dv_switching_pp > 2*pwm.dv_switching_rms
    @test pwm.pwm_modulation_margin > 0.0
    @test pwm.p_cap_loss ≈
          0.02*(2.0^2 + pwm.i_cap_switching_reserved^2) rtol=1e-4
    @test pwm.p_dc ≈ pwm.p_conv + pwm.p_loss + pwm.p_cap_loss rtol=1e-9

    finite_parameters = (id="i", i_sw=2.0, pwm_strategy=:CENTERED,
        f_sw=10e3, pwm_fundamental_samples=72, pwm_carrier_samples=128,
        pwm_dc_harmonics=64, pwm_dc_source_r=0.01, common...)
    finite = solve_advanced_inverter(net, AdvancedInverter(; finite_parameters...);
        solve_kw...)
    @test finite.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test finite.i_dc_bridge_switching_rms > 0.0
    @test finite.i_dc_source_switching_rms > 0.0
    @test finite.i_cap_switching < finite.i_dc_bridge_switching_rms
    @test finite.p_dc_source_switching_loss ≈
          0.01*finite.i_dc_source_switching_rms^2 rtol=1e-10
    @test finite.pwm_dc_network_margin > 0.0
    # The split link runs with DC-rail and bank headroom on purpose: this case
    # regresses the SPWM reserve closure, not the modulation or thermal limit
    # (those bind in their own testsets). Sitting on several active nonlinear
    # constraints at once is what made this instance solver-path-sensitive.
    split = solve_advanced_inverter(inv_grid3_unbal(), AdvancedInverter(; id="i",
        topology=:SPLIT_DC, v_dc=950.0, c_dc=3e-3,
        c_dc_upper=2.5e-3, c_dc_lower=3.5e-3, In_max=40.0,
        i_cap_max=30.0, pwm_strategy=:SPWM, f_sw=12e3,
        pwm_fundamental_samples=72, pwm_carrier_samples=128, _TOPO_COMMON...);
        solve_kw...)
    @test split.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test split.i_cap_switching > 0.0
    @test split.pwm_reserve_margin >= -2e-2
    @test split.dv_switching_pp > split.dv_switching_rms > 0.0
    @test split.pwm_modulation_margin > 0.0
end

@testset "Advanced inverter: split-link AC ripple matches published closed forms" begin
    a = cis(2pi/3)
    vdc, l, fsw = 100.0, 1.73e-3, 3.6e3
    scale = vdc/(2l*fsw)
    for m in (0.1, 0.3, 0.5)
        U = (m*vdc/sqrt(2)) .* ComplexF64[1, a^-1, a]
        inv = AdvancedInverter(id="acpwm", bus="poc",
            phase_terminals=["a","b","c"], neutral="n",
            topology=:SPLIT_DC, s_max=10e3, In_max=100.0,
            v_dc=vdc, c_dc=2e-3, i_cap_max=100.0,
            f=50.0, f_sw=fsw, pwm_strategy=:SPWM,
            pwm_ac_ripple=true, pwm_fundamental_samples=72,
            pwm_carrier_samples=256, pwm_ac_harmonics=128,
            x_filter=2pi*50*l)
        audit = PowerOptLab._pwm_ac_ripple_audit(
            inv, real.(U), imag.(U))
        phase_expected = scale/(4sqrt(3)) *
            sqrt(1 - 4m^2 + 6m^4) # Mandrioli et al. (2021), Eq. (15)
        neutral_expected = scale *
            sqrt(3/16 - 9m^2/8 + 2sqrt(3)/pi*m^3) # Eq. (27)
        @test audit.converter_rms ≈ fill(phase_expected, 3) rtol=3e-3
        @test audit.neutral_rms ≈ neutral_expected rtol=3e-3
        @test audit.grid_rms ≈ audit.converter_rms rtol=1e-12
        @test maximum(audit.shunt_rms) < 1e-12
        @test audit.converter_pp ≈ fill(vdc/(4l*fsw), 3) rtol=1e-2
    end

    # The ideal inductive result scales inversely with L and switching frequency.
    m = 0.3
    U = (m*vdc/sqrt(2)) .* ComplexF64[1, a^-1, a]
    base = AdvancedInverter(id="acpwm", bus="poc",
        phase_terminals=["a","b","c"], neutral="n",
        topology=:SPLIT_DC, s_max=10e3, In_max=100.0,
        v_dc=vdc, c_dc=2e-3, i_cap_max=100.0, f=50.0,
        f_sw=fsw, pwm_strategy=:SPWM, pwm_ac_ripple=true,
        pwm_fundamental_samples=36, pwm_carrier_samples=128,
        pwm_ac_harmonics=64, x_filter=2pi*50*l)
    doubled = AdvancedInverter(; merge(
        NamedTuple{fieldnames(AdvancedInverter)}(
            Tuple(getfield(base, n) for n in fieldnames(AdvancedInverter))),
        (f_sw=2fsw, x_filter=2pi*50*2l))...)
    r1 = PowerOptLab._pwm_ac_ripple_audit(base, real.(U), imag.(U))
    r2 = PowerOptLab._pwm_ac_ripple_audit(doubled, real.(U), imag.(U))
    @test r2.converter_rms ≈ r1.converter_rms ./ 4 rtol=2e-3
    @test r2.neutral_rms ≈ r1.neutral_rms / 4 rtol=2e-3
end

@testset "Advanced inverter: harmonic AC ripple traverses topology and LCL circuit" begin
    a = cis(2pi/3); vdc = 700.0; m = 0.3
    U = (m*vdc/sqrt(2)) .* ComplexF64[1, a^-1, a]
    common = (id="acpwm", bus="poc", phase_terminals=["a","b","c"],
              neutral="n", s_max=20e3, v_dc=vdc, c_dc=1e-3,
              i_cap_max=100.0, f_sw=10e3, pwm_strategy=:CENTERED,
              pwm_ac_ripple=true, pwm_fundamental_samples=36,
              pwm_carrier_samples=128, pwm_ac_harmonics=64)

    reduced = AdvancedInverter(; topology=:FOUR_LEG, In_max=100.0,
        x_filter=0.15, common...)
    series_lcl = AdvancedInverter(; topology=:FOUR_LEG, In_max=100.0,
        x_filter=0.05, x_filter_grid=0.10, common...)
    rr = PowerOptLab._pwm_ac_ripple_audit(reduced, real.(U), imag.(U))
    rs = PowerOptLab._pwm_ac_ripple_audit(series_lcl, real.(U), imag.(U))
    @test rs.converter_rms ≈ rr.converter_rms rtol=2e-10
    @test rs.grid_rms ≈ rr.grid_rms rtol=2e-10
    matrix_reduced = AdvancedInverter(; topology=:FOUR_LEG, In_max=100.0,
        x_filter_matrix=[0.15 0 0 0; 0 0.15 0 0;
                         0 0 0.15 0; 0 0 0 0.0], common...)
    rm = PowerOptLab._pwm_ac_ripple_audit(matrix_reduced, real.(U), imag.(U))
    @test rm.converter_rms ≈ rr.converter_rms rtol=2e-10
    @test rm.neutral_rms ≈ rr.neutral_rms rtol=2e-10

    lcl = AdvancedInverter(; topology=:FOUR_LEG, In_max=100.0,
        x_filter=0.05, x_filter_grid=0.10,
        c_filter_mid=20e-6, r_filter_damping=1.0, common...)
    rlcl = PowerOptLab._pwm_ac_ripple_audit(lcl, real.(U), imag.(U))
    @test maximum(rlcl.grid_rms) < 0.1 * maximum(rlcl.converter_rms)
    @test minimum(rlcl.shunt_rms) > 1.0

    no_ln = PowerOptLab._pwm_ac_ripple_audit(
        AdvancedInverter(; topology=:FOUR_LEG, In_max=100.0,
            x_filter=0.15, x_filter_neutral=0.0, common...), real.(U), imag.(U))
    with_ln = PowerOptLab._pwm_ac_ripple_audit(
        AdvancedInverter(; topology=:FOUR_LEG, In_max=100.0,
            x_filter=0.15, x_filter_neutral=0.30, common...), real.(U), imag.(U))
    @test with_ln.neutral_rms < 0.4 * no_ln.neutral_rms

    three = PowerOptLab._pwm_ac_ripple_audit(
        AdvancedInverter(; merge(common, (topology=:THREE_LEG,
            neutral=nothing, x_filter=0.15))...), real.(U), imag.(U))
    @test three.neutral_rms < 1e-12
    @test maximum(three.converter_rms) > 0.0
end

@testset "Advanced inverter: AC PWM reserve closes physical current ratings" begin
    common = (bus="poc", phase_terminals=["a","b","c"], neutral="n",
              topology=:FOUR_LEG, s_max=20e3, i_max=30.0, In_max=20.0,
              v_dc=700.0, c_dc=1.1e-3, i_cap_max=50.0,
              r_filter=0.05, x_filter=0.15, f_sw=10e3,
              pwm_strategy=:CENTERED, pwm_ac_ripple=true,
              pwm_fundamental_samples=24, pwm_carrier_samples=64,
              pwm_ac_harmonics=32)
    r = solve_advanced_inverter(inv_grid3_unbal(),
        AdvancedInverter(; id="i", common...); per_unit=true, s_base=100e3)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test all(r.i_ac_switching_reserved .+ 3e-2 .>= r.i_ac_switching_rms)
    @test r.i_neutral_switching_reserved + 3e-2 >= r.i_neutral_switching_rms
    @test maximum(r.i_ac_total_rms) <= 30.0 + 0.05
    @test r.i_neutral_total_rms <= 20.0 + 0.05
    @test r.pwm_iterations >= 2
    @test maximum(r.i_ac_switching_pp) > maximum(r.i_ac_switching_rms)

end

@testset "Advanced inverter: individual split half-bank ratings compose" begin
    cu, cl = 2.5e-3, 3.5e-3
    common = (bus="poc", phase_terminals=["a","b","c"], neutral="n",
              topology=:SPLIT_DC, s_max=20e3, i_max=50.0,
              v_dc=900.0, c_dc=3e-3, c_dc_upper=cu, c_dc_lower=cl,
              m_max=0.96, r_filter=0.05, x_filter=0.15)
    # A lower-half rating alone bounds neutral current and therefore satisfies
    # the topology's bounded-return-path requirement without In_max.
    lower = solve_advanced_inverter(inv_grid3_unbal(), AdvancedInverter(; id="i",
        i_cap_lower_max=3.0, common...))
    @test lower.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test lower.i_cap_thermal_lower <= 3.0 + 1e-3
    @test lower.i_cap_lower <= 3.0 + 1e-3 # equal weights by default
    @test lower.i_neutral > 0.1

    # A common rating constrains both; the larger-capacitance lower bank binds.
    both = solve_advanced_inverter(inv_grid3_unbal(), AdvancedInverter(; id="i",
        i_cap_max=3.0, common...))
    @test both.i_cap_thermal <= 3.0 + 1e-3
    @test both.i_cap_thermal_lower >= both.i_cap_thermal_upper
end

@testset "Advanced inverter: a larger bank rating cannot worsen the optimum" begin
    # A larger i_cap_max enlarges the feasible set, so the OBJECTIVE is
    # non-decreasing. That is the property actually guaranteed — the chosen
    # neutral current is NOT: the optimizer is free to re-mix zero- against
    # negative-sequence current and may serve more export with less |I_n|.
    # A zero-sequence-heavy grid makes the rating genuinely bite, so the
    # assertion is not vacuous.
    a = cis(2pi/3); vp = 230.0; v0 = 0.18*230.0
    Va = vp + v0; Vb = a^2*vp + v0; Vc = a*vp + v0
    net = inv_grid3_src(mags=[abs(Va), abs(Vb), abs(Vc)],
                        angs=[angle(Va), angle(Vb), angle(Vc)], vmax=300.0)
    split(icm) = solve_advanced_inverter(net, AdvancedInverter(; id="i",
        topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3, i_cap_max=icm, _TOPO_COMMON...))
    ps = Float64[]
    for icm in (1.5, 2.0, 3.0, 4.0)
        r = split(icm)
        @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
        @test r.i_cap <= icm + 1e-3            # the hard constraint always holds
        push!(ps, r.p_poc)
    end
    @test all(diff(ps) .> -1.0)                # non-decreasing (1 W solver slack)
    @test ps[end] > ps[1] + 100.0              # and the rating genuinely bites
    # Tightening it well below the sweep costs real export. It never renders the
    # problem infeasible — |I_n| ≤ 2·i_cap_max admits I_n = 0, so an idle inverter
    # is always feasible — the bank rating buys export, it is not a hard gate.
    tight = split(0.5)
    @test tight.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test tight.i_cap <= 0.5 + 1e-3
    @test tight.p_poc < ps[1] - 500.0
end

@testset "Advanced inverter: 4-leg caps carry the 2ω term but not the neutral" begin
    # The 4-leg's neutral current flows through its fourth leg, NOT the DC caps,
    # so a bank rating limits the ripple without costing neutral capability —
    # the physical distinction from the split link above.
    net = inv_grid3_unbal()
    four = (topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3, In_max=40.0)
    loose = solve_advanced_inverter(net, AdvancedInverter(; id="i", four..., _TOPO_COMMON...))
    tight = solve_advanced_inverter(net, AdvancedInverter(; id="i", four..., i_cap_max=1.0, _TOPO_COMMON...))
    @test tight.i_cap <= 1.0 + 1e-3               # rating binds
    @test tight.dv2 < loose.dv2 - 0.2             # ripple is what gives way …
    @test abs(tight.i_neutral - loose.i_neutral) < 0.2   # … neutral is untouched
    # Single-cap topology: the bank carries only the 2ω component, k = 2ωC/√2.
    k = 2 * 2pi * 50.0 * 1.1e-3 / sqrt(2)
    @test tight.i_cap ≈ k * tight.dv2  rtol=1e-4
end

@testset "Advanced inverter: switching polytope matches its closed form" begin
    # With a zero filter at zero power the internal EMF equals the POC voltage
    # (230 V), so the sampled polytope constrains V_dc alone and can be checked
    # against the textbook linear-modulation limits:
    #   3-leg (pairwise, balanced set): √2·√3·|U| ≤ m·V_dc  ⇒  |U| ≤ m·V_dc/√6
    #   split link (half bus, N = 0):   √2·|U| ≤ (m/2)·V_dc ⇒  |U| ≤ m·V_dc/(2√2)
    feas(topo, vdc; N=72) = begin
        kw = topo == :SPLIT_DC ? (; topology=topo, v_dc=vdc, c_dc=50e-3, In_max=40.0) :
                                 (; topology=topo, v_dc=vdc, c_dc=50e-3)
        r = solve_advanced_inverter(inv_grid3_bal(),
                AdvancedInverter(; id="i", kw..., bus="poc",
                    phase_terminals=["a","b","c"], neutral="n",
                    s_max=20e3, m_max=0.96, n_samples=N);
                objective=:min_loss, p_set=0.0)
        r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    end
    min_vdc(topo; N=72) = begin
        lo, hi = 400.0, 1200.0
        for _ in 1:16
            mid = (lo + hi)/2
            feas(topo, mid; N=N) ? (hi = mid) : (lo = mid)
        end
        hi
    end
    three, split = min_vdc(:THREE_LEG), min_vdc(:SPLIT_DC)
    @test three ≈ 230.0*sqrt(6)/0.96      rtol=2e-3   # ≈ 586.9 V
    @test split ≈ 230.0*2*sqrt(2)/0.96    rtol=2e-3   # ≈ 677.6 V
    # The split link's factor-of-two rail is the documented utilisation penalty:
    # it needs 2/√3 ≈ +15.5 % more DC voltage for the same per-phase output.
    @test split/three ≈ 2/sqrt(3)         rtol=2e-3
end

@testset "Advanced inverter: sampled feasibility is an outer approximation" begin
    # More samples ⇒ the polytope tightens (the true peak stops being missed), so
    # a bus that a coarse grid accepts must be rejected by a fine one. 587 V is
    # the converged 3-leg minimum; a 6-sample grid under-constrains down to ~508.
    coarse(vdc, N) = solve_advanced_inverter(inv_grid3_bal(),
            AdvancedInverter(; id="i", topology=:THREE_LEG, v_dc=vdc, c_dc=50e-3,
                bus="poc", phase_terminals=["a","b","c"], neutral="n",
                s_max=20e3, m_max=0.96, n_samples=N);
            objective=:min_loss, p_set=0.0)
    loose = coarse(550.0, 6)
    fine_bad = coarse(550.0, 36)
    fine_good = coarse(650.0, 36)
    @test loose.termination_status in ("LOCALLY_SOLVED", "OPTIMAL") # too loose
    @test loose.switching_margin < -1.0 # dense audit exposes the missed peak
    @test !(fine_bad.termination_status in ("LOCALLY_SOLVED", "OPTIMAL"))
    @test isnan(fine_bad.switching_margin)
    @test fine_good.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test fine_good.switching_margin > 1.0

    # The audit itself reproduces the closed-form boundary for a balanced set.
    a = cis(2pi/3); u = 230.0
    vthreshold = u*sqrt(6)/0.96
    audit_inv = AdvancedInverter(id="audit", bus="poc", topology=:THREE_LEG,
        phase_terminals=["a","b","c"], neutral="n", s_max=20e3,
        v_dc=vthreshold, c_dc=50e-3, m_max=0.96)
    margin = PowerOptLab._switching_margin(audit_inv,
        real.([u, u*a^-1, u*a^-2]), imag.([u, u*a^-1, u*a^-2]),
        0.0, 0.0, nothing, nothing)
    @test margin ≈ 0.0 atol=1e-8
end

@testset "Advanced inverter: 2ω ripple uses C_eq = C/2 for the split link" begin
    # |D| = |S̃|/(2ω·C_eq·V_dc). The split link's series caps halve C_eq, so its
    # denominator is ωC — not 2ωC. Getting this wrong is a factor-of-two error.
    net = inv_grid3_unbal(); w = 2pi*50.0
    com = (bus="poc", phase_terminals=["a","b","c"], neutral="n", s_max=20e3,
           i_max=40.0, r_filter=0.05, x_filter=0.15, m_max=0.96)
    four = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:FOUR_LEG,
               v_dc=700.0, c_dc=1.1e-3, In_max=40.0, com...))
    split = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:SPLIT_DC,
               v_dc=800.0, c_dc=2.8e-3, In_max=21.0, com...))
    @test four.dv2  ≈ four.ripple /(2w*1.1e-3*700.0)  rtol=1e-4   # single cap: 2ωC
    @test split.dv2 ≈ split.ripple/(  w*2.8e-3*800.0) rtol=1e-4   # series pair: ωC
end

@testset "Advanced inverter: 3-wire bank current is the 2ω term alone" begin
    # No zero-sequence path ⇒ no neutral share, so i_cap is purely the 2ω current.
    r = solve_advanced_inverter(inv_grid3_unbal(),
            AdvancedInverter(; id="i", topology=:THREE_LEG, v_dc=700.0, c_dc=1.1e-3,
                             _TOPO_COMMON...))
    @test r.i_neutral < 1e-3
    @test r.i_cap ≈ (2*2pi*50.0*1.1e-3/sqrt(2)) * r.dv2  rtol=1e-4
end

@testset "Advanced inverter: i_cap_max and In_max compose (whichever binds)" begin
    net = inv_grid3_src(mags=[250.0, 205.0, 230.0], angs=[0.1, -2.2, 1.95])
    sp = (topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3)
    nb = solve_advanced_inverter(net, AdvancedInverter(; id="i", sp..., In_max=4.0,  i_cap_max=20.0, _TOPO_COMMON...))
    cb = solve_advanced_inverter(net, AdvancedInverter(; id="i", sp..., In_max=40.0, i_cap_max=3.0,  _TOPO_COMMON...))
    @test nb.i_neutral ≈ 4.0  rtol=5e-2      # neutral rating governs …
    @test nb.i_cap < 20.0 - 1.0              # … bank rating slack
    @test cb.i_cap ≈ 3.0  rtol=1e-2          # bank rating governs …
    @test cb.i_neutral < 40.0 - 1.0          # … neutral rating slack
end

@testset "Advanced inverter: capacitor-rating argument validation" begin
    common = (bus="poc", phase_terminals=["a","b","c"], neutral="n", s_max=20e3)
    # i_cap_max needs a three-phase topology (it is a DC-bank rating).
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", i_cap_max=10.0, common...))
    # The switching allowance cannot exhaust (or exceed) the bank rating.
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3,
                         In_max=40.0, i_cap_max=5.0, i_sw=5.0, common...))
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3,
                         In_max=40.0, i_cap_max=-1.0, common...))
    # i_sw is an allowance carved out of i_cap_max — without it, it would be a
    # silent no-op, so it is rejected rather than ignored.
    @test_throws ArgumentError solve_advanced_inverter(inv_grid3_bal(),
        AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3,
                         In_max=40.0, i_sw=2.0, common...))
end

@testset "Advanced inverter: reactive setpoint q_set is met" begin
    # The converter delivers the requested reactive power at the POC, trading it
    # against active power on the apparent-power circle.
    r = solve_advanced_inverter(inv_grid(),
            AdvancedInverter(id="i", bus="poc", s_max=5000.0); q_set=1500.0)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.q_poc ≈ 1500.0  rtol=1e-2
    @test r.p_poc ≈ sqrt(5000.0^2 - 1500.0^2)  rtol=1e-2   # rest of the circle → P
end

@testset "Advanced inverter: per-conductor current limit i_max binds" begin
    base = solve_advanced_inverter(inv_grid(), AdvancedInverter(id="i", bus="poc", s_max=5000.0))
    lim  = solve_advanced_inverter(inv_grid(), AdvancedInverter(id="i", bus="poc", s_max=5000.0, i_max=15.0))
    @test lim.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test lim.i_mag[1] ≈ 15.0  rtol=1e-2          # current pinned at the limit
    @test lim.p_poc < base.p_poc - 100.0          # and export curtailed below nameplate
end

@testset "Advanced inverter: internal EMF cap v_int_max binds and curtails export" begin
    common = (bus="poc", s_max=5000.0, r_filter=0.2, x_filter=0.5)   # filter ⇒ V_int > V_poc
    loose = solve_advanced_inverter(inv_grid(), AdvancedInverter(; id="i", v_int_max=250.0, common...))
    tight = solve_advanced_inverter(inv_grid(), AdvancedInverter(; id="i", v_int_max=228.0, common...))
    @test tight.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test tight.v_int_mag[1] ≈ 228.0  atol=0.5    # EMF held at its cap …
    @test tight.p_poc < loose.p_poc - 50.0        # … curtailing export
end

@testset "Advanced inverter: legacy single-phase scalar cap binds (|V_int| ≤ m·Vdc/√3)" begin
    common = (bus="poc", s_max=5000.0, r_filter=0.2, x_filter=0.5)
    loose = solve_advanced_inverter(inv_grid(), AdvancedInverter(; id="i", common...))
    capd  = solve_advanced_inverter(inv_grid(), AdvancedInverter(; id="i", modulation_max=1.0, v_dc=400.0, common...))
    @test capd.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test capd.v_int_mag[1] ≈ 1.0 * 400.0 / sqrt(3)  atol=1.0   # pins the retained convention
    @test capd.p_poc < loose.p_poc - 20.0
end

@testset "Advanced inverter: POC Q includes the grid-side shunt" begin
    net = inv_grid_x()
    b = 0.05
    shun = solve_advanced_inverter(net,
        AdvancedInverter(id="i", bus="poc", s_max=5000.0, b_filter_shunt=b);
        objective=:min_loss, p_set=1000.0, q_set=0.0)
    @test shun.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    q_shunt = b * shun.bus["poc"]["1"]["vm"]^2
    @test shun.q_poc ≈ 0.0 atol=1e-3
    @test shun.q_poc - shun.q_conv ≈ q_shunt rtol=1e-5
end

@testset "Advanced inverter: iteration-limited points are not published" begin
    r = solve_advanced_inverter(inv_grid(),
        AdvancedInverter(id="i", bus="poc", s_max=5000.0);
        solver_options=(max_iter=0,))
    @test !(r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL"))
    @test isnan(r.p_poc)
    @test isnan(r.q_poc)
    @test isnan(r.bus["poc"]["1"]["vm"])
end

@testset "Advanced inverter: neutral-current limit In_max binds" begin
    net = inv_grid3_src(mags=[250.0, 205.0, 230.0], angs=[0.1, -2.2, 1.95])   # strong unbalance
    common = (bus="poc", phase_terminals=["a","b","c"], neutral="n", s_max=20e3,
              i_max=60.0, r_filter=0.05, x_filter=0.15, m_max=0.96, topology=:FOUR_LEG,
              v_dc=750.0, c_dc=1.1e-3)
    loose = solve_advanced_inverter(net, AdvancedInverter(; id="i", In_max=40.0, common...))
    tight = solve_advanced_inverter(net, AdvancedInverter(; id="i", In_max=5.0, common...))
    @test tight.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test loose.i_neutral > 8.0                   # unbalance wants a large neutral current …
    @test tight.i_neutral ≈ 5.0  rtol=5e-2        # … which the tight rating pins down
    @test tight.i_neutral < loose.i_neutral - 1.0
end

@testset "Advanced inverter: per-unit agrees with SI on well-scaled cases" begin
    # Three-phase topology with neutral current and ripple.
    net = inv_grid3_unbal()
    inv = AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3,
                            In_max=40.0, r_filter_neutral=0.03, _TOPO_COMMON...)
    si = solve_advanced_inverter(net, inv; per_unit=false)
    pu = solve_advanced_inverter(net, inv; per_unit=true)
    @test pu.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test pu.p_poc ≈ si.p_poc            rtol=1e-3
    @test pu.q_poc ≈ si.q_poc            rtol=1e-2
    @test pu.i_neutral ≈ si.i_neutral    rtol=1e-2
    @test pu.i_zero ≈ si.i_zero          rtol=1e-2
    @test pu.i_negative ≈ si.i_negative  rtol=1e-2
    @test pu.dv2 ≈ si.dv2                rtol=1e-2
    @test pu.i_cap ≈ si.i_cap            rtol=1e-2
    @test pu.switching_margin ≈ si.switching_margin rtol=1e-3
    @test pu.v_int_mag[1] ≈ si.v_int_mag[1]  rtol=1e-4

    # Split link with a binding bank rating: exercises the i_base scaling of the
    # neutral share inside the capacitor allocation, which the 4-leg above lacks.
    invs = AdvancedInverter(; id="i", topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3,
                            In_max=21.0, i_cap_max=2.5, _TOPO_COMMON...)
    csi = solve_advanced_inverter(net, invs; per_unit=false)
    cpu = solve_advanced_inverter(net, invs; per_unit=true)
    @test csi.i_cap ≈ 2.5                rtol=1e-3    # the rating binds in SI …
    @test cpu.i_cap ≈ csi.i_cap          rtol=1e-2    # … and per-unit agrees
    @test cpu.i_neutral ≈ csi.i_neutral  rtol=2e-2
    @test cpu.dv_mid ≈ csi.dv_mid        rtol=2e-2

    # Single-phase with filter + losses (guards per-unit scaling of the existing features).
    invsp = AdvancedInverter(id="i", bus="poc", s_max=5e3, r_filter=0.2, x_filter=0.5,
                             p_loss_fixed=20.0, a_loss=0.3, c_loss=0.02)
    ssi = solve_advanced_inverter(inv_grid(), invsp; objective=:min_loss, p_set=3000.0, per_unit=false)
    spu = solve_advanced_inverter(inv_grid(), invsp; objective=:min_loss, p_set=3000.0, per_unit=true)
    @test spu.p_poc ≈ ssi.p_poc          rtol=1e-4
    @test spu.p_loss ≈ ssi.p_loss        rtol=1e-3
    @test spu.p_dc ≈ ssi.p_dc            rtol=1e-4

    # Explicit LCL: guards admittance, second-arm impedance, midpoint voltage,
    # and both current bases.
    invlcl = AdvancedInverter(id="i", bus="poc", s_max=5e3,
        r_filter=0.05, x_filter=0.10,
        r_filter_grid=0.05, x_filter_grid=0.10,
        c_filter_mid=20e-6, r_filter_damping=1.0, i_grid_max=30.0)
    lsi = solve_advanced_inverter(inv_grid(), invlcl;
        objective=:min_loss, p_set=3000.0, q_set=0.0, per_unit=false)
    lpu = solve_advanced_inverter(inv_grid(), invlcl;
        objective=:min_loss, p_set=3000.0, q_set=0.0, per_unit=true)
    @test lpu.p_conv ≈ lsi.p_conv rtol=1e-4
    @test lpu.p_filter_loss ≈ lsi.p_filter_loss rtol=2e-3
    @test lpu.v_filter_mag ≈ lsi.v_filter_mag rtol=1e-4
    @test lpu.i_grid_mag ≈ lsi.i_grid_mag rtol=1e-4
    @test lpu.i_filter_shunt_mag ≈ lsi.i_filter_shunt_mag rtol=1e-4

    # Asymmetric split-link quantities stay in SI even though their constraints
    # mix AC per-unit phasors with SI capacitance, charge, ESR, and DC voltage.
    invasu = AdvancedInverter(; id="i", bus="poc",
        phase_terminals=["a","b","c"], neutral="n", topology=:SPLIT_DC,
        s_max=20e3, i_max=50.0, In_max=40.0, v_dc=900.0, c_dc=3e-3,
        c_dc_upper=2.5e-3, c_dc_lower=3.5e-3,
        q_mid_balance_max=0.25, v_mid_mean_max=1.0,
        i_cap_max=15.0, cap_thermal_weights=(2.0,1.0,0.2),
        esr_dc_upper=0.03, esr_dc_lower=0.05,
        m_max=0.96, r_filter=0.05, x_filter=0.15)
    asi = solve_advanced_inverter(net, invasu; per_unit=false)
    apu = solve_advanced_inverter(net, invasu; per_unit=true)
    @test apu.p_poc ≈ asi.p_poc rtol=1e-3
    @test apu.dv2 ≈ asi.dv2 rtol=1e-2
    @test apu.dv_mid ≈ asi.dv_mid rtol=1e-2
    @test apu.v_mid_mean ≈ asi.v_mid_mean atol=1e-3
    @test apu.q_mid_balance ≈ asi.q_mid_balance atol=1e-6
    @test apu.i_cap_thermal_upper ≈ asi.i_cap_thermal_upper rtol=1e-2
    @test apu.i_cap_thermal_lower ≈ asi.i_cap_thermal_lower rtol=1e-2
    @test apu.p_cap_loss ≈ asi.p_cap_loss rtol=1e-2
end
