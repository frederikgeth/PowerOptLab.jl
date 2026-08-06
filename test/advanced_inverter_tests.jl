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
    big = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=2.0e-3, In_max=40.0, _TOPO_COMMON...))
    small = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=0.3e-3, In_max=40.0, _TOPO_COMMON...))
    @test small.dv2 > big.dv2 + 1.0        # smaller capacitor ⇒ larger bus ripple
    capped = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=0.3e-3, In_max=40.0, dv2_max=3.0, _TOPO_COMMON...))
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
    # Zero sequence: the mirror image — neutral current, negligible ripple.
    @test zers.i_neutral > 5.0
    @test zers.ripple < negs.ripple/10
    # And the two superpose rather than trading against each other.
    @test both.i_neutral > 5.0
    @test both.ripple > 3000.0
    # The bank current tracks the ripple (4-leg: neutral bypasses the caps).
    @test negs.i_cap ≈ negs.ripple/(sqrt(2)*700.0)  rtol=1e-3
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
            objective=:min_loss, p_set=0.0).termination_status
    @test coarse(550.0, 6)  in ("LOCALLY_SOLVED", "OPTIMAL")        # too loose
    @test !(coarse(550.0, 36) in ("LOCALLY_SOLVED", "OPTIMAL"))     # correctly rejected
    @test coarse(650.0, 36) in ("LOCALLY_SOLVED", "OPTIMAL")        # genuinely feasible
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

@testset "Advanced inverter: single-phase modulation cap binds (|V_int| ≤ m·Vdc/√3)" begin
    common = (bus="poc", s_max=5000.0, r_filter=0.2, x_filter=0.5)
    loose = solve_advanced_inverter(inv_grid(), AdvancedInverter(; id="i", common...))
    capd  = solve_advanced_inverter(inv_grid(), AdvancedInverter(; id="i", modulation_max=1.0, v_dc=400.0, common...))
    @test capd.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test capd.v_int_mag[1] ≈ 1.0 * 400.0 / sqrt(3)  atol=1.0   # pinned at m·Vdc/√3 ≈ 231 V
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

@testset "Advanced inverter: per-unit matches SI" begin
    # Three-phase topology with neutral current and ripple.
    net = inv_grid3_unbal()
    inv = AdvancedInverter(; id="i", topology=:FOUR_LEG, v_dc=700.0, c_dc=1.1e-3, In_max=40.0, _TOPO_COMMON...)
    si = solve_advanced_inverter(net, inv; per_unit=false)
    pu = solve_advanced_inverter(net, inv; per_unit=true)
    @test pu.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test pu.p_poc ≈ si.p_poc            rtol=1e-3
    @test pu.q_poc ≈ si.q_poc            rtol=1e-2
    @test pu.i_neutral ≈ si.i_neutral    rtol=1e-2
    @test pu.dv2 ≈ si.dv2                rtol=1e-2
    @test pu.i_cap ≈ si.i_cap            rtol=1e-2
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

    # Single-phase with filter + losses (guards per-unit scaling of the existing features).
    invsp = AdvancedInverter(id="i", bus="poc", s_max=5e3, r_filter=0.2, x_filter=0.5,
                             p_loss_fixed=20.0, a_loss=0.3, c_loss=0.02)
    ssi = solve_advanced_inverter(inv_grid(), invsp; objective=:min_loss, p_set=3000.0, per_unit=false)
    spu = solve_advanced_inverter(inv_grid(), invsp; objective=:min_loss, p_set=3000.0, per_unit=true)
    @test spu.p_poc ≈ ssi.p_poc          rtol=1e-4
    @test spu.p_loss ≈ ssi.p_loss        rtol=1e-3
    @test spu.p_dc ≈ ssi.p_dc            rtol=1e-4
end
