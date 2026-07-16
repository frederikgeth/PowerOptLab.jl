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
    rs = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3, In_max=24.0, _TOPO_COMMON...))
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
    split_lo = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:SPLIT_DC, v_dc=650.0, c_dc=2.8e-3, In_max=24.0, _TOPO_COMMON...))
    split_hi = solve_advanced_inverter(net, AdvancedInverter(; id="i", topology=:SPLIT_DC, v_dc=800.0, c_dc=2.8e-3, In_max=24.0, _TOPO_COMMON...))
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
    @test pu.v_int_mag[1] ≈ si.v_int_mag[1]  rtol=1e-4

    # Single-phase with filter + losses (guards per-unit scaling of the existing features).
    invsp = AdvancedInverter(id="i", bus="poc", s_max=5e3, r_filter=0.2, x_filter=0.5,
                             p_loss_fixed=20.0, a_loss=0.3, c_loss=0.02)
    ssi = solve_advanced_inverter(inv_grid(), invsp; objective=:min_loss, p_set=3000.0, per_unit=false)
    spu = solve_advanced_inverter(inv_grid(), invsp; objective=:min_loss, p_set=3000.0, per_unit=true)
    @test spu.p_poc ≈ ssi.p_poc          rtol=1e-4
    @test spu.p_loss ≈ ssi.p_loss        rtol=1e-3
    @test spu.p_dc ≈ ssi.p_dc            rtol=1e-4
end
