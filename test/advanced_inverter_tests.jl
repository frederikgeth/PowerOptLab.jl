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
end
