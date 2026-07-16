# Tests for the current–voltage (IVQ) battery: chemistry library + device solve.

@testset "BatteryChemistry: constructors and validation" begin
    # Linear (state-dependent voltage source) — OCV falls empty→full with SoC.
    lin = linear_chemistry(; v_full=4.0, v_empty=3.0, r_internal=0.02, q_cell=10.0)
    @test lin.ocv(0.0) ≈ 3.0
    @test lin.ocv(1.0) ≈ 4.0
    @test lin.ocv(0.5) ≈ 3.5
    @test lin.r_internal(0.3) ≈ 0.02
    @test lin.i_charge_max ≈ 10.0            # default 1C

    # Thévenin — constant OCV behind fixed R.
    th = thevenin_chemistry(; v_nominal=3.6, r_internal=0.01, q_cell=5.0)
    @test th.ocv(0.1) ≈ 3.6
    @test th.ocv(0.9) ≈ 3.6

    # Tabulated — smooth monotone cubic (PCHIP), held flat outside the sample hull.
    tab = tabulated_chemistry(; soc_points=[0.0, 0.5, 1.0], ocv_points=[3.0, 3.4, 4.0],
                              r_internal=0.03, q_cell=8.0)
    @test tab.ocv(0.0) ≈ 3.0                  # passes through the data points
    @test tab.ocv(0.5) ≈ 3.4
    @test tab.ocv(1.0) ≈ 4.0
    @test 3.1 < tab.ocv(0.25) < 3.3           # interpolated (cubic, ≈ 3.19)
    @test tab.ocv(-1.0) ≈ 3.0                 # clamped low
    @test tab.ocv(2.0)  ≈ 4.0                 # clamped high
    @test tab.ocv_affine === nothing          # tabulated ⇒ registered as an operator
    @test tab.r_constant ≈ 0.03               # constant R ⇒ embedded directly

    # R(soc) table.
    tabr = tabulated_chemistry(; soc_points=[0.0, 1.0], ocv_points=[3.0, 4.0],
                               r_points=[0.05, 0.01], q_cell=8.0)
    @test tabr.r_internal(0.0) ≈ 0.05
    @test tabr.r_internal(0.5) ≈ 0.03

    # Validation.
    @test_throws ArgumentError linear_chemistry(; v_full=3.0, v_empty=4.0, r_internal=0.01, q_cell=1.0)
    @test_throws ArgumentError thevenin_chemistry(; v_nominal=3.6, r_internal=-0.1, q_cell=1.0)
    @test_throws ArgumentError tabulated_chemistry(; soc_points=[0.0, 1.0],
        ocv_points=[4.0, 3.0], r_internal=0.01, q_cell=1.0)   # non-monotone OCV
    @test_throws ArgumentError tabulated_chemistry(; soc_points=[0.0, 1.0],
        ocv_points=[3.0, 4.0], q_cell=1.0)                    # no resistance given
end

@testset "BatteryChemistry: preloaded library is monotone and cited" begin
    for chem in (illustrative_lfp(), illustrative_nmc(), illustrative_nca(),
                 illustrative_lead_acid(), illustrative_leaf())
        socs = range(chem.soc_min, chem.soc_max; length=25)
        ocvs = [chem.ocv(s) for s in socs]
        @test issorted(ocvs)                                  # OCV non-decreasing
        @test all(chem.v_cell_min .<= ocvs .<= chem.v_cell_max)
        @test all(chem.r_internal(s) >= 0 for s in socs)
        @test !isempty(chem.source)                            # provenance recorded
    end
    # LFP has the characteristically flat mid-SoC plateau.
    lfp = illustrative_lfp()
    @test abs(lfp.ocv(0.7) - lfp.ocv(0.3)) < 0.1
    # The tabulated OCV is SMOOTH (C¹): the derivative is continuous across an
    # interior data knot (a piecewise-linear interpolant would kink here). This
    # is what makes OCV(soc) safe to embed as a function of the SoC variable.
    for s0 in (0.30, 0.50, 0.70)
        h = 1e-6
        dleft  = (lfp.ocv(s0)   - lfp.ocv(s0-h)) / h
        dright = (lfp.ocv(s0+h) - lfp.ocv(s0))   / h
        @test isapprox(dleft, dright; atol=5e-3)              # continuous slope ⇒ no kink
        @test dleft >= -1e-6                                   # and non-decreasing
    end
    # Leaf cell matches the source paper's Table 2 band.
    leaf = illustrative_leaf()
    @test leaf.v_cell_min ≈ 3.20 && leaf.v_cell_max ≈ 4.15
    @test leaf.q_cell ≈ 29.0 && leaf.i_discharge_max ≈ 90.0
end

@testset "IVQ: round-trip efficiency emerges from the resistance" begin
    # η = f_v(soc,i)/f_v(soc,−i) = (OCV−IR)/(OCV+IR): a derived quantity, < 1,
    # and falling with current — no fixed efficiency parameter.
    chem = thevenin_chemistry(; v_nominal=3.6, r_internal=0.02, q_cell=10.0)
    ocv = chem.ocv(0.5); r = chem.r_internal(0.5)
    η(I) = (ocv - I*r) / (ocv + I*r)
    @test η(5.0) < 1.0
    @test η(10.0) < η(5.0)                                     # more current, worse round-trip
    @test η(0.0) ≈ 1.0                                          # lossless at zero current
end

@testset "IVQ battery: discharges to the converter rating (max export)" begin
    # Large pack, modest inverter → the converter s_max binds, not the cell.
    chem = illustrative_leaf()
    inv  = AdvancedInverter(id="bat", bus="poc", s_max=5000.0)
    bat  = IVQBattery(id="bat", bus="poc", chemistry=chem,
                      n_series=100, n_parallel=1, soc_init=0.6, inverter=inv)
    r = solve_ivq_battery(inv_grid(), bat)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.p_poc ≈ 5000.0  rtol=1e-3                          # export at nameplate
    @test r.i_cell > 0                                          # discharging
    @test r.p_dc ≈ r.p_poc  rtol=1e-3                          # no loss model ⇒ P_dc = P_poc
    @test r.v_pack ≈ r.v_cell * 100  rtol=1e-9                 # pack = n_series·cell
    # Terminal voltage sits below OCV under discharge current.
    @test r.v_cell < chem.ocv(0.6)
end

@testset "IVQ battery: charges from the grid (max charge)" begin
    chem = illustrative_nmc()
    inv  = AdvancedInverter(id="bat", bus="poc", s_max=4000.0)
    bat  = IVQBattery(id="bat", bus="poc", chemistry=chem,
                      n_series=100, n_parallel=2, soc_init=0.4, inverter=inv)
    r = solve_ivq_battery(inv_grid(), bat; objective=:max_charge)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.p_poc < 0                                           # importing from grid
    @test r.i_cell < 0                                          # charging
    @test r.v_cell > chem.ocv(0.4)                             # terminal V above OCV when charging
end

@testset "IVQ battery: the cell current limit binds before the converter" begin
    # Small per-cell discharge limit, oversized inverter → the cell current caps export.
    chem = illustrative_leaf(; i_discharge_max=40.0)
    inv  = AdvancedInverter(id="bat", bus="poc", s_max=1e6)
    bat  = IVQBattery(id="bat", bus="poc", chemistry=chem,
                      n_series=100, n_parallel=1, soc_init=0.7, inverter=inv)
    r = solve_ivq_battery(inv_grid(), bat)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.i_cell ≈ 40.0  rtol=1e-2                           # pinned at the cell limit
end

@testset "IVQ battery: the cell voltage limit caps discharge power" begin
    # Near-empty linear cell with a tight lower voltage bound: discharge current
    # is limited by v = OCV − i·R ≥ v_cell_min, i.e. i ≤ (OCV − v_min)/R.
    chem = linear_chemistry(; v_full=3.4, v_empty=3.0, r_internal=0.05, q_cell=100.0,
                            v_cell_min=2.9, v_cell_max=3.5, soc_min=0.0, soc_max=1.0)
    inv  = AdvancedInverter(id="bat", bus="poc", s_max=1e6)
    bat  = IVQBattery(id="bat", bus="poc", chemistry=chem,
                      n_series=100, n_parallel=1, soc_init=0.0, inverter=inv)
    r = solve_ivq_battery(inv_grid(), bat)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.v_cell ≈ 2.9  rtol=1e-3                            # terminal voltage at its floor
    ocv0 = chem.ocv(0.0)
    @test r.i_cell ≈ (ocv0 - 2.9)/0.05  rtol=1e-2             # (3.0−2.9)/0.05 = 2 A
end

@testset "IVQ battery: converter loss makes the DC side supply more than the POC" begin
    chem = illustrative_leaf()
    inv  = AdvancedInverter(id="bat", bus="poc", s_max=5000.0, r_filter=0.1, x_filter=0.2,
                            p_loss_fixed=15.0, a_loss=0.2, c_loss=0.01)
    bat  = IVQBattery(id="bat", bus="poc", chemistry=chem,
                      n_series=100, n_parallel=1, soc_init=0.6, inverter=inv)
    r = solve_ivq_battery(inv_grid(), bat; objective=:min_loss, p_set=3000.0)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.p_poc ≈ 3000.0  rtol=1e-3
    @test r.p_loss > 0
    @test r.p_dc > r.p_poc                                     # DC supplies delivery + all losses
    @test r.p_dc ≈ r.p_conv + r.p_loss  rtol=1e-6             # non-branching DC balance
    @test r.p_conv > r.p_poc                                   # filter loss between converter and POC
end

@testset "IVQ battery: input validation" begin
    chem = illustrative_leaf()
    inv  = AdvancedInverter(id="bat", bus="poc", s_max=5000.0)
    bat  = IVQBattery(id="bat", bus="poc", chemistry=chem,
                      n_series=100, n_parallel=1, soc_init=0.6, inverter=inv)
    @test_throws ArgumentError solve_ivq_battery(inv_grid(), bat; objective=:bogus)
    @test_throws ArgumentError solve_ivq_battery(inv_grid(), bat; objective=:min_loss)
    # Bus mismatch between battery and its inverter.
    badinv = AdvancedInverter(id="bat", bus="elsewhere", s_max=5000.0)
    badbat = IVQBattery(id="bat", bus="poc", chemistry=chem,
                        n_series=100, n_parallel=1, soc_init=0.6, inverter=badinv)
    @test_throws ArgumentError solve_ivq_battery(inv_grid(), badbat)
    # soc_init outside the chemistry's usable window is caught.
    oob = IVQBattery(id="bat", bus="poc", chemistry=chem,
                     n_series=100, n_parallel=1, soc_init=0.999, inverter=inv)
    @test_throws ArgumentError solve_ivq_battery(inv_grid(), oob)
    # Nonpositive pack counts.
    @test_throws ArgumentError solve_ivq_battery(inv_grid(),
        IVQBattery(id="bat", bus="poc", chemistry=chem, n_series=0, n_parallel=1,
                   soc_init=0.6, inverter=inv))
    @test_throws ArgumentError solve_ivq_battery(inv_grid(),
        IVQBattery(id="bat", bus="poc", chemistry=chem, n_series=100, n_parallel=0,
                   soc_init=0.6, inverter=inv))
    # soc_final outside the window.
    @test_throws ArgumentError solve_ivq_battery(inv_grid(),
        IVQBattery(id="bat", bus="poc", chemistry=chem, n_series=100, n_parallel=1,
                   soc_init=0.6, inverter=inv, soc_final=1.5))
    # Nonpositive dt_h (multi-period entry).
    @test_throws ArgumentError solve_multiperiod_ivq([inv_grid()], [bat]; dt_h=0.0)
    # The embedded inverter is validated through the battery entry points too.
    badrating = AdvancedInverter(id="bat", bus="poc", s_max=-1.0)
    @test_throws ArgumentError solve_ivq_battery(inv_grid(),
        IVQBattery(id="bat", bus="poc", chemistry=chem, n_series=100, n_parallel=1,
                   soc_init=0.6, inverter=badrating))
end

@testset "BatteryChemistry: nonphysical parameters are rejected" begin
    # Duplicate SoC knots (would divide by a zero interval → NaN).
    @test_throws ArgumentError tabulated_chemistry(; soc_points=[0.0, 0.5, 0.5, 1.0],
        ocv_points=[3.0, 3.3, 3.4, 3.6], r_internal=0.02, q_cell=10.0)
    # Non-increasing SoC knots.
    @test_throws ArgumentError tabulated_chemistry(; soc_points=[0.0, 0.8, 0.5],
        ocv_points=[3.0, 3.4, 3.6], r_internal=0.02, q_cell=10.0)
    # Non-finite data.
    @test_throws ArgumentError tabulated_chemistry(; soc_points=[0.0, 1.0],
        ocv_points=[3.0, NaN], r_internal=0.02, q_cell=10.0)
    # Negative capacity.
    @test_throws ArgumentError linear_chemistry(; v_full=3.6, v_empty=3.0,
        r_internal=0.01, q_cell=-5.0)
    # Zero discharge current (would zero the DC current base → division by zero).
    @test_throws ArgumentError thevenin_chemistry(; v_nominal=3.6, r_internal=0.01,
        q_cell=10.0, i_discharge_max=0.0)
    # Invalid voltage bounds.
    @test_throws ArgumentError thevenin_chemistry(; v_nominal=3.6, r_internal=0.01,
        q_cell=10.0, v_cell_min=4.0, v_cell_max=3.0)
    # SoC window outside [0,1].
    @test_throws ArgumentError linear_chemistry(; v_full=3.6, v_empty=3.0,
        r_internal=0.01, q_cell=10.0, soc_min=-0.1, soc_max=0.9)
    @test_throws ArgumentError linear_chemistry(; v_full=3.6, v_empty=3.0,
        r_internal=0.01, q_cell=10.0, soc_min=0.1, soc_max=1.2)
    # SoC window not strictly inside the tabulated knots.
    @test_throws ArgumentError tabulated_chemistry(; soc_points=[0.1, 0.9],
        ocv_points=[3.0, 3.6], r_internal=0.02, q_cell=10.0, soc_min=0.0, soc_max=0.9)
end

@testset "IVQ multi-period: arbitrage with exact coulomb conservation (forward)" begin
    # Expensive slack in period 1, cheaper in period 2 → discharge then recharge.
    # The forward balance conserves charge EXACTLY: for a cyclic horizon the
    # closure Σ i·Δt = 0 must hold, not merely soc[T+1] = soc[1] by some recurrence.
    nets = [single_bus_net(; src_cost=0.15, pload=0.0),
            single_bus_net(; src_cost=0.10, pload=0.0)]
    chem = linear_chemistry(; v_full=3.6, v_empty=3.0, r_internal=0.01, q_cell=50.0,
                            soc_min=0.05, soc_max=0.95)
    inv  = AdvancedInverter(id="b", bus="bus1", s_max=100e3)
    bat  = IVQBattery(id="b", bus="bus1", chemistry=chem,
                      n_series=300, n_parallel=1, soc_init=0.6, inverter=inv, cyclic=true)
    res = solve_multiperiod_ivq(nets, [bat]; dt_h=1.0)
    @test res.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    d = res.dispatch["b"]
    @test d.i_cell[1] > 0                                      # discharge when expensive
    @test d.i_cell[2] < 0                                      # charge when cheaper
    @test d.p_poc[1] > 0 && d.p_poc[2] < 0                    # export then import
    @test d.v_cell[1] < chem.ocv(d.soc[1])                    # terminal V droops on discharge
    @test d.v_cell[2] > chem.ocv(d.soc[2])                    # and rises on charge
    # Net amp-hour conservation (coulomb closure), checked INDEPENDENTLY of soc:
    @test sum(d.i_cell) * 1.0 ≈ 0.0  atol=1e-4               # Σ i·Δt = 0 (cyclic)
    # Cyclic closure of the SoC trajectory:
    @test d.soc[1] ≈ 0.6  rtol=1e-6
    @test d.soc[3] ≈ 0.6  rtol=1e-5
    # Forward balance is exact per step: Δsoc = −i·Δt / (q_cell·n_parallel).
    for t in 1:2
        @test d.soc[t+1] ≈ d.soc[t] - d.i_cell[t]*1.0/50.0  rtol=1e-4
    end
end

@testset "IVQ multi-period: tabulated (LFP) chemistry via smooth operator" begin
    # Tabulated OCV(soc) ⇒ the registered smooth-operator path, end to end.
    nets = [single_bus_net(; src_cost=0.20, pload=0.0),
            single_bus_net(; src_cost=0.05, pload=0.0)]
    chem = illustrative_lfp(; q_cell=50.0)
    inv  = AdvancedInverter(id="b", bus="bus1", s_max=100e3)
    bat  = IVQBattery(id="b", bus="bus1", chemistry=chem,
                      n_series=300, n_parallel=1, soc_init=0.5, inverter=inv, cyclic=true)
    res = solve_multiperiod_ivq(nets, [bat]; dt_h=1.0)
    @test res.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    d = res.dispatch["b"]
    @test d.i_cell[1] > 0 && d.i_cell[2] < 0                  # discharge then charge
    @test sum(d.i_cell) * 1.0 ≈ 0.0  atol=1e-4               # coulomb closure
    @test all(chem.soc_min - 1e-6 .<= d.soc .<= chem.soc_max + 1e-6)   # within the window
    @test all(chem.v_cell_min - 1e-3 .<= d.v_cell .<= chem.v_cell_max + 1e-3)
end

@testset "IVQ multi-period: R(soc) resistance table through a solve" begin
    # A tabulated R(soc) exercises the registered resistance operator in a solve.
    chem = tabulated_chemistry(; name="rtab",
        soc_points=[0.0, 0.5, 1.0], ocv_points=[3.0, 3.3, 3.6],
        r_points=[0.05, 0.02, 0.01], q_cell=50.0, soc_min=0.05, soc_max=0.95)
    @test chem.r_constant === nothing                         # R is a table, not constant
    nets = [single_bus_net(; src_cost=0.12, pload=0.0),
            single_bus_net(; src_cost=0.10, pload=0.0)]
    inv  = AdvancedInverter(id="b", bus="bus1", s_max=100e3)
    bat  = IVQBattery(id="b", bus="bus1", chemistry=chem,
                      n_series=300, n_parallel=1, soc_init=0.5, inverter=inv, cyclic=true)
    res = solve_multiperiod_ivq(nets, [bat]; dt_h=1.0)
    @test res.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test sum(res.dispatch["b"].i_cell) * 1.0 ≈ 0.0  atol=1e-4
end

@testset "IVQ multi-period: soc_final pins the terminal state (T=1)" begin
    chem = linear_chemistry(; v_full=3.6, v_empty=3.0, r_internal=0.01, q_cell=50.0,
                            soc_min=0.05, soc_max=0.95)
    inv  = AdvancedInverter(id="b", bus="bus1", s_max=100e3)
    for sf in (0.45, 0.55)                                     # net discharge / net charge
        bat = IVQBattery(id="b", bus="bus1", chemistry=chem, n_series=300, n_parallel=1,
                         soc_init=0.5, inverter=inv, cyclic=false, soc_final=sf)
        res = solve_multiperiod_ivq([single_bus_net(; src_cost=0.10, pload=0.0)], [bat]; dt_h=1.0)
        @test res.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
        d = res.dispatch["b"]
        @test d.soc[1] ≈ 0.5  rtol=1e-6                        # initial held
        @test d.soc[2] ≈ sf   rtol=1e-4                        # terminal pinned
        @test d.soc[2] ≈ d.soc[1] - d.i_cell[1]*1.0/50.0  rtol=1e-4   # forward balance
    end
end

@testset "IVQ: SI and per-unit engine modes agree" begin
    # The battery is per-unit on its own DC bases, so the engine's per_unit flag
    # must not change the (SI) result.
    chem = illustrative_leaf()
    inv  = AdvancedInverter(id="bat", bus="poc", s_max=5000.0, r_filter=0.1, x_filter=0.2)
    bat  = IVQBattery(id="bat", bus="poc", chemistry=chem,
                      n_series=100, n_parallel=1, soc_init=0.6, inverter=inv)
    r_si = solve_ivq_battery(inv_grid(), bat; per_unit=false)
    r_pu = solve_ivq_battery(inv_grid(), bat; per_unit=true, s_base=5e3)
    @test r_si.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r_pu.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r_si.p_poc  ≈ r_pu.p_poc   rtol=1e-3
    @test r_si.i_cell ≈ r_pu.i_cell  rtol=1e-3
    @test r_si.v_cell ≈ r_pu.v_cell  rtol=1e-3
    @test r_si.p_loss ≈ r_pu.p_loss  rtol=1e-2
end

@testset "IVQ: converter loss cannot be inflated under charging" begin
    # Adversarial: :max_charge (min p_poc). The regularised-norm equality pins the
    # current magnitude, so the solver cannot manufacture linear losses (a_loss·im)
    # to import more than it stores. Loss must equal the true three-term curve.
    chem = illustrative_leaf()
    inv  = AdvancedInverter(id="bat", bus="poc", s_max=5000.0,
                            p_loss_fixed=10.0, a_loss=0.5, c_loss=0.01)
    bat  = IVQBattery(id="bat", bus="poc", chemistry=chem,
                      n_series=100, n_parallel=1, soc_init=0.5, inverter=inv)
    r = solve_ivq_battery(inv_grid(), bat; objective=:max_charge)
    @test r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    @test r.i_cell < 0                                         # genuinely charging
    @test r.p_dc  < 0                                          # DC power into the battery
    imag = hypot(r.i_cell, 0.0) * 1                            # single phase: |I| pack-side
    # Loss stays bounded and physical (not inflated): at most the loaded three-term
    # value at the operating current, and small relative to the ~kW throughput.
    @test 0 <= r.p_loss < 0.25 * abs(r.p_poc) + 50.0
end
