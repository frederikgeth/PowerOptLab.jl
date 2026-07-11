@testset "Multi-period OPF: battery arbitrage" begin
    # Two periods: expensive slack import in period 1, cheap in period 2. A cyclic
    # battery must discharge into the expensive period and recharge in the cheap
    # one, returning to its initial charge.
    nets = [single_bus_net(; src_cost=0.20), single_bus_net(; src_cost=0.05)]
    bat = StorageDevice(id="bat", bus="bus1",
                        p_charge_max=40e3, p_discharge_max=40e3,
                        energy_max=100e3, energy_init=40e3, cyclic=true)
    res = solve_multiperiod_opf(nets, [bat]; dt_h=1.0)

    @test res.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    d = res.dispatch["bat"]
    @test d.p_net[1] ≈  40e3  rtol=1e-2      # discharge in the expensive period
    @test d.p_net[2] ≈ -40e3  rtol=1e-2      # recharge in the cheap period
    @test d.soc[1]  ≈ 40e3   rtol=1e-6       # start
    @test d.soc[3]  ≈ 40e3   rtol=1e-6       # cyclic closure
    @test d.soc[2]  ≈  0.0   atol=1e2        # drained after discharging
    # SOC conservation across each step (unit efficiency).
    for t in 1:2
        @test d.soc[t+1] ≈ d.soc[t] - d.p_net[t] * 1.0  rtol=1e-6
    end
    # Per-snapshot voltages are in band.
    for t in 1:2
        @test 900.0 <= res.snapshots[t]["bus"]["bus1"]["1"]["vm"] <= 1100.0
    end
end

@testset "Multi-period OPF: round-trip efficiency loss" begin
    # With <100% efficiency, arbitrage still happens but SOC dynamics reflect the
    # one-way losses: draining 40 kWh of throughput stores less than 40 kWh.
    nets = [single_bus_net(; src_cost=0.20), single_bus_net(; src_cost=0.05)]
    bat = StorageDevice(id="bat", bus="bus1",
                        p_charge_max=40e3, p_discharge_max=40e3,
                        energy_max=100e3, energy_init=40e3,
                        eff_charge=0.95, eff_discharge=0.95, cyclic=true)
    res = solve_multiperiod_opf(nets, [bat]; dt_h=1.0)
    @test res.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    d = res.dispatch["bat"]
    # Charge SOC gain uses eff_charge; discharge SOC drain divides by eff_discharge.
    @test d.soc[2] ≈ d.soc[1] + (0.95*d.p_charge[1] - d.p_discharge[1]/0.95)*1.0  rtol=1e-6
    @test d.soc[3] ≈ d.soc[2] + (0.95*d.p_charge[2] - d.p_discharge[2]/0.95)*1.0  rtol=1e-6
    @test d.soc[3] ≈ 40e3  rtol=1e-6         # still cyclic
end
