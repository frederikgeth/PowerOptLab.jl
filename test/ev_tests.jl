@testset "EV charging: availability + departure energy (V1G)" begin
    # Four periods. The EV is plugged in for periods 1-3 and leaves at period 4;
    # it must reach 30 kWh by departure. Cheapest import is period 2, so a
    # charge-only (V1G) charger should fill mostly then, respect the unplugged
    # period, and meet the departure target.
    prices = [0.15, 0.05, 0.12, 0.20]
    nets = [single_bus_net(; src_cost=p, pload=10000.0) for p in prices]
    ev = EVDevice(id="ev1", bus="bus1",
                  p_charge_max=20e3,               # 20 kW charger, V1G (no discharge)
                  energy_max=40e3, energy_init=10e3,
                  available=[true, true, true, false],
                  departure_energy=30e3, departure_period=3)
    res = solve_multiperiod_opf(nets, [ev]; dt_h=1.0)

    @test res.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    d = res.dispatch["ev1"]
    @test all(d.p_discharge .<= 1.0)                 # V1G: never discharges (≈0 W)
    @test all(d.p_net .<= 1.0)                        # charger only draws ⇒ injection ≤ 0
    @test d.p_net[4] ≈ 0.0  atol=1.0                  # unplugged period is idle
    @test d.soc[4] >= 30e3 - 1.0                     # departure target met by period 3
    @test d.soc[1] ≈ 10e3  rtol=1e-6
    # Total energy charged over periods 1-3 = required rise (unit efficiency).
    @test sum(d.p_charge[1:3]) ≈ (d.soc[4] - d.soc[1])  rtol=1e-6
    # Cheapest available period (period 2) should be used at/near full power.
    @test d.p_charge[2] ≈ 20e3  rtol=1e-2
end


@testset "EV charging: horizon validation" begin
    nets = [single_bus_net(), single_bus_net()]
    short = EVDevice(id="ev", bus="bus1", p_charge_max=10e3,
        energy_max=40e3, energy_init=10e3, available=[true],
        departure_energy=20e3)
    @test_throws ArgumentError solve_multiperiod_opf(nets, [short])
    impossible = EVDevice(id="ev", bus="bus1", p_charge_max=10e3,
        energy_max=40e3, energy_init=10e3, available=[true, true],
        departure_energy=50e3)
    @test_throws ArgumentError solve_multiperiod_opf(nets, [impossible])
end

@testset "EV charging: V2G discharge into an expensive peak" begin
    # Plugged in the whole horizon, bidirectional. A cheap period then an
    # expensive one; with a modest departure target the EV can arbitrage by
    # discharging into the peak, unlike the V1G case.
    prices = [0.05, 0.25]
    nets = [single_bus_net(; src_cost=p, pload=10000.0) for p in prices]
    ev = EVDevice(id="ev2", bus="bus1",
                  p_charge_max=15e3, p_discharge_max=15e3,
                  energy_max=40e3, energy_init=20e3,
                  available=[true, true],
                  departure_energy=5e3, departure_period=2)   # low target frees energy to sell
    res = solve_multiperiod_opf(nets, [ev]; dt_h=1.0)
    @test res.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    d = res.dispatch["ev2"]
    @test d.p_net[2] > 0.0                            # discharges into the expensive period
    @test d.soc[3] >= 5e3 - 1.0                       # still meets the (low) departure floor
end
