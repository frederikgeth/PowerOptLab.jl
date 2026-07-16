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

@testset "Multi-period OPF: duration and device contracts" begin
    nets = [single_bus_net(; src_cost=0.20), single_bus_net(; src_cost=0.05)]
    bat = StorageDevice(id="bat", bus="bus1", p_charge_max=40e3,
        p_discharge_max=40e3, energy_max=100e3, energy_init=40e3)
    @test_throws ArgumentError solve_multiperiod_opf(nets, [bat]; dt_h=0.0)
    @test_throws ArgumentError solve_multiperiod_opf(nets, [bat]; dt_h=Inf)
    @test_throws ArgumentError TimeGrid(Float64[])
    @test_throws ArgumentError TimeGrid([1.0, 0.0])
    @test_throws ArgumentError solve_multiperiod_opf(nets, [bat];
        time_grid=TimeGrid([1.0]))
    @test_throws ArgumentError solve_multiperiod_opf(nets,
        [StorageDevice(id="bad", bus="bus1", p_charge_max=-1.0,
            p_discharge_max=1.0, energy_max=10.0, energy_init=5.0)])
    @test_throws ArgumentError solve_multiperiod_opf(nets,
        [StorageDevice(id="bad", bus="missing", p_charge_max=1.0,
            p_discharge_max=1.0, energy_max=10.0, energy_init=5.0)])
    @test_throws ArgumentError solve_multiperiod_opf(nets,
        [StorageDevice(id="bad", bus="bus1", p_charge_max=1.0,
            p_discharge_max=1.0, energy_max=10.0, energy_init=11.0)])
    @test_throws ArgumentError solve_multiperiod_opf(nets,
        [StorageDevice(id="bad", bus="bus1", p_charge_max=1.0,
            p_discharge_max=1.0, energy_max=10.0, energy_init=5.0,
            eff_charge=0.0)])

    # With no inter-temporal device, duration only converts each cost rate to
    # an interval cost and must scale the reported objective exactly.
    hourly = solve_multiperiod_opf(nets, StorageDevice[]; dt_h=1.0)
    half_hourly = solve_multiperiod_opf(nets, StorageDevice[]; dt_h=0.5)
    @test half_hourly.objective ≈ 0.5 * hourly.objective rtol=1e-8

    # Nonuniform periods weight each rate separately rather than applying one
    # scalar to the whole horizon.
    p1 = solve_multiperiod_opf([nets[1]], StorageDevice[]).objective
    p2 = solve_multiperiod_opf([nets[2]], StorageDevice[]).objective
    nonuniform = solve_multiperiod_opf(nets, StorageDevice[];
        time_grid=TimeGrid([0.25, 0.75]))
    @test nonuniform.objective ≈ 0.25p1 + 0.75p2 rtol=1e-8

    multi = build_multi_context([nets[1]]; per_unit=false)
    @test length(multi) == 1
    @test multi[1].model === multi.model

    @test bat isa AbstractDevice
    @test validate_device(bat, nets; periods=2) === nothing
    @test solve_status(hourly).publishable
    @test solve_diagnostics(hourly).periods == 2
end

@testset "Multi-period OPF: nonuniform SOC integration" begin
    nets = [single_bus_net(; src_cost=0.20), single_bus_net(; src_cost=0.05)]
    grid = TimeGrid([0.5, 1.5])
    bat = StorageDevice(id="bat", bus="bus1",
        p_charge_max=40e3, p_discharge_max=40e3,
        energy_max=100e3, energy_init=40e3, cyclic=true)
    result = solve_multiperiod_opf(nets, [bat]; time_grid=grid)
    dispatch = result.dispatch["bat"]
    for t in eachindex(grid.durations_h)
        @test dispatch.soc[t+1] ≈
              dispatch.soc[t] - dispatch.p_net[t] * grid[t] rtol=1e-6
    end
    @test dispatch.soc[end] ≈ dispatch.soc[1] rtol=1e-6
end

@testset "Multi-period OPF: iteration-limited points are not published" begin
    nets = [single_bus_net(; src_cost=0.20), single_bus_net(; src_cost=0.05)]
    r = solve_multiperiod_opf(nets, StorageDevice[];
        solver_options=(max_iter=0,))
    @test !(r.termination_status in ("LOCALLY_SOLVED", "OPTIMAL"))
    @test !solve_status(r).publishable
    @test isnan(r.objective)
    @test all(isnan(r.snapshots[t]["bus"]["bus1"]["1"]["vm"]) for t in eachindex(nets))
end
