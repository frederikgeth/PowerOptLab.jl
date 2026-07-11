using BMOPFTools: solve_pf

@testset "State estimation: exact recovery from noiseless measurements" begin
    # Ground truth from a determined power flow, then estimate from perfect
    # measurements — the WLS solution must reproduce the true state.
    pf = solve_pf(two_bus_net(; load=true); per_unit=false)
    true_vm = Dict(b => hypot(pf["bus"][b]["1"]["vr"], pf["bus"][b]["1"]["vi"])
                   for b in ("bus1","bus2"))
    true_pinj = Dict("bus1" => -20000.0, "bus2" => -20000.0)   # loads draw nominal

    meas = Measurement[]
    for b in ("bus1","bus2")
        push!(meas, Measurement(kind=:vmag, bus=b, value=true_vm[b], sigma=2.0))
        push!(meas, Measurement(kind=:pinj, bus=b, value=true_pinj[b], sigma=400.0))
        push!(meas, Measurement(kind=:qinj, bus=b, value=0.0, sigma=400.0))
    end

    se = solve_state_estimation(two_bus_net(; load=false), meas)
    @test se.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    for b in ("bus1","bus2")
        @test se.bus[b]["1"]["vm"] ≈ true_vm[b]  atol=1e-2
    end
    # Residuals essentially zero for perfect data.
    @test all(abs(r.residual) < 1e-1 for r in se.residuals if r.kind == :vmag)
    @test se.objective < 1e-3
end

@testset "State estimation: noise reduction vs raw measurements" begin
    # A fixed, deterministic perturbation. Redundant measurements (6 for 4 voltage
    # unknowns) let the fused estimate beat the raw voltage readings.
    pf = solve_pf(two_bus_net(; load=true); per_unit=false)
    true_vm = Dict(b => hypot(pf["bus"][b]["1"]["vr"], pf["bus"][b]["1"]["vi"])
                   for b in ("bus1","bus2"))
    # Deterministic perturbations (volts / watts).
    dv = Dict("bus1" => 1.5, "bus2" => -2.5)
    dp = Dict("bus1" => -300.0, "bus2" => 350.0)

    meas = Measurement[]
    for b in ("bus1","bus2")
        push!(meas, Measurement(kind=:vmag, bus=b, value=true_vm[b] + dv[b], sigma=2.0))
        push!(meas, Measurement(kind=:pinj, bus=b, value=-20000.0 + dp[b], sigma=400.0))
        push!(meas, Measurement(kind=:qinj, bus=b, value=0.0, sigma=400.0))
    end

    se = solve_state_estimation(two_bus_net(; load=false), meas)
    @test se.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")

    raw_rms = sqrt(sum(dv[b]^2 for b in ("bus1","bus2")) / 2)
    est_rms = sqrt(sum((se.bus[b]["1"]["vm"] - true_vm[b])^2 for b in ("bus1","bus2")) / 2)
    @test est_rms < raw_rms                                   # fusion filtered noise
    for b in ("bus1","bus2")
        @test abs(se.bus[b]["1"]["vm"] - true_vm[b]) <= 3*2.0 # within 3σ of truth
    end
    # Residual bookkeeping is consistent: measured − estimated, normalized by σ.
    for r in se.residuals
        @test r.residual ≈ r.measured - r.estimated  atol=1e-6
        @test r.normalized ≈ r.residual / (r.kind == :vmag ? 2.0 : 400.0)  rtol=1e-6
    end
end
