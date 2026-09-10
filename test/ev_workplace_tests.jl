include(joinpath(@__DIR__,"..","examples","ev_workplace.jl"))

@testset "Workplace EVSE current limits, costs and deadlines" begin
    base = EVWorkplaceExample.study()
    upgrade = EVWorkplaceExample.study(weak_current=32.0)
    si = EVWorkplaceExample.study(per_unit=false)
    larger = EVWorkplaceExample.study(weak_current=32.0,visitor_gain=8e3)
    for case in (base,upgrade,si,larger)
        @test solve_status(case.result).publishable
        for s in case.sessions
            d = case.result.dispatch[s.id]
            @test maximum(abs.(d.diagnostics.energy_balance_wh)) < 0.01
            @test d.diagnostics.departure_shortfall_wh < 0.01
            @test maximum(d.diagnostics.current_limit_violation_a) < 1e-4
            @test maximum(d.diagnostics.apparent_power_violation_va) < 0.01
            @test maximum(d.diagnostics.disconnected_current_a) == 0.0
            @test all(d.p_charge[t] == 0 for t in eachindex(s.available) if !s.available[t])
        end
    end
    @test base.result.dispatch["morning"].p_charge[1] > 1000.0
    @test upgrade.result.dispatch["morning"].p_charge[1] < 1.0
    @test base.result.dispatch["afternoon"].p_charge[4] > 1000.0
    @test upgrade.result.dispatch["afternoon"].p_charge[4] < 1.0
    @test upgrade.result.objective < base.result.objective
    @test si.result.objective ≈ base.result.objective atol=1e-4
    morning=base.result.dispatch["morning"]
    @test only(morning.current_magnitude_a[2]) ≈ 16.0 atol=1e-4
    @test 200.0 < morning.p_charge[2]/only(morning.current_magnitude_a[2]) < 215.0
    # This necessary bound proves infeasibility without treating a local solver
    # status as a global infeasibility certificate.
    @test EVWorkplaceExample.visitor_energy_bound() == 7200.0 < 8000.0
    impossible = EVWorkplaceExample.study(visitor_gain=8e3)
    @test !solve_status(impossible.result).publishable
    @test all(isnan,impossible.result.dispatch["morning"].p_charge)
end
