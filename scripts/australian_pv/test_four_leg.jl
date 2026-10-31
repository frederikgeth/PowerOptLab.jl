include("four_leg.jl")
using Test

@testset "Australian four-leg counterfactual" begin
    plan=four_leg_cases()
    @test length(plan)==11 && length(unique(c.case_id for c in plan))==11
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    original=deepcopy(base)
    case=only(c for c in plan if c.case_id=="matched_s2_p3")
    net=voltage_aware_network(base,case);spec=four_leg_fleet(net,case.pv_fraction)
    for (id,c) in spec.devices
        inv=c.device
        @test inv.topology===:FOUR_LEG && inv.neutral=="n" && inv.In_max==10.
        @test inv.i_max==HARDWARE.i_max && inv.s_max==HARDWARE.s_max
        @test inv.c_dc==HARDWARE.c_dc && inv.v_dc==HARDWARE.v_dc && inv.i_cap_max==HARDWARE.i_cap_max
        @test inv.r_filter_neutral==inv.r_filter && inv.x_filter_neutral==inv.x_filter
        @test net["ibr"][id]["terminal_map"]==["a","b","c","n"]
    end
    r=solve_voltage_aware_inverter_fleet(net,spec;solver_options=VA_OPTIONS)
    @test r.solve.publishable
    replay=four_leg_pq_replay(net,r)
    @test replay.accepted && replay.max_difference_V<1e-5
    @test all(p.passed for p in values(physical_feasibility(r)))
    for (id,p) in r.replayed
        inv=spec.devices[id].device
        power=p.state.grid_voltage.*conj.(p.state.grid_current)
        @test maximum(real.(power))-minimum(real.(power))<.01
        @test abs(sum(p.state.converter_current))<=inv.In_max+1e-5
        @test p.p_dc<=5000.01
    end
    @test base==original
    @test !four_leg_pq_replay(net,(solve=(publishable=false,),)).accepted
    night=only(c for c in plan if c.case_id=="matched_s2_p1")
    nn=voltage_aware_network(base,night);ns=four_leg_fleet(nn,0.)
    on=solve_voltage_aware_inverter_fleet(nn,ns;solver_options=VA_OPTIONS)
    @test on.solve.publishable && four_leg_pq_replay(nn,on).accepted
    @test sum(p.p_poc for p in values(on.replayed))<0.
    off=solve_voltage_aware_inverter_fleet(nn,ns;operating_modes=:disconnected,solver_options=VA_OPTIONS)
    @test off.solve.publishable && four_leg_pq_replay(nn,off).accepted
    @test isempty(off.devices) && isempty(off.physical_reports)
end
