include("thermal_energy.jl")
using Test
@testset "Mixed Australian controller/thermal/energy setup" begin
    plan=mixed_thermal_cases()
    @test length(plan)==8 && length(unique(c.case_id for c in plan))==8
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    original=deepcopy(base)
    for case in plan
        net=voltage_aware_network(base,case);mixed=mixed_thermal_fleet(net,case)
        @test length(mixed.spec.devices)==12
        @test Set(values(mixed.labels))==Set(MIXED_CONTROLLERS)
        for (id,c) in mixed.spec.devices
            inv=c.device
            @test inv.s_max==6000. && inv.i_max==10. && inv.v_dc==750. && inv.c_dc==300e-6 && inv.i_cap_max==1.
            @test inv.pwm_strategy== (case.model===:fundamental ? :NONE : :CENTERED)
            @test inv.dc_energy_model== (case.model===:energy ? :nonlinear : :linearized)
            @test net["ibr"][id]["topology"]==String(inv.topology)
        end
        @test count(==(:disconnected),values(mixed.modes))==(case.operation===:mixed_disconnection ? 4 : 0)
        @test sum(r.p_available for r in values(mixed.spec.requests)) ≈ 37000case.pv_fraction
    end
    @test base==original
end
