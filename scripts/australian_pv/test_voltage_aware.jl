# Standalone scientific setup/accounting regressions; no solver runs required.
include("voltage_aware.jl")
using Test

@testset "Voltage-aware Australian experiment plan and energy accounting" begin
    snapshots=voltage_aware_cases(:snapshots)
    chronology=voltage_aware_cases(:chronology)
    @test length(snapshots)==41
    @test length(chronology)==36
    @test all(sum(x.duration_h for x in chronology if x.group==g && x.controller==m)==24
        for g in ("energized_night","disconnected_night") for m in TERMINAL_MODES)
    @test all(x.source_negative_angle_rad in (pi/2,pi,3pi/2) for x in snapshots if x.group=="source_angle")
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    original=deepcopy(base)
    for case in snapshots
        net=voltage_aware_network(base,case)
        for source in values(net["voltage_source"])
            v=source["v_magnitude"].*cis.(source["v_angle"]);a=cis(2pi/3)
            u1=(v[1]+a*v[2]+a^2*v[3])/3;u2=(v[1]+a^2*v[2]+a*v[3])/3
            @test u2/u1 ≈ case.source_vuf*cis(case.source_negative_angle_rad) atol=1e-12
            @test abs(sum(v))<1e-8
        end
        @test sum(sum(ld["p_nom"]) for ld in values(net["load"])) ≈ 18000case.load_scale
    end
    @test base==original
    @test all(!r.accepted for r in voltage_aware_energy([],[],chronology))
    # Independent analytic fixture: one 5 kW PV, 80% DC capture, 100 W active
    # AC-model losses (70 W converter + 30 W filter), signed Q=-300 var.
    solved=[merge(c,(accepted=true,)) for c in chronology]
    devices=[let active=c.operating_mode===:energized,available=5000c.pv_fraction
        (case_id=c.case_id,model_present=active,p_available_W=available,
            p_poc_W=active ? .8available-100 : 0.,p_dc_W=active ? .8available : nothing,
            converter_loss_W=active ? 70. : nothing,filter_loss_W=active ? 30. : nothing,
            q_poc_var=active ? -300. : 0.)
        end for c in chronology]
    energy=voltage_aware_energy(solved,devices,chronology)
    @test length(energy)==6 && all(r.accepted for r in energy)
    for row in energy
        on=row.group=="energized_night"
        @test row.available_kWh ≈ 33.75
        @test row.unused_dc_kWh ≈ 6.75
        @test row.export_kWh ≈ 25.8
        @test row.import_kWh ≈ (on ? 1.2 : 0.)
        @test row.net_poc_kWh ≈ (on ? 24.6 : 25.8)
        @test row.converter_loss_kWh ≈ (on ? 1.68 : .84)
        @test row.filter_loss_kWh ≈ (on ? .72 : .36)
        @test row.reactive_support_kvarh ≈ (on ? -7.2 : -3.6)
        @test abs(row.balance_residual_Wh)<1e-9
    end
    failed=copy(solved);failed[1]=merge(failed[1],(accepted=false,))
    withheld=voltage_aware_energy(failed,devices,chronology)
    @test !withheld[1].accepted && withheld[1].net_poc_kWh===nothing
    @test count(r->r.accepted,withheld)==5
end
