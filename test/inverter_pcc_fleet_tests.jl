# Reproduce the reported terminal-profile failures on the actual Australian
# feeder. Include its builders in an isolated namespace without running or
# rewriting any study campaign. Positive-only laws need no I₂ gain design.
module AustralianPVTerminalFixture
    include(joinpath(@__DIR__,"..","scripts","australian_pv","run.jl"))
end

@testset "PCC fleet: matched reactive-priority voltage controls" begin
    study=AustralianPVTerminalFixture
    base=study.parse_bmopf(joinpath(study.ROOT,"networks","LV3_55bus.bmopf.json"))
    options=study.TERMINAL_SOLVER_OPTIONS
    modes=study.TERMINAL_MODES
    cases=[(mode,source,level,1e5,:loss) for mode in modes
        for source in (1.,1.03) for level in (0.,.5,1.)]
    append!(cases,[(mode,1.03,1.,1e5,:zero) for mode in modes])
    reference=Dict{Tuple{String,Float64,Float64},Any}()
    for (mode,source,level,s_base,objective) in cases
        @testset "$mode source=$source PV=$level base=$s_base objective=$objective" begin
            net=study.study_network(base,level,source)
            spec=study.fleet(net,mode,level;settings=study.TERMINAL_CONTROL_SETTINGS,
                design=(kappa=0.,phi2=0.))
            r=solve_controlled_inverter_fleet(net,spec;s_base=s_base,
                selection_objective=objective,solver_options=options)
            @test solve_status(r).publishable
            if solve_status(r).publishable
                rows=controlled_inverter_rows(r)
                @test length(rows)==12
                @test study.accepted_control_rows(rows,study.HARDWARE)
                @test maximum(x.exact_smooth_current_residual_A for x in rows)<.005
                @test all(abs(d.plant.q_poc-d.control.q_request)<.01 for d in values(r.devices))
                # The numerical selection objective must not select a different
                # dispatch for this fully controlled power-flow comparison.
                key=(mode,source,level)
                if !haskey(reference,key)
                    reference[key]=r
                else
                    ref=reference[key]
                    @test maximum(abs(d.plant.p_poc-ref.devices[id].plant.p_poc)
                        for (id,d) in r.devices)<.01
                    @test maximum(abs(d.plant.q_poc-ref.devices[id].plant.q_poc)
                        for (id,d) in r.devices)<.01
                    @test maximum(abs(d.control.phase_voltage[k]-
                        ref.devices[id].control.phase_voltage[k])
                        for (id,d) in r.devices for k in 1:3)<1e-4
                end
            end
        end
    end
end
