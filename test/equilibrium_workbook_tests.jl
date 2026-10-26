include(joinpath(@__DIR__,"..","scripts","formulations","equilibrium_workbook.jl"))

@testset "Equilibrium workbook: fixed points and guarded branch graphs" begin
    p = wb_parameters()
    @test wb_equilibrium_oracle(p) ≈ 5/11
    direct = wb_iteration(p;steps=20)
    @test direct.power == [iseven(i) ? 0. : 1. for i in 0:20]
    @test abs(last(direct.residual)) > .9
    relaxed = wb_iteration(p;relaxation=.5,steps=20)
    @test last(relaxed.power) ≈ 5/11 atol=1e-9
    @test abs(last(relaxed.residual)) < 1e-9
    for pp in (p,wb_parameters(source_v=.8),wb_parameters(source_v=1.2),
               wb_parameters(gain=0.),wb_parameters(available=0.),
               wb_parameters(resistance=0.),wb_parameters(source_v=.9),
               wb_parameters(source_v=1+1/12))
        result = wb_solve(pp)
        @test result.selected !== nothing
        @test length(result.solutions) == 1
        s = result.selected
        @test s.power ≈ wb_equilibrium_oracle(pp) atol=1e-7
        @test s.voltage ≈ pp.source_v+pp.resistance*s.power atol=1e-8
        @test s.power ≈ clamp(pp.available-pp.gain*(s.voltage-pp.reference_v),0.,pp.available) atol=1e-7
        report = wb_report(pp,result)
        @test maximum(abs(x.residual) for x in report.equalities) < 1e-7
        @test minimum(x.margin for x in report.margins) >= -1e-7
        for c in result.candidates
            if !c.solved
                @test isnan(c.power) && isnan(c.voltage)
            end
        end
    end
    controller = wb_solve(p)
    @test only(filter(c -> c.solved,controller.candidates)).branch === :slope
    # An unguarded full-output assignment looks executable but is not an
    # equilibrium of this controller: its assumed branch contradicts its input.
    false_voltage = p.source_v+p.resistance*p.available
    @test wb_raw(p,false_voltage) < p.available
    @test wb_command(p,false_voltage) != p.available
end

@testset "Equilibrium workbook: a ceiling and a command have different promises" begin
    for p in (wb_parameters(),wb_parameters(source_v=.8),wb_parameters(source_v=1.2),
              wb_parameters(gain=0.),wb_parameters(available=0.))
        oracle = wb_equilibrium_oracle(p)
        for (objective,target,expected) in ((:max_p,.2,oracle),(:min_p,.2,0.),
                                           (:tracking,.2,clamp(.2,0.,oracle)))
            result = wb_solve(p;relationship=:capability,objective,target)
            @test result.selected !== nothing
            @test result.selected.power ≈ expected atol=1e-6
            report = wb_report(p,result)
            @test only(report.equalities).name === :network
            @test minimum(x.margin for x in report.margins) >= -1e-7
            @test result.selected.power <= wb_command(p,result.selected.voltage)+1e-7
        end
    end
    p = wb_parameters()
    strict = wb_solve(p;voltage_max=1.03)
    flexible = wb_solve(p;relationship=:capability,objective=:max_p,voltage_max=1.03)
    @test strict.selected === nothing
    @test wb_report(p,strict) === nothing
    @test flexible.selected.power ≈ .3 atol=1e-7
    @test flexible.selected.voltage ≈ 1.03 atol=1e-8
    report = wb_report(p,flexible)
    @test report.command_p ≈ .64 atol=1e-7
    @test report.controller_gap ≈ -.34 atol=1e-7
    @test :voltage_upper in report.binding
    @test !(:command_ceiling in report.binding)
    @test only(report.equalities).residual ≈ 0. atol=1e-8
    # Even a feasible-set model fails when the ceiling is below the source V.
    @test wb_solve(p;relationship=:capability,voltage_max=.99).selected === nothing
end

@testset "Equilibrium workbook: physical time and memory survive optimisation" begin
    for (limit,power,energy,revenue) in ((1.03,[-1.,1.,0.],[1.,2.,1.,1.],4.),
                                       (1.015,[-1.,.5,.5],[1.,2.,1.5,1.],3.))
        r = wb_battery_schedule(voltage_max=limit)
        @test r.solved
        @test r.power ≈ power atol=1e-7
        @test r.energy ≈ energy atol=1e-7
        @test r.replay_energy ≈ r.energy atol=1e-8
        @test r.revenue ≈ revenue atol=1e-7
        @test r.voltage ≈ 1 .+ .03r.power atol=1e-8
        @test maximum(r.voltage) <= limit+1e-8
        @test sum(r.duration .* r.power) ≈ 0. atol=1e-8
        @test r.energy[1] ≈ 1. atol=1e-8
        @test r.energy[end] ≈ 1. atol=1e-8
    end
    # A ceiling below the source voltage would require charging in every
    # period, contradicting the same initial/final energy. Do not publish it.
    failed = wb_battery_schedule(voltage_max=.99)
    @test !failed.solved
    @test all(isnan,failed.power)
    @test all(isnan,failed.energy)
    @test isnan(failed.revenue)
    measurements = [1.,1.03,1.,.97,1.]
    @test wb_hysteresis(measurements;initial=false) == [false,true,true,false,false]
    @test wb_hysteresis(measurements;initial=true) == [true,true,true,false,false]
    @test wb_hysteresis([.98,1.02];initial=true) == [false,true]
    @test wb_hysteresis([1.,1.];initial=false) != wb_hysteresis([1.,1.];initial=true)
    @test isempty(wb_hysteresis(Float64[]))
end

@testset "Equilibrium workbook: malformed teaching inputs fail clearly" begin
    for bad in (NaN,Inf,-Inf)
        @test_throws ArgumentError wb_parameters(gain=bad)
        @test_throws ArgumentError wb_iteration(wb_parameters();relaxation=bad)
        @test_throws ArgumentError wb_solve(wb_parameters();target=bad)
        @test_throws ArgumentError wb_battery_schedule(voltage_max=bad)
        @test_throws ArgumentError wb_hysteresis([bad])
    end
    @test_throws ArgumentError wb_parameters(resistance=-1.)
    @test_throws ArgumentError wb_parameters(available=-1.)
    @test_throws ArgumentError wb_parameters(source_v=0.)
    @test_throws ArgumentError wb_iteration(wb_parameters();relaxation=0.)
    @test_throws ArgumentError wb_iteration(wb_parameters();steps=-1)
    @test_throws ArgumentError wb_iteration(wb_parameters();initial=2.)
    @test_throws ArgumentError wb_solve(wb_parameters();relationship=:unknown)
    @test_throws ArgumentError wb_solve(wb_parameters();objective=:unknown)
    @test_throws ArgumentError wb_hysteresis([1.];lower=1.02,upper=.98)
    @test_throws ArgumentError wb_report(wb_parameters(),wb_solve(wb_parameters());tolerance=0.)
end
