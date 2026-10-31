# Independent acceptance/scaling checks for the Australian capability reference.
include("coordinated.jl")
using Test

@testset "Australian coordinated capability reference" begin
    @test length(coordinated_cases())==24
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    net=study_network(base,1.,1.03);original=deepcopy(net)
    spec=fleet(net,"mean_vv_vw",1.;settings=TERMINAL_CONTROL_SETTINGS,design=(kappa=0.,phi2=0.))
    results=Dict{Tuple,Any}()
    for allocation in (:balanced_converter,:free),cap in (nothing,253.)
        r=solve_australian_coordinated(net,spec;current_allocation=allocation,voltage_max=cap)
        @test r.accepted && r.solve.publishable
        @test r.replay.accepted && r.replay.max_difference_V<1e-5
        @test all(c.passed for c in r.network_checks)
        @test all(d.plant.physical_feasibility.passed for d in values(r.devices))
        @test isempty(r.network["ibr"])
        if cap!==nothing
            @test length(r.voltage_checks)==36 && all(c.passed for c in r.voltage_checks)
            @test maximum(c.value for c in r.voltage_checks) ≈ 253. atol=1e-4
        end
        for d in values(r.devices)
            @test d.plant.p_dc<=5000.01
            @test d.plant.p_dc-d.plant.p_poc ≈ d.plant.p_loss+d.plant.p_cap_loss+d.plant.p_filter_loss atol=.001
            allocation===:balanced_converter && @test d.plant.i_negative<1e-5
        end
        results[(allocation,cap)]=r
    end
    total(r)=sum(d.plant.p_poc for d in values(r.devices))
    @test total(results[(:free,253.)])>total(results[(:balanced_converter,253.)])+2000.
    @test total(results[(:free,nothing)])>=total(results[(:free,253.)])-.01
    @test total(results[(:balanced_converter,nothing)])>=total(results[(:balanced_converter,253.)])-.01
    for base in (1e4,1e6),allocation in (:balanced_converter,:free)
        r=solve_australian_coordinated(net,spec;current_allocation=allocation,s_base=base)
        @test r.accepted
        @test total(r) ≈ total(results[(allocation,253.)]) atol=.1
    end
    coarse=("tol"=>1e3,"constr_viol_tol"=>1e3,"dual_inf_tol"=>1e3,"compl_inf_tol"=>1e3)
    rejected=solve_australian_coordinated(net,spec;solver_options=coarse)
    @test rejected.raw.optimal && !rejected.accepted && !rejected.solve.publishable
    @test any(!c.passed for c in rejected.network_checks)
    @test rejected.replay.termination=="NOT_RUN"
    @test_throws ArgumentError solve_australian_coordinated(net,spec;current_allocation=:balanced_grid)
    @test_throws ArgumentError solve_australian_coordinated(net,spec;voltage_max=-1.)
    @test_throws ArgumentError solve_australian_coordinated(net,spec;s_base=0.)
    @test net==original
end
