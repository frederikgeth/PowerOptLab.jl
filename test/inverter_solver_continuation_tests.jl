@testset "Inverter continuation: rejected stages cannot seed final solve" begin
    m=Model(Ipopt.Optimizer);set_silent(m)
    @variable(m,0<=x<=1,start=.1)
    @variable(m,z);@variable(m,w)
    a=@constraint(m,z==0);b=@constraint(m,w==0)
    starts=Dict{Symbol,Any}()
    function stage!(name,weight)
        starts[name]=start_value(x)
        target=name===:neutral_initialization ? .2 : name===:neutral_release_one ? .8 : .9
        @objective(m,Max,-(x-target)^2-z^2-w^2)
        optimize!(m)
        (stage=name,raw=solve_status(m))
    end
    result=PowerOptLab._inverter_solver_continuation!(m,["one"=>[a],"two"=>[b]],(0.,);
        solve_stage! = stage!,accepted=r->r.raw.publishable && r.stage!==:neutral_release_one)
    @test all(r.raw.publishable for r in result.attempts)
    @test length(result.attempts)==3
    @test starts[:full_capability]≈.2 atol=1e-6
    @test value(x)≈.9 atol=1e-6
    @test result.temporary_constraints_remaining==0
    @test !is_valid(m,a) && !is_valid(m,b)
    @test_throws ArgumentError PowerOptLab._inverter_solver_continuation!(m,[],(.1,);
        solve_stage! = stage!,accepted=r->true)
end

@testset "Inverter continuation: exception releases temporary capabilities" begin
    m=Model(Ipopt.Optimizer);set_silent(m);@variable(m,x)
    cut=@constraint(m,x==0)
    @test_throws ErrorException PowerOptLab._inverter_solver_continuation!(m,["pv"=>[cut]],(0.,);
        solve_stage! = (name,w)->error("retained numerical failure"),accepted=r->true)
    @test !is_valid(m,cut)
end

@testset "Controller numerical policy: bounded cold starts and caller options" begin
    limit=SolveStatus("ITERATION_LIMIT","INFEASIBLE_POINT",true,false,false,false)
    timed=SolveStatus("TIME_LIMIT","INFEASIBLE_POINT",true,false,false,false)
    @test PowerOptLab._inverter_numerical_restart(limit,("max_iter"=>300,))
    @test !PowerOptLab._inverter_numerical_restart(limit,("max_iter"=>0,))
    @test !PowerOptLab._inverter_numerical_restart(timed,())
    @test_throws ArgumentError PowerOptLab._validate_inverter_initialization(:unknown)
    m=Model(Ipopt.Optimizer)
    PowerOptLab._set_inverter_solver_options!(m,();controller_graph=true)
    @test get_attribute(m,"nlp_scaling_method")=="none"
    PowerOptLab._set_inverter_solver_options!(m,("nlp_scaling_method"=>"gradient-based",);controller_graph=true)
    @test get_attribute(m,"nlp_scaling_method")=="gradient-based"
end
