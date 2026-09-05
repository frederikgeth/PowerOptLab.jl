# Run in the environment produced by setup.jl. Each method gets one solve;
# any relaxation/homotopy iterations belong to CCOpt, not PowerOptLab.
using Test, PowerOptLab, JuMP, HiGHS, PiecewiseLinearOpt, MadNLP
using MathOptComplements, NLPModelsJuMP, CCOpt
include("control_lowering.jl")
methods = ControlLoweringExample.smooth_methods(MadNLP.Optimizer;solver_name="MadNLP")
append!(methods,[
    FormulationMethod("exact / HiGHS",ExactPWLGraph(),HiGHS.Optimizer;configure! = set_silent),
    FormulationMethod("hull / HiGHS",PWLConvexHull(),HiGHS.Optimizer;configure! = set_silent),
    FormulationMethod("MPCC / CCOpt",ComplementarityGraph(scale=1.),CCOpt.Optimizer;
        configure! = m -> (MathOptComplements.Bridges.add_all_bridges(m);set_silent(m)),
        options=(tol=1e-9,))])
rows = run_formulation_experiment([ControlLoweringExample.curve_case()],methods;on_error=:throw)
if haskey(ENV,"POL_FORMULATION_RESULTS")
    write_formulation_results(ENV["POL_FORMULATION_RESULTS"],rows;
        sources=[@__FILE__,joinpath(@__DIR__,"control_lowering.jl")],
        metadata=(purpose="control lowering worked example",stationarity="unassessed"))
end
@testset "Worked control lowering through external solvers" begin
    @test all(r -> r["run_status"]=="finished",rows)
    for r in rows[1:5]
        @test r["strict_solver_success"]
        @test r["candidate_available"]
    end
    for r in rows[1:3]
        @test abs(only(r["observations"]).surrogate_equation_error)<1e-7
        @test abs(only(r["observations"]).exact_graph_error)<=1.001e-3
    end
    @test rows[4]["metrics"].fraction ≈ .5 atol=1e-7
    @test rows[5]["metrics"].fraction ≈ .875 atol=1e-7
    # CCOpt is characterized: a candidate or strict success is not guaranteed.
    @test haskey(rows[6],"candidate_available")
end
for r in rows
    println(r["method"],": ",r["termination_status"],"; candidate=",r["candidate_available"],
        "; graph error=",r["candidate_available"] ? only(r["observations"]).exact_graph_error : missing)
end

include("bounded_relations.jl")
bounded_rows=BoundedRelationExample.run(HiGHS.Optimizer;options=[],graphs=true)
if haskey(ENV,"POL_FORMULATION_RESULTS")
    write_formulation_results(splitext(ENV["POL_FORMULATION_RESULTS"])[1]*"-bounded.toml",
        bounded_rows;sources=[@__FILE__,joinpath(@__DIR__,"bounded_relations.jl")])
end
@testset "Bound-dependent LP and MILP relation pathways" begin
    @test length(bounded_rows)==30
    for r in bounded_rows
        if r["run_status"]=="unsupported"
            @test r["method"]=="auto"
        else
            @test r["strict_solver_success"]
            @test r["metrics"].domain_violation<1e-6
            if r["metrics"].semantics==:exact_relation
                @test r["metrics"].canonical_violation<1e-6
            elseif r["configuration"].upper_V==270.
                @test r["metrics"].output ≈ 5/6 atol=1e-6
            end
        end
    end
    f=BoundedRelationExample.INTENT.volt_watt
    for (relation,status) in ((:upper,MOI.OPTIMAL),(:equal,MOI.INFEASIBLE))
        m=Model(HiGHS.Optimizer); set_silent(m)
        @variable(m,220<=v<=250); @constraint(m,v==245.)
        @variable(m,p); @constraint(m,p==.1)
        formulate_pwl_relation!(m,f,v,p;relation,formulation=PWLConvexHull(),specialize=false)
        optimize!(m)
        @test termination_status(m)==status
    end
end

mpcc_method=FormulationMethod("bounded MPCC / CCOpt",ComplementarityGraph(scale=1.),CCOpt.Optimizer;
    configure! = m -> (MathOptComplements.Bridges.add_all_bridges(m);set_silent(m)),options=(tol=1e-9,))
mpcc_rows=run_formulation_experiment([BoundedRelationExample.case()],[mpcc_method];
    configurations=[(lower_V=220.,upper_V=270.,voltage_V=245.,relation=:upper,specialize=true)],on_error=:throw)
if haskey(ENV,"POL_FORMULATION_RESULTS")
    write_formulation_results(splitext(ENV["POL_FORMULATION_RESULTS"])[1]*"-bounded-mpcc.toml",mpcc_rows;sources=[@__FILE__])
end
@testset "Bounded MPCC relation characterization" begin
    r=only(mpcc_rows)
    @test r["run_status"]=="finished"
    @test haskey(r,"candidate_available")
    if r["candidate_available"]
        @test haskey(r["metrics"].graph_audit,:complementarity_product)
    end
end
