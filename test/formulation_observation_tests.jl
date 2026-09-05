using TOML

@testset "Mixed graph and relation experiment evidence" begin
    f=PWLFunction([0.,1.,2.],[0.,1.,0.];input_unit=:V,output_unit=:A)
    case=FormulationCase("slack concave limit",(rep,c)->begin
        m=Model()
        @variable(m,0<=x<=1,start=.5)
        @constraint(m,x==.5) # Physical input 1 V; retain the declared interval.
        @variable(m,y,start=.01)
        @constraint(m,y==.01) # Physical output .1 A, below the 1 A limit.
        h=formulate_pwl_relation!(m,f,x,y;relation=:upper,formulation=rep,
            input_scale=2.,output_scale=10.,specialize=false)
        if h.graph!==nothing
            # Choose a feasible lifted hull point that is off the canonical graph.
            @constraint(m,h.graph.output==.08)
        end
        g=formulate_pwl!(m,f,x,LocalC2Formulation(.1);input_scale=2.,output_scale=10.)
        @objective(m,Min,0.)
        (model=m,observations=[g,h])
    end)
    methods=[FormulationMethod(string(rep),rep,Ipopt.Optimizer;
        options=(tol=1e-10,),configure! = set_silent)
        for rep in (PWLConvexHull(),LocalC2Formulation(.1))]
    rows=run_formulation_experiment([case],methods;on_error=:throw)
    for row in rows
        @test row["strict_solver_success"]
        g,h=row["observations"]
        @test g.observation_kind==:graph
        @test abs(g.surrogate_equation_error)<1e-7
        @test h.observation_kind==:relation && h.relation==:upper
        @test h.canonical_residual ≈ -.9 atol=1e-7
        @test h.canonical_violation==0.
        @test h.domain==(0.,2.) && h.input_scale==2. && h.output_scale==10.
        @test h.curve.input_unit==:V && h.curve.output_unit==:A
    end
    hull=rows[1]["observations"][2]
    @test hull.semantics==:exact_relation
    @test abs(hull.graph_audit.exact_graph_error) ≈ .2 atol=1e-7
    @test hull.approximation_contract===nothing
    smooth=rows[2]["observations"][2]
    @test smooth.semantics==:smooth_surrogate
    @test smooth.approximation_contract.width==.1
    @test smooth.surrogate_violation==0.
    mktempdir() do dir
        file=write_formulation_results(joinpath(dir,"mixed.toml"),rows)
        saved=TOML.parsefile(file)["runs"][1]["observations"][2]
        @test saved["semantics"]=="exact_relation"
        @test saved["curve"]["input_unit"]=="V"
        @test saved["canonical_violation"]==0.
        @test saved["graph_audit"]["exact_graph_error"] ≈ -.2 atol=1e-7
    end

    # An old candidate from a different model must not become evidence for this run.
    foreign=Model(); @variable(foreign,z)
    h=formulate_pwl!(foreign,f,z,LocalC2Formulation(.1))
    wrong=FormulationCase("foreign observation",(rep,c)->begin
        m=Model(); @variable(m,w,start=1.); @objective(m,Min,(w-1)^2)
        (model=m,observations=[h])
    end)
    row=only(run_formulation_experiment([wrong],methods[1:1]))
    @test row["run_status"]=="error"
    @test occursin("belongs to a different model",row["error"])
end
