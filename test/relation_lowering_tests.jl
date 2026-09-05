@testset "Bounds and relation-aware lowering" begin
    f = PWLFunction([220.,240.,250.,270.],[1.,1.,0.,0.])
    p = plan_pwl_relation(f,(220.,250.);relation=:upper)
    @test p.shape == :concave && p.strategy == :supporting_lines
    @test p.semantics == :exact_relation
    @test plan_pwl_relation(f,(220.,270.);relation=:upper).strategy == :unresolved
    @test plan_pwl_relation(f,(220.,250.);relation=:equal).strategy == :unresolved
    @test plan_pwl_relation(f,(240.,270.);relation=:lower).strategy == :supporting_lines
    @test plan_pwl_relation(f,(242.,248.)).strategy == :affine
    @test plan_pwl_relation(f,(220.,250.);relation=:upper,formulation=PWLConvexHull()).semantics == :exact_relation
    @test plan_pwl_relation(f,(220.,270.);relation=:upper,formulation=PWLConvexHull()).semantics == :outer_relaxation
    @test plan_pwl_relation(f,(220.,250.);relation=:upper,formulation=PWLConvexHull(),specialize=false).strategy == :graph
    @test plan_pwl_relation(f,(220.,250.);relation=:upper,formulation=PWLConvexHull(),specialize=false).semantics == :exact_relation
    @test_throws ArgumentError plan_pwl_relation(f,(1.,0.))
    @test_throws ArgumentError plan_pwl_relation(f,(1.,Inf))
    @test_throws ArgumentError plan_pwl_relation(f,(220.,250.);relation=:band)
    # No near-convexity tolerance that silently changes a subtly nonconvex curve.
    subtle = PWLFunction([0.,1.,2.],[1.,.5,eps(.5)])
    @test plan_pwl_relation(subtle,(0.,2.);relation=:upper).strategy == :unresolved
    m=Model(); @variable(m,220<=v<=250); @variable(m,y)
    h=formulate_pwl_relation!(m,f,v,y;relation=:upper)
    n=num_constraints(m;count_variable_in_set_constraints=true)
    @test length(h.rows)==3 && h.graph===nothing
    @test formulate_pwl_relation!(m,f,v,y;relation=:upper) === h
    @test num_constraints(m;count_variable_in_set_constraints=true)==n
    @variable(m,220<=v_other<=250); @variable(m,y_other)
    another=formulate_pwl_relation!(m,f,v_other,y_other;relation=:upper)
    @test another !== h
    @test another.plan === h.plan
    @test length(m.ext[:PowerOptLabPWLRelationPlans])==1
    set_upper_bound(v,270.)
    @test_throws UnsupportedFormulation formulate_pwl_relation!(m,f,v,y;relation=:upper)
    @test h.plan.domain == (220.,250.) # retained interval row survives loosening
    @test constraint_object(first(h.rows)).set == MOI.Interval(220.,250.)
    set_lower_bound(v,242.); set_upper_bound(v,248.)
    h2=formulate_pwl_relation!(m,f,v,y;relation=:upper)
    @test h2 !== h && h2.plan.strategy == :affine
    @test h2.plan.domain == (242.,248.)
    @test_throws ArgumentError formulate_pwl_relation!(m,f,1.0v,y;relation=:upper)
    fixed=Model(); @variable(fixed,x); fix(x,245.)
    @variable(fixed,z)
    @test formulate_pwl_relation!(fixed,f,x,z).plan.domain == (245.,245.)
    # Direct line inequalities are exact over the WHOLE declared domain,
    # including slack output values, not only at an optimizer's maximizer.
    for (domain,relation) in (((220.,250.),:upper),((240.,270.),:lower)),
            voltage in range(domain...;length=21), output in range(-.2,1.2;length=11)
        plan=plan_pwl_relation(f,domain;relation)
        exact=primitive_value(f,voltage;domain_policy=:flat_extension)
        violation=relation==:upper ? output-exact : exact-output
        encoded=maximum(relation==:upper ? output-(a*voltage+b) : a*voltage+b-output for (a,b) in plan.lines)
        @test encoded ≈ violation atol=1e-13
    end
end

@testset "Faithful specialized smoothing and residual direction" begin
    f=PWLFunction([220.,240.,250.,270.],[1.,1.,0.,0.])
    @test plan_pwl_relation(f,(242.,248.);formulation=LocalC2Formulation(.2)).strategy == :affine
    @test plan_pwl_relation(f,(242.,248.);formulation=SoftplusFormulation(.2)).strategy == :smooth
    @test plan_pwl_relation(f,(239.,242.);formulation=LocalC2Formulation(.2)).strategy == :local_c2
    for rep in (LocalC2Formulation(.2),SoftplusFormulation(.2)),
            domain in ((242.,248.),(239.,242.),(239.9,240.1),(251.,275.))
        m=Model(Ipopt.Optimizer); set_silent(m)
        voltage=sum(domain)/2
        @variable(m,domain[1]/230<=x<=domain[2]/230,start=voltage/230)
        @constraint(m,x==voltage/230)
        @variable(m,y)
        h=formulate_pwl_relation!(m,f,x,y;formulation=rep,input_scale=230.,output_scale=.5)
        @objective(m,Min,y)
        optimize!(m)
        @test is_solved_and_feasible(m)
        @test .5value(y) ≈ primitive_value(f,voltage,rep;domain_policy=:flat_extension) atol=1e-7
        @test audit_pwl_relation(h).surrogate_violation<1e-7
    end
    m=Model(Ipopt.Optimizer); set_silent(m)
    @variable(m,242<=v<=248,start=245.); @constraint(m,v==245.)
    @variable(m,y); @constraint(m,y==.2)
    h=formulate_pwl_relation!(m,f,v,y;relation=:upper)
    optimize!(m)
    @test audit_pwl_relation(h).canonical_residual ≈ -.3 atol=1e-7
    @test audit_pwl_relation(h).canonical_violation==0.
end
