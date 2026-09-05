@testset "Magnitude direction and reproducible upstream provenance" begin
    B=PowerOptLab.BMOPFTools
    for (rating,ib) in ((40.,40/0.0092),(37.,1234.),(nothing,230.))
        old=PowerOptLab.build_opf_model(inv_grid3_bal();add_objective=false)
        new=PowerOptLab.build_opf_model(inv_grid3_bal();add_objective=false)
        scale=something(rating,1.)/ib
        B.smooth_norm(old,.1,.2;scale,eps_rel=PowerOptLab._magnitude_epsilon(rating,ib)/scale,name="loss")
        PowerOptLab._push_magnitude(new,nothing,.1,.2,rating,ib;name="loss")
        a,b=only(values(B.opf_differentiability_annotations(old))),only(values(B.opf_differentiability_annotations(new)))
        @test a.name==b.name
        @test a.metadata==b.metadata
        @test a.description==b.description
        @test B.opf_research_hashes(old)["differentiability_annotations_sha256"] ==
            B.opf_research_hashes(new)["differentiability_annotations_sha256"]
    end
    for n in (2,3)
        ctx=PowerOptLab.build_opf_model(inv_grid3_bal();add_objective=false)
        m=PowerOptLab._opf_model(ctx)
        @variable(m,x,start=.3); @constraint(m,x==.3)
        rep=MagnitudeApproximation(.5;direction=:upper,unit=:A,scale=20.)
        expr=magnitude_expression(ctx,[x;zeros(n-1)],rep;component_scale=10.,output_scale=2.,name="upper")
        a=only(values(B.opf_differentiability_annotations(ctx)))
        @test a.owner==:PowerOptLab && a.metadata["form"]=="upper_2norm"
        @test occursin("overestimates",a.description) && !occursin("underestimates",a.description)
        @test a.metadata["scale"]==2. && a.metadata["eps_rel"]==.025
        @objective(m,Min,expr); PowerOptLab.enforce_kcl!(ctx); set_silent(m); optimize!(m)
        @test value(expr)*2 ≈ hypot(3.,.5) atol=1e-7
    end
    # A relative override is honored even on a plain model (no annotation).
    # This admissible one-ulp difference detects silently ignoring the keyword.
    plain=Model(); upper=MagnitudeApproximation(1.;direction=:upper)
    @test magnitude_expression(plain,(0.,0.),upper)==1.
    @test magnitude_expression(plain,(0.,0.),upper;eps_rel=nextfloat(1.))==nextfloat(1.)
    @test_throws ArgumentError magnitude_expression(plain,(0.,0.),upper;eps_rel=2.)
    @test_throws ArgumentError MagnitudeApproximation(.1;scale=0.)
    ctx=PowerOptLab.build_opf_model(inv_grid3_bal();add_objective=false)
    @test_throws ArgumentError magnitude_expression(ctx,(1.,2.),MagnitudeApproximation(.1);eps_rel=1.)
end

@testset "Prepared immutable controller curves and curvature contracts" begin
    xs,ys=[220.,240.,250.,270.],[1.,1.,0.,0.]
    law=PiecewiseLinearLaw(xs,ys;smoothing_epsilon=.1)
    xs[1]=0.;ys[1]=0.
    @test law.breakpoints[1]==220. && law.values[1]==1.
    @test_throws Exception setindex!(law.values,.5,1)
    @test PWLFunction(law) === law.curve
    for family in (SoftplusFormulation,LocalC2Formulation,AlgebraicFormulation)
        r=smoothing_for_error(law.curve,family,.01)
        c=formulation_contract(law.curve,r)
        @test c.second_derivative_bound ≈ .2hinge_contract(r).second_derivative_bound
        for x in (219.,240.,245.,250.,275.)
            l=PiecewiseLinearLaw(collect(law.breakpoints),collect(law.values);formulation=r)
            @test PowerOptLab._pwl_numeric_smooth(l,x) ≈ primitive_value(l.curve,x,r;domain_policy=:flat_extension)
            @test abs(primitive_derivatives(l.curve,x,r;domain_policy=:flat_extension)[2])<=c.second_derivative_bound+1e-12
        end
    end
    @test formulation_contract(PWLFunction([0.,1.],[2.,2.]),SoftplusFormulation(1e-320)).second_derivative_bound==0.
    # A rejected foreign variable must not stamp domain rows, even without reuse.
    m=Model(); other=Model(); @variable(other,x)
    intent=VoltVarWattIntent(volt_watt=law.curve)
    for reuse in (false,true)
        @test_throws ArgumentError formulate_control_curve!(m,intent,:volt_watt,x,
            VoltVarWattEncoding(volt_watt=LocalC2Formulation(.1));domain=(220.,270.),reuse)
        @test num_variables(m)==0 && num_constraints(m;count_variable_in_set_constraints=true)==0
    end
end

@testset "Graph construction coverage and degenerate references" begin
    f=PWLFunction([220.,240.,250.,270.],[1.,1.,0.,0.])
    @test Base.get_extension(PowerOptLab,:PowerOptLabPiecewiseLinearOptExt)!==nothing
    for r in (ExactPWLGraph(),ComplementarityGraph())
        m=Model(); @variable(m,242<=x<=260); @variable(m,y)
        h=formulate_pwl_relation!(m,f,x,y;formulation=r,specialize=false)
        @test h.plan.strategy==:graph
        @test h.graph.curve.breakpoints==(242.,250.,260.)
        @test collect(h.graph.curve.values) ≈ [.8,0.,0.]
        @test_throws ArgumentError audit_pwl_relation(h)
        r isa ComplementarityGraph && (@test length(h.graph.pairs)==1)
    end
    f=PWLFunction([0.,1.,2.,3.],[0.,1.,2.,1.])
    reference=resistive_equilibria(f,0.,1.)
    @test isempty(reference.points) # The root at 2 belongs to a continuum.
    @test length(reference.intervals)==2
end

struct ZeroAtKinkHinge <: AbstractPWLSmoothing end
PowerOptLab.hinge_contract(::ZeroAtKinkHinge)=(error_lower=0.,error_upper=1.,width=1.,regularity=:C2,second_derivative_bound=1.)
PowerOptLab.hinge_value(x,::ZeroAtKinkHinge)=max(x,0.)
@testset "Reject invalid nonnegative minimum contract" begin
    @test_throws ArgumentError selector_value(0.,1.,ZeroAtKinkHinge();kind=:nonnegative_min)
end

@testset "Inspectable automatic policy and conservative relations" begin
    for (ys,shape) in (([0.,1.,2.],:affine),([0.,1.,1.],:concave),([0.,0.,1.],:convex),([0.,1.,0.,1.],:neither)),
            relation in (:equal,:upper,:lower)
        f=PWLFunction(collect(0.:length(ys)-1),ys)
        p=plan_pwl_relation(f,(first(f.breakpoints),last(f.breakpoints));relation)
        expected=shape==:affine ? :affine :
            (shape==:concave && relation==:upper || shape==:convex && relation==:lower) ? :supporting_lines : :unresolved
        @test p.strategy==expected
        @test p.reason_code==(expected==:unresolved ? :explicit_encoding_required : :exact_polyhedral)
        @test occursin(string(expected),sprint(show,p))
        @test !occursin("PWLFunction",sprint(show,p))
    end
    f=PWLFunction([220.,240.,250.,270.],[1.,1.,0.,0.])
    for rep in (LocalC2Formulation(.2),SoftplusFormulation(.2),AlgebraicFormulation(.2)),
            relation in (:upper,:lower), domain in ((239.,242.),(242.,248.))
        m=Model(Ipopt.Optimizer); set_silent(m); set_optimizer_attribute(m,"tol",1e-10); set_optimizer_attribute(m,"bound_relax_factor",0.)
        voltage=sum(domain)/2
        @variable(m,domain[1]/230<=x<=domain[2]/230,start=voltage/230)
        @constraint(m,x==voltage/230); @variable(m,y,start=.5)
        h=formulate_pwl_relation!(m,f,x,y;formulation=rep,relation,conservative=true,input_scale=230.,output_scale=2.)
        @objective(m,Min,relation==:upper ? -y : y); optimize!(m)
        @test is_solved_and_feasible(m)
        a=audit_pwl_relation(h)
        @test a.canonical_violation<1e-7 && a.surrogate_violation<1e-7
        @test a.semantics==:inner_approximation && a.conservative
        contract=h.plan.approximation_contract
        @test PowerOptLab._formulation_observation(h).approximation_contract==contract
        @test a.output_shift==-(relation==:upper ? contract.error_upper : contract.error_lower)
        @test a.output ≈ primitive_value(f,voltage,rep)+a.output_shift atol=1e-7
        @test formulate_pwl_relation!(m,f,x,y;formulation=rep,relation,conservative=true,input_scale=230.,output_scale=2.)===h
    end
    # Compact patches incur no error on the zero tail; global softplus bounds
    # can still make that same nonnegative-dispatch model infeasible.
    for rep in (LocalC2Formulation(.2),SoftplusFormulation(.2))
        m=Model(Ipopt.Optimizer); set_silent(m)
        @variable(m,252<=x<=258,start=255.); @constraint(m,x==255.)
        @variable(m,y>=0,start=0.)
        h=formulate_pwl_relation!(m,f,x,y;relation=:upper,formulation=rep,conservative=true)
        optimize!(m)
        if rep isa LocalC2Formulation
            @test h.plan.output_shift==0.
            @test is_solved_and_feasible(m)
            @test abs(value(y))<1e-7
        else
            @test h.plan.output_shift<0.
            @test termination_status(m)==MOI.LOCALLY_INFEASIBLE
        end
    end
    for (relation,rep) in ((:equal,SoftplusFormulation(.1)),(:upper,PWLConvexHull()),(:upper,:auto))
        @test_throws ArgumentError plan_pwl_relation(f,(220.,270.);relation,formulation=rep,conservative=true)
    end
    # Requested hulls and the actual supporting rows are distinct export fields.
    case=FormulationCase("specialized",(rep,c)->begin
        m=Model(); @variable(m,220<=x<=250,start=245.); @constraint(m,x==245.)
        @variable(m,y); @objective(m,Max,y)
        h=formulate_pwl_relation!(m,f,x,y;formulation=rep,relation=:upper)
        (model=m,observations=[h])
    end)
    row=only(run_formulation_experiment([case],[FormulationMethod("hull",PWLConvexHull(),Ipopt.Optimizer;configure! = set_silent)];on_error=:throw))
    a=only(row["observations"])
    @test a.built_formulation==:linear_inequalities
    @test endswith(a.requested_formulation,"PWLConvexHull")
    @test a.reason_code==:exact_polyhedral
    @test !hasproperty(a,:formulation_type)
end

@testset "Solver-free physical graph hull gap bounds" begin
    f=PWLFunction([220.,240.,250.,270.],[20.,20.,0.,0.];input_unit=:V,output_unit=:A)
    g=hull_gap_bound(f)
    @test g.upper_gap>=40/3 && g.lower_gap>=40/3
    @test g.upper_gap ≈ 40/3 && g.lower_gap ≈ 40/3
    @test Rational{BigInt}(g.upper_gap)>=40//3 && Rational{BigInt}(g.lower_gap)>=40//3
    @test g.upper_witness.input==250. && g.lower_witness.input==240.
    @test g.input_unit==:V && g.output_unit==:A
    @test hull_gap_bound(f,(242.,248.)).upper_gap==0.
    @test hull_gap_bound(f,(245.,245.)).lower_gap==0.
    @test hull_gap_bound(f,(200.,210.)).upper_gap==0.
    @test_throws ArgumentError hull_gap_bound(f,(250.,220.))
    # Independent oracle: enumerate every chord at each restricted knot. The
    # extreme vertical slice of a 2D graph hull uses at most two vertices.
    for ys in ([0.,2.,0.,1.],[2.,2.,1.,0.],[0.,.1,.9,1.],[0.,1e-12,0.,0.]), domain in ((.2,2.7),(-1.,4.))
        f=PWLFunction([0.,1.,2.,3.],ys)
        g=hull_gap_bound(f,domain)
        xs=[domain[1];[x for x in f.breakpoints if domain[1]<x<domain[2]];domain[2]]
        vals=[primitive_value(f,x;domain_policy=:flat_extension) for x in xs]
        up,lo=0.,0.
        for (x,y) in zip(xs,vals)
            slice=[vals[a]+(vals[b]-vals[a])*(x-xs[a])/(xs[b]-xs[a])
                for a in eachindex(xs) for b in a+1:length(xs) if xs[a]<=x<=xs[b]]
            up=max(up,maximum(slice)-y);lo=max(lo,y-minimum(slice))
        end
        @test g.upper_gap ≈ up atol=1e-15
        @test g.lower_gap ≈ lo atol=1e-15
    end
end


@testset "Domain-certified compact-patch conservative shifts" begin
    f=PWLFunction([220.,240.,250.,270.],[1.,1.,0.,0.])
    r=LocalC2Formulation(.2)
    # Independent analytic bounds: the slope changes at 240 and 250 are
    # -0.1 and +0.1; each active quartic hinge contributes at most 3δ/16.
    B=.1*3*.2/16
    domains=(((242.,248.),0.,0.), ((239.,242.),-B,0.),
        ((249.,252.),0.,B), ((220.,270.),-B,B),
        ((240.2,249.8),0.,0.), ((200.,210.),0.,0.),
        ((252.,258.),0.,0.), ((240.,240.),-B,0.), ((250.,250.),0.,B))
    for (domain,lower,upper) in domains, specialize in (false,true), relation in (:upper,:lower)
        p=plan_pwl_relation(f,domain;relation,formulation=r,specialize,conservative=true)
        c=p.approximation_contract
        @test c.error_scope==:domain && c.domain==domain
        @test c.error_lower ≈ lower atol=1e-16
        @test c.error_upper ≈ upper atol=1e-16
        @test p.output_shift ≈ -(relation==:upper ? upper : lower) atol=1e-16
        @test c.second_derivative_bound ≈ (upper-lower)/(.2^2/4) atol=1e-14
        # Compare actual functions, including fixed inputs and flat extensions.
        for x in range(domain...;length=31)
            error=primitive_value(f,x,r;domain_policy=:flat_extension)-primitive_value(f,x;domain_policy=:flat_extension)
            @test lower-1e-13<=error<=upper+1e-13
            @test relation==:upper ? error+p.output_shift<=1e-13 : error+p.output_shift>=-1e-13
        end
    end
    p=plan_pwl_relation(f,(242.,248.);relation=:upper,formulation=r,conservative=true)
    @test p.strategy==:affine && isempty(p.active_hinges) && p.output_shift==0.
    # Affine strategy alone is insufficient: a fixed point inside a patch
    # evaluates the surrogate, whose constant differs from the canonical value.
    p=plan_pwl_relation(f,(250.,250.);relation=:upper,formulation=r,conservative=true)
    @test p.strategy==:affine && isempty(p.active_hinges) && p.output_shift<0.
    for rep in (SoftplusFormulation(.2),AlgebraicFormulation(.2))
        p=plan_pwl_relation(f,(242.,248.);relation=:upper,formulation=rep,conservative=true)
        @test p.approximation_contract.error_scope==:global
        @test p.output_shift==-formulation_contract(f,rep).error_upper
    end
end
