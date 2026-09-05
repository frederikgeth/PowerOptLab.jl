using ForwardDiff

@testset "PWL primitive contracts and analytic reference" begin
    xs, ys = [-2., -1., 1., 2.], [1., 0., 0., 1.]
    f = PWLFunction(xs,ys; input_unit=:V,output_unit=:A)
    xs[1] = -99. # Input arrays are copied into immutable canonical data.
    @test first(f.breakpoints) == -2.
    @test primitive_value(f,-1.5) == .5
    @test primitive_value(f,0.) == 0.
    @test primitive_value(f,1.5) == .5
    @test_throws DomainError primitive_value(f,3.)
    @test_throws DomainError primitive_value(f,NaN)
    @test_throws ArgumentError PWLFunction([0.,0.],[1.,2.])
    @test_throws ArgumentError PWLFunction([0.,1.],[1.,Inf])
    @test_throws ArgumentError PWLFunction([0.,1.],[1.])
    for T in (SoftplusFormulation,LocalC2Formulation), width in (0.,-1.,Inf,NaN)
        @test_throws ArgumentError T(width)
    end
    for T in (SoftplusFormulation,LocalC2Formulation), width in (.01,.1,2.)
        r = T(width)
        c = formulation_contract(f,r)
        @test c.semantics == :surrogate_graph
        @test c.input_unit == :V && c.output_unit == :A
        for x in range(-2.,2.;length=51)
            error = primitive_value(f,x,r)-primitive_value(f,x)
            @test c.error_lower-1e-13 <= error <= c.error_upper+1e-13
            # Numerically differentiate the value oracle independently of the
            # supplied derivative formulas (softplus oracle is Float64-only).
            d1,d2 = primitive_derivatives(f,x,r)
            if -1.99 < x < 1.99
                step = 1e-5*width
                fd = (primitive_value(f,x+step,r)-primitive_value(f,x-step,r))/(2step)
                @test isapprox(d1,fd; atol=2e-8)
            end
            if r isa LocalC2Formulation
                @test isapprox(d1,ForwardDiff.derivative(z -> primitive_value(f,z,r),x); atol=1e-12)
                @test isapprox(d2,ForwardDiff.derivative(z -> ForwardDiff.derivative(
                    w -> primitive_value(f,w,r),z),x); atol=1e-10)
            end
        end
    end
    # The hinge is linear to the right until the distant endpoint clamp.
    hinge = PWLFunction([-2.,0.,2.],[0.,0.,2.])
    r = LocalC2Formulation(.2)
    @test primitive_value(hinge,0.,r) == 3*.2/16
    for x in (-.2,.2)
        d1,d2 = primitive_derivatives(hinge,x,r)
        @test d1 == (x < 0 ? 0. : 1.)
        @test d2 == 0.
        @test primitive_value(hinge,x,r) ≈ max(0.,x)
    end
    @test primitive_value(f,0.,LocalC2Formulation(.1)) == 0.
    # Derivatives and widths transform with physical input/output units.
    g = PWLFunction(collect(f.breakpoints).*1000,collect(f.values).*2)
    for T in (SoftplusFormulation,LocalC2Formulation), x in (-1.01,0.,.99,1.8)
        @test primitive_value(g,1000x,T(100.)) ≈ 2primitive_value(f,x,T(.1)) atol=1e-12
        d1,d2 = primitive_derivatives(f,x,T(.1))
        e1,e2 = primitive_derivatives(g,1000x,T(100.))
        @test e1 ≈ 2d1/1000 atol=1e-14
        @test e2 ≈ 2d2/1e6 atol=1e-14
    end
    @test formulation_contract(f,PWLConvexHull()).error_upper === nothing
    @test formulation_contract(f,ComplementarityGraph()).semantics == :exact_graph
    @test formulation_contract(f,ExactPWLGraph()).error_upper == 0.
    constant = PWLFunction([0.,1.],[2.,2.])
    @test formulation_contract(constant,SoftplusFormulation(.1)).error_upper == 0.
    @test primitive_value(constant,.5,LocalC2Formulation(.1)) == 2.
end

@testset "JuMP PWL representations and coordinate scaling" begin
    f = PWLFunction([-2.,-1.,1.,2.],[1.,0.,0.,1.]; input_unit=:V,output_unit=:A)
    for r in (SoftplusFormulation(.1),LocalC2Formulation(.1)),
        (si,so) in ((1.,1.),(230.,40.),(1000.,100.)), target in (-1.5,-1.,0.,1.,1.5)
        m = Model(Ipopt.Optimizer)
        set_silent(m)
        set_optimizer_attribute(m,"tol",1e-10)
        @variable(m,x,start=target/si)
        h = formulate_pwl!(m,f,x,r; input_scale=si,output_scale=so)
        @constraint(m,x == target/si)
        @objective(m,Min,h.output^2)
        optimize!(m)
        @test termination_status(m) == JuMP.MOI.LOCALLY_SOLVED
        a = audit_pwl(h)
        @test a.input ≈ target atol=1e-7
        @test a.output ≈ primitive_value(f,target,r) atol=1e-7
        @test abs(a.surrogate_equation_error) <= 1e-7
        @test a.semantics == :surrogate_graph
    end
    # Hull admits a convex combination that is not a point on the controller.
    m = Model(Ipopt.Optimizer)
    set_silent(m)
    @variable(m,x)
    h = formulate_pwl!(m,f,x,PWLConvexHull())
    @constraint(m,x == 0.)
    @objective(m,Max,h.output)
    optimize!(m)
    a = audit_pwl(h)
    @test a.semantics == :outer_relaxation
    @test a.output ≈ 1. atol=1e-7
    @test a.exact_graph_error ≈ 1. atol=1e-7
    # Endpoint hinge elimination: an affine graph needs no complementarity.
    m = Model()
    @variable(m,x)
    affine = PWLFunction([0.,1.],[2.,3.])
    h = formulate_pwl!(m,affine,x,ComplementarityGraph())
    @test isempty(h.pairs)
    h2 = formulate_pwl!(m,f,x,ComplementarityGraph())
    @test length(h2.pairs) == 2
    @test_throws ArgumentError audit_pwl(h2)
    n = num_variables(m)
    @test_throws ArgumentError formulate_pwl!(m,f,x,SoftplusFormulation(.1); input_scale=0)
    @test num_variables(m) == n
    if Base.get_extension(PowerOptLab,:PowerOptLabPiecewiseLinearOptExt) === nothing
        @test_throws ArgumentError formulate_pwl!(m,f,x,ExactPWLGraph())
        @test num_variables(m) == n
    end
end

@testset "BMOPFTools staged softplus reuse" begin
    ctx = PowerOptLab.build_opf_model(single_bus_net(); add_objective=false)
    PowerOptLab.enforce_kcl!(ctx)
    m = PowerOptLab.BMOPFTools.opf_model(ctx)
    @variable(m,x,start=.5)
    f = PWLFunction([0.,1.,2.],[0.,1.,1.])
    h = formulate_pwl!(ctx,f,x,SoftplusFormulation(.1))
    @constraint(m,x == .5)
    @objective(m,Min,h.output^2)
    optimize!(m)
    @test termination_status(m) == JuMP.MOI.LOCALLY_SOLVED
    @test audit_pwl(h).output ≈ primitive_value(f,.5,SoftplusFormulation(.1)) atol=1e-7
    @test !haskey(m.ext,:PowerOptLabPWLOperators) # Upstream builder owns registration.
end
