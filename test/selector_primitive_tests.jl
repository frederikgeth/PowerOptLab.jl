using ForwardDiff

@testset "Shared selectors: contracts and independent derivatives" begin
    r = AlgebraicFormulation(.02)
    @test_throws ArgumentError AlgebraicFormulation(0.)
    @test hinge_value(-1e8,r)>0 # stable negative tail
    for x in (-.02,0.,.02,1.)
        d1,d2 = hinge_derivatives(x,r)
        @test d1 ≈ ForwardDiff.derivative(t -> hinge_value(t,r),x) atol=1e-12
        @test d2 ≈ ForwardDiff.derivative(t -> ForwardDiff.derivative(
            u -> hinge_value(u,r),t),x) atol=1e-10
    end
    f = PWLFunction([0.,1.],[0.,1.])
    @test formulation_contract(f,smoothing_for_error(f,AlgebraicFormulation,.01)).error_upper ≈ .01
    for family in (AlgebraicFormulation,SoftplusFormulation,LocalC2Formulation),
            (a,b) in ((0.,0.),(.02,0.),(1.,1.),(-1.,2.)), kind in (:min,:max)
        rep = family(.02)
        contract = selector_contract(rep;kind)
        exact = kind == :min ? min(a,b) : max(a,b)
        err = selector_value(a,b,rep;kind)-exact
        @test contract.error_lower-1e-14 <= err <= contract.error_upper+1e-14
    end
    for (a,b) in ((0.,0.),(0.,1.),(.1,.1),(.2,.9))
        value = selector_value(a,b,r;kind=:nonnegative_min)
        @test 0. <= value <= min(a,b)
        @test min(a,b)-value <= .01
    end
    @test_throws DomainError selector_value(-1.,1.,r;kind=:nonnegative_min)
    @test selector_value(0.,0.,LocalC2Formulation(.1);kind=:nonnegative_min)==0.
    @test symmetric_clip_value(0.,1.,r)==0.
    @test abs(symmetric_clip_value(2.,1.,r)-1.)<=.01
end

@testset "Shared selectors: physical numeric and stamped agreement" begin
    for family in (AlgebraicFormulation,SoftplusFormulation,LocalC2Formulation), si in (1.,230.)
        rep = family(.1)
        m = Model(Ipopt.Optimizer); set_silent(m)
        set_optimizer_attribute(m,"tol",1e-10)
        @variable(m,x,start=.03/si); @variable(m,y,start=.04/si)
        @constraint(m,x==.03/si); @constraint(m,y==.04/si)
        expressions = [selector_expression(m,x,y,rep;kind=k,input_scale=si,output_scale=20.)
            for k in (:min,:max,:nonnegative_min)]
        push!(expressions,symmetric_clip_expression(m,x,y,rep;input_scale=si,output_scale=20.))
        push!(expressions,hinge_expression(m,x,rep;input_scale=si,output_scale=20.))
        @objective(m,Min,sum(expressions))
        optimize!(m)
        @test termination_status(m)==MOI.LOCALLY_SOLVED
        expected = [selector_value(.03,.04,rep;kind=k) for k in (:min,:max,:nonnegative_min)]
        append!(expected,[symmetric_clip_value(.03,.04,rep),hinge_value(.03,rep)])
        @test value.(expressions).*20. ≈ expected atol=1e-8
    end
    # The staged lower norm uses BMOPFTools with physical scales and provenance.
    using_ctx = PowerOptLab.build_opf_model(inv_grid3_bal();per_unit=true,add_objective=false)
    m = PowerOptLab._opf_model(using_ctx)
    @variable(m,a,start=.03); @constraint(m,a==.03)
    rep = MagnitudeApproximation(.01;unit=:A)
    expr = magnitude_expression(using_ctx,(a,0.),rep;component_scale=100.,output_scale=10.,
        name="shared_norm_test")
    @objective(m,Min,expr)
    PowerOptLab.enforce_kcl!(using_ctx)
    set_optimizer(m,Ipopt.Optimizer); set_silent(m); optimize!(m)
    @test value(expr)*10. ≈ magnitude_value((3.,0.),rep) atol=1e-7
end
