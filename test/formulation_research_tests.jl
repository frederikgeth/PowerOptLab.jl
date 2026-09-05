using LinearAlgebra, ForwardDiff, TOML

# A third family staged entirely through the public extension interface.
struct ResearchRootHinge <: AbstractPWLSmoothing
    width::Float64
end
PowerOptLab.hinge_value(x,r::ResearchRootHinge) = (x+hypot(x,r.width))/2
PowerOptLab.hinge_derivatives(x,r::ResearchRootHinge) =
    ((1+x/hypot(x,r.width))/2,r.width^2/(2hypot(x,r.width)^3))
PowerOptLab.hinge_contract(r::ResearchRootHinge) =
    (error_lower=0.,error_upper=r.width/2,regularity=:C_infinity,
     width=r.width,second_derivative_bound=1/(2r.width))

@testset "Research primitives, extension and physical budgets" begin
    f = PWLFunction([220.,240.,250.,270.],[20.,20.,0.,0.])
    for family in (SoftplusFormulation,LocalC2Formulation)
        rep = smoothing_for_error(f,family,.01)
        c = formulation_contract(f,rep)
        @test c.error_lower ≈ -.01
        @test c.error_upper ≈ .01
    end
    @test_throws ArgumentError smoothing_for_error(f,ResearchRootHinge,.01)
    @test primitive_value(f,300.;domain_policy=:flat_extension) == 0.
    rep = ResearchRootHinge(.1)
    @test primitive_derivatives(f,240.,rep)[1] ≈
        ForwardDiff.derivative(x -> primitive_value(f,x,rep),240.)
    @test primitive_derivatives(f,240.,rep)[2] ≈
        ForwardDiff.derivative(x -> ForwardDiff.derivative(
            y -> primitive_value(f,y,rep),x),240.)
    low,up = MagnitudeApproximation(.01),MagnitudeApproximation(.01;direction=:upper)
    @test magnitude_value([0.,0.],low) == 0.
    @test magnitude_value([0.,0.],up) == .01
    @test norm(ForwardDiff.gradient(x -> magnitude_value(x,low),[0.,0.])) == 0.
    for x in ([3.,4.],[1e-8,0.])
        @test -.01 <= magnitude_value(x,low)-norm(x) <= 0.
        @test 0. <= magnitude_value(x,up)-norm(x) <= .01
    end
    @test magnitude_contract(low).hessian_bound == 100.
    @test affine_error_bound([2.,-3.],[(-.1,.2),(0.,.1)]) == (lower=-.5,upper=.4)
    @test local_error_response([3.;;],[.03]).delta ≈ [-.01]
    @test local_error_response([0.;;],[.03]).status == :ill_conditioned
    m = Model(Ipopt.Optimizer); set_silent(m)
    @variable(m,x,start=3.); fix(x,3.)
    mag = magnitude_expression(m,[x,4.],low;component_scale=2.,output_scale=10.)
    root = positive_root_expression(m,x;lower_bound=1.)
    @objective(m,Min,mag+root)
    optimize!(m)
    @test value(mag)*10 ≈ magnitude_value([6.,8.],low) atol=1e-8
    @test value(root) ≈ sqrt(3.)
    @test_throws ArgumentError ComplementarityGraph(scale=0.)
    c = Model(); @variable(c,v)
    h = formulate_pwl!(c,f,2v,ComplementarityGraph(scale=.1);input_scale=230.)
    @test h.complementarity_scale == .1
    # Complementarity difference equals physical hinge input divided by its own scale.
    aff = [constraint_object(ref).func for ref in all_constraints(c,AffExpr,MOI.EqualTo{Float64})]
    @test any(a -> abs(coefficient(a,v)) ≈ 4600.,aff)
end

@testset "Configurable cases, diagnostics and export" begin
    f = PWLFunction([220.,240.,250.,270.],[20.,20.,0.,0.])
    reference = resistive_equilibria(f,230.,1.)
    @test only(reference.points).current ≈ 40/3
    @test isempty(reference.intervals)
    @test length(resistive_equilibria(PWLFunction([0.,1.,2.],[0.,2.,0.]),0.,1.).points)==2
    @test length(resistive_equilibria(PWLFunction([0.,1.],[0.,1.]),0.,1.).intervals)==1
    @test only(resistive_equilibria(f,245.,0.).points).current==10.
    case = resistive_control_case(f;source_voltage=230.,resistance=1.)
    methods = [FormulationMethod("custom",ResearchRootHinge(.01),Ipopt.Optimizer;
        configure! = set_silent,options=(tol=1e-9,)),
        FormulationMethod("C2",(f,c)->smoothing_for_error(f,LocalC2Formulation,c.budget),
            Ipopt.Optimizer;configure! = set_silent,options=(tol=1e-9,))]
    configs = [(input_scale=230.,output_scale=20.,budget=.01)]
    rows = run_formulation_experiment([case],methods;configurations=configs,
        assess=row -> (investigator_accepted=row["strict_solver_success"],),on_error=:throw)
    @test all(r -> r["run_status"]=="finished" && r["strict_solver_success"],rows)
    @test all(r -> abs(r["metrics"].electrical_residual_V)<1e-7,rows)
    @test all(r -> abs(only(r["observations"]).surrogate_equation_error)<1e-7,rows)
    stopped = FormulationMethod("stopped",LocalC2Formulation(.1),Ipopt.Optimizer;
        configure! = set_silent,options=(max_iter=0,))
    r = only(run_formulation_experiment([case],[stopped];on_error=:throw))
    @test r["termination_status"] == "ITERATION_LIMIT"
    @test !r["strict_solver_success"]
    @test r["candidate_available"] && haskey(r,"metrics")
    # A primal ray/certificate is not an operating point to feed into metrics.
    mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    MOI.Utilities.set_mock_optimize!(mock,opt -> MOI.Utilities.mock_optimize!(opt,
        MOI.DUAL_INFEASIBLE,(MOI.INFEASIBILITY_CERTIFICATE,[1.])))
    ray_case = FormulationCase("ray",(rep,c) -> begin
        m = Model(); @variable(m,x); @objective(m,Min,-x)
        (model=m,metrics=()->error("certificate must not be evaluated as a point"))
    end)
    ray = only(run_formulation_experiment([ray_case],
        [FormulationMethod("mock ray",nothing,()->mock)];on_error=:throw))
    @test ray["run_status"]=="finished" && ray["primal_is_certificate"]
    @test !ray["candidate_available"] && !haskey(ray,"metrics")
    unsupported = FormulationCase("unsupported",(r,c)->throw(UnsupportedFormulation("by design")))
    broken = FormulationCase("broken",(r,c)->error("intentional failure"))
    statuses = run_formulation_experiment([unsupported,broken],methods[1:1])
    @test [r["run_status"] for r in statuses] == ["unsupported","error"]
    @test_throws ErrorException run_formulation_experiment([broken],methods[1:1];on_error=:throw)
    mktempdir() do dir
        path = write_formulation_results(joinpath(dir,"runs.toml"),[rows;r;statuses];
            sources=[@__FILE__],metadata=(purpose="smoke",))
        data = TOML.parsefile(path)
        @test data["schema_version"]==1 && length(data["runs"])==5
        @test length(only(values(data["source_sha256"])))==64
    end
end

@testset "Selected controller family and inverter adapter" begin
    build_device = (r,c) -> begin
        law = PiecewiseLinearLaw([220.,240.,250.,270.],[1.,1.,.2,.2];formulation=r)
        controller = SequenceController(AverageVoltageVoltVarWatt(volt_watt=law))
        inverter = AdvancedInverter(id="research",bus="poc",phase_terminals=["a","b","c"],
            neutral="n",topology=:THREE_LEG,s_max=20e3,i_max=40.,v_dc=700.,
            c_dc=1.1e-3,r_filter=.05,x_filter=.15)
        ControlledDevice(inverter,controller)
    end
    case = controlled_inverter_case(inv_grid3_bal(),build_device,
        InverterControlRequest(p_available=2e3,q_scale=4e3))
    methods = [FormulationMethod(string(T),T(.05),Ipopt.Optimizer;
        options=(tol=1e-8,max_iter=400),configure! = set_silent)
        for T in (SoftplusFormulation,LocalC2Formulation)]
    rows = run_formulation_experiment([case],methods;on_error=:throw)
    @test all(r -> r["strict_solver_success"],rows)
    @test all(r -> length(r["metrics"].phase_voltage_V)==3,rows)
    @test all(r -> r["metrics"].pre_capability_target_gap_A < .01,rows)
    for (row,method) in zip(rows,methods)
        device = build_device(method.representation,NamedTuple())
        numeric = evaluate_smooth(device.controller,
            InverterControlMeasurement(row["metrics"].phase_voltage_V),
            InverterControlRequest(p_available=2e3,q_scale=4e3),
            InverterControlRatings(device.device))
        @test row["metrics"].requested_active_power_W ≈ numeric.p_request atol=1e-4
    end
    law = PiecewiseLinearLaw([0.,1.],[0.,1.];formulation=LocalC2Formulation(.1))
    @test law.smoothing_epsilon == .1
    @test primitive_value(PWLFunction(law),2.,law.formulation;domain_policy=:flat_extension)==1.
    @test_throws ArgumentError PiecewiseLinearLaw([0.,1.],[0.,1.];
        formulation=LocalC2Formulation(.1),smoothing_epsilon=.2)
end
