@testset "Control intent and faithful lowering" begin
    watt = PWLFunction([220.,240.,250.],[1.,1.,0.])
    var = PWLFunction([210.,220.,240.,250.],[.5,0.,0.,-.5])
    intent = VoltVarWattIntent(volt_watt=watt,volt_var=var,volt_watt_basis=:rated,
        conflict_policy=:net,conflict_epsilon=.02)
    @test intent.volt_watt.input_unit == :V
    @test intent.volt_var.output_unit == :pu
    @test intent.volt_watt.hinges == ((-.1,240.),(.1,250.))
    @test_throws ArgumentError VoltVarWattIntent(volt_watt=PWLFunction([0.,1.],[0.,1.]))
    @test_throws ArgumentError VoltVarWattIntent(volt_watt=PWLFunction([0.,1.],[2.,0.]))
    @test_throws ArgumentError VoltVarWattIntent(volt_var=PWLFunction([0.,1.],[1.,0.];input_unit=:kV))
    @test_throws ArgumentError VoltVarWattIntent(sensing=:unknown)
    @test_throws ArgumentError VoltVarWattEncoding(extrema_epsilon=0.)
    @test_throws UnsupportedFormulation lower_positive_policy(intent,VoltVarWattEncoding())
    @test_throws UnsupportedFormulation lower_positive_policy(intent,
        VoltVarWattEncoding(volt_watt=ExactPWLGraph(),volt_var=LocalC2Formulation(.1)))
    ratings = InverterControlRatings(s_max=20e3,i_max=40.)
    request = InverterControlRequest(p_available=10e3,q_scale=8e3)
    voltage = InverterControlMeasurement([215.0 + 0im,245cis(-2pi/3),252cis(2pi/3)])
    results = []
    for family in (SoftplusFormulation,LocalC2Formulation,AlgebraicFormulation)
        encoding = VoltVarWattEncoding(volt_watt=family(.1),volt_var=family(.2))
        policy = lower_positive_policy(intent,encoding)
        @test policy isa WorstPhaseVoltVarWatt
        @test policy.conflict_policy == :net && policy.conflict_epsilon == .02
        @test policy.volt_watt_basis == :rated
        push!(results,evaluate_exact(SequenceController(policy),voltage,request,ratings))
    end
    @test all(r -> all(isapprox.(r.phase_current,results[1].phase_current)),results)
    for (sensing,T) in ((:average_voltage,AverageVoltageVoltVarWatt),
                       (:positive_sequence,PositiveSequenceVoltVarWatt))
        i = VoltVarWattIntent(volt_watt=watt,sensing=sensing,worst_phase_watt_guard=false)
        p = lower_positive_policy(i,VoltVarWattEncoding(volt_watt=LocalC2Formulation(.1)))
        @test p isa T
        sensing == :positive_sequence && (@test !p.worst_phase_watt_guard)
    end
end

@testset "Curve block reuse and explicit domain" begin
    intent = VoltVarWattIntent(volt_watt=PWLFunction([220.,240.,250.],[1.,1.,0.]))
    encoding = VoltVarWattEncoding(volt_watt=LocalC2Formulation(.2))
    m = Model(); @variable(m,v)
    h = formulate_control_curve!(m,intent,:volt_watt,v,encoding;domain=(200.,270.))
    before = lowering_statistics(m)
    @test formulate_control_curve!(m,intent,:volt_watt,v,encoding;domain=(200.,270.)) === h
    @test lowering_statistics(m) == before
    fresh = formulate_control_curve!(m,intent,:volt_watt,v,encoding;domain=(200.,270.),reuse=false)
    @test fresh !== h
    @test lowering_statistics(m).variables == before.variables+1
    @test lowering_statistics(m).pwl_operators == 1
    @test lowering_statistics(m).shared_curve_blocks == 1
    for kwargs in ((domain=(230.,245.),),(domain=(200.,270.),input_scale=230.))
        @test formulate_control_curve!(m,intent,:volt_watt,v,encoding;kwargs...) !== h
    end
    @test lowering_statistics(m).shared_curve_blocks == 3
    n = Model(); @variable(n,w)
    @test_throws ArgumentError formulate_control_curve!(n,intent,:volt_watt,v,encoding;domain=(200.,270.))
    @test_throws ArgumentError formulate_control_curve!(m,intent,:missing,v,encoding;domain=(200.,270.))
    @test_throws ArgumentError formulate_control_curve!(m,intent,:volt_var,v,encoding;domain=(200.,270.))
    @test_throws ArgumentError formulate_control_curve!(m,intent,:volt_watt,v,encoding;domain=(270.,200.))
    @test_throws ArgumentError formulate_control_curve!(m,intent,:volt_watt,v,encoding;domain=(200.,Inf))
    @test_throws ArgumentError formulate_control_curve!(m,intent,:volt_watt,v,encoding;domain=(200.,270.),input_scale=0.)
    other = formulate_control_curve!(n,intent,:volt_watt,w,encoding;domain=(200.,270.))
    @test other !== h && lowering_statistics(n).shared_curve_blocks == 1
    expression = 1.0v
    a = formulate_control_curve!(m,intent,:volt_watt,expression,encoding;domain=(200.,270.))
    add_to_expression!(expression,1.)
    b = formulate_control_curve!(m,intent,:volt_watt,expression,encoding;domain=(200.,270.))
    @test a !== b && lowering_statistics(m).shared_curve_blocks == 3
    @test a.input.constant == 0.
    @test b.input.constant == 1.
    # Truncating a segment before smoothing would change the law at the new
    # endpoint. Check faithful smoothing there and in both physical flat tails.
    for family in (SoftplusFormulation,LocalC2Formulation,AlgebraicFormulation),
            (domain,voltage) in (((242.,248.),242.),((200.,270.),200.),((200.,270.),270.))
        model = Model(Ipopt.Optimizer); set_silent(model)
        @variable(model,x,start=voltage/230.)
        @constraint(model,x == voltage/230.)
        rep = family(1.)
        handle = formulate_control_curve!(model,intent,:volt_watt,x,
            VoltVarWattEncoding(volt_watt=rep);domain,input_scale=230.,output_scale=.5)
        @objective(model,Min,handle.output)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test value(handle.output)*.5 ≈ primitive_value(intent.volt_watt,voltage,rep;
            domain_policy=:flat_extension) atol=1e-7
        audit = audit_pwl(handle)
        @test abs(audit.surrogate_equation_error)<1e-7
        @test audit.domain_violation<1e-7
        @test handle.domain == domain
    end
    hull = formulate_control_curve!(m,intent,:volt_watt,v,
        VoltVarWattEncoding(volt_watt=PWLConvexHull());domain=(242.,260.))
    @test hull.curve.breakpoints == (242.,250.,260.)
    @test collect(hull.curve.values) ≈ [.8,0.,0.]
end
