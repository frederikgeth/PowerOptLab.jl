# Check the exact computation graph against the independent firmware oracle at
# prescribed voltages. Starts supply a witness, not a claim of solver success.
@testset "Controller complementarity encoding" begin
    controller = _study_controller()
    @test controller.encoding === :smooth
    @test_throws ArgumentError SequenceController(controller.positive;encoding=:unknown)
    @test_throws ArgumentError SequenceController(controller.positive;encoding=ExactPWLGraph())
    for policy in (
            controller.positive,
            WorstPhaseVoltVarWatt(volt_watt=controller.positive.volt_watt,
                volt_var=controller.positive.volt_var,conflict_policy=:dominant),
            AverageVoltageVoltVarWatt(volt_watt=controller.positive.volt_watt,
                volt_var=controller.positive.volt_var),
            PositiveSequenceVoltVarWatt(volt_watt=controller.positive.volt_watt,
                volt_var=controller.positive.volt_var)),
        voltage in (190.,220.,230.,245.,280.),
        negative in (0.,5.),
        sb in (1e4,1e6)
        m = Model()
        b = PowerOptLab._ComplementarityControlModel(m,ComplementarityGraph(scale=1.))
        ctx = PowerOptLab._ComplementarityControlContext(nothing,b)
        vb = 230.; ib = sb/vb
        phasors = _ctrl_phases(complex(voltage),complex(negative))
        @variable(m, vre[k=1:3],start=real(phasors[k])/vb)
        @variable(m, vim[k=1:3],start=imag(phasors[k])/vb)
        magnitudes = [PowerOptLab._implicit_sqrt!(b,vre[k]^2+vim[k]^2;scale=1.) for k in 1:3]
        _,u1,u2 = PowerOptLab._sequence_pair(vre,vim)
        p,q,lo,hi = PowerOptLab._positive_smooth(ctx,policy,magnitudes,u1,_CTRL_REQUEST,sb,vb)
        startval(x) = JuMP.value(v -> something(start_value(v),0.),x)
        exact = evaluate_exact(SequenceController(policy),InverterControlMeasurement(phasors),
            _CTRL_REQUEST,InverterControlRatings(s_max=1e9,i_max=1e6))
        @test startval(p)*sb ≈ exact.p_request atol=1e-7
        @test startval(q)*sb ≈ exact.q_request atol=1e-7
        @test startval(lo)*vb ≈ minimum(abs,phasors) atol=1e-9
        @test startval(hi)*vb ≈ maximum(abs,phasors) atol=1e-9
        @test !has_lower_bound(vre[1]) && !has_upper_bound(vre[1])
        i1 = (real(exact.i1_request)/ib,imag(exact.i1_request)/ib)
        _,_,i2,eta = PowerOptLab._unbalance_smooth(ctx,controller.unbalance,u1,u2,i1,vb,ib,1.)
        _,_,expected_i2,expected_eta = PowerOptLab._unbalance_exact(
            controller.unbalance,exact.voltage_sequence[2],exact.voltage_sequence[3],exact.i1_request)
        @test complex(startval.(i2)...)*ib ≈ expected_i2 atol=1e-8
        @test startval(eta) ≈ expected_eta atol=1e-10
        for (F,S) in list_of_constraint_types(m)
            if S <: MOI.EqualTo
                for ref in all_constraints(m,F,S)
                    constraint = constraint_object(ref)
                    @test startval(constraint.func) ≈ constraint.set.value atol=1e-8
                end
            end
            S <: MOI.Complements || continue
            for ref in all_constraints(m,F,S)
                values = startval(constraint_object(ref).func)
                @test minimum(values) >= 0
                @test abs(prod(values)) <= 1e-12
            end
        end
    end
end

@testset "Exact MPCC limiter and allocator witnesses" begin
    for priority in (:proportional,:watt,:var), smax in (.02,2.),
        (pf,qf) in ((0.,0.),(.2,.3),(1.5,.4),(.7,-1.5))
        m = Model()
        b = PowerOptLab._ComplementarityControlModel(m,ComplementarityGraph())
        @variable(m,p,start=pf*smax)
        @variable(m,q,start=qf*smax)
        actual = PowerOptLab._limit_positive_power_smooth!(b,p,q,smax,1e-6*smax,priority,.001)
        expected = PowerOptLab._limit_positive_power_exact(pf*smax,qf*smax,smax,priority,.001)
        startval(x) = JuMP.value(v -> something(start_value(v),0.),x)
        @test collect(startval.(actual)) ≈ collect(expected) atol=1e-10
    end
    for offset in ((0.,0.),(.3,.4),(1.,0.)), direction in ((0.,0.),(.1,0.),(2.,1.)), limit in (0.,.5,1.)
        m = Model()
        b = PowerOptLab._ComplementarityControlModel(m,ComplementarityGraph())
        @variable(m,d[k=1:2],start=direction[k])
        actual = PowerOptLab._safe_direction_scale_implicit!(b,offset,d,limit,1e-6;
            magnitude_start=.5,scale=1.)
        expected = PowerOptLab._safe_direction_scale_exact(complex(offset...),complex(direction...),limit)
        startval(x) = JuMP.value(v -> something(start_value(v),0.),x)
        @test startval(actual) ≈ expected atol=1e-12
    end
end

@testset "MPCC publication requires firmware replay" begin
    # Deliberately misreport a stale controller witness as a successful solve.
    # This tests our acceptance gate, not the reliability of any real optimizer.
    c = _study_controller()
    controller = SequenceController(c.positive,c.unbalance,c.limiter;
        encoding=ComplementarityGraph())
    inverter = AdvancedInverter(id="audit",bus="poc",phase_terminals=["a","b","c"],
        neutral="n",topology=:THREE_LEG,s_max=20e3,i_max=40.,v_dc=700.,
        c_dc=1.1e-3,r_filter=.05,x_filter=.15)
    candidate = Ref{Vector{Float64}}()
    mock = MOI.Utilities.MockOptimizer(MOI.Utilities.UniversalFallback(
        MOI.Utilities.Model{Float64}());eval_objective_value=false)
    MOI.Utilities.set_mock_optimize!(mock, m -> begin
        MOI.Utilities.mock_optimize!(m,MOI.LOCALLY_SOLVED,(MOI.FEASIBLE_POINT,candidate[]))
        MOI.set(m,MOI.ObjectiveValue(),0.)
        MOI.set(m,MOI.SolveTimeSec(),0.)
    end)
    configure! = m -> begin
        variables = all_variables(m)
        candidate[] = [something(start_value(v),0.) for v in variables]
        index = findfirst(v -> name(v)=="vr_poc_a",variables)
        @test index !== nothing
        candidate[][index] = 1.5
    end
    result = solve_controlled_inverter(inv_grid3_bal(),ControlledDevice(inverter,controller),
        _CTRL_REQUEST;optimizer=()->mock,configure!,verbose=true)
    @test result.termination_status == "LOCALLY_SOLVED"
    @test !solve_status(result).publishable
    @test !solve_status(result).feasible
    @test !solve_status(result).optimal
    @test result.encoding_audit.status == :FAILED
    @test result.exact_smooth_current_residual > result.encoding_audit.current_tolerance_A
    @test all(z -> isnan(real(z)),result.control.phase_current)
end
