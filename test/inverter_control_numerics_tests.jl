# Independent, single-solve checks of the stamped equations and physical units.
# Each start/base is a fresh solve of the same target formulation; no warm starts,
# continuation, controller-law replay or solver retry is performed.
function _controller_numerics_case(net, device, request; s_base=1e6, start_factor=1.0)
    handles = Ref{Any}()
    hook! = ctx -> begin
        handles[] = stamp_device!(ctx, device; request=request)
        @objective(PowerOptLab._opf_model(ctx), Min,
            handles[].plant.p_loss + handles[].plant.p_cap_loss)
    end
    ctx = PowerOptLab.build_opf_model(net; per_unit=true, s_base=s_base, add_objective=false,
        model_hook! = hook!, optimizer=Ipopt.Optimizer)
    PowerOptLab.enforce_kcl!(ctx)
    model = PowerOptLab._opf_model(ctx)
    set_optimizer_attribute(model, "tol", 1e-8)
    set_optimizer_attribute(model, "max_iter", 500)
    for v in all_variables(model)
        start = start_value(v)
        if start !== nothing && !is_fixed(v)
            candidate = start_factor*start
            has_lower_bound(v) && (candidate = max(candidate, lower_bound(v)))
            has_upper_bound(v) && (candidate = min(candidate, upper_bound(v)))
            set_start_value(v, candidate)
        end
    end
    optimize!(model)
    outcome = PowerOptLab._solve_outcome(model)
    result = PowerOptLab._controlled_inverter_result(device, request, handles[],
        PowerOptLab.SolveStatus(outcome), outcome, Dict{String,Any}())
    return model, result, handles[], ctx
end

@testset "Controller numerics: norm and selector error budgets" begin
    # Upper norms protect capability denominators. The bound holds at zero,
    # around the smoothing width and far outside it.
    for epsilon in (1e-6, 1e-3, 0.1), radius in (0., epsilon/10, epsilon, 1., 40.)
        upper = magnitude_value((radius,),MagnitudeApproximation(epsilon;direction=:upper))
        @test radius <= upper <= radius + epsilon + 1e-14
    end
    # A nonnegative conservative scale selector must remain nonnegative even
    # when one protection reduces the command to zero.
    for (a,b) in ((0.,0.), (0.,1.), (1.,0.), (1e-8,1e-8), (.4,.4), (.2,.9)), epsilon in (1e-6, 1e-3)
        m = Model()
        selected = PowerOptLab._combine_safe_scale_implicit!(m,a,b,epsilon)
        @test 0 <= selected <= min(a,b)
        @test min(a,b)-selected <= epsilon/2 + 1e-15
    end
    # Symbolic zero coefficients must not introduce a squared-root equality.
    m = Model()
    @variable(m, x)
    count_before = num_variables(m)
    @test PowerOptLab._implicit_sqrt!(m, 0.0*x^2) == 0.0
    @test num_variables(m) == count_before
end

@testset "Controller numerics: physical accuracy across starts and bases" begin
    common = (id="numerics", bus="poc", phase_terminals=["a","b","c"],
        neutral="n", topology=:THREE_LEG, s_max=20e3, i_max=40.,
        v_dc=700., c_dc=1.1e-3, r_filter=.05, x_filter=.15)
    ripple = ControlledDevice(AdvancedInverter(; common..., dv2_max=.1),
        SequenceController(AverageVoltageVoltVarWatt()))
    zero = ControlledDevice(AdvancedInverter(; common...), _study_controller())
    law = PiecewiseLinearLaw([240.,250.],[1.,.2]; smoothing_epsilon=.01)
    volt_watt = ControlledDevice(AdvancedInverter(; common...),
        SequenceController(WorstPhaseVoltVarWatt(volt_watt=law)))
    lcl = ControlledDevice(AdvancedInverter(id="lcl", bus="poc",
        phase_terminals=["a","b","c"], neutral="n", topology=:THREE_LEG,
        s_max=20e3, i_max=40., i_grid_max=35., v_dc=750., c_dc=1.1e-3,
        r_filter=.02, x_filter=.06, r_filter_grid=.03, x_filter_grid=.09,
        c_filter_mid=30e-6, r_filter_damping=.5),
        SequenceController(PositiveSequenceVoltVarWatt(); current_target=:grid))
    balanced = inv_grid3_src(mags=fill(230.,3), angs=[0.,-2pi/3,2pi/3])
    cases = (("tight ripple", inv_grid3_unbal(), ripple, _CTRL_REQUEST),
        ("balanced zero dispatch", balanced, zero, InverterControlRequest(p_available=0.,q_scale=0.)),
        ("Volt-watt plateau", inv_grid3_bal(255.), volt_watt, InverterControlRequest(p_available=12e3,q_scale=0.)),
        ("LCL grid target", inv_grid3_bal(), lcl, InverterControlRequest(p_available=9e3,q_scale=0.)))
    for (name,net,device,request) in cases
        references = []
        # A decade in either direction changes numerical units, not physics.
        for base in (1e6,1e5,1e7), start in (1.,.98,1.02)
            @testset "$name base=$base start=$start" begin
                model,r,h,ctx = _controller_numerics_case(net,device,request;
                    s_base=base,start_factor=start)
                @test r.solve.publishable
                if r.solve.publishable
                    @test maximum(values(primal_feasibility_report(model; atol=0.)); init=0.) <= 1e-7
                    target = device.controller.current_target isa PowerOptLab.ConverterCurrentTarget ?
                        r.converter_terminal.phase_current : r.grid_phase_current
                    @test maximum(abs,target .- r.control.phase_current) <= 1e-4 # A
                    @test maximum(r.plant.i_mag) <= device.device.i_max + 1e-4 # A
                    if device.device.i_grid_max !== nothing
                        @test maximum(abs,r.grid_phase_current) <= device.device.i_grid_max + 1e-4
                    end
                    @test abs(r.converter_terminal.total_power) <= device.device.s_max + .1 # VA
                    @test 0 <= r.control.current_scale <= 1.0
                    @test r.exact_smooth_current_residual <= 1e-3 # A at this operating point
                    if device.device.dv2_max !== nothing
                        @test r.plant.dv2 <= device.device.dv2_max + 1e-6 # V
                    end
                    if name == "balanced zero dispatch"
                        @test maximum(r.plant.i_mag) <= 1e-4 # A
                        @test abs(r.control.i2) <= 1e-6 # A
                    end
                    push!(references,r)
                end
            end
        end
        if length(references)==9
            for r in references[2:end]
                @test collect(r.control.phase_current) ≈ collect(references[1].control.phase_current) atol=1e-4
                @test r.plant.p_poc ≈ references[1].plant.p_poc atol=.1
            end
        end
    end
end
