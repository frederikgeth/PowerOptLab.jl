function _variant_external_solve(device, configure!; base=1e6, net=_constraint_grid(230.))
    handles = Ref{Any}()
    hook! = ctx -> begin
        handles[] = stamp_device!(ctx, device; request=_CONSTRAINT_REQUEST)
        configure!(PowerOptLab._opf_model(ctx), inverter_handles(device, handles[]))
    end
    ctx = PowerOptLab.build_opf_model(net; per_unit=true, s_base=base,
        add_objective=false, model_hook! = hook!, optimizer=Ipopt.Optimizer)
    m = PowerOptLab._opf_model(ctx)
    set_optimizer_attribute(m, "tol", 1e-9)
    PowerOptLab.BMOPFTools.optimize_opf!(ctx)
    status = SolveStatus(PowerOptLab._solve_outcome(m))
    return (; status, result=extract_device(device, handles[], status))
end

@testset "Operating variants: input domains and contradictory laws" begin
    for pf in (0., -0., Inf, -Inf, NaN, 1.01, -1.01)
        @test_throws ArgumentError FixedPowerFactor(pf)
    end
    for pf in (1., -1., .8, -.8, nextfloat(0.), prevfloat(1.))
        @test all(isfinite, PowerOptLab._pf_coefficients(FixedPowerFactor(pf)))
    end
    for bad in (NaN, Inf, -Inf)
        @test_throws ArgumentError FixedReactivePower(bad)
        @test_throws ArgumentError DCSourceBounds(p_min=bad)
        @test_throws ArgumentError DCSourceBounds(p_max=bad)
        @test_throws ArgumentError InverterOperatingConstraints(p_min=bad)
        @test_throws ArgumentError InverterOperatingConstraints(q_max=bad)
    end
    @test_throws ArgumentError DCSourceBounds(p_min=2., p_max=1.)
    @test_throws ArgumentError InverterOperatingConstraints(p_min=2., p_max=1.)
    @test_throws ArgumentError InverterOperatingConstraints(q_min=2., q_max=1.)
    @test_throws ArgumentError InverterOperatingConstraints(_constraint_policy().positive)
    @test_throws ArgumentError ConstrainedInverter(_constraint_plant(),
        InverterOperatingConstraints(); current_allocation=:independent_phases)
end

@testset "Operating variants: missing-curve defaults remain explicit" begin
    request = InverterControlRequest(p_available=3000.,p_rated=1000.,q_scale=6000.)
    for (positive,expected) in ((AverageVoltageVoltVarWatt(),3000.),
                               (AverageVoltageVoltVarWatt(volt_watt_basis=:rated),1000.))
        legacy = solve_constrained_inverter(_constraint_grid(230.),
            ConstrainedInverter(_constraint_plant(),VoltVarWattConstraints(positive)),
            request;solver_options=("tol"=>1e-9,))
        @test legacy.solve.publishable
        @test legacy.plant.q_poc ≈ 0. atol=.03
        @test legacy.operation.q_target_var ≈ 0. atol=.03
        @test legacy.plant.p_poc ≈ expected atol=.03
        @test legacy.operation.volt_watt_cap_W ≈ expected atol=.03
    end
    general = solve_constrained_inverter(_constraint_grid(230.),
        ConstrainedInverter(_constraint_plant(),InverterOperatingConstraints()),
        request;solver_options=("tol"=>1e-9,))
    @test general.solve.publishable
    @test general.plant.p_poc ≈ 3000. atol=.03
    @test general.operation.q_target_var === nothing
    @test general.operation.volt_watt_cap_W === nothing
end

@testset "Operating variants: rating-scaled contradictory dispatch" begin
    for base in (6000., 1e6, 1e8)
        device = ConstrainedInverter(_constraint_plant(p_loss_fixed=100.),
            InverterOperatingConstraints(reactive=FixedReactivePower(600.));
            dc_source=DCSourceBounds(p_min=-850.))
        failed = solve_constrained_inverter(_constraint_grid(230.), device,
            InverterControlRequest(p_available=0.,p_rated=6000.,q_scale=6000.);
            objective=:feasibility,p_set=-1000.,s_base=base,
            solver_options=("tol"=>1e-9,))
        @test !failed.solve.publishable
        @test isnan(failed.plant.p_poc)
        @test isnan(failed.operation.dc_lower_margin_W)
    end
end

@testset "Operating variants: internal voltage and split midpoint bind" begin
    for base in (6000.,1e6)
        plant = _constraint_plant(r_filter=1.,v_int_max=232.)
        result = solve_constrained_inverter(_constraint_grid(230.),
            ConstrainedInverter(plant,InverterOperatingConstraints(reactive=FixedReactivePower(0.))),
            _CONSTRAINT_REQUEST;s_base=base,solver_options=("tol"=>1e-9,))
        @test result.solve.publishable
        @test result.plant.p_poc ≈ 3*230*2 atol=.03
        @test maximum(result.plant.v_int_mag) ≈ 232. atol=1e-5
        # Equal split banks: dv_mid = In/(2ωC), with phase B/C unloaded.
        split = _constraint_plant(topology=:SPLIT_DC,neutral="n",In_max=20.,dv_mid_max=.1)
        device = ConstrainedInverter(split,InverterOperatingConstraints();current_allocation=:free)
        solved = _variant_external_solve(device, function (m,h)
            for k in 2:3
                @constraint(m,h.cri[k]==0)
                @constraint(m,h.cii[k]==0)
            end
            @constraint(m,h.cii[1]==0)
            @objective(m,Max,h.p_poc/(split.s_max/h.sb))
        end;base=base)
        @test solved.status.publishable
        @test solved.result.plant.i_neutral ≈ 2*2pi*50*.002*.1 atol=1e-6
        @test solved.result.plant.dv_mid ≈ .1 atol=1e-6
        @test solved.result.operation.midpoint_ripple_margin_V >= -1e-6
    end
    minimum_voltage = ConstrainedInverter(_constraint_plant(r_filter=1.,v_int_min=231.),
        InverterOperatingConstraints(reactive=FixedReactivePower(0.)))
    impossible = solve_constrained_inverter(_constraint_grid(230.),minimum_voltage,
        _CONSTRAINT_REQUEST;objective=:feasibility,p_set=0.,solver_options=("tol"=>1e-9,))
    @test !impossible.solve.publishable
end

@testset "Operating variants: fixed Q and signed PF including import and zero" begin
    for pf in (1., -1., .8, -.8, .95, -.95), p in (-1000., 0., 1000.)
        device = ConstrainedInverter(_constraint_plant(),
            InverterOperatingConstraints(reactive=FixedPowerFactor(pf));
            dc_source=DCSourceBounds(p_min=-2000.))
        result = solve_constrained_inverter(_constraint_grid(230.), device,
            _CONSTRAINT_REQUEST; objective=:feasibility, p_set=p,
            solver_options=("tol"=>1e-9,))
        expected_q = -sign(pf)*sqrt((1-abs(pf))*(1+abs(pf)))/abs(pf)*p
        @test result.solve.publishable
        @test result.plant.p_poc ≈ p atol=.03
        @test result.plant.q_poc ≈ expected_q atol=.03
        @test abs(result.operation.fixed_pf_residual_VA) < .03
        @test result.operation.q_target_var === nothing
        @test result.operation.volt_watt_cap_W === nothing
        @test result.operation.reactive_mode === :fixed_pf
    end
    for q in (-1200., 0., 1200.)
        result = solve_constrained_inverter(_constraint_grid(),
            ConstrainedInverter(_constraint_plant(), InverterOperatingConstraints(
                reactive=FixedReactivePower(q), p_max=2000., q_min=-1200., q_max=1200.)),
            _CONSTRAINT_REQUEST; solver_options=("tol"=>1e-9,))
        @test result.solve.publishable
        @test result.plant.p_poc ≈ 2000. atol=.03
        @test result.plant.q_poc ≈ q atol=.03
        @test result.operation.q_exact_target_var == q
        @test result.operation.reactive_lower_margin_var >= -.03
        @test result.operation.reactive_upper_margin_var >= -.03
    end
    # Hard-coded operating points, independent of the coefficient formula:
    # tan(acos(.9)) = 0.4843221, tan(acos(.95)) = 0.3286841.
    for (law, p, q) in ((FixedPowerFactor(.9, :absorbing), 1000., -484.3221),
                        (FixedPowerFactor(.9, :injecting), 1000., 484.3221),
                        (FixedPowerFactor(.95, :absorbing), 2000., -657.3682),
                        (FixedPowerFactor(.95, :injecting), -1000., -328.6841),
                        (FixedPowerFactor(1., :injecting), 1500., 0.))
        result = solve_constrained_inverter(_constraint_grid(230.),
            ConstrainedInverter(_constraint_plant(), InverterOperatingConstraints(reactive=law);
                dc_source=DCSourceBounds(p_min=-2000.)),
            _CONSTRAINT_REQUEST; objective=:feasibility, p_set=p,
            solver_options=("tol"=>1e-9,))
        @test result.solve.publishable
        @test result.plant.q_poc ≈ q atol=.03
    end
    @test FixedPowerFactor(.9, :absorbing).pf == .9
    @test FixedPowerFactor(.9, :injecting).pf == -.9
    for (pf, direction) in ((.9, :leading), (-.9, :absorbing), (0., :injecting),
                            (1.1, :absorbing), (NaN, :injecting))
        @test_throws ArgumentError FixedPowerFactor(pf, direction)
    end
    impossible = solve_constrained_inverter(_constraint_grid(),
        ConstrainedInverter(_constraint_plant(),
            InverterOperatingConstraints(reactive=FixedPowerFactor(1e-8))),
        _CONSTRAINT_REQUEST; objective=:feasibility, p_set=1000.,
        solver_options=("tol"=>1e-9,))
    @test !impossible.solve.publishable
    @test isnan(impossible.operation.fixed_pf_residual_VA)
end

@testset "Operating variants: explicit source direction and effective availability" begin
    plant = _constraint_plant(p_loss_fixed=100.)
    policy = InverterOperatingConstraints(reactive=FixedReactivePower(600.))
    night = InverterControlRequest(p_available=0., p_rated=6000., q_scale=6000.)
    pv = ConstrainedInverter(plant, policy; dc_source=DCSourceBounds(p_min=0.))
    result = solve_constrained_inverter(_constraint_grid(230.), pv, night;
        solver_options=("tol"=>1e-9,))
    @test result.solve.publishable
    @test result.plant.p_poc ≈ -100. atol=.03
    @test result.plant.p_dc ≈ 0. atol=.03
    @test result.operation.dc_lower_margin_W >= -.03
    disallowed = solve_constrained_inverter(_constraint_grid(230.), pv, night;
        objective=:feasibility, p_set=-1000., solver_options=("tol"=>1e-9,))
    @test !disallowed.solve.publishable
    @test isnan(disallowed.plant.p_dc)
    # The same import is admissible for explicitly bidirectional DC equipment.
    for lower in (-900., -850.)
        battery = ConstrainedInverter(plant, policy; dc_source=DCSourceBounds(p_min=lower))
        charged = solve_constrained_inverter(_constraint_grid(230.), battery, night;
            objective=:feasibility, p_set=-1000., solver_options=("tol"=>1e-9,))
        @test charged.solve.publishable == (lower == -900.)
        if charged.solve.publishable
            @test charged.plant.p_dc ≈ -900. atol=.03
        end
    end
    watt = PiecewiseLinearLaw([200.,270.], [.5,.5]; formulation=LocalC2Formulation(.01))
    for (basis, expected) in ((:rated,2900.), (:available,1450.))
        policy = InverterOperatingConstraints(AverageVoltageVoltVarWatt(
            volt_watt=watt, volt_watt_basis=basis); reactive=FixedPowerFactor(.8))
        device = ConstrainedInverter(plant, policy;
            dc_source=DCSourceBounds(p_min=0., p_max=3000.))
        result = solve_constrained_inverter(_constraint_grid(230.), device,
            _CONSTRAINT_REQUEST; solver_options=("tol"=>1e-9,))
        @test result.solve.publishable
        @test result.plant.p_poc ≈ expected atol=.03
        @test result.plant.q_poc ≈ -.75expected atol=.03
        @test result.operation.effective_dc_availability_W == 3000.
        @test result.operation.ac_available_W ≈ 2900. atol=.03
    end
    inconsistent = ConstrainedInverter(plant, policy; dc_source=DCSourceBounds(p_min=100.))
    failed = solve_constrained_inverter(_constraint_grid(), inconsistent, night)
    @test !failed.solve.publishable
end



@testset "Operating variants: free Q is selected by the external objective" begin
    device = ConstrainedInverter(_constraint_plant(), InverterOperatingConstraints(
        p_min=1000., p_max=1000., q_min=-1000., q_max=1000.))
    result = _variant_external_solve(device, (m,h) -> @objective(m, Min,
        ((h.q_poc-900/h.sb)/(6000/h.sb))^2))
    @test result.status.publishable
    @test result.result.plant.q_poc ≈ 900. atol=.1
    @test result.result.operation.q_target_var === nothing
    @test result.result.operation.q_tracking_residual_var === nothing
    @test result.result.operation.reactive_mode === :free_q
end

@testset "Operating variants: free negative sequence and LCL allocation location" begin
    a = cis(2pi/3)
    prescribed = (3.4+.2im, a^2*3+a*(.4+.2im), a*3+a^2*(.4+.2im))
    device = ConstrainedInverter(_constraint_plant(), InverterOperatingConstraints();
        current_allocation=:free)
    result = _variant_external_solve(device, function (m,h)
        for k in 1:2
            @constraint(m, h.cri[k] == real(prescribed[k])/h.ib)
            @constraint(m, h.cii[k] == imag(prescribed[k])/h.ib)
        end
        @objective(m, Min, 0.)
    end)
    @test result.status.publishable
    @test result.result.plant.i_negative ≈ abs(.4+.2im) atol=1e-6
    @test collect(result.result.converter_terminal.phase_current) ≈ collect(prescribed) atol=1e-6
    for allocation in (:balanced_converter, :balanced_grid)
        plant = _constraint_plant(topology=:FOUR_LEG, neutral="n", In_max=20.,
            r_filter=.05, x_filter=.15, r_filter_grid=.03, x_filter_grid=.10,
            c_filter_mid=30e-6, r_filter_damping=1.)
        device = ConstrainedInverter(plant,
            InverterOperatingConstraints(reactive=FixedReactivePower(600.), p_max=2000.);
            current_allocation=allocation)
        r = solve_constrained_inverter(_constraint_grid(magnitudes=[245.,215.,230.]),
            device, _CONSTRAINT_REQUEST; solver_options=("tol"=>1e-9,))
        @test r.solve.publishable
        conv = r.converter_terminal.phase_current
        grid = r.grid_phase_current
        i0(x) = sum(x)/3
        i2(x) = (x[1]+a^2*x[2]+a*x[3])/3
        balanced, other = allocation === :balanced_converter ? (conv,grid) : (grid,conv)
        @test abs(i0(balanced)) < 1e-6
        @test abs(i2(balanced)) < 1e-6
        @test abs(i0(other)) > 1e-3
        @test abs(i2(other)) > 1e-3
    end
end

@testset "Operating variants: four-leg neutral and split capacitor bounds" begin
    for base in (6000.,1e6), topo in (:FOUR_LEG,:SPLIT_DC)
        plant = topo === :FOUR_LEG ? _constraint_plant(topology=topo, neutral="n", In_max=1.5) :
            _constraint_plant(topology=topo, neutral="n", In_max=20.,
                i_cap_max=1., i_cap_upper_max=.8, i_cap_lower_max=1.2)
        device = ConstrainedInverter(plant, InverterOperatingConstraints(); current_allocation=:free)
        solved = _variant_external_solve(device, function (m,h)
            for k in 2:3
                @constraint(m, h.cri[k] == 0)
                @constraint(m, h.cii[k] == 0)
            end
            @constraint(m, h.cii[1] == 0)
            @objective(m, Max, h.p_poc)
        end; base=base)
        @test solved.status.publishable
        r = solved.result
        expected_current = topo === :FOUR_LEG ? 1.5 : .8/sqrt(.25+(230/(sqrt(2)*800))^2)
        @test r.plant.i_neutral ≈ expected_current atol=1e-5
        @test r.plant.p_poc ≈ 230expected_current atol=.03
        @test r.operation.neutral_current_margin_A >= -1e-5
        if topo === :SPLIT_DC
            @test r.plant.i_cap_thermal_upper ≈ .8 atol=1e-6
            @test r.operation.upper_capacitor_margin_A >= -1e-6
            @test r.operation.lower_capacitor_margin_A >= -1e-6
        end
    end
end
