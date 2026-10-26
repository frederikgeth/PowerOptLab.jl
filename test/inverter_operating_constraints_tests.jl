# Fixed POC voltage keeps the capability oracle independent of feeder feedback.
_constraint_grid(v=207.0; magnitudes=fill(v, 3)) = parse_bmopf("""
{"bus":{"poc":{"terminal_names":["a","b","c","n"],
 "perfectly_grounded_terminals":["n"],"v_min":[180,180,180],"v_max":[280,280,280]}},
 "voltage_source":{"vs":{"bus":"poc","terminal_map":["a","b","c"],
 "v_magnitude":$(Float64.(magnitudes)),"v_angle":[0,-2.0943951023931953,2.0943951023931953]}}}
"""; from_string=true)

_constraint_plant(; kwargs...) = AdvancedInverter(;
    id="pv", bus="poc", topology=:THREE_LEG,
    phase_terminals=["a", "b", "c"], neutral=nothing,
    s_max=6000.0, i_max=1.1*6000/(3*230), v_dc=800.0, c_dc=2e-3,
    r_filter=0.0, x_filter=0.0, kwargs...)

function _constraint_policy(; q_fraction=0.44, watt_fraction=nothing,
        basis=:rated, sensing=:average_voltage)
    vv = PWLFunction([180.0, 280.0], [q_fraction, q_fraction])
    vw = watt_fraction === nothing ? PWLFunction([253., 260.], [1., .2]) :
        PWLFunction([180., 280.], [watt_fraction, watt_fraction])
    intent = VoltVarWattIntent(volt_var=vv, volt_watt=vw,
        sensing=sensing, volt_watt_basis=basis, conflict_policy=:net)
    encoding = VoltVarWattEncoding(volt_var=LocalC2Formulation(.01),
        volt_watt=LocalC2Formulation(.01))
    VoltVarWattConstraints(intent, encoding)
end

const _CONSTRAINT_REQUEST = InverterControlRequest(
    p_available=6000.0, p_rated=6000.0, q_scale=6000.0)

@testset "Inverter operating constraints: validation" begin
    policy = _constraint_policy()
    device = ConstrainedInverter(_constraint_plant(), policy)
    @test device_id(device) == "pv"
    @test inverter_spec(device) === device.device
    @test validate_device(device, (_constraint_grid(),)) === device
    @test sensing_terminal(policy, device.device) === nothing
    @test_throws ArgumentError VoltVarWattConstraints(policy.positive; sensing_reference=:bad)
    @test_throws ArgumentError validate_device(ConstrainedInverter(
        _constraint_plant(topology=:SINGLE_PHASE, phase_terminals=["a"], neutral="n"), policy))
    @test_throws ArgumentError validate_device(ConstrainedInverter(
        _constraint_plant(grid_forming=true), policy))
    @test_throws ArgumentError validate_device(ConstrainedInverter(
        _constraint_plant(pwm_strategy=:CENTERED, f_sw=10000.0), policy))
    for ref in ("a", "missing")
        bad = ConstrainedInverter(device.device,
            VoltVarWattConstraints(policy.positive; sensing_reference=ref))
        @test_throws ArgumentError validate_device(bad, (_constraint_grid(),))
    end
    @test_throws ArgumentError solve_constrained_inverter(
        _constraint_grid(), device, _CONSTRAINT_REQUEST; objective=:loss)
    @test_throws ArgumentError solve_constrained_inverter(
        _constraint_grid(), device, _CONSTRAINT_REQUEST; objective=:feasibility)
    @test_throws ArgumentError solve_constrained_inverter(
        _constraint_grid(), device, _CONSTRAINT_REQUEST; p_set=1000.0)
    @test_throws ArgumentError solve_constrained_inverter(
        _constraint_grid(), device, _CONSTRAINT_REQUEST; s_base=NaN)
    @test_throws ArgumentError solve_constrained_inverter(
        _constraint_grid(), device, _CONSTRAINT_REQUEST; per_unit=false)
end

@testset "Inverter operating constraints: current and apparent capability oracles" begin
    plant = _constraint_plant()
    device = ConstrainedInverter(plant, _constraint_policy())
    for v in (200., 207., 210., 230.)
        result = solve_constrained_inverter(_constraint_grid(v), device,
            _CONSTRAINT_REQUEST; solver_options=("tol" => 1e-9,))
        @test result.solve.publishable
        capacity = min(plant.s_max, 3v*plant.i_max)
        @test result.plant.q_poc ≈ 2640.0 atol=.02
        @test result.plant.p_poc ≈ sqrt(capacity^2-2640.0^2) atol=.05
        @test maximum(result.plant.i_mag) <= plant.i_max + 1e-6
        @test abs(sum(result.grid_phase_current)) < 1e-6
        @test result.plant.i_negative < 1e-7
        @test abs(result.operation.q_constraint_residual_var) < .01
        @test abs(result.operation.q_tracking_residual_var) < .01
        @test result.operation.dc_availability_margin_W >= -.02
        @test result.operation.converter_apparent_power_margin_VA >= -.02
        @test minimum(result.operation.converter_current_margin_A) >= -1e-6
        @test solve_status(result) === result.solve
        @test solve_diagnostics(result) === result.operation
        @test !haskey(result.operation, :current_scale)
    end

    # The unchanged forward controller deliberately has a different saturation law.
    old = solve_controlled_inverter(_constraint_grid(),
        ControlledDevice(plant, SequenceController(device.policy.positive;
            limiter=CommonScaleLimiter(pq_priority=:var))), _CONSTRAINT_REQUEST;
        solver_options=("tol" => 1e-9,))
    new = solve_constrained_inverter(_constraint_grid(), device, _CONSTRAINT_REQUEST;
        solver_options=("tol" => 1e-9,))
    @test old.solve.publishable && new.solve.publishable
    @test old.plant.q_poc ≈ 2613.6 atol=1.0
    @test old.control.current_scale < .991
    @test new.plant.q_poc > old.plant.q_poc + 25.0
    @test new.plant.p_poc < old.plant.p_poc - 10.0

    reserved = ConstrainedInverter(_constraint_plant(
        pwm_ac_converter_reserve=(3., 3., 3.)), _constraint_policy())
    rms = solve_constrained_inverter(_constraint_grid(), reserved, _CONSTRAINT_REQUEST;
        solver_options=("tol" => 1e-9,))
    @test rms.solve.publishable
    fundamental = sqrt(plant.i_max^2-9)
    @test rms.plant.p_poc ≈ sqrt((3*207*fundamental)^2-2640^2) atol=.05
    @test minimum(rms.operation.converter_current_margin_A) >= -1e-6
end

@testset "Inverter operating constraints: availability, losses and volt-watt basis" begin
    request = InverterControlRequest(p_available=3000., p_rated=6000., q_scale=6000.)
    for (basis, expected) in ((:rated, 2900.), (:available, 1450.))
        device = ConstrainedInverter(_constraint_plant(p_loss_fixed=100.),
            _constraint_policy(q_fraction=.1, watt_fraction=.5, basis=basis))
        result = solve_constrained_inverter(_constraint_grid(230.), device, request;
            solver_options=("tol" => 1e-9,))
        @test result.solve.publishable
        @test result.plant.p_poc ≈ expected atol=.03
        @test result.plant.p_dc ≈ expected+100 atol=.03
        @test result.operation.ac_available_W ≈ 2900 atol=.03
        @test result.operation.dc_availability_margin_W >= -.03
        @test result.operation.volt_watt_margin_W >= -.03

        night = solve_constrained_inverter(_constraint_grid(230.), device,
            InverterControlRequest(p_available=0., p_rated=6000., q_scale=6000.);
            solver_options=("tol" => 1e-9,))
        @test night.solve.publishable
        @test night.plant.p_poc ≈ -100.0 atol=.03
        @test night.plant.p_dc ≈ 0.0 atol=.03
        @test night.plant.q_poc ≈ 600.0 atol=.03
    end
    # Exercise the sloped curve rather than only constant fractions.
    result = solve_constrained_inverter(_constraint_grid(256.5),
        ConstrainedInverter(_constraint_plant(), _constraint_policy(q_fraction=0.)),
        _CONSTRAINT_REQUEST; solver_options=("tol" => 1e-9,))
    @test result.solve.publishable
    @test result.plant.p_poc ≈ 3600.0 atol=.03
end

@testset "Inverter operating constraints: POC Q, filters and voltage policies" begin
    inv = _constraint_plant(r_filter=.05, x_filter=.15,
        r_filter_grid=.03, x_filter_grid=.10, c_filter_mid=10e-6,
        r_filter_damping=1., b_filter_shunt=1e-4, i_grid_max=8.0)
    for sensing in (:average_voltage, :worst_phase, :positive_sequence)
        device = ConstrainedInverter(inv, _constraint_policy(q_fraction=.2, sensing=sensing))
        result = solve_constrained_inverter(_constraint_grid(
            magnitudes=[245., 215., 230.]), device, _CONSTRAINT_REQUEST;
            solver_options=("tol" => 1e-9,))
        @test result.solve.publishable
        @test result.plant.q_poc ≈ 1200. atol=.03
        @test abs(result.plant.q_conv-result.plant.q_poc) > 10
        @test result.plant.i_negative < 1e-6
        @test maximum(result.plant.i_grid_mag) <= 8.0+1e-6
        @test minimum(result.operation.grid_current_margin_A) >= -1e-6
        # Physical POC exchange includes the passive POC shunt as well as grid-arm current.
        voltages = result.operation.phase_voltage
        q_arm = imag(sum(voltages .* conj.(result.grid_phase_current)))
        q_shunt = inv.b_filter_shunt*sum(abs2, voltages)
        @test result.plant.q_poc ≈ q_arm+q_shunt atol=.03
        @test result.plant.p_dc ≈ result.plant.p_poc + result.plant.p_filter_loss +
              result.plant.p_loss + result.plant.p_cap_loss atol=.03
    end
end

@testset "Inverter operating constraints: voltage-dependent Q and ripple limits" begin
    vv = PiecewiseLinearLaw([207.,220.,240.,258.], [.44,0.,0.,-.60];
        formulation=LocalC2Formulation(.01))
    for (positive, expected) in (
            (AverageVoltageVoltVarWatt(volt_var=vv), 6000*.44*12/13),
            (PositiveSequenceVoltVarWatt(volt_var=vv), 6000*.44*12/13),
            (WorstPhaseVoltVarWatt(volt_var=vv, conflict_policy=:low_voltage), 2640.))
        device = ConstrainedInverter(_constraint_plant(), VoltVarWattConstraints(positive))
        result = solve_constrained_inverter(_constraint_grid(magnitudes=[218.,208.,198.]),
            device, _CONSTRAINT_REQUEST; solver_options=("tol" => 1e-9,))
        @test result.solve.publishable
        @test result.plant.q_poc ≈ expected atol=2.0
        @test abs(result.operation.q_constraint_residual_var) < .03
        @test result.plant.i_negative < 1e-6
    end

    # Balanced current under unbalanced voltage creates real 2ω link stress.
    # Check that the physical constraints curtail P while Q still tracks.
    for inv in (_constraint_plant(dv2_max=.15),
                _constraint_plant(i_cap_max=.10),
                _constraint_plant(i_positive_max=5.))
        result = solve_constrained_inverter(_constraint_grid(magnitudes=[245.,215.,230.]),
            ConstrainedInverter(inv, _constraint_policy(q_fraction=.2)),
            _CONSTRAINT_REQUEST; solver_options=("tol" => 1e-9,))
        @test result.solve.publishable
        @test result.plant.q_poc ≈ 1200. atol=.03
        @test result.plant.p_poc < 4000.
        if inv.dv2_max !== nothing
            @test result.plant.dv2 ≈ inv.dv2_max atol=1e-6
            @test result.operation.dc_ripple_margin_V >= -1e-6
        elseif inv.i_cap_max !== nothing
            @test result.plant.i_cap_thermal ≈ inv.i_cap_max atol=1e-6
        else
            @test result.plant.i_positive ≈ inv.i_positive_max atol=1e-6
        end
    end
end

@testset "Inverter operating constraints: sensing reference" begin
    vv = PiecewiseLinearLaw([200.,260.], [.3,-.3]; formulation=LocalC2Formulation(.01))
    plant = _constraint_plant(i_max=15., r_filter=.05, x_filter=.15)
    request = InverterControlRequest(p_available=3000., p_rated=6000., q_scale=6000.)
    for reference in (:ground, "n")
        device = ConstrainedInverter(plant, VoltVarWattConstraints(
            AverageVoltageVoltVarWatt(volt_var=vv); sensing_reference=reference))
        result = solve_constrained_inverter(inv_grid3_floating_neutral(), device, request;
            solver_options=("tol" => 1e-9,))
        @test result.solve.publishable
        bus = result.bus["poc"]
        phasor(t) = complex(bus[t]["vr"], bus[t]["vi"])
        offset = reference === :ground ? 0im : phasor("n")
        sensed = [phasor(p)-offset for p in ("a", "b", "c")]
        @test collect(result.operation.phase_voltage) ≈ sensed atol=1e-6
        @test abs(phasor("n")) > .5
        mean_voltage = sum(abs, sensed)/3
        @test result.plant.q_poc ≈ 6000*clamp(.3-.01*(mean_voltage-200),-.3,.3) atol=.03
    end
end

@testset "Inverter operating constraints: fixed dispatch and infeasibility" begin
    device = ConstrainedInverter(_constraint_plant(), _constraint_policy())
    for p in (1000., 4000.)
        result = solve_constrained_inverter(_constraint_grid(), device,
            _CONSTRAINT_REQUEST; objective=:feasibility, p_set=p,
            solver_options=("tol" => 1e-9,))
        @test result.solve.publishable
        @test result.plant.p_poc ≈ p atol=.03
        @test result.plant.q_poc ≈ 2640 atol=.03
    end
    impossible = solve_constrained_inverter(_constraint_grid(),
        ConstrainedInverter(_constraint_plant(), _constraint_policy(q_fraction=1.2)),
        _CONSTRAINT_REQUEST; solver_options=("tol" => 1e-9,))
    @test !impossible.solve.publishable
    @test occursin("INFEASIBLE", impossible.termination_status)
    @test isnan(impossible.plant.p_poc)
    @test isnan(impossible.operation.q_target_var)
    @test all(isnan, real.(impossible.grid_phase_current))
end

@testset "Inverter operating constraints: external objective and scale invariance" begin
    device = ConstrainedInverter(_constraint_plant(), _constraint_policy())
    for base in (1e4, 1e6)
        maximum_export = solve_constrained_inverter(_constraint_grid(), device,
            _CONSTRAINT_REQUEST; s_base=base, solver_options=("tol" => 1e-9,))
        @test maximum_export.solve.publishable
        @test maximum_export.plant.p_poc ≈ sqrt(5940.0^2-2640.0^2) atol=.05
        @test minimum(maximum_export.operation.converter_current_margin_A) >= -1e-6
        handles = Ref{Any}()
        hook! = ctx -> begin
            m = PowerOptLab._opf_model(ctx)
            @objective(m, Min, 123.0)
            handles[] = stamp_device!(ctx, device; request=_CONSTRAINT_REQUEST)
            @test objective_sense(m) == JuMP.MOI.MIN_SENSE
            @test objective_function(m) == 123.0
            h = inverter_handles(device, handles[])
            # The caller selects an interior P rather than maximum feasible export.
            @objective(m, Min, ((h.p_poc-2000/h.sb)/(device.device.s_max/h.sb))^2)
        end
        ctx = PowerOptLab.build_opf_model(_constraint_grid(); per_unit=true, s_base=base,
            add_objective=false, model_hook! = hook!, optimizer=Ipopt.Optimizer)
        PowerOptLab.BMOPFTools.optimize_opf!(ctx)
        status = SolveStatus(PowerOptLab._solve_outcome(PowerOptLab._opf_model(ctx)))
        @test status.publishable
        result = extract_device(device, handles[], status)
        @test result.plant.p_poc ≈ 2000.0 atol=.1
        @test result.plant.q_poc ≈ 2640.0 atol=.03
    end
end
