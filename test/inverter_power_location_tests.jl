# A fixed voltage source isolates device limits from feeder constraints.
_power_location_grid() = parse_bmopf("""
{"bus":{"poc":{"terminal_names":["a","b","c","n"],
 "perfectly_grounded_terminals":["n"],"v_min":[180,180,180],"v_max":[280,280,280]}},
 "voltage_source":{"vs":{"bus":"poc","terminal_map":["a","b","c"],
 "v_magnitude":[230,230,230],"v_angle":[0,-2.0943951023931953,2.0943951023931953]}}}
"""; from_string=true)
_power_location_plant(; kwargs...) = AdvancedInverter(;
    id="pv", bus="poc", phase_terminals=["a", "b", "c"], neutral="n",
    topology=:THREE_LEG, s_max=10e3, v_dc=800., c_dc=2e-3,
    r_filter=.2, x_filter=.5, kwargs...)
const _POWER_LOCATION_OPTIONS = ("tol"=>1e-9, "bound_relax_factor"=>0., "max_iter"=>1000)
const _POWER_LOCATION_REQUEST = InverterControlRequest(
    p_available=20e3, p_rated=20e3, q_scale=10e3)

@testset "Inverter power boundary: configuration and optional current allocation" begin
    @test _power_location_plant().power_limit_location === :converter
    @test_throws ArgumentError validate_device(_power_location_plant(power_limit_location=:grid))
    for topology in (:FOUR_LEG, :SPLIT_DC)
        inv = _power_location_plant(; topology)
        @test validate_device(inv) === nothing
    end
    @test_throws ArgumentError InverterControlRatings(s_max=1., i_scale=0.)
    @test_throws ArgumentError InverterControlRatings(s_max=1., i_max=Inf)
    @test InverterControlRatings(_power_location_plant()).i_max === nothing
    for target in (ConverterCurrentTarget(), GridCurrentTarget())
        @test InverterControlRatings(_power_location_plant(), target).i_max === nothing
    end
    # No hidden current cap: a negative-sequence request far above the numerical
    # reference is retained, for both controller limiter policies and evaluators.
    voltage = ComplexF64[250, 215cis(-2pi/3), 225cis(2pi/3)]
    measurement = InverterControlMeasurement(voltage)
    unbalance = NegativeSequenceAdmittanceDroop(
        PiecewiseLinearLaw([0., 1.], [10., 10.]; smoothing_epsilon=1e-4))
    for limiter in (CommonScaleLimiter(), SequencePriorityLimiter(priority=:power))
        controller = SequenceController(AverageVoltageVoltVarWatt(), unbalance, limiter)
        ratings = InverterControlRatings(s_max=10e3, i_scale=.01)
        for evaluate in (evaluate_exact, evaluate_smooth)
            r = evaluate(controller, measurement, _POWER_LOCATION_REQUEST, ratings)
            @test r.current_scale == 1.
            @test r.i2 ≈ r.i2_request
            @test maximum(abs, r.phase_current) > 10.
        end
    end
end

@testset "Inverter power boundary: total PCC circle, losses and shunts" begin
    for location in (:converter, :poc), lcl in (false, true), base in (1e4, 1e8)
        inv = _power_location_plant(power_limit_location=location,
            r_filter_grid=lcl ? .1 : 0., x_filter_grid=lcl ? .2 : 0.,
            c_filter_mid=lcl ? 30e-6 : 0., b_filter_shunt=lcl ? .002 : 0.)
        device = ConstrainedInverter(inv,
            InverterOperatingConstraints(reactive=FixedReactivePower(0.));
            current_allocation=:balanced_grid, dc_source=DCSourceBounds(p_min=0.))
        r = solve_constrained_inverter(_power_location_grid(), device, _POWER_LOCATION_REQUEST;
            s_base=base, solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        p, op = r.plant, r.operation
        spoc, sconv = hypot(p.p_poc, p.q_poc), hypot(p.p_conv, p.q_conv)
        @test (location === :poc ? spoc : sconv) ≈ inv.s_max atol=.02
        @test op.power_limit_location === location
        @test op.apparent_power_margin_VA ≈ 0. atol=.02
        @test op.converter_current_margin_A === nothing
        @test op.grid_current_margin_A === nothing
        @test p.p_dc ≈ p.p_poc + p.p_filter_loss + p.p_loss + p.p_cap_loss atol=1e-6
        @test p.p_filter_loss > 100.
        # Independently reconstruct net boundary exchange, including the shunt.
        snet = sum(p.poc_phase_voltage .* conj.(r.grid_phase_current)) +
            im*inv.b_filter_shunt*sum(abs2, p.poc_phase_voltage)
        @test snet ≈ complex(p.p_poc, p.q_poc) atol=1e-6
        if location === :poc
            @test sconv > inv.s_max + 100.
            @test op.converter_apparent_power_margin_VA === nothing
            @test op.poc_apparent_power_margin_VA ≈ 0. atol=.02
        else
            @test op.poc_apparent_power_margin_VA === nothing
            @test op.converter_apparent_power_margin_VA ≈ 0. atol=.02
        end
    end
end

@testset "Inverter power boundary: declared current limits remain hard" begin
    for (ic, ig) in ((8., nothing), (nothing, 7.))
        inv = _power_location_plant(power_limit_location=:poc, i_max=ic, i_grid_max=ig,
            r_filter_grid=.1, x_filter_grid=.2, c_filter_mid=30e-6)
        device = ConstrainedInverter(inv,
            InverterOperatingConstraints(reactive=FixedReactivePower(0.));
            current_allocation=:balanced_grid)
        r = solve_constrained_inverter(_power_location_grid(), device, _POWER_LOCATION_REQUEST;
            solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        @test r.plant.p_poc < 6000.
        if ic !== nothing
            @test maximum(r.plant.i_mag) ≈ ic atol=1e-5
        else
            @test maximum(r.plant.i_grid_mag) ≈ ig atol=1e-5
        end
    end
end

@testset "Inverter power boundary: controller exact and smooth physical backoff" begin
    for target in (:converter, :grid), priority in (:common, :power), lcl in (false, true)
        inv = _power_location_plant(power_limit_location=:poc,
            r_filter_grid=lcl ? .1 : 0., x_filter_grid=lcl ? .2 : 0.,
            c_filter_mid=lcl ? 30e-6 : 0., b_filter_shunt=lcl ? .003 : 0.)
        limiter = priority === :common ? CommonScaleLimiter() : SequencePriorityLimiter(priority=:power)
        controller = SequenceController(AverageVoltageVoltVarWatt(), NoUnbalanceControl(), limiter;
            current_target=target)
        r = solve_controlled_inverter(_power_location_grid(), ControlledDevice(inv, controller),
            _POWER_LOCATION_REQUEST; solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        @test hypot(r.plant.p_poc, r.plant.q_poc) <= inv.s_max + .02
        if !lcl
            @test hypot(r.plant.p_conv, r.plant.q_conv) > inv.s_max + 50.
        end
        @test r.plant.power_limit_location === :poc
        @test r.plant.converter_apparent_power_margin_VA === nothing
        @test r.exact_smooth_current_residual < .01
        @test solve_diagnostics(r).power_limit_location === :poc
    end

end

@testset "Inverter power boundary: standalone plant reporting and SI parity" begin
    net = parse_bmopf("""
    {"bus":{"poc":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]}},
     "voltage_source":{"vs":{"bus":"poc","terminal_map":["1"],"v_magnitude":[230],"v_angle":[0]}}}
    """; from_string=true)
    inv = AdvancedInverter(id="pv", bus="poc", s_max=5000., power_limit_location=:poc,
        r_filter=.2, x_filter=.5, b_filter_shunt=.003)
    for pu in (false, true)
        r = solve_advanced_inverter(net, inv; per_unit=pu, q_set=0.,
            solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        @test r.p_poc ≈ 5000. atol=.02
        @test hypot(r.p_conv, r.q_conv) > 5050.
        @test r.power_limit_location === :poc
        @test r.apparent_power_rating_VA == 5000.
        @test r.apparent_power_margin_VA ≈ 0. atol=.02
        @test solve_diagnostics(r).converter_apparent_power_margin_VA === nothing
        @test solve_diagnostics(r).poc_apparent_power_margin_VA ≈ 0. atol=.02
    end
end

@testset "Inverter power boundary: fleet and sizing report the selected rating" begin
    net = _power_location_grid()
    net["ibr"] = Dict("pv" => Dict{String,Any}(
        "bus"=>"poc", "terminal_map"=>["a","b","c","n"],
        "topology"=>"FOUR_LEG", "prime_mover"=>"PV", "s_max"=>ones(3),
        "p_min"=>zeros(3), "p_max"=>zeros(3), "q_min"=>zeros(3), "q_max"=>zeros(3)))
    inv = _power_location_plant(power_limit_location=:poc)
    device = ControlledDevice(inv, SequenceController(AverageVoltageVoltVarWatt()))
    fleet = ControlledInverterFleetSpec(Dict("pv"=>device), Dict("pv"=>_POWER_LOCATION_REQUEST))
    study = run_inverter_control_study([InverterControlStudyCase(
        scenario_id="pcc", variant_id="pcc", network=net, fleet=fleet)];
        solver_options=_POWER_LOCATION_OPTIONS)
    @test solve_status(study).publishable
    row = only(controlled_inverter_rows(only(study.cases).result))
    @test row.power_limit_location === :poc
    @test row.control_power_basis === :net_poc
    @test row.pq_priority === :proportional
    @test_throws ArgumentError controlled_inverter_fleet_loop_gain(
        fleet, Dict(), zeros(4,4))
    @test row.converter_apparent_power_VA > inv.s_max + 100.
    @test row.poc_apparent_power_VA <= inv.s_max + .02
    @test row.apparent_power_margin_VA == inv.s_max - row.poc_apparent_power_VA
    requirements = only(inverter_control_hardware_requirement_rows(study;
        allowed_dc_ripple_fraction=.02))
    @test requirements.power_limit_location === :poc
    @test requirements.installed_apparent_power_rating_VA == inv.s_max
    @test requirements.installed_poc_apparent_power_rating_VA == inv.s_max
    @test isnan(requirements.installed_converter_apparent_power_rating_VA)
    @test requirements.apparent_power_binding
    @test requirements.achieved_converter_apparent_power_VA > inv.s_max + 100.
    @test requirements.achieved_poc_apparent_power_VA <= inv.s_max + .02
    resized = resize_controlled_inverter_fleet(fleet, InverterHardwareSweepPoint(id="identity"))
    @test inverter_spec(resized.devices["pv"]).power_limit_location === :poc
    @test inverter_spec(resized.devices["pv"]).i_max === nothing
end

@testset "Inverter power boundary: optional arm ratings in controlled LCL plants" begin
    for target in (:converter, :grid), rating in (:converter, :grid)
        inv = _power_location_plant(power_limit_location=:poc,
            i_max=rating === :converter ? 8. : nothing,
            i_grid_max=rating === :grid ? 7. : nothing,
            c_filter_mid=30e-6, r_filter_grid=.1, x_filter_grid=.2)
        controlled = ControlledDevice(inv, SequenceController(AverageVoltageVoltVarWatt();
            current_target=target))
        r = solve_controlled_inverter(_power_location_grid(), controlled, _POWER_LOCATION_REQUEST;
            solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        @test r.exact_smooth_current_residual < .01
        current = rating === :converter ? r.plant.i_mag : r.plant.i_grid_mag
        limit = rating === :converter ? inv.i_max : inv.i_grid_max
        @test maximum(current) <= limit + 1e-4
        @test 1000. < r.plant.p_poc < 6000.
        @test r.control.current_scale < 1.
    end
end

@testset "Inverter power boundary: unbalance and ripple still use converter physics" begin
    for location in (:converter, :poc)
        inv = _power_location_plant(power_limit_location=location, dv2_max=.5)
        negative = NegativeSequenceAdmittanceDroop(
            PiecewiseLinearLaw([0., 1.], [.2, .2]; smoothing_epsilon=1e-4))
        controlled = ControlledDevice(inv, SequenceController(AverageVoltageVoltVarWatt(), negative,
            SequencePriorityLimiter(priority=:balance)))
        r = solve_controlled_inverter(inv_grid3_unbal(), controlled, _POWER_LOCATION_REQUEST;
            solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        @test r.plant.dv2 <= inv.dv2_max + 1e-5
        @test r.plant.apparent_power_margin_VA >= -.02
        @test r.exact_smooth_current_residual < .01
        @test r.plant.i_negative > .01
    end
end

@testset "PCC review: physical reference, signed shunt and independent backoff oracle" begin
    for (target, base) in ((:converter, 1e4), (:grid, 1e8)), b in (-.003, .003)
        inv = _power_location_plant(power_limit_location=:poc, b_filter_shunt=b)
        controller = SequenceController(AverageVoltageVoltVarWatt();
            current_target=target, sensing_reference=:ground)
        r = solve_controlled_inverter(inv_grid3_floating_neutral(),
            ControlledDevice(inv, controller), _POWER_LOCATION_REQUEST;
            s_base=base, solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        bus = r.bus["poc"]
        ground_voltage = Tuple(complex(bus[t]["vr"], bus[t]["vi"]) for t in ("a", "b", "c"))
        vn = complex(bus["n"]["vr"], bus["n"]["vi"])
        vpcc = ground_voltage .- vn
        @test abs(vn) > .5
        @test collect(r.control.phase_voltage) ≈ collect(ground_voltage) atol=1e-6
        @test collect(r.plant.poc_phase_voltage) ≈ collect(vpcc) atol=1e-6
        # Reconstruct physical power from all four external conductor injections.
        # I_shunt = -j B V_pn; its neutral return is included independently here.
        ipcc = r.grid_phase_current .- im*b .* vpcc
        sneutral = vn*conj(-sum(ipcc))
        snet = sum(ground_voltage .* conj.(ipcc)) + sneutral
        sarm = sum(vpcc .* conj.(r.grid_phase_current))
        @test snet ≈ complex(r.plant.p_poc, r.plant.q_poc) atol=1e-5
        @test imag(snet-sarm) ≈ b*sum(abs2, vpcc) atol=1e-5
        @test abs(snet) <= inv.s_max + .02
        # The terminal law compensates the signed shunt instead of reserving
        # its magnitude against output kVA. At zero requested Q, all capacity
        # is available for terminal active power.
        @test r.plant.p_poc ≈ inv.s_max atol=.03
        @test r.plant.q_poc ≈ 0. atol=.03
        @test r.exact_smooth_current_residual < .01
    end

end

@testset "PCC review: omitted current limits preserve loss imports at night" begin
    for target in (:converter, :grid)
        inv = _power_location_plant(power_limit_location=:poc, b_filter_shunt=.003,
            r_filter_grid=.1, x_filter_grid=.2, c_filter_mid=30e-6,
            p_loss_fixed=20., a_loss=.3, c_loss=.02)
        controlled = ControlledDevice(inv, SequenceController(AverageVoltageVoltVarWatt();
            current_target=target))
        request = InverterControlRequest(p_available=0., p_rated=10e3, q_scale=0.)
        r = solve_controlled_inverter(_power_location_grid(), controlled, request;
            solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        @test r.plant.p_dc ≈ 0. atol=.01
        @test r.plant.p_poc < -20.
        @test r.plant.p_dc ≈ r.plant.p_poc+r.plant.p_filter_loss+r.plant.p_loss+r.plant.p_cap_loss atol=1e-6
        @test r.exact_smooth_current_residual < .01
        @test r.plant.apparent_power_margin_VA >= -.02
    end
end

@testset "PCC review: an omitted current limiter has no numerical scale dependence" begin
    measurement = InverterControlMeasurement(ComplexF64[230, 230cis(-2pi/3), 230cis(2pi/3)])
    request = InverterControlRequest(p_available=0., p_rated=10e3, q_scale=0.)
    for limiter in (CommonScaleLimiter(), SequencePriorityLimiter(priority=:balance))
        controller = SequenceController(AverageVoltageVoltVarWatt(), NoUnbalanceControl(), limiter)
        baseline = evaluate_smooth(controller, measurement, request,
            InverterControlRatings(s_max=10e3); dc_power_offset=25.)
        for scale in (nextfloat(0.), 1e-3, 1e3, floatmax(Float64))
            result = evaluate_smooth(controller, measurement, request,
                InverterControlRatings(s_max=10e3, i_scale=scale); dc_power_offset=25.)
            @test result.current_scale == 1.
            @test result.phase_current == baseline.phase_current
            @test real(result.total_power) ≈ -25. atol=.001
        end
    end
end

@testset "PCC terminal controller: exact service allocation" begin
    voltage = Tuple(ComplexF64[230,230cis(-2pi/3),230cis(2pi/3)])
    # Deliberately large internal and external capacitive offsets. Neither is
    # an extra terminal VA budget; both must be compensated algebraically.
    shunt = -.02im .* voltage
    for target in (:converter,:grid), priority in (:var,:watt,:proportional), b in (-.08,.08), sign in (-1.,1.)
        inv = _power_location_plant(power_limit_location=:poc,b_filter_shunt=b)
        control = SequenceController(AverageVoltageVoltVarWatt(
            volt_var=PiecewiseLinearLaw([0.,400.],sign .* [.6,.6];smoothing_epsilon=.01)), NoUnbalanceControl(),
            CommonScaleLimiter(pq_priority=priority); current_target=target)
        request = InverterControlRequest(p_available=12e3,p_rated=12e3,q_scale=10e3)
        exact = PowerOptLab._evaluate_pcc_control_exact(control,request,inv,
            voltage,voltage,1.2 .* voltage,shunt,0.)
        grid = target === :converter ? exact.phase_current .- shunt : exact.phase_current
        net = sum(voltage .* conj.(grid)) + im*b*sum(abs2,voltage)
        expected = priority === :var ? 8000+sign*6000im : priority === :watt ?
            10000+0im : 10000*(12000+sign*6000im)/hypot(12000,6000)
        @test net ≈ expected atol=1e-7
        @test abs(net) ≈ 10000. atol=1e-7
    end
    # Exact segment geometry admits an inward direction even at the boundary.
    @test PowerOptLab._pcc_segment_exact((1.,0.),(-2.,0.),1.,0.) == 1.
    @test PowerOptLab._pcc_segment_exact((0.,.6),(.9,0.),1.,0.) ≈ 8/9
    @test_throws ArgumentError PowerOptLab._pcc_segment_exact((2.,0.),(0.,1.),1.,0.)
end

@testset "PCC terminal controller: priority through physical arm limits" begin
    for priority in (:var,:watt,:proportional), target in (:converter,:grid)
        inv = _power_location_plant(power_limit_location=:poc,
            i_grid_max=10., b_filter_shunt=.003,
            r_filter_grid=.1,x_filter_grid=.2,c_filter_mid=30e-6)
        controller = SequenceController(AverageVoltageVoltVarWatt(
            volt_var=PiecewiseLinearLaw([0.,400.],[.4,.4];smoothing_epsilon=.01)), NoUnbalanceControl(),
            CommonScaleLimiter(pq_priority=priority);current_target=target)
        request = InverterControlRequest(p_available=20e3,p_rated=20e3,q_scale=10e3)
        r = solve_controlled_inverter(_power_location_grid(),ControlledDevice(inv,controller),request;
            solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        if r.solve.publishable
            qsh = inv.b_filter_shunt*3*230^2
            arm_rating = 3*230*inv.i_grid_max
            if priority === :var
                @test r.plant.q_poc ≈ 4000. atol=.15
                @test r.plant.p_poc ≈ sqrt(arm_rating^2-(4000-qsh)^2) atol=.15
            elseif priority === :watt
                @test r.plant.p_poc ≈ sqrt(arm_rating^2-qsh^2) atol=.15
                # The remaining Q can move inward across the current circle:
                # preserve P and retain as much of the requested Q as possible.
                @test 0. <= r.plant.q_poc <= 4000.1
            else
                @test r.plant.q_poc/r.plant.p_poc ≈ r.control.q_request/r.control.p_request atol=1e-5
            end
            @test maximum(r.plant.i_grid_mag) ≈ 10. atol=1e-4
            @test r.exact_smooth_current_residual < .02
        end
    end
end

@testset "PCC terminal controller: full terminal circle despite filter power" begin
    for target in (:converter,:grid), priority in (:var,:watt,:proportional)
        inv = _power_location_plant(power_limit_location=:poc,b_filter_shunt=.08,
            r_filter_grid=.1,x_filter_grid=.2,c_filter_mid=30e-6)
        controller = SequenceController(AverageVoltageVoltVarWatt(
            volt_var=PiecewiseLinearLaw([0.,400.],[.6,.6];smoothing_epsilon=.01)),
            NoUnbalanceControl(),CommonScaleLimiter(pq_priority=priority);current_target=target)
        request = InverterControlRequest(p_available=20e3,p_rated=12e3,q_scale=10e3)
        # Rated watt basis makes the unconstrained terminal request 12 kW.
        controller = SequenceController(AverageVoltageVoltVarWatt(
            volt_var=controller.positive.volt_var,volt_watt_basis=:rated),
            NoUnbalanceControl(),controller.limiter;current_target=target)
        r = solve_controlled_inverter(_power_location_grid(),ControlledDevice(inv,controller),request;
            s_base=target === :converter ? 1e4 : 1e8,
            # At the P-axis, Q error is proportional to sqrt(P headroom).
            # Resolve that endpoint tightly rather than loosening P/Q assertions.
            solver_options=("tol"=>1e-11,"bound_relax_factor"=>0.,"max_iter"=>1000))
        @test r.solve.publishable
        if r.solve.publishable
            expected = priority === :var ? 8000+6000im : priority === :watt ?
                10000+0im : 10000*(12000+6000im)/hypot(12000,6000)
            @test complex(r.plant.p_poc,r.plant.q_poc) ≈ expected atol=.2
            @test r.exact_smooth_current_residual < .02
            @test solve_diagnostics(r).control_power_basis === :net_poc
            @test maximum(abs, r.converter_terminal.phase_current .- r.grid_phase_current) > 1.
        end
    end
end

@testset "PCC terminal controller: converter current preserves var service" begin
    v = Tuple(ComplexF64[230,230cis(-2pi/3),230cis(2pi/3)])
    shunt = -.005im .* v
    for target in (:converter,:grid)
        inv = _power_location_plant(power_limit_location=:poc,b_filter_shunt=.003,i_max=10.)
        controller = SequenceController(AverageVoltageVoltVarWatt(
            volt_var=PiecewiseLinearLaw([0.,400.],[.4,.4];smoothing_epsilon=.01)),
            NoUnbalanceControl(),CommonScaleLimiter(pq_priority=:var);current_target=target)
        r = PowerOptLab._evaluate_pcc_control_exact(controller,_POWER_LOCATION_REQUEST,inv,v,v,v,shunt,0.)
        grid = target === :converter ? r.phase_current .- shunt : r.phase_current
        conv = grid .+ shunt
        net = sum(v .* conj.(grid))+im*.003*sum(abs2,v)
        @test imag(net) ≈ 4000. atol=1e-7
        @test real(net) ≈ sqrt(6900^2-(4000-(.003-.005)*sum(abs2,v))^2) atol=1e-7
        @test maximum(abs,conv) ≈ 10. atol=1e-7
    end
end

@testset "PCC terminal controller: independent plant replay" begin
    inv = _power_location_plant(power_limit_location=:poc,b_filter_shunt=.003,
        r_filter_grid=.1,x_filter_grid=.2,c_filter_mid=30e-6)
    controller = SequenceController(AverageVoltageVoltVarWatt(
        volt_var=PiecewiseLinearLaw([0.,400.],[.4,.4];smoothing_epsilon=.01)),
        NoUnbalanceControl(),CommonScaleLimiter(pq_priority=:var))
    device = ControlledDevice(inv,controller)
    request = InverterControlRequest(p_available=5e3,p_rated=5e3,q_scale=10e3)
    result = solve_inverter_control_network_fixed_point(_power_location_grid(),device,request;
        max_iterations=20,solver_options=_POWER_LOCATION_OPTIONS)
    @test result.solve.publishable
    @test result.converged
    @test result.plant.q_poc ≈ 4000. atol=.1
    net = _power_location_grid()
    net["ibr"] = Dict("pv" => Dict{String,Any}(
        "bus"=>"poc", "terminal_map"=>["a","b","c","n"],
        "topology"=>"FOUR_LEG", "prime_mover"=>"PV", "s_max"=>ones(3),
        "p_min"=>zeros(3), "p_max"=>zeros(3), "q_min"=>zeros(3), "q_max"=>zeros(3)))
    fleet = ControlledInverterFleetSpec(Dict("pv"=>device),Dict("pv"=>request))
    replay = solve_controlled_inverter_fleet_network_fixed_point(net,fleet;
        max_iterations=20,solver_options=_POWER_LOCATION_OPTIONS)
    @test replay.solve.publishable
    @test replay.converged
    @test replay.plants["pv"].q_poc ≈ 4000. atol=.1

end

@testset "PCC terminal controller: LCL zero-sequence compensation" begin
    inv = _power_location_plant(power_limit_location=:poc,b_filter_shunt=-.003,
        r_filter_grid=.1,x_filter_grid=.2,c_filter_mid=30e-6)
    request = InverterControlRequest(p_available=5e3,p_rated=5e3,q_scale=0.)
    results = []
    for target in (:converter,:grid)
        controller = SequenceController(AverageVoltageVoltVarWatt();
            current_target=target,sensing_reference=:ground)
        r = solve_controlled_inverter(inv_grid3_floating_neutral(),ControlledDevice(inv,controller),request;
            solver_options=_POWER_LOCATION_OPTIONS)
        @test r.solve.publishable
        if r.solve.publishable
            @test abs(sum(r.converter_terminal.phase_current)) < 1e-5
            @test abs(sum(r.grid_phase_current)) > .01
            @test r.plant.q_poc ≈ 0. atol=.02
            @test r.plant.p_dc <= request.p_available+.02
            @test r.exact_smooth_current_residual < .01
            push!(results,complex(r.plant.p_poc,r.plant.q_poc))
        end
    end
    @test length(results) == 2
    length(results) == 2 && @test results[1] ≈ results[2] atol=.02
end

@testset "PCC terminal controller: operating-domain failures are explicit" begin
    v = Tuple(ComplexF64[230,230cis(-2pi/3),230cis(2pi/3)])
    z = (0.0im,0.0im,0.0im)
    controller = SequenceController(AverageVoltageVoltVarWatt())
    inv = _power_location_plant(power_limit_location=:poc,b_filter_shunt=.01,i_max=.1)
    @test_throws ArgumentError PowerOptLab._evaluate_pcc_control_exact(controller,
        _POWER_LOCATION_REQUEST,inv,v,v,v,z,0.)
    inv = _power_location_plant(power_limit_location=:poc)
    @test_throws ArgumentError PowerOptLab._evaluate_pcc_control_exact(controller,
        _POWER_LOCATION_REQUEST,inv,z,z,z,z,0.)
    night = InverterControlRequest(p_available=0.,p_rated=10e3,q_scale=0.)
    @test_throws ArgumentError PowerOptLab._evaluate_pcc_control_exact(controller,
        night,inv,v,v,v,z,11e3)
    @test PowerOptLab._pcc_segment_exact((0.,0.),(1.,0.),0.,0.) == 0.
end
