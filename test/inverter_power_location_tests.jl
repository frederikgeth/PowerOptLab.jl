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
    # Exact boundary replay must use physical POC voltage, not converter voltage
    # or the potentially different controller sensing reference. With a fixed
    # shunt offset the conservative allocator has an independent scalar oracle.
    inv = _power_location_plant(power_limit_location=:poc, b_filter_shunt=.01)
    v = Tuple(ComplexF64[230, 230cis(-2pi/3), 230cis(2pi/3)])
    vint = 1.3 .* v
    direction = v ./ 10
    zero = (0.0im, 0.0im, 0.0im)
    scale = PowerOptLab._plant_scale_exact(vint, zero, zero, direction, inv; poc_voltage=v)
    expected = (inv.s_max - inv.b_filter_shunt*sum(abs2, v)) / abs(sum(v .* conj.(direction)))
    @test scale ≈ expected atol=1e-12
    @test_throws ArgumentError PowerOptLab._plant_scale_exact(vint, zero, zero, direction, inv)
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
        # No target-current limiter: independently replay the conservative
        # physical PCC circle, using the unscaled firmware command.
        offset = dc_power_offset(r.plant.p_dc, r.control.phase_voltage,
            target === :converter ? r.converter_terminal.phase_current : r.grid_phase_current)
        raw = evaluate_exact(controller, InverterControlMeasurement(collect(ground_voltage)),
            _POWER_LOCATION_REQUEST, InverterControlRatings(inv, controller.current_target);
            dc_power_offset=offset)
        @test raw.current_scale == 1.
        direction = sum(vpcc .* conj.(raw.phase_current))
        expected = min(1., (inv.s_max - abs(b*sum(abs2, vpcc)))/abs(direction))
        @test r.exact_control.current_scale ≈ expected atol=1e-10
        @test r.exact_smooth_current_residual < .01
    end
    # A reactive arm-current offset makes the sign of B matter to backoff,
    # not just to reported Q. An unrelated internal voltage must not enter it.
    vpcc = Tuple(ComplexF64[240, 220cis(-2pi/3), 230cis(2pi/3)])
    offset_g = .01im .* vpcc
    direction_i = vpcc ./ 10
    for b in (-.003, .003)
        inv = _power_location_plant(power_limit_location=:poc, b_filter_shunt=b)
        offset_s = sum(vpcc .* conj.(offset_g .- im*b .* vpcc))
        direction_s = sum(vpcc .* conj.(direction_i))
        expected = (inv.s_max-abs(offset_s))/abs(direction_s)
        actual = PowerOptLab._plant_scale_exact(1.5 .* vpcc, 2 .* offset_g,
            offset_g, direction_i, inv; poc_voltage=vpcc)
        @test actual ≈ expected atol=1e-12
        @test abs(offset_s + actual*direction_s) <= inv.s_max
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
