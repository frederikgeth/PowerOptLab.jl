# Regressions for the October 2026 advanced-IBR review
# (dev/notes/ibr_review_resolution.md). Each testset names its finding.

const _REVIEW_A = cis(2pi / 3)
_review_phases(u1::Complex, u2::Complex=0.0im) = ComplexF64[
    u1 + u2, _REVIEW_A^2*u1 + _REVIEW_A*u2, _REVIEW_A*u1 + _REVIEW_A^2*u2]

# Exact flat-clamped PWL evaluation, independent of the package helpers.
function _review_pwl(law, x)
    xs, ys = collect(law.breakpoints), collect(law.values)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    k = searchsortedlast(xs, x)
    return ys[k] + (ys[k+1] - ys[k])*(x - xs[k])/(xs[k+1] - xs[k])
end

_review_lossy_inverter(; kwargs...) = AdvancedInverter(;
    id="pv", bus="poc", phase_terminals=["a", "b", "c"], neutral=nothing,
    topology=:THREE_LEG, s_max=6000.0, i_max=10.0, v_dc=750.0, c_dc=2e-3,
    m_max=0.96, r_filter=0.05, x_filter=0.15,
    a_loss=1.0, c_loss=0.05, p_loss_fixed=20.0, kwargs...)

@testset "Review A1: exact DC-availability cap" begin
    measurement = InverterControlMeasurement(_review_phases(230.0 + 0im))
    ratings = InverterControlRatings(s_max=6e3, i_max=10.0)
    volt_watt = PiecewiseLinearLaw(
        [250.0, 260.0], [1.0, 0.0]; smoothing_epsilon=0.05)
    available = SequenceController(AverageVoltageVoltVarWatt(
        volt_watt=volt_watt, volt_watt_basis=:available))
    rated = SequenceController(AverageVoltageVoltVarWatt(
        volt_watt=volt_watt, volt_watt_basis=:rated))
    day = InverterControlRequest(p_available=5e3, p_rated=6e3, q_scale=0.0)
    night = InverterControlRequest(p_available=0.0, p_rated=6e3, q_scale=0.0)

    # The default offset is the lossless plant and leaves the command unchanged.
    @test evaluate_exact(available, measurement, day, ratings).p_request ≈ 5e3
    @test evaluate_exact(rated, measurement, day, ratings).p_request ≈ 5e3
    for controller in (available, rated)
        capped = evaluate_exact(controller, measurement, day, ratings;
                                dc_power_offset=57.0)
        @test capped.p_request ≈ 5e3 - 57.0
        importing = evaluate_exact(controller, measurement, night, ratings;
                                   dc_power_offset=20.0)
        @test importing.p_request ≈ -20.0
        smooth = evaluate_smooth(controller, measurement, day, ratings;
                                 dc_power_offset=57.0)
        @test smooth.p_request <= capped.p_request + 1e-9
        @test smooth.p_request ≈ capped.p_request atol=1e-3
    end

    # Volt-watt curtails positive availability but never the loss import.
    high = InverterControlMeasurement(_review_phases(255.0 + 0im))
    @test evaluate_exact(available, high, day, ratings;
                         dc_power_offset=57.0).p_request ≈ 0.5(5e3 - 57.0)
    @test evaluate_exact(available, high, night, ratings;
                         dc_power_offset=20.0).p_request ≈ -20.0
    @test_throws ArgumentError evaluate_exact(
        available, measurement, day, ratings; dc_power_offset=NaN)
end

@testset "Review A1: plant never overdraws DC availability" begin
    inverter = _review_lossy_inverter()
    for target in (:converter, :grid), p_available in (5e3, 0.0)
        controller = SequenceController(
            AverageVoltageVoltVarWatt(); current_target=target)
        request = InverterControlRequest(
            p_available=p_available, p_rated=6e3, q_scale=6e3)
        result = solve_controlled_inverter(
            inv_grid3_bal(), ControlledDevice(inverter, controller), request;
            solver_options=("tol" => 1e-8,))
        @test result.solve.publishable
        @test result.plant.p_dc <= p_available + 1e-3
        @test result.exact_smooth_current_residual < 0.01
        if iszero(p_available)
            # A zero-irradiance converter imports exactly its losses.
            @test result.plant.p_dc ≈ 0.0 atol=1e-3
            @test result.plant.p_poc < -20.0
        else
            @test result.plant.p_dc ≈ p_available atol=1.0
        end
    end

    # A filter-only offset equals the series-resistance loss of balanced current.
    lossless = _review_lossy_inverter(a_loss=0.0, c_loss=0.0, p_loss_fixed=0.0)
    result = solve_controlled_inverter(
        inv_grid3_bal(),
        ControlledDevice(lossless, SequenceController(AverageVoltageVoltVarWatt())),
        InverterControlRequest(p_available=5e3, p_rated=6e3, q_scale=6e3);
        solver_options=("tol" => 1e-8,))
    offset = dc_power_offset(
        result.plant.p_dc, result.control.phase_voltage,
        result.converter_terminal.phase_current)
    @test offset ≈ 3*0.05*abs2(result.converter_terminal.phase_current[1]) rtol=1e-4
end

@testset "Review A2: exact sequence-priority fill" begin
    limit = 10.0
    for (a, b) in ((3.0 + 4.0im, 8.0 + 0.0im), (6.0im, -9.0 + 2.0im),
                   (2.0 + 0.0im, 0.0 + 0.0im), (9.99 + 0.0im, 1.0im))
        t = PowerOptLab._priority_fill_exact(a, b, limit)
        @test 0.0 <= t <= 1.0
        @test abs(a + t*b) <= limit + 1e-10
        t < 1.0 && @test abs(a + t*b) ≈ limit rtol=1e-10
    end
    @test PowerOptLab._priority_fill_exact(10.0 + 0im, 1.0 + 0im, 10.0) == 0.0
    @test_throws ArgumentError SequencePriorityLimiter(priority=:neutral)
end

@testset "Review A2: power priority preserves the positive-sequence command" begin
    gain = PiecewiseLinearLaw([0.0, 0.002, 0.02], [0.0, 0.0, 3.0];
                              smoothing_epsilon=1e-4)
    unbalance = NegativeSequenceAdmittanceDroop(gain)
    positive = AverageVoltageVoltVarWatt()
    common = SequenceController(positive, unbalance, CommonScaleLimiter())
    power = SequenceController(
        positive, unbalance, SequencePriorityLimiter(priority=:power))
    balance = SequenceController(
        positive, unbalance, SequencePriorityLimiter(priority=:balance))
    measurement = InverterControlMeasurement(
        _review_phases(230.0 + 0im, 3.0cis(0.7)))
    request = InverterControlRequest(p_available=5e3, p_rated=6e3, q_scale=0.0)
    ratings = InverterControlRatings(s_max=6e3, i_max=10.0)

    exact_common = evaluate_exact(common, measurement, request, ratings)
    exact_power = evaluate_exact(power, measurement, request, ratings)
    exact_balance = evaluate_exact(balance, measurement, request, ratings)
    # The negative-sequence request alone would breach the leg rating.
    @test exact_common.current_scale < 1.0
    @test exact_power.current_scale == 1.0
    @test exact_power.i1 ≈ exact_power.i1_request
    @test abs(exact_power.i2) < abs(exact_power.i2_request)
    @test real(exact_power.sequence_power) > real(exact_common.sequence_power)
    @test abs(exact_balance.i2) ≈ abs(exact_balance.i2_request)
    @test abs(exact_balance.i1) < abs(exact_power.i1)
    for result in (exact_common, exact_power, exact_balance)
        @test maximum(abs, result.phase_current) <= ratings.i_max + 1e-9
    end
    # At least one phase is exactly at the rating when the fill binds.
    @test maximum(abs, exact_power.phase_current) ≈ ratings.i_max rtol=1e-9

    for controller in (power, balance)
        exact = evaluate_exact(controller, measurement, request, ratings)
        smooth = evaluate_smooth(controller, measurement, request, ratings)
        @test maximum(abs, smooth.phase_current) <= ratings.i_max + 1e-6
        @test maximum(abs, smooth.phase_current .- exact.phase_current) < 5e-3
    end

    # Balanced voltage: no negative-sequence request, so every limiter agrees.
    balanced = InverterControlMeasurement(_review_phases(230.0 + 0im))
    for controller in (common, power, balance)
        @test evaluate_exact(controller, balanced, request, ratings).i1 ≈
              evaluate_exact(common, balanced, request, ratings).i1
        @test evaluate_smooth(controller, balanced, request, ratings).i1 ≈
              evaluate_exact(common, balanced, request, ratings).i1 rtol=1e-3
    end
end

@testset "Review A2: power priority in the network solve" begin
    gain = PiecewiseLinearLaw([0.0, 0.002, 0.02], [0.0, 0.0, 3.0];
                              smoothing_epsilon=1e-4)
    unbalance = NegativeSequenceAdmittanceDroop(gain)
    inverter = _review_lossy_inverter(a_loss=0.0, c_loss=0.0, p_loss_fixed=0.0)
    request = InverterControlRequest(p_available=5e3, p_rated=6e3, q_scale=0.0)
    results = Dict{Symbol,Any}()
    for (name, limiter) in ((:common, CommonScaleLimiter()),
                            (:power, SequencePriorityLimiter(priority=:power)))
        controller = SequenceController(
            AverageVoltageVoltVarWatt(), unbalance, limiter)
        result = solve_controlled_inverter(
            inv_grid3_unbal(), ControlledDevice(inverter, controller), request;
            solver_options=("tol" => 1e-8,))
        @test result.solve.publishable
        @test result.exact_smooth_current_residual < 0.02
        @test maximum(result.plant.i_mag) <= inverter.i_max + 1e-5
        results[name] = result
    end
    @test results[:common].control.current_scale < 1.0
    @test results[:power].plant.p_poc > results[:common].plant.p_poc + 100.0
end

@testset "Review A5: dominant rule is the standard curve without conflict" begin
    volt_var = PiecewiseLinearLaw(
        [207.0, 220.0, 240.0, 258.0], [0.44, 0.0, 0.0, -0.60];
        smoothing_epsilon=0.05)
    request = InverterControlRequest(p_available=0.0, p_rated=6e3, q_scale=6e3)
    ratings = InverterControlRatings(s_max=6e3, i_max=10.0)
    for epsilon in (0.01, 0.05)
        dominant = SequenceController(WorstPhaseVoltVarWatt(
            volt_var=volt_var, conflict_policy=:dominant,
            conflict_epsilon=epsilon))
        net = SequenceController(WorstPhaseVoltVarWatt(
            volt_var=volt_var, conflict_policy=:net))
        # Single-branch requests at, near and beyond the deadband edges match
        # the plain curve exactly; the former rule lost up to ~0.15ε here.
        for voltage in (219.99, 219.9, 219.0, 215.0, 240.01, 240.2, 245.0, 230.0)
            measurement = InverterControlMeasurement(_review_phases(voltage + 0im))
            expected = request.q_scale*_review_pwl(volt_var, voltage)
            @test evaluate_exact(dominant, measurement, request, ratings).q_request ≈
                  expected atol=1e-9
            @test evaluate_exact(net, measurement, request, ratings).q_request ≈
                  expected atol=1e-9
        end
    end

    # Conflict between a low and a high phase. Sweeping the unbalance angle
    # crosses the equal-severity tie: the response changes sign there,
    # continuously, and each side tends to its own branch.
    dominant = SequenceController(WorstPhaseVoltVarWatt(
        volt_var=volt_var, conflict_policy=:dominant, conflict_epsilon=0.01))
    q_at(theta) = evaluate_exact(dominant,
        InverterControlMeasurement(_review_phases(230.0 + 0im, 22.0cis(theta))),
        request, ratings).q_request
    values = [q_at(deg2rad(theta)) for theta in range(15.0, 30.0; length=301)]
    @test values[1] < -0.5*0.60*request.q_scale    # high branch dominates
    @test values[end] > 0.5*0.44*request.q_scale   # low branch dominates
    # Continuity, not a jump: refining the angle step near the tie shrinks the
    # largest increment proportionally (a discontinuity would not shrink).
    crossing = findfirst(k -> values[k] <= 0 < values[k+1], 1:length(values)-1)
    @test crossing !== nothing
    tie = deg2rad(range(15.0, 30.0; length=301)[crossing])
    jump(step) = maximum(abs, diff([q_at(tie + x)
        for x in range(-50step, 50step; length=101)]))
    @test jump(deg2rad(1e-3)) < 0.6*jump(deg2rad(2e-3))
    @test jump(deg2rad(1e-4)) < 0.01*request.q_scale
end

@testset "Review A4: sensing reference is part of the law, not the wiring" begin
    net = inv_grid3_floating_neutral()
    three_wire = AdvancedInverter(
        id="pv", bus="poc", phase_terminals=["a", "b", "c"], neutral=nothing,
        topology=:THREE_LEG, s_max=6e3, i_max=10.0, v_dc=750.0, c_dc=2e-3,
        r_filter=0.05, x_filter=0.15)
    volt_var = PiecewiseLinearLaw(
        [207.0, 220.0, 240.0, 258.0], [0.44, 0.0, 0.0, -0.60];
        smoothing_epsilon=0.05)
    request = InverterControlRequest(p_available=2e3, p_rated=6e3, q_scale=6e3)
    @test_throws ArgumentError SequenceController(
        AverageVoltageVoltVarWatt(); sensing_reference=:neutral)
    @test_throws ArgumentError SequenceController(
        AverageVoltageVoltVarWatt(); sensing_reference="")
    @test_throws ArgumentError solve_controlled_inverter(net,
        ControlledDevice(three_wire, SequenceController(
            AverageVoltageVoltVarWatt(); sensing_reference="x")), request)
    @test_throws ArgumentError solve_controlled_inverter(net,
        ControlledDevice(three_wire, SequenceController(
            AverageVoltageVoltVarWatt(); sensing_reference="a")), request)

    results = Dict{Any,Any}()
    for reference in (:device, :ground, "n")
        controller = SequenceController(
            AverageVoltageVoltVarWatt(volt_var=volt_var);
            sensing_reference=reference)
        @test sensing_terminal(controller, three_wire) ==
              (reference === :device || reference === :ground ? nothing : "n")
        result = solve_controlled_inverter(
            net, ControlledDevice(three_wire, controller), request;
            solver_options=("tol" => 1e-8,))
        @test result.solve.publishable
        @test result.exact_smooth_current_residual < 0.01
        results[reference] = result
    end
    # The sensed phasors are exactly the bus voltages referred to the declared
    # terminal; the neutral is floating, so the two references differ.
    bus = results["n"].bus["poc"]
    phasor(t) = ComplexF64(bus[t]["vr"], bus[t]["vi"])
    @test collect(results["n"].control.phase_voltage) ≈
          [phasor(p) - phasor("n") for p in ("a", "b", "c")] atol=1e-6
    @test abs(phasor("n")) > 0.5
    @test results[:device].control.phase_voltage ==
          results[:ground].control.phase_voltage
    # Metamorphic contract in the active Volt-var region: a 5 V common-mode
    # offset leaves the sequence law unchanged but moves the worst-phase law
    # to first order. (The mean-phase law moves only at second order: the
    # offset's projections onto the three phases nearly cancel.)
    pn = _review_phases(214.0 + 0im, 4.0cis(0.3))
    pg = pn .+ 5.0cis(0.2)
    ratings = InverterControlRatings(s_max=6e3, i_max=10.0)
    q(law, phases) = evaluate_exact(law, InverterControlMeasurement(phases),
                                    request, ratings).q_request
    sequence_law = SequenceController(PositiveSequenceVoltVarWatt(
        volt_var=volt_var, worst_phase_watt_guard=false))
    worst_law = SequenceController(WorstPhaseVoltVarWatt(volt_var=volt_var))
    mean_law = SequenceController(AverageVoltageVoltVarWatt(volt_var=volt_var))
    @test q(sequence_law, pn) > 100.0
    @test q(sequence_law, pg) ≈ q(sequence_law, pn) rtol=1e-12
    @test abs(q(worst_law, pg) - q(worst_law, pn)) > 300.0
    @test abs(q(mean_law, pg) - q(mean_law, pn)) <
          0.1abs(q(worst_law, pg) - q(worst_law, pn))
end

@testset "Review A3: one DC-source model for the split link" begin
    # Independent oracle: nodal analysis of the averaged three-node DC circuit
    # (rails P, M, N) at the fundamental. The bridge draws J/2 from P and
    # returns J/2 to N... equivalently each rail supplies -J/2 and J returns
    # into M; the source Z_s connects P and N.
    function circuit(cu, cl, zs, j; omega=2pi*50.0)
        yu, yl = im*omega*cu, im*omega*cl
        ys = zs === nothing ? 0.0im : 1/zs
        # Unknowns v_P, v_M with v_N = 0. KCL (currents leaving each node):
        #   P: yu(vP-vM) + ys vP + J/2 = 0
        #   M: yu(vM-vP) + yl vM - J = 0
        A = [yu+ys -yu; -yu yu+yl]
        vP, vM = A \ [-j/2, j]
        midpoint = vM - vP/2          # displacement from the bus centre
        return midpoint, vP           # vP is the total-bus ripple
    end
    j = 7.0cis(0.4)
    for (cu, cl) in ((2.5e-3, 3.5e-3), (3e-3, 3e-3)),
            source in (nothing, 0.0, 0.2, complex(0.05, 2pi*50*1e-3))
        inv = AdvancedInverter(id="s", bus="poc", s_max=20e3,
            phase_terminals=["a", "b", "c"], neutral="n", topology=:SPLIT_DC,
            v_dc=900.0, c_dc=3e-3, c_dc_upper=cu, c_dc_lower=cl, In_max=40.0,
            dc_source_r=source === nothing ? nothing : real(source),
            dc_source_l=source === nothing ? 0.0 : imag(source)/(2pi*50))
        k = PowerOptLab._split_fundamental_coefficients(inv)
        zs = source === nothing ? nothing :
            (iszero(source) ? 1e-12 + 0im : complex(source))
        midpoint, bus = circuit(cu, cl, zs, j)
        @test k.midpoint*j ≈ midpoint rtol=1e-8
        @test abs(k.bus*j - bus) <= 1e-8*abs(midpoint)
        if cu == cl
            # Equal banks: unchanged by the source model.
            @test k.midpoint ≈ 1/(im*2pi*50*(cu + cl)) rtol=1e-12
            @test abs(k.bus) <= 1e-12*abs(k.midpoint)
            @test k.upper_weight ≈ 0.25 && k.lower_weight ≈ 0.25
        end
    end
    @test_throws ArgumentError AdvancedInverter(id="s", bus="poc", s_max=1e3,
        dc_source_l=1e-3) |> PowerOptLab._validate_inverter
end

@testset "Review C2: grid-forming virtual-impedance comparator" begin
    r, x = 0.1, 0.3
    policy = negative_sequence_virtual_impedance(r, x)
    @test policy.impedance_angle ≈ atan(x, r)
    @test_throws ArgumentError negative_sequence_virtual_impedance(0.0, 0.0)
    @test_throws ArgumentError negative_sequence_virtual_impedance(-0.1, 0.3)

    # Exact law: I2 = -U2/Z_v whenever the current limit does not bind.
    controller = SequenceController(AverageVoltageVoltVarWatt(), policy,
        SequencePriorityLimiter(priority=:power))
    u2 = 2.0cis(-0.7)
    exact = evaluate_exact(controller,
        InverterControlMeasurement(_review_phases(230.0 + 0im, u2)),
        InverterControlRequest(p_available=0.0, p_rated=6e3, q_scale=0.0),
        InverterControlRatings(s_max=20e3, i_max=40.0))
    @test exact.i2 ≈ -u2/complex(r, x) rtol=1e-9

    # Plant oracle: a three-leg grid-forming plant with a balanced internal EMF
    # behind the same series filter draws the same negative-sequence current.
    # A mild source unbalance keeps |I2| = |U2|/|Z_v| inside generous ratings;
    # the strongly unbalanced fixture would demand more than the plant can
    # carry, which is exactly why this comparator needs current limiting.
    plant = AdvancedInverter(id="gfm", bus="poc", phase_terminals=["a", "b", "c"],
        neutral=nothing, topology=:THREE_LEG, s_max=50e3, i_max=100.0,
        v_dc=800.0, c_dc=2e-3, r_filter=r, x_filter=x, grid_forming=true,
        v_int_min=200.0, v_int_max=260.0)
    mild = inv_grid3_src(mags=[232.0, 228.0, 230.0], angs=[0.0, -2.0944, 2.0944])
    result = solve_advanced_inverter(mild, plant; objective=:max_export)
    @test result.termination_status in ("LOCALLY_SOLVED", "OPTIMAL")
    bus = result.bus["poc"]
    v = [complex(bus[p]["vr"], bus[p]["vi"]) for p in ("a", "b", "c")]
    poc_u2 = (v[1] + _REVIEW_A^2*v[2] + _REVIEW_A*v[3])/3
    @test result.i_negative ≈ abs(poc_u2)/abs(complex(r, x)) rtol=1e-4
end
