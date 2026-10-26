# Run with julia --project=. scripts/formulations/ibr_reproductions.jl
# These are comparisons of algebraic steady-state models, not dynamic tests.
using PowerOptLab, JuMP, Ipopt
import BMOPFTools

_repro_phases(i1, i2=0im) = (i1+i2, cis(-2pi/3)*i1+cis(2pi/3)*i2,
    cis(2pi/3)*i1+cis(-2pi/3)*i2)

function _repro_grid(voltage; feeder=false)
    u = voltage isa Real ? _repro_phases(complex(voltage)) : voltage
    # Explicit phasors avoid source-angle truncation introducing spurious U0.
    bus = feeder ? "grid" : "poc"
    extra = feeder ? """,
      "grid":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"]}""" : ""
    lines = feeder ? """,
      "linecode":{"lc":{"R_series_1_1":0.05,"R_series_2_2":0.05,"R_series_3_3":0.05}},
      "line":{"l":{"bus_from":"grid","bus_to":"poc","terminal_map_from":["a","b","c"],
        "terminal_map_to":["a","b","c"],"linecode":"lc","length":1.0}}""" : ""
    BMOPFTools.parse_bmopf("""
    {"bus":{"poc":{"terminal_names":["a","b","c","n"],"perfectly_grounded_terminals":["n"],
      "v_min":[180,180,180],"v_max":[280,280,280]}$extra},
     "voltage_source":{"vs":{"bus":"$bus","terminal_map":["a","b","c"],
       "v_magnitude":$(collect(abs.(u))),"v_angle":$(collect(angle.(u)))}}$lines}
    """; from_string=true)
end

function _repro_plant(; kwargs...)
    AdvancedInverter(; id="reproduction", bus="poc", topology=:THREE_LEG,
        phase_terminals=["a","b","c"], neutral=nothing, s_max=6000.,
        i_max=1.1*6000/(3*230), v_dc=800., c_dc=.002, kwargs...)
end

function _repro_policy(; constant_q=false)
    vv = constant_q ? PWLFunction([180.,280.],[.44,.44]) :
        PWLFunction([207.,220.,240.,258.],[.44,0.,0.,-.60])
    lower_positive_policy(VoltVarWattIntent(volt_var=vv,
        volt_watt=PWLFunction([253.,260.],[1.,.2]), sensing=:average_voltage,
        volt_watt_basis=:rated), VoltVarWattEncoding(
        volt_var=LocalC2Formulation(.01), volt_watt=LocalC2Formulation(.01)))
end

# A ray is an extra controller-specific restriction. The general feasible-set
# model does not impose it. Fix two phases only: THREE_LEG supplies zero sum.
function _repro_ray_solve(net, plant, request, reserved, direction)
    device = ConstrainedInverter(plant, InverterOperatingConstraints();
        current_allocation=:free, dc_source=DCSourceBounds(p_min=0.))
    handles = Ref{Any}()
    scale = Ref{VariableRef}()
    hook! = ctx -> begin
        handles[] = stamp_device!(ctx, device; request)
        h = inverter_handles(device, handles[])
        m = BMOPFTools.opf_model(ctx)
        scale[] = @variable(m, lower_bound=0., upper_bound=1., start=.5)
        for k in 1:2
            @constraint(m, h.cri[k]*h.sb/h.vb == real(reserved[k])+scale[]*real(direction[k]))
            @constraint(m, h.cii[k]*h.sb/h.vb == imag(reserved[k])+scale[]*imag(direction[k]))
        end
        @objective(m, Max, scale[])
    end
    ctx = BMOPFTools.build_opf_model(net; per_unit=true, s_base=1e6,
        add_objective=false, model_hook! = hook!, optimizer=Ipopt.Optimizer)
    m = BMOPFTools.opf_model(ctx)
    set_optimizer_attribute(m,"tol",1e-9)
    BMOPFTools.optimize_opf!(ctx)
    status = solve_status(m)
    result = extract_device(device, handles[], status)
    return (; status, result, scale=status.publishable ? value(scale[]) : NaN)
end

"""
    ibr_formulation_reproductions()

Independent feeder solves, proportional-current rays, and fixed-measurement
sequence-priority reproductions. Each row retains both results and physical
residuals. A successful solver status alone does not establish equivalence.
"""
function ibr_formulation_reproductions()
    feeder = map((215.,230.,259.)) do voltage
        net = _repro_grid(voltage; feeder=true)
        plant = _repro_plant(i_max=15., r_filter=.05, x_filter=.15, p_loss_fixed=20.)
        positive = _repro_policy()
        request = InverterControlRequest(p_available=3000.,p_rated=6000.,q_scale=6000.)
        controller = SequenceController(positive; limiter=CommonScaleLimiter(pq_priority=:var),
            power_voltage_floor=1e-3)
        forward = solve_controlled_inverter(net,ControlledDevice(plant,controller),request;
            solver_options=("tol"=>1e-9,))
        constrained = solve_constrained_inverter(net,ConstrainedInverter(plant,
            VoltVarWattConstraints(positive);dc_source=DCSourceBounds(p_min=0.)),request;
            solver_options=("tol"=>1e-9,))
        (; voltage, forward, constrained,
            p_error_W=abs(forward.plant.p_poc-constrained.plant.p_poc),
            q_error_var=abs(forward.plant.q_poc-constrained.plant.q_poc),
            current_error_A=maximum(abs.(collect(forward.converter_terminal.phase_current).-collect(constrained.converter_terminal.phase_current))),
            voltage_error_V=maximum(abs.(collect(forward.control.phase_voltage).-collect(constrained.operation.phase_voltage))))
    end
    proportional = map((230.,207.)) do voltage
        plant, positive = _repro_plant(), _repro_policy(constant_q=true)
        request = InverterControlRequest(p_available=6000.,p_rated=6000.,q_scale=6000.)
        controller = SequenceController(positive;limiter=CommonScaleLimiter(pq_priority=:proportional),
            power_voltage_floor=1e-3)
        forward = solve_controlled_inverter(_repro_grid(voltage),ControlledDevice(plant,controller),request;
            solver_options=("tol"=>1e-9,))
        direction = _repro_phases(complex(6000.,-2640.)/(3voltage))
        optimisation = _repro_ray_solve(_repro_grid(voltage),plant,request,(0im,0im,0im),direction)
        analytic_scale = min(1.,plant.s_max/hypot(6000.,2640.),
            3voltage*plant.i_max/hypot(6000.,2640.))
        (;voltage,forward,optimisation,analytic_scale,
            p_error_W=abs(forward.plant.p_poc-optimisation.result.plant.p_poc),
            q_error_var=abs(forward.plant.q_poc-optimisation.result.plant.q_poc))
    end
    sequence_priority = map((:power,:balance)) do priority
        voltage = _repro_phases(230.0+0im,10cis(.2))
        measurement = InverterControlMeasurement(collect(voltage))
        plant = _repro_plant(s_max=20000.,i_max=6.)
        gain = PiecewiseLinearLaw([0.,1.],[.3,.3];formulation=LocalC2Formulation(.001))
        negative = NegativeSequenceAdmittanceDroop(gain;impedance_angle=pi/2,ripple_blend=0.)
        controller = SequenceController(AverageVoltageVoltVarWatt(),negative,
            SequencePriorityLimiter(;priority);power_voltage_floor=1e-3)
        request = InverterControlRequest(p_available=3000.,p_rated=6000.,q_scale=0.)
        exact = evaluate_exact(controller,measurement,request,InverterControlRatings(plant))
        reserved, direction = priority === :power ?
            (_repro_phases(exact.i1_request),_repro_phases(0im,exact.i2_request)) :
            (_repro_phases(0im,exact.i2_request),_repro_phases(exact.i1_request))
        optimisation = _repro_ray_solve(_repro_grid(voltage),plant,request,reserved,direction)
        # Independent intersection of the three phase-current disks along a ray.
        phase_scales = map(reserved,direction) do r,d
            c = real(r*conj(d))
            (-c+sqrt(c^2+abs2(d)*(plant.i_max^2-abs2(r))))/abs2(d)
        end
        analytic_scale = min(1.,phase_scales...)
        (;priority,exact,optimisation,analytic_scale,
            current_error_A=maximum(abs.(collect(exact.phase_current).-collect(optimisation.result.converter_terminal.phase_current))))
    end
    return (;feeder,proportional,sequence_priority)
end

function _repro_errors(forward, constrained)
    (; p_error_W=abs(forward.plant.p_poc-constrained.plant.p_poc),
       q_error_var=abs(forward.plant.q_poc-constrained.plant.q_poc),
       current_error_A=maximum(abs.(collect(forward.converter_terminal.phase_current).-collect(constrained.converter_terminal.phase_current))),
       voltage_error_V=maximum(abs.(collect(forward.control.phase_voltage).-collect(constrained.operation.phase_voltage))))
end

"""Balanced lossless tutorial: four operating choices and closed-form oracles.
The same 2640 var request is held at both voltages to isolate capability geometry.
"""
function ibr_pq_tutorial()
    plant = _repro_plant()
    request = InverterControlRequest(p_available=6000.,p_rated=6000.,q_scale=6000.)
    rows = NamedTuple[]
    for voltage in (230.,207.)
        net = _repro_grid(voltage)
        capacity = min(plant.s_max,3voltage*plant.i_max)
        scale = capacity/hypot(6000.,2640.)
        direction = _repro_phases(complex(6000.,-2640.)/(3voltage))
        ray = _repro_ray_solve(net,plant,request,(0im,0im,0im),direction)
        push!(rows,(;voltage,mode=:proportional,capacity_VA=capacity,
            expected_p_W=6000scale,expected_q_var=2640scale,
            result=(;ray.result...,solve=ray.status)))
        choices = (
            (:strict_q,VoltVarWattConstraints(_repro_policy(constant_q=true)),
                sqrt(capacity^2-2640^2),2640.),
            (:fixed_pf,InverterOperatingConstraints(reactive=FixedPowerFactor(-.95)),
                .95capacity,sqrt(1-.95^2)*capacity),
            (:free_q,InverterOperatingConstraints(),capacity,0.))
        for (mode,policy,expected_p_W,expected_q_var) in choices
            result = solve_constrained_inverter(net,ConstrainedInverter(plant,policy;
                dc_source=DCSourceBounds(p_min=0.)),request;solver_options=("tol"=>1e-9,))
            push!(rows,(;voltage,mode,capacity_VA=capacity,expected_p_W,expected_q_var,result))
        end
    end
    return rows
end

# Reproduce a converter-current target's Q definition using the general free-Q
# policy. This is POC measurement voltage times CONVERTER current, not internal
# converter power and not total POC exchange. The source budget includes losses.
function _repro_converter_q_solve(net, plant, request, q_command)
    device = ConstrainedInverter(plant,InverterOperatingConstraints();
        dc_source=DCSourceBounds(p_min=0.))
    handles = Ref{Any}()
    hook! = ctx -> begin
        handles[] = stamp_device!(ctx,device;request)
        h = inverter_handles(device,handles[])
        u = handles[].phase_voltage
        m = BMOPFTools.opf_model(ctx)
        q_at_target = @expression(m,sum(u[k][2]*h.cri[k]-u[k][1]*h.cii[k] for k in 1:3))
        @constraint(m,q_at_target/(plant.s_max/h.sb) == q_command/plant.s_max)
        @objective(m,Max,h.p_poc/(plant.s_max/h.sb))
    end
    ctx = BMOPFTools.build_opf_model(net;per_unit=true,s_base=1e6,
        add_objective=false,model_hook! = hook!,optimizer=Ipopt.Optimizer)
    m = BMOPFTools.opf_model(ctx)
    set_optimizer_attribute(m,"tol",1e-9)
    BMOPFTools.optimize_opf!(ctx)
    status = solve_status(m)
    extracted = extract_device(device,handles[],status)
    return (;extracted...,solve=status)
end

"""
    ibr_network_reproductions()

Unbalanced independent feeder comparisons, LCL target-location matching, binding
current limits and two kinds of infeasibility. `expected` declares the scientific
comparison; tests also check physical powers, phasors and analytical obstructions.
"""
function ibr_network_reproductions()
    request = InverterControlRequest(p_available=3000.,p_rated=6000.,q_scale=6000.)
    q_policy = AverageVoltageVoltVarWatt(volt_var=PiecewiseLinearLaw(
        [180.,280.],[.1,.1];formulation=LocalC2Formulation(.01)),volt_watt_basis=:rated)
    unbalanced_voltage = _repro_phases(230.0+0im,8cis(.2))
    net = _repro_grid(unbalanced_voltage;feeder=true)
    pairs = NamedTuple[]
    for target in (:converter,:grid)
        plant = _repro_plant(i_max=15.,r_filter=.05,x_filter=.15,
            r_filter_grid=.03,x_filter_grid=.10,c_filter_mid=30e-6,r_filter_damping=.5,
            p_loss_fixed=20.)
        controller = SequenceController(q_policy;current_target=target,
            limiter=CommonScaleLimiter(pq_priority=:var),power_voltage_floor=1e-3)
        forward = solve_controlled_inverter(net,ControlledDevice(plant,controller),request;
            solver_options=("tol"=>1e-9,))
        allocation = target === :converter ? :balanced_converter : :balanced_grid
        strict = solve_constrained_inverter(net,ConstrainedInverter(plant,
            VoltVarWattConstraints(q_policy);current_allocation=allocation,
            dc_source=DCSourceBounds(p_min=0.)),request;solver_options=("tol"=>1e-9,))
        matched = target === :converter ? _repro_converter_q_solve(net,plant,request,600.) : strict
        push!(pairs,(;id=Symbol(:lcl_,target),expected=:agreement,target,plant,request,
            forward,constrained=matched,_repro_errors(forward,matched)...))
        if target === :converter
            push!(pairs,(;id=:lcl_poc_q_mismatch,expected=:different_q_location,target,plant,request,
                forward,constrained=strict,_repro_errors(forward,strict)...))
        end
    end
    # Reduced filter with an unbalanced feeder: no midpoint current to change
    # the Q location, so a balanced converter command and strict POC Q agree.
    plant = _repro_plant(i_max=15.,r_filter=.05,x_filter=.15,p_loss_fixed=20.)
    controller = SequenceController(q_policy;power_voltage_floor=1e-3)
    forward = solve_controlled_inverter(net,ControlledDevice(plant,controller),request;
        solver_options=("tol"=>1e-9,))
    constrained = solve_constrained_inverter(net,ConstrainedInverter(plant,
        VoltVarWattConstraints(q_policy);dc_source=DCSourceBounds(p_min=0.)),request;
        solver_options=("tol"=>1e-9,))
    push!(pairs,(;id=:unbalanced_reduced,expected=:agreement,target=:converter,plant,request,
        forward,constrained,_repro_errors(forward,constrained)...))

    # Current saturation in a separately solved unbalanced feeder. The two
    # models deliberately allocate P/Q differently; do not assert equivalence.
    saturated_net = _repro_grid(_repro_phases(207.0+0im,5cis(.2));feeder=true)
    plant = _repro_plant()
    full = InverterControlRequest(p_available=6000.,p_rated=6000.,q_scale=6000.)
    positive = _repro_policy(constant_q=true)
    forward = solve_controlled_inverter(saturated_net,ControlledDevice(plant,
        SequenceController(positive;limiter=CommonScaleLimiter(pq_priority=:var),
            power_voltage_floor=1e-3)),full;solver_options=("tol"=>1e-9,))
    constrained = solve_constrained_inverter(saturated_net,ConstrainedInverter(plant,
        VoltVarWattConstraints(positive);dc_source=DCSourceBounds(p_min=0.)),full;
        solver_options=("tol"=>1e-9,))
    saturated = (;expected=:different_saturation,plant,request=full,forward,constrained,
        _repro_errors(forward,constrained)...)

    impossible_q = AverageVoltageVoltVarWatt(volt_var=PiecewiseLinearLaw(
        [180.,280.],[1.2,1.2];formulation=LocalC2Formulation(.01)))
    forward = solve_controlled_inverter(net,ControlledDevice(plant,
        SequenceController(impossible_q;limiter=CommonScaleLimiter(pq_priority=:var),
            power_voltage_floor=1e-3)),full;solver_options=("tol"=>1e-9,))
    constrained = solve_constrained_inverter(net,ConstrainedInverter(plant,
        VoltVarWattConstraints(impossible_q)),full;solver_options=("tol"=>1e-9,))
    impossible_service = (;expected=:infeasible_service,plant,forward,constrained,
        requested_q_var=7200.,minimum_required_S_VA=7200.)

    # At fixed POC U0 != 0, balanced grid current implies midpoint U0=POC U0.
    # The grounded capacitor then requires I0_shunt=Y_mid*U0, incompatible with
    # BOTH I0_grid=0 and the three-leg I0_converter=0. No objective resolves KCL.
    u0 = 3.0+0im
    zero_net = _repro_grid(unbalanced_voltage .+ u0)
    plant = _repro_plant(i_max=15.,r_filter=.05,x_filter=.15,
        r_filter_grid=.03,x_filter_grid=.10,c_filter_mid=30e-6,r_filter_damping=.5)
    balanced = solve_constrained_inverter(zero_net,ConstrainedInverter(plant,
        VoltVarWattConstraints(q_policy);current_allocation=:balanced_grid),request;
        solver_options=("tol"=>1e-9,))
    released = solve_constrained_inverter(zero_net,ConstrainedInverter(plant,
        VoltVarWattConstraints(q_policy);current_allocation=:balanced_converter),request;
        solver_options=("tol"=>1e-9,))
    y_mid = inv(complex(plant.r_filter_damping,-1/(2pi*plant.f*plant.c_filter_mid)))
    zero_sequence = (;expected=:incompatible_allocation,plant,u0,balanced,released,
        required_shunt_i0_A=y_mid*u0)
    return (;pairs,saturated,impossible_service,zero_sequence)
end

if abspath(PROGRAM_FILE) == @__FILE__
    rows = ibr_formulation_reproductions()
    for row in rows.feeder
        println("feeder $(row.voltage) V: ΔP=$(row.p_error_W) W, ΔQ=$(row.q_error_var) var, ΔI=$(row.current_error_A) A")
    end
    for row in rows.proportional
        println("proportional $(row.voltage) V: ΔP=$(row.p_error_W) W, ΔQ=$(row.q_error_var) var")
    end
    for row in rows.sequence_priority
        println("sequence $(row.priority): ΔI=$(row.current_error_A) A")
    end
    network = ibr_network_reproductions()
    for row in network.pairs
        println("$(row.id) [$(row.expected)]: ΔP=$(row.p_error_W) W, ΔQ=$(row.q_error_var) var")
    end
    println("saturated feeder [$(network.saturated.expected)]: ΔP=$(network.saturated.p_error_W) W, ΔQ=$(network.saturated.q_error_var) var")
    println("impossible Q: $(network.impossible_service.constrained.termination_status)")
    println("incompatible grid balance: $(network.zero_sequence.balanced.termination_status)")
end
