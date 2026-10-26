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
    status = SolveStatus(PowerOptLab._solve_outcome(m))
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
end
