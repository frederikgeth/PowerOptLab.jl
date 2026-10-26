# Operating relationships on the existing physical plant. Unlike SequenceController,
# this formulation defines a feasible set and leaves dispatch to the caller.

"""
    VoltVarWattConstraints(intent, encoding; sensing_reference=:device)
    VoltVarWattConstraints(positive_policy; sensing_reference=:device)

Volt-var equality and volt-watt upper bound on total **POC** power. Reuses the
voltage sensing, flat-extended curves and smoothing of the positive-sequence
policies, but does not allocate apparent power or scale a current command.
`sensing_reference` has the same meaning as in [`SequenceController`](@ref).

With `volt_watt_basis=:rated`, the watt ceiling is `p_rated * fraction`.
With `:available`, it is `fraction * (p_available - (p_dc - p_poc))`.
An independent `p_dc <= p_available` bound accounts for all plant losses and
preserves required loss imports at zero availability. No Q curve means Q = 0;
no watt curve means fraction = 1. Curve encodings must currently be smooth.
These relationships do not certify compliance with a regional standard.
"""
struct VoltVarWattConstraints{P<:AbstractPositiveSequencePolicy}
    positive::P
    sensing_reference::Union{Symbol,String}
    function VoltVarWattConstraints(positive::P;
            sensing_reference::Union{Symbol,AbstractString}=:device) where
            {P<:AbstractPositiveSequencePolicy}
        new{P}(positive, _sensing_reference(sensing_reference))
    end
end

VoltVarWattConstraints(intent::VoltVarWattIntent, encoding::VoltVarWattEncoding;
                      kwargs...) =
    VoltVarWattConstraints(lower_positive_policy(intent, encoding); kwargs...)

"""
    ConstrainedInverter(plant::AdvancedInverter, policy::VoltVarWattConstraints)

Three-leg physical plant with voltage-dependent operating constraints. Converter
current is balanced positive sequence (`I2_converter = 0`); the plant supplies
`I0_converter = 0`. Under unbalance or with an LCL filter, grid current need not
be balanced. Q is constrained at the POC, including the declared POC shunt.
All hardware inequalities remain in the plant. An impossible Q request makes
the model infeasible rather than activating a limiter or a tracking slack.

Use [`stamp_device!`](@ref) with an `InverterControlRequest` to compose this
device with an external OPF objective, or [`solve_constrained_inverter`](@ref)
for a standalone dispatch solve. There is no implicit dispatch objective.
The initial formulation supports `topology=:THREE_LEG`, `grid_forming=false`
and `pwm_strategy=:NONE`; fixed RMS current reserves remain supported.
"""
struct ConstrainedInverter{P<:VoltVarWattConstraints} <: AbstractDevice
    device::AdvancedInverter
    policy::P
end

device_id(device::ConstrainedInverter) = device_id(device.device)
inverter_spec(device::ConstrainedInverter) = device.device
sensing_terminal(policy::VoltVarWattConstraints, inverter::AdvancedInverter) =
    policy.sensing_reference === :device ? inverter.neutral :
    policy.sensing_reference === :ground ? nothing : policy.sensing_reference

function validate_device(device::ConstrainedInverter, nets=();
                         periods::Integer=length(nets))
    inv = device.device
    validate_device(inv, nets; periods)
    inv.topology === :THREE_LEG && length(inv.phase_terminals) == 3 ||
        throw(ArgumentError("ConstrainedInverter requires a three-leg, three-phase plant"))
    inv.grid_forming && throw(ArgumentError(
        "ConstrainedInverter specifies balanced current, not grid-forming voltage"))
    inv.pwm_strategy === :NONE || throw(ArgumentError(
        "ConstrainedInverter requires pwm_strategy=:NONE; carrier-reserve iteration is not implemented"))
    reference = sensing_terminal(device.policy, inv)
    if reference !== nothing
        reference in inv.phase_terminals && throw(ArgumentError(
            "sensing reference '$reference' must not be a phase terminal"))
        for net in nets
            bus = get(get(net, "bus", Dict()), inv.bus, nothing)
            bus === nothing && continue
            reference in String.(get(bus, "terminal_names", String[])) ||
                throw(ArgumentError("sensing reference '$reference' is not a terminal of bus '$(inv.bus)'"))
        end
    end
    return device
end

struct _ConstrainedInverterHandles
    plant::_InvHandles
    phase_voltage
    p_fraction
    q_target
    p_cap
    request::InverterControlRequest
end

inverter_handles(::ConstrainedInverter, handles::_ConstrainedInverterHandles) =
    handles.plant

function stamp_device!(ctx, device::ConstrainedInverter;
                       request::InverterControlRequest)
    validate_device(device, (_opf_network(ctx),); periods=1)
    inv, policy = device.device, device.policy
    # At large network s_base, an unnormalised squared per-unit bound can be
    # smaller than the solver's constraint-bound relaxation. Normalise to the
    # device rating rather than relying on a current-command backoff for safety.
    h = stamp_device!(ctx, inv; retain_lossless_magnitude_epigraph=false,
        normalize_capability_constraints=true)
    m = _opf_model(ctx)
    phase_voltage = ntuple(k -> _dv(ctx, inv.bus, inv.phase_terminals[k],
        sensing_terminal(policy, inv)), 3)
    vre, vim = [v[1] for v in phase_voltage], [v[2] for v in phase_voltage]
    voltage_start = _control_voltage_start(policy.positive, h.vb)
    magnitudes = [_implicit_sqrt!(m, vre[k]^2 + vim[k]^2;
        start=voltage_start, scale=voltage_start) for k in 1:3]
    _, u1, _ = _sequence_pair(vre, vim)
    fraction, q, _, _ = _positive_smooth(
        ctx, policy.positive, magnitudes, u1, request, h.sb, h.vb)

    # Only the negative sequence needs additional rows: the three-leg plant
    # already stamps the two real I0 equalities. Do not repeat them.
    _, _, i2 = _sequence_pair(h.cri, h.cii)
    JuMP.@constraint(m, i2[1] == 0)
    JuMP.@constraint(m, i2[2] == 0)
    JuMP.@constraint(m, h.q_poc == q)
    JuMP.@constraint(m, h.p_dc <= request.p_available/h.sb)
    cap = if _volt_watt_basis(policy.positive) === :rated
        request.p_rated/h.sb * fraction
    else
        # With the separate DC bound, P <= f*A also implements P <= min(A,f*A).
        # For A < 0 the DC bound retains the full loss import, without a max()
        # selector or a smoothed approximation of the availability transition.
        fraction * (request.p_available/h.sb - (h.p_dc - h.p_poc))
    end
    JuMP.@constraint(m, h.p_poc <= cap)
    return _ConstrainedInverterHandles(h, phase_voltage, fraction, q, cap, request)
end

link_device!(model, device::ConstrainedInverter, handles, sb, grid::TimeGrid) = nothing

function extract_device(device::ConstrainedInverter,
                        handles::_ConstrainedInverterHandles, status::SolveStatus)
    h, request = handles.plant, handles.request
    inv = device.device
    scalar(x, scale=1.0) = status.publishable ? Float64(JuMP.value(x))*scale : NaN
    phase_voltage = ntuple(k -> ComplexF64(
        scalar(handles.phase_voltage[k][1], h.vb),
        scalar(handles.phase_voltage[k][2], h.vb)), 3)
    plant = extract_device(inv, h, status)
    fraction_exact, q_exact = if status.publishable
        u1 = _phasor_sequences(phase_voltage)[2]
        f, q, _, _ = _positive_exact(
            device.policy.positive, u1, abs.(phase_voltage), request)
        f, q
    else
        NaN, NaN
    end
    ac_available = request.p_available - (plant.p_dc - plant.p_poc)
    exact_cap = _volt_watt_basis(device.policy.positive) === :rated ?
        request.p_rated*fraction_exact : ac_available*fraction_exact
    current_margin(limit, current, reserve) = limit === nothing ? nothing :
        ntuple(k -> limit - hypot(current[k], reserve[k]), 3)
    operation = (
        phase_voltage=phase_voltage,
        volt_watt_fraction=scalar(handles.p_fraction),
        q_target_var=scalar(handles.q_target, h.sb),
        q_exact_target_var=q_exact,
        q_constraint_residual_var=plant.q_poc-scalar(handles.q_target, h.sb),
        q_tracking_residual_var=plant.q_poc-q_exact,
        volt_watt_cap_W=scalar(handles.p_cap, h.sb),
        volt_watt_margin_W=scalar(handles.p_cap, h.sb)-plant.p_poc,
        exact_volt_watt_margin_W=exact_cap-plant.p_poc,
        ac_available_W=ac_available,
        dc_availability_margin_W=request.p_available-plant.p_dc,
        converter_apparent_power_margin_VA=inv.s_max-hypot(plant.p_conv, plant.q_conv),
        converter_current_margin_A=current_margin(
            inv.i_max, plant.i_mag, inv.pwm_ac_converter_reserve),
        grid_current_margin_A=current_margin(
            inv.i_grid_max, plant.i_grid_mag, inv.pwm_ac_grid_reserve),
        dc_ripple_margin_V=inv.dv2_max === nothing ? nothing : inv.dv2_max-plant.dv2,
        switching_margin_V=plant.switching_margin,
    )
    return (; plant, operation,
        converter_terminal=_extract_converter_terminal(h, status),
        grid_phase_current=_extract_grid_current(h, status))
end

"""
Result of [`solve_constrained_inverter`](@ref). `plant` contains physical SI
quantities; `operation` contains curve targets, POC tracking residuals and signed
capability margins (positive is feasible). `q_constraint_residual_var` checks
the encoded equality, while `q_tracking_residual_var` compares with the exact
PWL policy at the solved sensing voltage. Neither is a current-command replay.
Unpublishable solves retain their status and return NaN for numerical outputs.
"""
struct ConstrainedInverterResult <: AbstractSolveResult
    termination_status::String
    plant::NamedTuple
    operation::NamedTuple
    converter_terminal::ConverterTerminalResult
    grid_phase_current::NTuple{3,ComplexF64}
    bus::Dict{String,Any}
    solve::SolveStatus
end

solve_status(result::ConstrainedInverterResult) = result.solve
solve_diagnostics(result::ConstrainedInverterResult) = result.operation

"""
    solve_constrained_inverter(net, device, request;
        objective=:max_export, p_set=nothing, per_unit=true,
        s_base=1e6, optimizer=Ipopt.Optimizer, verbose=false, solver_options=())

Solve a three-leg inverter subject to [`VoltVarWattConstraints`](@ref).
`:max_export` maximizes POC active injection. `:feasibility` requires a finite
`p_set` (W, injection positive) and checks that fixed dispatch. Both retain the
same Q equality and hard plant limits. Active power can import to supply losses
when availability is zero. No controller-side current scaling is introduced.

This is a nonconvex OPF/dispatch formulation, not an autonomous inverter power
flow or a proof of global optimality. Use `stamp_device!(ctx, device; request)`
inside an external model hook to select another objective; stamping does not
replace the objective. The standalone solver requires `per_unit=true`.
"""
function solve_constrained_inverter(net::Dict{String,Any},
        device::ConstrainedInverter, request::InverterControlRequest;
        objective::Symbol=:max_export, p_set::Union{Nothing,Real}=nothing,
        per_unit::Bool=true, s_base::Real=1e6, optimizer=Ipopt.Optimizer,
        verbose::Bool=false, solver_options=())
    objective in (:max_export, :feasibility) || throw(ArgumentError(
        "objective must be :max_export or :feasibility"))
    if objective === :feasibility
        p_set !== nothing && isfinite(p_set) || throw(ArgumentError(
            ":feasibility requires a finite p_set in watts"))
    else
        p_set === nothing || throw(ArgumentError("p_set requires objective=:feasibility"))
    end
    per_unit || throw(ArgumentError("solve_constrained_inverter requires per_unit=true"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be finite and positive"))
    validate_device(device, (net,); periods=1)
    handles = Ref{_ConstrainedInverterHandles}()
    hook! = ctx -> begin
        handles[] = stamp_device!(ctx, device; request)
        h = handles[].plant
        if objective === :max_export
            JuMP.@objective(_opf_model(ctx), Max, h.p_poc)
        else
            JuMP.@constraint(_opf_model(ctx), h.p_poc == p_set/h.sb)
            JuMP.@objective(_opf_model(ctx), Min, 0.0)
        end
    end
    ctx = build_opf_model(net; per_unit=true, s_base=Float64(s_base),
        add_objective=false, model_hook! = hook!, optimizer, verbose)
    _set_solver_options!(_opf_model(ctx), solver_options)
    BMOPFTools.optimize_opf!(ctx)
    outcome = _solve_outcome(_opf_model(ctx))
    status = SolveStatus(outcome)
    extracted = extract_device(device, handles[], status)
    network = _extract_result(ctx, outcome)
    return ConstrainedInverterResult(string(outcome.termination_status),
        extracted.plant, extracted.operation, extracted.converter_terminal,
        extracted.grid_phase_current, network["bus"], status)
end
