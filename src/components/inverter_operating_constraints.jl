# Operating relationships on the existing physical plant. Unlike SequenceController,
# this formulation defines a feasible set and leaves dispatch to the caller.

abstract type _AbstractInverterOperatingPolicy end
abstract type _AbstractReactivePowerRelation end

_operating_bound(x::Nothing, label) = nothing
function _operating_bound(x::Real, label)
    value = Float64(x)
    isfinite(value) || throw(ArgumentError("$label must be finite when supplied"))
    return value
end
function _operating_bounds(lo, hi, label)
    lower, upper = _operating_bound(lo, "$label minimum"), _operating_bound(hi, "$label maximum")
    lower === nothing || upper === nothing || lower <= upper ||
        throw(ArgumentError("$label minimum must not exceed its maximum"))
    return lower, upper
end

"""Leave total POC reactive power dispatchable within the physical and declared bounds."""
struct FreeReactivePower <: _AbstractReactivePowerRelation end

"""`FixedReactivePower(q)` imposes total POC Q = q (var, injection positive)."""
struct FixedReactivePower <: _AbstractReactivePowerRelation
    q::Float64
    FixedReactivePower(q::Real) = new(_operating_bound(q, "fixed reactive power"))
end

"""
    FixedPowerFactor(pf, direction)
    FixedPowerFactor(pf)

Fixed power-factor relation on total POC powers. Prefer the explicit form:
`direction` is `:absorbing` (Q < 0 while exporting P, under-excited) or
`:injecting` (Q > 0 while exporting P, over-excited), with `0 < pf <= 1`.
Standards and datasheets attach different signs and leading/lagging labels to
the same operating point, so the direction is named rather than signed.

The single-argument form is the stored signed value: for positive P, positive
pf absorbs Q and negative pf supplies Q. Requires `0 < abs(pf) <= 1`.
The relation uses bounded sine/cosine coefficients, including near zero PF,
instead of an overflowing tangent. It also applies to negative P with the same
algebraic ratio; it does not impose a generation or charging direction. At
P = Q = 0 the equation is satisfied but a measured PF is undefined.
"""
struct FixedPowerFactor <: _AbstractReactivePowerRelation
    pf::Float64
    function FixedPowerFactor(pf::Real)
        value = Float64(pf)
        isfinite(value) && 0 < abs(value) <= 1 || throw(ArgumentError(
            "fixed power factor must be finite with 0 < abs(pf) <= 1"))
        new(value)
    end
end

function FixedPowerFactor(pf::Real, direction::Symbol)
    direction in (:absorbing, :injecting) || throw(ArgumentError(
        "power-factor direction must be :absorbing or :injecting"))
    value = Float64(pf)
    isfinite(value) && 0 < value <= 1 || throw(ArgumentError(
        "with an explicit direction, pf must satisfy 0 < pf <= 1"))
    return FixedPowerFactor(direction === :absorbing ? value : -value)
end

"""
    DCSourceBounds(; p_min=nothing, p_max=nothing)

Additional signed bounds on physical DC-link power (W, positive from the source
into the converter). Request availability always remains an upper bound.
`p_min=0` models a non-charging PV source while permitting AC loss imports.
Signed bounds permit bidirectional sources. Omitted bounds retain the legacy
availability-only source model; they do not assert that the source is PV.
"""
struct DCSourceBounds
    p_min::Union{Nothing,Float64}
    p_max::Union{Nothing,Float64}
    function DCSourceBounds(; p_min::Union{Nothing,Real}=nothing,
                             p_max::Union{Nothing,Real}=nothing)
        lo, hi = _operating_bounds(p_min, p_max, "DC source power")
        new(lo, hi)
    end
end

"""
    InverterOperatingConstraints(positive_policy=AverageVoltageVoltVarWatt();
        reactive=FreeReactivePower(), p_min=nothing, p_max=nothing,
        q_min=nothing, q_max=nothing, sensing_reference=:device)

POC capability/dispatch relationships with an explicit reactive operating mode:
[`FreeReactivePower`](@ref), [`FixedReactivePower`](@ref), or
[`FixedPowerFactor`](@ref). Optional finite P/Q bounds are in W/var and compose
with every mode. The positive policy may supply a volt-watt ceiling. **No watt
curve means no watt ceiling at all**, even with `volt_watt_basis=:rated`; add
`p_max=p_rated` for a rated cap. This differs deliberately from
[`VoltVarWattConstraints`](@ref), which mirrors the forward controller. A policy
containing a volt-var curve is rejected: use `VoltVarWattConstraints` for that
Q equality.

An objective selects the remaining freedom. Neither a PF relation nor balanced
current is an implication of the hardware. DC source bounds and current
allocation are specified on [`ConstrainedInverter`](@ref).
"""
struct InverterOperatingConstraints{P<:AbstractPositiveSequencePolicy,
                                    R<:_AbstractReactivePowerRelation} <: _AbstractInverterOperatingPolicy
    positive::P
    reactive::R
    p_min::Union{Nothing,Float64}
    p_max::Union{Nothing,Float64}
    q_min::Union{Nothing,Float64}
    q_max::Union{Nothing,Float64}
    sensing_reference::Union{Symbol,String}
    function InverterOperatingConstraints(positive::P=AverageVoltageVoltVarWatt();
            reactive::_AbstractReactivePowerRelation=FreeReactivePower(), p_min::Union{Nothing,Real}=nothing,
            p_max::Union{Nothing,Real}=nothing, q_min::Union{Nothing,Real}=nothing,
            q_max::Union{Nothing,Real}=nothing,
            sensing_reference::Union{Symbol,AbstractString}=:device) where
            {P<:AbstractPositiveSequencePolicy}
        positive.volt_var === nothing || throw(ArgumentError(
            "use VoltVarWattConstraints for a volt-var equality; do not combine two Q laws"))
        pl, pu = _operating_bounds(p_min, p_max, "POC active power")
        ql, qu = _operating_bounds(q_min, q_max, "POC reactive power")
        new{P,typeof(reactive)}(positive, reactive, pl, pu, ql, qu, _sensing_reference(sensing_reference))
    end
end

"""
    VoltVarWattConstraints(intent, encoding; sensing_reference=:device)
    VoltVarWattConstraints(positive_policy; sensing_reference=:device)

Volt-var equality and volt-watt upper bound on total **POC** power. Reuses the
voltage sensing, flat-extended curves and smoothing of the positive-sequence
policies, but does not allocate apparent power or scale a current command.
`sensing_reference` has the same meaning as in [`SequenceController`](@ref).

With `volt_watt_basis=:rated`, the watt ceiling is `p_rated * fraction`.
With `:available`, it is `fraction * (available_dc - (p_dc - p_poc))`, where
`available_dc` is the tighter of request availability and a DC source upper bound.
An independent `p_dc <= p_available` bound accounts for all plant losses and
preserves required loss imports at zero availability. No Q curve means Q = 0.

No watt curve means fraction = 1, and the ceiling is **still imposed**, matching
the forward controller's `min(p_rated*fraction, available)` command: with
`:rated` this caps P at `p_rated`; with `:available` it coincides with the DC
availability bound. [`InverterOperatingConstraints`](@ref) instead imposes no
ceiling without a watt curve.

The `:available` ceiling subtracts the **solved** plant losses `p_dc - p_poc`,
so it depends on the operating point: more reactive current raises filter and
conduction losses and lowers the ceiling. Firmware that scales a separately
measured or estimated available power can differ by those losses. The product
of `fraction` and the losses is an additional nonconvex term.

Curve encodings must currently be smooth. These relationships do not certify
compliance with a regional standard.
"""
struct VoltVarWattConstraints{P<:AbstractPositiveSequencePolicy} <: _AbstractInverterOperatingPolicy
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
    ConstrainedInverter(plant, policy;
        current_allocation=:balanced_converter, dc_source=DCSourceBounds())

Three-phase physical plant with explicit operating constraints. The default
`:balanced_converter` imposes zero/negative sequence at the converter;
`:balanced_grid` imposes balance on the grid arm; `:free` leaves sequence
allocation to the optimisation within topology and capability constraints.
The plant alone determines the three-leg zero-sum restriction. LCL grid and
converter currents can differ. POC powers include the declared POC shunt.
All hardware inequalities remain in the plant. An impossible Q request makes
the model infeasible rather than activating a limiter or a tracking slack.

Use [`stamp_device!`](@ref) with an `InverterControlRequest` to compose this
device with an external OPF objective, or [`solve_constrained_inverter`](@ref)
for a standalone dispatch solve. There is no implicit dispatch objective.
Supports THREE_LEG, FOUR_LEG and SPLIT_DC with `grid_forming=false` and
`pwm_strategy=:NONE`; fixed RMS current reserves remain supported. The default
DC source retains the request availability upper bound. Supply
`dc_source=DCSourceBounds(p_min=0)` for non-charging PV.
"""
struct ConstrainedInverter{P<:_AbstractInverterOperatingPolicy} <: AbstractDevice
    device::AdvancedInverter
    policy::P
    current_allocation::Symbol
    dc_source::DCSourceBounds
    function ConstrainedInverter(device::AdvancedInverter, policy::P;
            current_allocation::Symbol=:balanced_converter,
            dc_source::DCSourceBounds=DCSourceBounds()) where {P<:_AbstractInverterOperatingPolicy}
        current_allocation in (:balanced_converter, :balanced_grid, :free) ||
            throw(ArgumentError("current_allocation must be :balanced_converter, :balanced_grid or :free"))
        new{P}(device, policy, current_allocation, dc_source)
    end
end

device_id(device::ConstrainedInverter) = device_id(device.device)
inverter_spec(device::ConstrainedInverter) = device.device
sensing_terminal(policy::_AbstractInverterOperatingPolicy, inverter::AdvancedInverter) =
    policy.sensing_reference === :device ? inverter.neutral :
    policy.sensing_reference === :ground ? nothing : policy.sensing_reference

function validate_device(device::ConstrainedInverter, nets=();
                         periods::Integer=length(nets))
    inv = device.device
    validate_device(inv, nets; periods)
    inv.topology in (:THREE_LEG, :FOUR_LEG, :SPLIT_DC) && length(inv.phase_terminals) == 3 ||
        throw(ArgumentError("ConstrainedInverter requires a supported three-phase plant"))
    inv.grid_forming && throw(ArgumentError(
        "ConstrainedInverter does not yet support grid-forming voltage laws"))
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

_poc_operating_bounds(::VoltVarWattConstraints) =
    (p_min=nothing, p_max=nothing, q_min=nothing, q_max=nothing)
_poc_operating_bounds(p::InverterOperatingConstraints) =
    (p_min=p.p_min, p_max=p.p_max, q_min=p.q_min, q_max=p.q_max)
_has_watt_cap(::VoltVarWattConstraints) = true
_has_watt_cap(p::InverterOperatingConstraints) = p.positive.volt_watt !== nothing
_reactive_mode(p::VoltVarWattConstraints) = p.positive.volt_var === nothing ? :zero_q : :volt_var
_reactive_mode(p::InverterOperatingConstraints) =
    p.reactive isa FreeReactivePower ? :free_q :
    p.reactive isa FixedReactivePower ? :fixed_q : :fixed_pf
_pf_coefficients(law::FixedPowerFactor) =
    (abs(law.pf), sign(law.pf)*sqrt((1-abs(law.pf))*(1+abs(law.pf))))
_available_dc(device, request) = device.dc_source.p_max === nothing ? request.p_available :
    min(request.p_available, device.dc_source.p_max)

function _stamp_operating_q!(m, policy::VoltVarWattConstraints, h, curve_q, power_scale)
    JuMP.@constraint(m, (h.q_poc - curve_q)/power_scale == 0)
    return curve_q
end
function _stamp_operating_q!(m, policy::InverterOperatingConstraints, h, curve_q, power_scale)
    law = policy.reactive
    if law isa FixedReactivePower
        target = law.q/h.sb
        JuMP.@constraint(m, (h.q_poc - target)/power_scale == 0)
        return target
    elseif law isa FixedPowerFactor
        cosine, signed_sine = _pf_coefficients(law)
        JuMP.@constraint(m, (cosine*h.q_poc + signed_sine*h.p_poc)/power_scale == 0)
    end
    return nothing
end

function _stamp_current_allocation!(m, device, h)
    allocation = device.current_allocation
    allocation === :free && return
    converter = allocation === :balanced_converter
    ir, ii = converter ? (h.cri, h.cii) : (h.gri, h.gii)
    i0, _, i2 = _sequence_pair(ir, ii)
    # With no midpoint capacitor, grid and converter I0 coincide even if the
    # LCL series arms use different variables. Avoid redundant topology rows.
    plant_zero_sum = device.device.topology === :THREE_LEG &&
        (converter || device.device.c_filter_mid == 0)
    if !plant_zero_sum
        JuMP.@constraint(m, i0[1] == 0)
        JuMP.@constraint(m, i0[2] == 0)
    end
    JuMP.@constraint(m, i2[1] == 0)
    JuMP.@constraint(m, i2[2] == 0)
end

# Device validation uses SI parameters, including the voltage bounds used to
# diagnose negligible filters. A staged context owns a per-unit working net.
# Restore only this bus's diagnostic bounds in a shallow copy; never change
# the network data or coordinates used by the electrical model.
function _constrained_validation_network(ctx, inv)
    net, bases = _opf_network(ctx), _opf_bases(ctx)
    bases === nothing && return net
    bus = get(get(net, "bus", Dict()), inv.bus, nothing)
    bus === nothing && return net
    physical_bus = copy(bus)
    vb = get(bases.v_base, inv.bus, 1.0)
    for field in ("v_min", "v_max")
        bound = get(bus, field, nothing)
        bound === nothing || (physical_bus[field] = bound .* vb)
    end
    validation_net = copy(net)
    validation_net["bus"] = copy(net["bus"])
    validation_net["bus"][inv.bus] = physical_bus
    return validation_net
end

function stamp_device!(ctx, device::ConstrainedInverter;
                       request::InverterControlRequest)
    inv, policy = device.device, device.policy
    validate_device(device, (_constrained_validation_network(ctx, inv),); periods=1)
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

    _stamp_current_allocation!(m, device, h)
    power_scale = inv.s_max/h.sb
    q_target = _stamp_operating_q!(m, policy, h, q, power_scale)
    bounds = _poc_operating_bounds(policy)
    for (quantity, lo, hi) in ((h.p_poc, bounds.p_min, bounds.p_max),
                                (h.q_poc, bounds.q_min, bounds.q_max))
        lo === nothing || JuMP.@constraint(m, quantity/power_scale >= lo/inv.s_max)
        hi === nothing || JuMP.@constraint(m, quantity/power_scale <= hi/inv.s_max)
    end
    source = device.dc_source
    upper = _available_dc(device, request)
    JuMP.@constraint(m, h.p_dc/power_scale <= upper/inv.s_max)
    source.p_min === nothing || JuMP.@constraint(m, h.p_dc/power_scale >= source.p_min/inv.s_max)
    cap = if !_has_watt_cap(policy)
        nothing
    elseif _volt_watt_basis(policy.positive) === :rated
        request.p_rated/h.sb * fraction
    else
        # With the separate DC bound, P <= f*A also implements P <= min(A,f*A).
        # For A < 0 the DC bound retains the full loss import, without a max()
        # selector or a smoothed approximation of the availability transition.
        fraction * (upper/h.sb - (h.p_dc - h.p_poc))
    end
    cap === nothing || JuMP.@constraint(m, (h.p_poc-cap)/power_scale <= 0)
    return _ConstrainedInverterHandles(h, phase_voltage, fraction, q_target, cap, request)
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
    fraction_exact, curve_q_exact = if status.publishable
        u1 = _phasor_sequences(phase_voltage)[2]
        f, q, _, _ = _positive_exact(
            device.policy.positive, u1, abs.(phase_voltage), request)
        f, q
    else
        NaN, NaN
    end
    policy = device.policy
    q_exact = policy isa VoltVarWattConstraints ? curve_q_exact :
        policy.reactive isa FixedReactivePower ? (status.publishable ? policy.reactive.q : NaN) : nothing
    q_target = handles.q_target === nothing ? nothing : scalar(handles.q_target, h.sb)
    ac_available = _available_dc(device, request) - (plant.p_dc - plant.p_poc)
    exact_cap = !_has_watt_cap(policy) ? nothing :
        _volt_watt_basis(policy.positive) === :rated ?
            request.p_rated*fraction_exact : ac_available*fraction_exact
    cap = handles.p_cap === nothing ? nothing : scalar(handles.p_cap, h.sb)
    fixed_pf_residual = if policy isa InverterOperatingConstraints && policy.reactive isa FixedPowerFactor
        c, s = _pf_coefficients(policy.reactive)
        c*plant.q_poc+s*plant.p_poc
    else
        nothing
    end
    bounds = _poc_operating_bounds(policy)
    current_margin(limit, current, reserve) = limit === nothing ? nothing :
        ntuple(k -> limit - hypot(current[k], reserve[k]), 3)
    operation = (
        reactive_mode=_reactive_mode(policy),
        current_allocation=device.current_allocation,
        dc_source_bounds_W=(p_min=device.dc_source.p_min, p_max=device.dc_source.p_max),
        phase_voltage=phase_voltage,
        volt_watt_fraction=scalar(handles.p_fraction),
        q_target_var=q_target,
        q_exact_target_var=q_exact,
        q_constraint_residual_var=q_target === nothing ? nothing : plant.q_poc-q_target,
        q_tracking_residual_var=q_exact === nothing ? nothing : plant.q_poc-q_exact,
        fixed_pf_residual_VA=fixed_pf_residual,
        volt_watt_cap_W=cap,
        volt_watt_margin_W=cap === nothing ? nothing : cap-plant.p_poc,
        exact_volt_watt_margin_W=exact_cap === nothing ? nothing : exact_cap-plant.p_poc,
        active_lower_margin_W=bounds.p_min === nothing ? nothing : plant.p_poc-bounds.p_min,
        active_upper_margin_W=bounds.p_max === nothing ? nothing : bounds.p_max-plant.p_poc,
        reactive_lower_margin_var=bounds.q_min === nothing ? nothing : plant.q_poc-bounds.q_min,
        reactive_upper_margin_var=bounds.q_max === nothing ? nothing : bounds.q_max-plant.q_poc,
        ac_available_W=ac_available,
        effective_dc_availability_W=_available_dc(device, request),
        dc_availability_margin_W=request.p_available-plant.p_dc,
        dc_lower_margin_W=device.dc_source.p_min === nothing ? nothing : plant.p_dc-device.dc_source.p_min,
        dc_upper_margin_W=device.dc_source.p_max === nothing ? nothing : device.dc_source.p_max-plant.p_dc,
        converter_apparent_power_margin_VA=inv.s_max-hypot(plant.p_conv, plant.q_conv),
        converter_current_margin_A=current_margin(
            inv.i_max, plant.i_mag, inv.pwm_ac_converter_reserve),
        grid_current_margin_A=current_margin(
            inv.i_grid_max, plant.i_grid_mag, inv.pwm_ac_grid_reserve),
        dc_ripple_margin_V=inv.dv2_max === nothing ? nothing : inv.dv2_max-plant.dv2,
        neutral_current_margin_A=inv.In_max === nothing ? nothing :
            inv.In_max-hypot(plant.i_neutral, inv.pwm_ac_neutral_reserve),
        capacitor_current_margin_A=inv.i_cap_max === nothing ? nothing : inv.i_cap_max-plant.i_cap_thermal,
        upper_capacitor_margin_A=inv.i_cap_upper_max === nothing ? nothing : inv.i_cap_upper_max-plant.i_cap_thermal_upper,
        lower_capacitor_margin_A=inv.i_cap_lower_max === nothing ? nothing : inv.i_cap_lower_max-plant.i_cap_thermal_lower,
        midpoint_ripple_margin_V=inv.dv_mid_max === nothing ? nothing : inv.dv_mid_max-plant.dv_mid,
        midpoint_mean_margin_V=inv.v_mid_mean_max === nothing ? nothing : inv.v_mid_mean_max-abs(plant.v_mid_mean),
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
Unpublishable solves retain their status and configuration, with NaN for solved
numerical quantities and `nothing` for absent optional targets.
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

Solve a three-phase [`ConstrainedInverter`](@ref) with its declared operating policy.
`:max_export` maximizes POC active injection. `:feasibility` requires a finite
`p_set` (W, injection positive) and checks that fixed dispatch. Both retain the
same reactive relationship and hard plant limits. Active power can import to supply losses
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
            JuMP.@objective(_opf_model(ctx), Max, h.p_poc/(device.device.s_max/h.sb))
        else
            JuMP.@constraint(_opf_model(ctx), h.p_poc/(device.device.s_max/h.sb) == p_set/device.device.s_max)
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
