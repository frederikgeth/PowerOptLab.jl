# Bilevel distribution-network proof of concept.
#
# The lower-level model is deliberately kept as a live BMOPFTools staged model:
# the utility tap is a JuMP Parameter linked to BMOPFTools' native tap variable,
# and DiffOpt differentiates the solved KKT system through that link. This is a
# local sensitivity tool, not a global bilevel optimality certificate.

"""
    BilevelPVResult

Result of [`solve_bilevel_pv_tap`](@ref). `tap` is the user-facing tap
multiplier, `exported_power_W` is the lower-level aggregate PV export, and
`voltages_V` contains the monitored phase-to-neutral voltage magnitudes.
`converged` and `termination_reason` describe the upper-level search; a
lower-level `LOCALLY_SOLVED` status alone is not an upper-level convergence
certificate.

The lower-level consumers are represented by one BMOPFTools native IBR each.
Their Volt-var and Volt-watt curves are therefore part of the differentiated
network model. `lower_level=:aggregate` uses a shared aggregate-export
objective; `lower_level=:local_controller` instead ties each PV's active
output to its own Volt-watt response, with Volt-var priority enforced through
a smooth apparent-power cap.
"""
struct BilevelPVResult <: AbstractSolveResult
    termination_status::String
    lower_level::Symbol
    tap::Float64
    exported_power_W::Float64
    voltages_V::Dict{String,Float64}
    utility_objective::Float64
    history::Vector{NamedTuple}
    lower_result::Dict{String,Any}
    differentiability_report::Any
    converged::Bool
    termination_reason::Symbol
end

"""
    BilevelPVResponse

Fixed-tap lower-level response used for sensitivity validation and experiment
scripts. The response contains the solved PV export, monitored voltages, and
the local DiffOpt derivative of each monitored voltage with respect to the
utility tap multiplier.
"""
struct BilevelPVResponse <: AbstractSolveResult
    termination_status::String
    lower_level::Symbol
    tap::Float64
    exported_power_W::Float64
    voltages_V::Dict{String,Float64}
    voltage_sensitivity_V_per_tap::Dict{String,Float64}
    lower_result::Dict{String,Any}
    differentiability_report::Any
end

"""Result of the centralized single-level comparison solve."""
struct SingleLevelPVResult <: AbstractSolveResult
    termination_status::String
    tap::Float64
    exported_power_W::Float64
    voltages_V::Dict{String,Float64}
    objective::Float64
    result::Dict{String,Any}
end

solve_status(result::BilevelPVResult) = _result_solve_status(
    result.termination_status,
    result.converged && result.termination_status in ("LOCALLY_SOLVED", "OPTIMAL"))
solve_status(result::BilevelPVResponse) = _result_solve_status(
    result.termination_status, result.termination_status in ("LOCALLY_SOLVED", "OPTIMAL"))
solve_status(result::SingleLevelPVResult) = _result_solve_status(
    result.termination_status, result.termination_status in ("LOCALLY_SOLVED", "OPTIMAL"))

function _bilevel_transformer(net::Dict{String,Any}, transformer_id::String)
    transformers = get(net, "transformer", Dict{String,Any}())
    for (family, entries) in transformers
        haskey(entries, transformer_id) && return String(family), entries[transformer_id]
    end
    throw(ArgumentError("transformer '$transformer_id' was not found"))
end

"""
    solve_bilevel_pv_response(net; transformer_id, pv_ids, monitored_buses,
        tap=1.0, tap_bounds=(0.95, 1.05), lower_level=:aggregate, ...)

Solve one fixed-tap lower-level response and return its local voltage
sensitivity with respect to the tap. This is the public diagnostic counterpart
to the differentiated response used by [`solve_bilevel_pv_tap`](@ref): it is
useful for comparing DiffOpt derivatives against finite differences and for
checking sensitivity quality near Volt-watt transitions.
"""
function solve_bilevel_pv_response(net::Dict{String,Any};
                                   transformer_id::AbstractString,
                                   pv_ids::AbstractVector{<:AbstractString},
                                   monitored_buses::AbstractVector{<:AbstractString},
                                   tap::Real=1.0,
                                   tap_bounds=(0.95, 1.05),
                                   lower_level::Symbol=:aggregate,
                                   voltage_reference::Real=230.0,
                                   voltage_band::Real=15.0,
                                   volt_var_watt_eps::Real=2e-3,
                                   optimizer=Ipopt.Optimizer,
                                   verbose::Bool=false)
    tid = String(transformer_id)
    pids = String.(pv_ids)
    buses = String.(monitored_buses)
    tap_value = Float64(tap)
    _bilevel_validation(net, tid, pids, buses)
    lower_level in (:aggregate, :local_controller) || throw(ArgumentError(
        "lower_level must be :aggregate or :local_controller"))
    _bilevel_voltage_data(net, buses, Float64(voltage_reference),
                          Float64(voltage_band))
    volt_var_watt_eps > 0 || throw(ArgumentError("volt_var_watt_eps must be > 0"))
    lo, hi = Float64.(tap_bounds)
    state = _bilevel_lower_state(net, tid, pids, buses, tap_value, (lo, hi);
                                 lower_level, relu_eps=Float64(volt_var_watt_eps),
                                 optimizer, verbose)
    snap = _bilevel_lower_solve!(state, tap_value; differentiate=true)
    lower_result = BMOPFTools.extract_result(state.context)
    BilevelPVResponse(snap.status, lower_level, tap_value, snap.exported,
        snap.voltages, snap.derivatives, lower_result,
        BMOPFTools.opf_differentiability_report(state.context))
end

function _bilevel_copy_with_tap_bounds(net, transformer_id, tap_bounds, tap_initial)
    lo, hi = Float64.(tap_bounds)
    isfinite(lo) && isfinite(hi) && lo < hi || throw(ArgumentError(
        "tap_bounds must be finite and strictly ordered"))
    lo <= tap_initial <= hi || throw(ArgumentError(
        "tap_initial must lie inside tap_bounds"))
    family, transformer = _bilevel_transformer(net, transformer_id)
    family == "single_phase" || throw(ArgumentError(
        "the proof of concept currently supports single_phase transformers; " *
        "extend _bilevel_tap_scale for other transformer families"))
    n = deepcopy(net)
    xfmr = n["transformer"][family][transformer_id]
    xfmr["tap"] = Float64(tap_initial)
    xfmr["tap_min"] = lo
    xfmr["tap_max"] = hi
    return n
end

function _bilevel_phase_terminals(net, bus::String)
    data = get(net["bus"], bus, nothing)
    data === nothing && throw(ArgumentError("monitored bus '$bus' was not found"))
    terminals = String.(get(data, "terminal_names", String[]))
    isempty(terminals) && throw(ArgumentError(
        "monitored bus '$bus' must define terminal_names"))
    neutral = findfirst(t -> lowercase(t) in ("n", "neutral"), terminals)
    neutral === nothing && throw(ArgumentError(
        "monitored bus '$bus' must define a neutral terminal for phase-to-neutral voltage"))
    phases = [t for (i, t) in enumerate(terminals) if i != neutral]
    isempty(phases) && throw(ArgumentError(
        "monitored bus '$bus' must have at least one phase terminal"))
    phases, terminals[neutral]
end

function _bilevel_voltage_handles(ctx, net, buses::Vector{String})
    model = _opf_model(ctx)
    vr, vi = _opf_voltage_maps(ctx)
    bases = _opf_bases(ctx)
    voltage_start = bases === nothing ? 230.0 : 1.0
    handles = Dict{String,Any}()
    for bus in buses
        phases, neutral = _bilevel_phase_terminals(net, bus)
        for phase in phases
            key = (bus, phase)
            haskey(vr, key) || throw(ArgumentError(
                "monitored bus '$bus' must have phase terminal '$phase'"))
            pair_key = length(phases) == 1 ? bus : "$(bus):$(phase)"
            vm = JuMP.@variable(model, lower_bound=0.0,
                                base_name="bilevel_vm_$(replace(pair_key, ':' => '_'))")
            JuMP.set_start_value(vm, voltage_start)
            JuMP.@constraint(model, vm^2 ==
                (vr[(bus, phase)] - vr[(bus, neutral)])^2 +
                (vi[(bus, phase)] - vi[(bus, neutral)])^2)
            handles[pair_key] = vm
        end
    end
    handles
end

function _bilevel_pv_handles(ctx, pv_ids::Vector{String})
    Dict{String,Any}(id => BMOPFTools.opf_object(ctx,
        BMOPFTools.opf_ibr_power_key(id, 1)) for id in pv_ids)
end

# These small curve helpers mirror BMOPFTools' built-in smooth ReLU
# representation. They are kept local so the public POC does not depend on
# extension-module internals while still constructing the same expressions.
struct _BilevelBuiltinSoftplus
    eps::Float64
end

(_op::_BilevelBuiltinSoftplus)(x::Real) = _op.eps * log1p(exp(x / _op.eps))
(_op::_BilevelBuiltinSoftplus)(x) = _op.eps * log1p(exp(x / _op.eps))

function _bilevel_breakpoints_to_triples(xs, ys)
    length(xs) == length(ys) || throw(ArgumentError(
        "curve breakpoints and values must have the same length"))
    length(xs) >= 2 || throw(ArgumentError("curve requires at least two breakpoints"))
    triples = Tuple{Float64,Float64}[]
    baseline = Float64(ys[1])
    for i in 1:(length(xs) - 1)
        dx = Float64(xs[i + 1]) - Float64(xs[i])
        dx > 0 || throw(ArgumentError("curve breakpoints must be strictly increasing"))
        slope = (Float64(ys[i + 1]) - Float64(ys[i])) / dx
        slope == 0.0 && continue
        # Turn each segment slope on and off. Without the negative term at the
        # segment end, slopes accumulate instead of reproducing a clamped
        # piecewise-linear curve.
        push!(triples, (slope, Float64(xs[i])))
        push!(triples, (-slope, Float64(xs[i + 1])))
    end
    baseline, triples
end

function _bilevel_curve_expression(op, voltage, baseline, triples)
    baseline + sum(slope * op(voltage - breakpoint)
                   for (slope, breakpoint) in triples)
end

function _bilevel_single_phase_voltage_reference(curve, law::AbstractString)
    raw = uppercase(String(get(curve, "voltage_reference", "PN_PER_PHASE")))
    raw in ("PN_PER_PHASE", "PG_PER_PHASE", "PP_PER_PHASE",
            "PN_AVERAGED", "PG_AVERAGED", "PP_AVERAGED") || throw(ArgumentError(
        "$law voltage_reference '$raw' is not supported"))
    startswith(raw, "PG") ? :single_pg : :single_diff
end

function _bilevel_local_controller_constraints!(ctx, net, pv_ids, relu_eps)
    model = _opf_model(ctx)
    bases = _opf_bases(ctx)
    bases === nothing && throw(ArgumentError(
        "local-controller mode requires per-unit OPF bases"))
    for id in pv_ids
        inverter = net["ibr"][id]
        topology = uppercase(String(get(inverter, "topology", "")))
        topology == "SINGLE_PHASE" || throw(ArgumentError(
            "local-controller mode currently supports SINGLE_PHASE PV '$id' only"))
        profile_id = get(inverter, "control_profile", nothing)
        profile_id === nothing && throw(ArgumentError(
            "PV '$id' must specify control_profile for local-controller mode"))
        profile = get(net["control_profile"], String(profile_id), nothing)
        profile === nothing && throw(ArgumentError(
            "control profile '$profile_id' for PV '$id' was not found"))
        haskey(profile, "volt_var") && haskey(profile, "volt_watt") ||
            throw(ArgumentError("PV '$id' local-controller mode requires both volt_var and volt_watt curves"))
        vv = profile["volt_var"]
        vw = profile["volt_watt"]
        bus = String(inverter["bus"])
        haskey(bases.v_base, bus) || throw(ArgumentError("no voltage base for PV bus '$bus'"))
        vbase = Float64(bases.v_base[bus])
        sbase = Float64(bases.s_base)

        raw_vbreaks = get(vv, "breakpoints", nothing)
        raw_pbreaks = get(vw, "breakpoints", nothing)
        raw_vbreaks isa AbstractVector && length(raw_vbreaks) == 4 ||
            throw(ArgumentError("PV '$id' Volt-var requires four breakpoints"))
        raw_pbreaks isa AbstractVector && length(raw_pbreaks) == 2 ||
            throw(ArgumentError("PV '$id' Volt-watt requires two breakpoints"))
        vbreaks = Float64.(raw_vbreaks) ./ vbase
        raw_q_limits = get(vv, "q_limits", nothing)
        raw_q_limits isa AbstractVector && length(raw_q_limits) == 2 ||
            throw(ArgumentError("PV '$id' Volt-var requires two q_limits"))
        q_limits = Float64.(raw_q_limits)
        length(q_limits) == 2 || throw(ArgumentError(
            "PV '$id' q_limits must have two values"))
        abs(q_limits[1]) <= 1.0 && abs(q_limits[2]) <= 1.0 ||
            throw(ArgumentError("PV '$id' q_limits must be normalized fractions in [-1, 1]"))
        q_limits[1] <= q_limits[2] || throw(ArgumentError(
            "PV '$id' q_limits must be ordered as [absorb, inject]"))

        pbreaks = Float64.(raw_pbreaks) ./ vbase
        raw_p_limits = get(vw, "p_limits", nothing)
        raw_p_limits isa AbstractVector && length(raw_p_limits) == 2 ||
            throw(ArgumentError("PV '$id' Volt-watt requires two p_limits"))
        p_limits = Float64.(raw_p_limits)
        length(p_limits) == 2 || throw(ArgumentError(
            "PV '$id' p_limits must have two values"))
        p_limits[1] <= p_limits[2] || throw(ArgumentError(
            "PV '$id' p_limits must be ordered as [p_low, p_high]"))
        p_values = [p_limits[2], p_limits[1]]
        p_baseline, p_triples = _bilevel_breakpoints_to_triples(pbreaks, p_values)
        _bilevel_single_phase_voltage_reference(vv, "Volt-var")
        p_voltage_key = BMOPFTools.opf_ibr_voltage_magnitude_key(
            id, 1; reference=_bilevel_single_phase_voltage_reference(vw, "Volt-watt"),
            controller=:single)
        p_voltage = BMOPFTools.opf_object(ctx, p_voltage_key)
        curve_eps = Float64(relu_eps) * sum(pbreaks) / length(pbreaks)
        curve_eps > 0 || throw(ArgumentError("local controller curve smoothing must be > 0"))
        op = _BilevelBuiltinSoftplus(curve_eps)

        s_max = Float64(inverter["s_max"][1]) / sbase
        s_max > 0 || throw(ArgumentError("PV '$id' s_max must be > 0"))
        p_max = Float64(inverter["p_max"][1]) / sbase
        p_ref = Symbol(uppercase(String(get(vw, "p_ref", "S_MAX"))))
        p_base = if p_ref in (:S_MAX, :VAR_MAX)
            s_max
        elseif p_ref == :P_MAX
            Float64(inverter["p_max"][1]) / sbase
        elseif p_ref == :P_AVAILABLE
            Float64(get(inverter, "p_available", inverter["p_max"])[1]) / sbase
        else
            throw(ArgumentError("unsupported local Volt-watt p_ref '$p_ref' for PV '$id'"))
        end
        # The native Volt-var equality is authoritative. Reusing its registered
        # reactive-power auxiliary prevents the local priority allocator from
        # reserving headroom against a second, drifting curve evaluation.
        q = BMOPFTools.opf_object(ctx,
            BMOPFTools.opf_ibr_power_key(id, 1; component=:reactive))
        p_curve = p_base * _bilevel_curve_expression(
            op, p_voltage, p_baseline, p_triples)
        # Volt-var priority reserves apparent power for reactive support. The
        # smooth minimum is below both arguments, preserving native PV bounds
        # while remaining differentiable for DiffOpt.
        q_cap = JuMP.@expression(model, sqrt(s_max^2 - q^2 + 1e-10))
        p_cap = JuMP.@expression(model, 0.5 * (p_curve + p_max -
            sqrt((p_curve - p_max)^2 + 1e-12)))
        p_target = JuMP.@expression(model, 0.5 * (p_cap + q_cap -
            sqrt((p_cap - q_cap)^2 + 1e-12)))
        p = BMOPFTools.opf_object(ctx, BMOPFTools.opf_ibr_power_key(id, 1))
        JuMP.@constraint(model, p == p_target)
    end
    nothing
end

function _bilevel_validation(net, transformer_id, pv_ids, monitored_buses)
    family, _ = _bilevel_transformer(net, transformer_id)
    family == "single_phase" || throw(ArgumentError(
        "the proof of concept currently supports single_phase transformers"))
    isempty(pv_ids) && throw(ArgumentError("pv_ids must not be empty"))
    isempty(monitored_buses) && throw(ArgumentError("monitored_buses must not be empty"))
    haskey(net, "ibr") || throw(ArgumentError("net must contain native IBR PV systems"))
    for id in pv_ids
        haskey(net["ibr"], id) || throw(ArgumentError("IBR PV '$id' was not found"))
    end
    nothing
end

function _bilevel_voltage_data(net, buses, voltage_reference, voltage_band)
    isfinite(voltage_reference) || throw(ArgumentError(
        "voltage_reference must be finite"))
    isfinite(voltage_band) || throw(ArgumentError("voltage_band must be finite"))
    voltage_band > 0 || throw(ArgumentError("voltage_band must be > 0"))
    for bus in buses
        data = get(net["bus"], bus, nothing)
        data === nothing && throw(ArgumentError("monitored bus '$bus' was not found"))
        vmin = get(data, "v_min", [voltage_reference - voltage_band])[1]
        vmax = get(data, "v_max", [voltage_reference + voltage_band])[1]
        isfinite(vmin) && isfinite(vmax) && vmin < vmax || throw(ArgumentError(
            "bus '$bus' must have a finite voltage interval"))
    end
    nothing
end

function _bilevel_status(model; throw_on_failure=true)
    status = JuMP.termination_status(model)
    solved = status in (JuMP.MOI.LOCALLY_SOLVED, JuMP.MOI.OPTIMAL)
    !solved && throw_on_failure && throw(ErrorException(
        "bilevel lower/centralized solve failed with $status"))
    (status=string(status), solved=solved)
end

function _bilevel_snapshot(model, p_handles, vm_handles, voltage_scales, power_scale)
    exported = power_scale * sum(JuMP.value(p) for p in values(p_handles))
    voltages = Dict(bus => voltage_scales[bus] * JuMP.value(vm)
                    for (bus, vm) in vm_handles)
    exported, voltages
end

function _bilevel_utility_value_gradient(tap, voltages, dvoltages;
                                         voltage_reference, voltage_band,
                                         voltage_weight, tap_penalty)
    stress = 0.0
    dstress = 0.0
    for (bus, voltage) in voltages
        x = (voltage - voltage_reference) / voltage_band
        stress += x^4
        dstress += 4x^3 / voltage_band * get(dvoltages, bus, 0.0)
    end
    value = voltage_weight * stress + tap_penalty * (tap - 1.0)^2
    gradient = voltage_weight * dstress + 2tap_penalty * (tap - 1.0)
    value, gradient
end

function _bilevel_tap_evaluation(net, transformer_id, pv_ids, monitored_buses,
                                 tap, tap_bounds; lower_level, relu_eps,
                                 voltage_reference, voltage_band, voltage_weight,
                                 tap_penalty, optimizer, verbose)
    state = _bilevel_lower_state(net, transformer_id, pv_ids, monitored_buses,
                                 tap, tap_bounds; lower_level,
                                 relu_eps, optimizer, verbose)
    snap = _bilevel_lower_solve!(state, tap; differentiate=true,
                                 throw_on_failure=false)
    snap.solved || return (solved=false, tap=Float64(tap), state=state,
                           snap=snap, objective=Inf, gradient=NaN)
    objective, gradient = _bilevel_utility_value_gradient(
        tap, snap.voltages, snap.derivatives;
        voltage_reference, voltage_band, voltage_weight, tap_penalty)
    (solved=true, tap=Float64(tap), state=state, snap=snap,
     objective=objective, gradient=gradient)
end

"""Mutable state for the differentiated lower-level response."""
struct _BilevelLowerState
    model::JuMP.Model
    context::Any
    tap_parameter::JuMP.VariableRef
    p_handles::Dict{String,Any}
    vm_handles::Dict{String,Any}
    voltage_scales::Dict{String,Float64}
    power_scale::Float64
end

function _bilevel_lower_state(net, transformer_id, pv_ids, monitored_buses,
                              tap_initial, tap_bounds; lower_level,
                              relu_eps, optimizer, verbose)
    local_net = _bilevel_copy_with_tap_bounds(
        net, transformer_id, tap_bounds, tap_initial)
    model = DiffOpt.nonlinear_diff_model(optimizer)
    verbose || JuMP.set_silent(model)
    ctx = build_opf_model(local_net; model=model, per_unit=true,
        add_objective=false, softplus=:builtin, volt_var_watt_eps=relu_eps)
    parameter = JuMP.@variable(model, utility_tap in JuMP.Parameter(tap_initial))
    tap_key = BMOPFTools.opf_transformer_tap_key(transformer_id)
    BMOPFTools.bind_opf_parameter!(ctx,
        BMOPFTools.OpfModelKey(:parameter, :utility_tap, transformer_id),
        parameter, tap_key; aliases=[:utility_tap],
        input_unit=:tap_multiplier, working_unit=:effective_turns_ratio,
        to_working_scale=1.0, owner=:PowerOptLabBilevel)
    p_handles = _bilevel_pv_handles(ctx, pv_ids)
    vm_handles = _bilevel_voltage_handles(ctx, local_net, monitored_buses)
    if lower_level == :aggregate
        JuMP.@objective(model, Max, sum(values(p_handles)))
    elseif lower_level == :local_controller
        _bilevel_local_controller_constraints!(ctx, local_net, pv_ids, relu_eps)
        JuMP.@objective(model, Min, 0.0)
    else
        throw(ArgumentError("lower_level must be :aggregate or :local_controller"))
    end
    BMOPFTools.enforce_kcl!(ctx)
    bases = _opf_bases(ctx)
    bases === nothing && throw(ArgumentError(
        "bilevel lower-level model requires per-unit OPF bases"))
    _, transformer = _bilevel_transformer(net, transformer_id)
    n0 = Float64(transformer["v_nom_from"]) / Float64(transformer["v_nom_to"])
    base_ratio = bases.v_base[String(transformer["bus_from"])] /
                 bases.v_base[String(transformer["bus_to"])]
    isapprox(base_ratio, n0; rtol=1e-10, atol=1e-12) || throw(ArgumentError(
        "per-unit voltage-base propagation is inconsistent with transformer " *
        "nominal turns ratio; the DiffOpt tap binding requires per_unit=true"))
    voltage_scales = Dict(key => bases.v_base[first(split(key, ":"))]
                          for key in keys(vm_handles))
    power_scale = bases.s_base
    _BilevelLowerState(model, ctx, parameter, p_handles, vm_handles,
                       voltage_scales, power_scale)
end

function _bilevel_lower_solve!(state::_BilevelLowerState, tap;
                               differentiate=false, throw_on_failure=true)
    JuMP.set_parameter_value(state.tap_parameter, tap)
    JuMP.optimize!(state.model)
    solve = _bilevel_status(state.model; throw_on_failure=throw_on_failure)
    solve.solved || return (status=solve.status, solved=false, exported=NaN,
                            voltages=Dict{String,Float64}(),
                            derivatives=Dict{String,Float64}(), error=nothing)
    exported, voltages = _bilevel_snapshot(state.model, state.p_handles,
        state.vm_handles, state.voltage_scales, state.power_scale)
    derivatives = Dict{String,Float64}()
    if differentiate
        try
            JuMP.MOI.set(state.model, DiffOpt.NonLinearKKTJacobianFactorization(),
                BMOPFTools.opf_checked_kkt_factorization(state.context))
            DiffOpt.empty_input_sensitivities!(state.model)
            DiffOpt.set_forward_parameter(state.model, state.tap_parameter, 1.0)
            DiffOpt.forward_differentiate!(state.model)
            for (bus, vm) in state.vm_handles
                # The live OPF is per-unit, while the public response and outer
                # utility objective are expressed in volts per tap multiplier.
                derivatives[bus] = state.voltage_scales[bus] *
                                   DiffOpt.get_forward_variable(state.model, vm)
            end
        catch err
            throw_on_failure && rethrow()
            return (status=solve.status, solved=false, exported=exported,
                    voltages=voltages, derivatives=Dict{String,Float64}(),
                    error=err)
        end
        DiffOpt.empty_input_sensitivities!(state.model)
    end
    (status=solve.status, solved=true, exported=exported, voltages=voltages,
     derivatives=derivatives, error=nothing)
end

"""
    solve_bilevel_pv_tap(net; transformer_id, pv_ids, monitored_buses,
        tap_initial=1.0, tap_bounds=(0.95, 1.05), ...)

Solve the POC hierarchy. With `lower_level=:aggregate`, the lower level
maximises aggregate PV active export with native BMOPFTools Volt-var/Volt-watt
control curves active. With `lower_level=:local_controller`, each PV follows
its own Volt-var/Volt-watt response and the lower-level objective is constant;
the extra equations are differentiated as part of the network equilibrium.
The upper level minimises a smooth fourth-power voltage-stress metric plus tap
movement, using DiffOpt's implicit KKT derivative and a safeguarded
one-dimensional derivative search.

This is intentionally a local, smooth-region experiment: it is not a global
solution method for a nonconvex bilevel problem. Each upper-level trial builds
and solves a fresh lower-level model, which avoids making the reported response
path depend on a rejected warm start. At a Volt-watt kink or active-set
transition, inspect `differentiability_report` and validate the final point by
perturbing the tap and resolving the lower level.
"""
function solve_bilevel_pv_tap(net::Dict{String,Any};
                              transformer_id::AbstractString,
                              pv_ids::AbstractVector{<:AbstractString},
                              monitored_buses::AbstractVector{<:AbstractString},
                              lower_level::Symbol=:aggregate,
                              tap_initial::Real=1.0,
                              tap_bounds=(0.95, 1.05),
                              voltage_reference::Real=230.0,
                              voltage_band::Real=15.0,
                              voltage_weight::Real=1.0,
                              tap_penalty::Real=0.05,
                              gradient_tolerance::Real=1e-6,
                              tap_tolerance::Real=1e-5,
                              max_iterations::Integer=16,
                              volt_var_watt_eps::Real=2e-3,
                              optimizer=Ipopt.Optimizer,
                              verbose::Bool=false)
    tid = String(transformer_id)
    pids = String.(pv_ids)
    buses = String.(monitored_buses)
    _bilevel_validation(net, tid, pids, buses)
    lower_level in (:aggregate, :local_controller) || throw(ArgumentError(
        "lower_level must be :aggregate or :local_controller"))
    _bilevel_voltage_data(net, buses, Float64(voltage_reference),
                          Float64(voltage_band))
    voltage_weight >= 0 || throw(ArgumentError("voltage_weight must be >= 0"))
    tap_penalty >= 0 || throw(ArgumentError("tap_penalty must be >= 0"))
    gradient_tolerance > 0 || throw(ArgumentError(
        "gradient_tolerance must be > 0"))
    tap_tolerance > 0 || throw(ArgumentError("tap_tolerance must be > 0"))
    max_iterations >= 1 || throw(ArgumentError("max_iterations must be >= 1"))
    volt_var_watt_eps > 0 || throw(ArgumentError("volt_var_watt_eps must be > 0"))
    lo, hi = Float64.(tap_bounds)
    lo < hi || throw(ArgumentError("tap_bounds must be finite and strictly ordered"))
    isfinite(lo) && isfinite(hi) || throw(ArgumentError(
        "tap_bounds must be finite and strictly ordered"))
    tap_margin = min(1e-4, 0.25 * (hi - lo))
    search_lo, search_hi = lo + tap_margin, hi - tap_margin
    history = NamedTuple[]

    evaluate(t) = _bilevel_tap_evaluation(
        net, tid, pids, buses, t, (lo, hi); lower_level,
        relu_eps=Float64(volt_var_watt_eps),
        voltage_reference=Float64(voltage_reference),
        voltage_band=Float64(voltage_band), voltage_weight=Float64(voltage_weight),
        tap_penalty=Float64(tap_penalty), optimizer, verbose)
    record!(evaluation, iteration) = evaluation.solved && push!(history,
        (iteration=iteration, tap=evaluation.tap, objective=evaluation.objective,
         gradient=evaluation.gradient, exported_power_W=evaluation.snap.exported))

    tap = clamp(Float64(tap_initial), search_lo, search_hi)
    current = evaluate(tap)
    record!(current, 1)
    if !current.solved
        # An infeasible requested starting tap is a normal operating possibility
        # (for example, a low tap can push a PV bus above its voltage limit).
        # Try the interior midpoint before declaring the whole upper problem
        # infeasible.
        for (i, fallback_tap) in zip(2:4, (0.5 * (search_lo + search_hi),
                                           search_lo, search_hi))
            fallback = evaluate(fallback_tap)
            if fallback.solved
                current = fallback
                tap = fallback_tap
                record!(current, i)
                break
            end
        end
    end
    current.solved || throw(ErrorException(
        "no feasible lower-level response was found inside tap_bounds"))

    points = Any[current]
    for endpoint in (search_lo, search_hi)
        endpoint == current.tap && continue
        evaluation = evaluate(endpoint)
        evaluation.solved && (push!(points, evaluation); record!(evaluation, length(history) + 1))
    end
    sort!(points, by=evaluation -> evaluation.tap)

    converged = abs(current.gradient) <= gradient_tolerance
    termination_reason = converged ? :gradient_tolerance : :iteration_limit
    selected = current
    if !converged
        bracket = nothing
        for i in 1:(length(points) - 1)
            left, right = points[i], points[i + 1]
            left.gradient == 0.0 && (selected = left; bracket = nothing; converged = true;
                                     termination_reason = :gradient_tolerance; break)
            right.gradient == 0.0 && (selected = right; bracket = nothing; converged = true;
                                      termination_reason = :gradient_tolerance; break)
            left.gradient * right.gradient < 0.0 && (bracket = (left, right); break)
        end
        if !converged && bracket !== nothing
            left, right = bracket
            for iteration in 1:max_iterations
                width = right.tap - left.tap
                width <= tap_tolerance && break
                # Bisection is deliberately used as the safeguarded core. The
                # lower-level response is nonconvex, so a fast open Newton step
                # can jump across an active-set transition; bisection preserves
                # the derivative sign bracket and has a predictable stopping
                # guarantee.
                candidate_tap = 0.5 * (left.tap + right.tap)
                candidate = evaluate(candidate_tap)
                if !candidate.solved
                    candidate_tap = 0.5 * (left.tap + candidate_tap)
                    candidate = evaluate(candidate_tap)
                end
                if !candidate.solved
                    termination_reason = :lower_level_failure
                    break
                end
                push!(points, candidate)
                record!(candidate, length(history) + 1)
                selected = candidate
                if abs(candidate.gradient) <= gradient_tolerance
                    converged = true
                    termination_reason = :gradient_tolerance
                    break
                elseif candidate.gradient * left.gradient < 0.0
                    right = candidate
                else
                    left = candidate
                end
                selected = abs(left.gradient) < abs(right.gradient) ? left : right
                if right.tap - left.tap <= tap_tolerance
                    converged = true
                    termination_reason = :tap_tolerance
                    break
                end
                iteration == max_iterations && (termination_reason = :iteration_limit)
            end
        elseif !converged
            selected = argmin(evaluation -> evaluation.objective, points)
            at_lower = selected.tap == search_lo && selected.gradient >= 0.0
            at_upper = selected.tap == search_hi && selected.gradient <= 0.0
            if at_lower || at_upper
                converged = true
                termination_reason = :at_bound
            else
                termination_reason = :no_bracket
            end
        end
    end

    snap = selected.snap
    lower_result = BMOPFTools.extract_result(selected.state.context)
    BilevelPVResult(snap.status, lower_level, selected.tap, snap.exported,
        snap.voltages, selected.objective, history, lower_result,
        BMOPFTools.opf_differentiability_report(selected.state.context),
        converged, termination_reason)
end

"""
    solve_single_level_pv_tap(net; transformer_id, pv_ids, monitored_buses, ...)

Solve the centralized comparison: utility tap and all PV operating points are
chosen in one BMOPFTools/JuMP OPF. The objective is the same voltage stress and
tap movement used by the bilevel upper level. An optional normalized PV export
reward can be enabled explicitly with `export_weight`. This is a benchmark for
the coordination gap, not the consumer hierarchy.
"""
function solve_single_level_pv_tap(net::Dict{String,Any};
                                   transformer_id::AbstractString,
                                   pv_ids::AbstractVector{<:AbstractString},
                                   monitored_buses::AbstractVector{<:AbstractString},
                                   tap_initial::Real=1.0,
                                   tap_bounds=(0.95, 1.05),
                                   voltage_reference::Real=230.0,
                                   voltage_band::Real=15.0,
                                   voltage_weight::Real=1.0,
                                   tap_penalty::Real=0.05,
                                   export_weight::Real=0.0,
                                   export_scale_W::Real=10_000.0,
                                   volt_var_watt_eps::Real=2e-3,
                                   optimizer=Ipopt.Optimizer,
                                   verbose::Bool=false)
    tid = String(transformer_id); pids = String.(pv_ids); buses = String.(monitored_buses)
    _bilevel_validation(net, tid, pids, buses)
    _bilevel_voltage_data(net, buses, Float64(voltage_reference),
                          Float64(voltage_band))
    export_scale_W > 0 || throw(ArgumentError("export_scale_W must be > 0"))
    voltage_weight >= 0 || throw(ArgumentError("voltage_weight must be >= 0"))
    tap_penalty >= 0 || throw(ArgumentError("tap_penalty must be >= 0"))
    export_weight >= 0 || throw(ArgumentError("export_weight must be >= 0"))
    volt_var_watt_eps > 0 || throw(ArgumentError("volt_var_watt_eps must be > 0"))
    local_net = _bilevel_copy_with_tap_bounds(net, tid, tap_bounds,
                                               Float64(tap_initial))
    model = JuMP.Model(optimizer)
    verbose || JuMP.set_silent(model)
    ctx = BMOPFTools.build_opf_model(local_net; model=model, per_unit=true,
        add_objective=false, softplus=:builtin,
        volt_var_watt_eps=Float64(volt_var_watt_eps))
    p_handles = _bilevel_pv_handles(ctx, pids)
    vm_handles = _bilevel_voltage_handles(ctx, local_net, buses)
    tap_multiplier = BMOPFTools.opf_object(ctx,
        BMOPFTools.opf_transformer_tap_key(tid))
    bases = _opf_bases(ctx)
    bases === nothing && throw(ArgumentError(
        "centralized bilevel benchmark requires per-unit OPF bases"))
    voltage_scales = Dict(key => bases.v_base[first(split(key, ":"))]
                          for key in keys(vm_handles))
    power_scale = bases.s_base
    voltage_expr(key) = voltage_scales[key] * vm_handles[key]
    stress = sum(((voltage_expr(key) - Float64(voltage_reference)) /
                  Float64(voltage_band))^4 for key in keys(vm_handles))
    export_expr = power_scale * sum(values(p_handles)) / Float64(export_scale_W)
    JuMP.@objective(model, Min,
        Float64(voltage_weight) * stress +
        Float64(tap_penalty) * (tap_multiplier - 1.0)^2 -
        Float64(export_weight) * export_expr)
    BMOPFTools.enforce_kcl!(ctx)
    JuMP.optimize!(model)
    status = _bilevel_status(model).status
    exported, voltages = _bilevel_snapshot(model, p_handles, vm_handles,
                                            voltage_scales, power_scale)
    result = BMOPFTools.extract_result(ctx)
    SingleLevelPVResult(status, JuMP.value(tap_multiplier), exported, voltages,
        JuMP.objective_value(model), result)
end

"""Return a small reproducible single-phase feeder for the bilevel tutorial."""
function bilevel_demo_network(; source_voltage::Real=11_000.0)
    v = Float64(source_voltage)
    parse_bmopf("""
    {"bus":{
        "hv":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
        "lv0":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"]},
        "lv1":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],"v_min":[216.0],"v_max":[244.0]},
        "lv2":{"terminal_names":["1","n"],"perfectly_grounded_terminals":["n"],"v_min":[216.0],"v_max":[244.0]}},
     "voltage_source":{"grid":{"bus":"hv","terminal_map":["1"],"v_magnitude":[$v],"v_angle":[0.0]}},
     "transformer":{"single_phase":{"reg":{"bus_from":"hv","bus_to":"lv0",
         "terminal_map_from":["1","n"],"terminal_map_to":["1","n"],
         "v_nom_from":$v,"v_nom_to":240.0,"s_rating":100000.0,"tap":1.0}}},
     "linecode":{"lc":{"R_series_1_1":0.08,"X_series_1_1":0.03}},
     "line":{
       "l1":{"bus_from":"lv0","bus_to":"lv1","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0},
       "l2":{"bus_from":"lv1","bus_to":"lv2","terminal_map_from":["1"],"terminal_map_to":["1"],"linecode":"lc","length":1.0}},
     "load":{
       "load1":{"bus":"lv1","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[3000.0],"q_nom":[900.0]},
       "load2":{"bus":"lv2","terminal_map":["1","n"],"configuration":"SINGLE_PHASE","p_nom":[2500.0],"q_nom":[700.0]}},
     "control_profile":{"vvw":{"volt_var":{"voltage_reference":"PN_PER_PHASE",
         "breakpoints":[216.0,225.0,235.0,244.0],"q_limits":[-0.44,0.44],"q_unit":"VA_FRACTION","q_ref":"VAR_MAX"},
       "volt_watt":{"voltage_reference":"PN_PER_PHASE","breakpoints":[236.0,244.0],
         "p_limits":[0.0,1.0],"p_unit":"VA_FRACTION","p_ref":"S_MAX"}}},
     "ibr":{
       "pv1":{"bus":"lv1","terminal_map":["1","n"],"topology":"SINGLE_PHASE","prime_mover":"PV",
         "s_max":[7000.0],"p_min":[0.0],"p_max":[6000.0],"q_min":[-7000.0],"q_max":[7000.0],"control_profile":"vvw"},
       "pv2":{"bus":"lv2","terminal_map":["1","n"],"topology":"SINGLE_PHASE","prime_mover":"PV",
         "s_max":[7000.0],"p_min":[0.0],"p_max":[6000.0],"q_min":[-7000.0],"q_max":[7000.0],"control_profile":"vvw"}}}
    """; from_string=true)
end
