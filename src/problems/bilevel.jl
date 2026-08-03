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
    return n, transformer, Float64(xfmr["v_nom_from"]) / Float64(xfmr["v_nom_to"])
end

function _bilevel_voltage_handles(ctx, buses::Vector{String})
    model = _opf_model(ctx)
    vr, vi = _opf_voltage_maps(ctx)
    bases = _opf_bases(ctx)
    voltage_start = bases === nothing ? 230.0 : 1.0
    handles = Dict{String,Any}()
    for bus in buses
        key = (bus, "1")
        haskey(vr, key) || throw(ArgumentError(
            "monitored bus '$bus' must have terminal '1'"))
        vm = JuMP.@variable(model, base_name="bilevel_vm_$(bus)", lower_bound=0.0)
        JuMP.set_start_value(vm, voltage_start)
        JuMP.@constraint(model, vm^2 == vr[key]^2 + vi[key]^2)
        handles[bus] = vm
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
        push!(triples, ((Float64(ys[i + 1]) - Float64(ys[i])) / dx,
                        Float64(xs[i])))
    end
    baseline, triples
end

function _bilevel_curve_expression(op, voltage, baseline, triples)
    baseline + sum(slope * op(voltage - breakpoint)
                   for (slope, breakpoint) in triples)
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

        vbreaks = Float64.(vv["breakpoints"]) ./ vbase
        q_limits = Float64.(vv["q_limits"])
        length(q_limits) == 2 || throw(ArgumentError(
            "PV '$id' q_limits must have two values"))
        abs(q_limits[1]) <= 1.0 && abs(q_limits[2]) <= 1.0 ||
            throw(ArgumentError("PV '$id' q_limits must be normalized fractions in [-1, 1]"))
        # BMOPFTools uses q injection at the low-voltage end and q absorption
        # at the high-voltage end, while q_limits is [absorb, inject].
        q_values = [q_limits[2], 0.0, 0.0, q_limits[1]]
        q_baseline, q_triples = _bilevel_breakpoints_to_triples(vbreaks, q_values)

        pbreaks = Float64.(vw["breakpoints"]) ./ vbase
        p_limits = Float64.(vw["p_limits"])
        length(p_limits) == 2 || throw(ArgumentError(
            "PV '$id' p_limits must have two values"))
        p_values = [p_limits[2], p_limits[1]]
        p_baseline, p_triples = _bilevel_breakpoints_to_triples(pbreaks, p_values)
        curve_eps = Float64(relu_eps) * sum(vbreaks) / length(vbreaks)
        curve_eps > 0 || throw(ArgumentError("local controller curve smoothing must be > 0"))
        op = _BilevelBuiltinSoftplus(curve_eps)
        voltage_key = BMOPFTools.opf_ibr_voltage_magnitude_key(
            id, 1; reference=:single_diff, controller=:single)
        voltage = BMOPFTools.opf_object(ctx, voltage_key)

        s_max = Float64(inverter["s_max"][1]) / sbase
        s_max > 0 || throw(ArgumentError("PV '$id' s_max must be > 0"))
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
        q_curve = s_max * _bilevel_curve_expression(
            op, voltage, q_baseline, q_triples)
        p_curve = p_base * _bilevel_curve_expression(
            op, voltage, p_baseline, p_triples)
        # Volt-var priority reserves apparent power for reactive support. The
        # smooth minimum is below both arguments, preserving native PV bounds
        # while remaining differentiable for DiffOpt.
        q_cap = JuMP.@expression(model, sqrt(s_max^2 - q_curve^2 + 1e-10))
        p_target = JuMP.@expression(model,
            0.5 * (p_curve + q_cap - sqrt((p_curve - q_cap)^2 + 1e-12)))
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
    refs = Dict{String,Float64}()
    bands = Dict{String,Float64}()
    for bus in buses
        data = get(net["bus"], bus, nothing)
        data === nothing && throw(ArgumentError("monitored bus '$bus' was not found"))
        vmin = get(data, "v_min", [voltage_reference - voltage_band])[1]
        vmax = get(data, "v_max", [voltage_reference + voltage_band])[1]
        refs[bus] = Float64(voltage_reference)
        bands[bus] = Float64(voltage_band)
        isfinite(vmin) && isfinite(vmax) && vmin < vmax || throw(ArgumentError(
            "bus '$bus' must have a finite voltage interval"))
        bands[bus] > 0 || throw(ArgumentError("voltage_band must be > 0"))
    end
    refs, bands
end

function _bilevel_status(model)
    status = JuMP.termination_status(model)
    status in (JuMP.MOI.LOCALLY_SOLVED, JuMP.MOI.OPTIMAL) ||
        throw(ErrorException("bilevel lower/centralized solve failed with $status"))
    string(status)
end

function _bilevel_snapshot(model, p_handles, vm_handles, voltage_scales, power_scale)
    exported = power_scale * sum(JuMP.value(p) for p in values(p_handles))
    voltages = Dict(bus => voltage_scales[bus] * JuMP.value(vm)
                    for (bus, vm) in vm_handles)
    exported, voltages
end

function _bilevel_utility_value_gradient(tap, _exported, voltages, dvoltages;
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
    local_net, _, _ = _bilevel_copy_with_tap_bounds(
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
    vm_handles = _bilevel_voltage_handles(ctx, monitored_buses)
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
    voltage_scales = Dict(bus => bases === nothing ? 1.0 : bases.v_base[bus]
                          for bus in monitored_buses)
    power_scale = bases === nothing ? 1.0 : bases.s_base
    _BilevelLowerState(model, ctx, parameter, p_handles, vm_handles,
                       voltage_scales, power_scale)
end

function _bilevel_lower_solve!(state::_BilevelLowerState, tap; differentiate=false)
    JuMP.set_parameter_value(state.tap_parameter, tap)
    JuMP.optimize!(state.model)
    status = _bilevel_status(state.model)
    exported, voltages = _bilevel_snapshot(state.model, state.p_handles,
        state.vm_handles, state.voltage_scales, state.power_scale)
    derivatives = Dict{String,Float64}()
    if differentiate
        JuMP.MOI.set(state.model, DiffOpt.NonLinearKKTJacobianFactorization(),
            BMOPFTools.opf_checked_kkt_factorization(state.context))
        DiffOpt.empty_input_sensitivities!(state.model)
        DiffOpt.set_forward_parameter(state.model, state.tap_parameter, 1.0)
        DiffOpt.forward_differentiate!(state.model)
        for (bus, vm) in state.vm_handles
            derivatives[bus] = DiffOpt.get_forward_variable(state.model, vm)
        end
        DiffOpt.empty_input_sensitivities!(state.model)
    end
    (status=status, exported=exported, voltages=voltages, derivatives=derivatives)
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
movement, using DiffOpt's implicit KKT derivative and a projected/backtracked
gradient search.

This is intentionally a local, smooth-region experiment: it is not a global
solution method for a nonconvex bilevel problem. At a Volt-watt kink or active-set
transition, rerun with a smaller step, inspect `differentiability_report`, and
validate the final point by perturbing the tap and resolving the lower level.
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
                              step_size::Real=0.01,
                              backtrack::Real=0.5,
                              armijo::Real=1e-4,
                              gradient_tolerance::Real=1e-6,
                              max_iterations::Integer=20,
                              volt_var_watt_eps::Real=2e-3,
                              optimizer=Ipopt.Optimizer,
                              verbose::Bool=false)
    tid = String(transformer_id)
    pids = String.(pv_ids)
    buses = String.(monitored_buses)
    _bilevel_validation(net, tid, pids, buses)
    lower_level in (:aggregate, :local_controller) || throw(ArgumentError(
        "lower_level must be :aggregate or :local_controller"))
    refs, bands = _bilevel_voltage_data(net, buses, Float64(voltage_reference),
                                        Float64(voltage_band))
    0 < step_size && step_size < 1 || throw(ArgumentError("step_size must be in (0,1)"))
    0 < backtrack < 1 || throw(ArgumentError("backtrack must be in (0,1)"))
    max_iterations >= 1 || throw(ArgumentError("max_iterations must be >= 1"))
    volt_var_watt_eps > 0 || throw(ArgumentError("volt_var_watt_eps must be > 0"))
    lo, hi = Float64.(tap_bounds)
    tap = Float64(tap_initial)
    state = _bilevel_lower_state(net, tid, pids, buses, tap, (lo, hi);
                                 lower_level, relu_eps=Float64(volt_var_watt_eps),
                                 optimizer, verbose)
    snap = _bilevel_lower_solve!(state, tap; differentiate=true)
    history = NamedTuple[]
    for iteration in 1:max_iterations
        value, gradient = _bilevel_utility_value_gradient(tap, snap.exported,
            snap.voltages, snap.derivatives; voltage_reference=Float64(voltage_reference),
            voltage_band=Float64(voltage_band), voltage_weight=Float64(voltage_weight),
            tap_penalty=Float64(tap_penalty))
        push!(history, (iteration=iteration, tap=tap, objective=value,
                        gradient=gradient, exported_power_W=snap.exported))
        abs(gradient) <= gradient_tolerance && break
        candidate_step = Float64(step_size)
        accepted = false
        while candidate_step > 1e-10
            candidate_tap = clamp(tap - candidate_step * gradient, lo, hi)
            candidate_tap == tap && break
            candidate = _bilevel_lower_solve!(state, candidate_tap)
            candidate_value, _ = _bilevel_utility_value_gradient(candidate_tap,
                candidate.exported, candidate.voltages, Dict{String,Float64}();
                voltage_reference=Float64(voltage_reference), voltage_band=Float64(voltage_band),
                voltage_weight=Float64(voltage_weight), tap_penalty=Float64(tap_penalty))
            candidate_value <= value - armijo * candidate_step * gradient^2 ||
                candidate_value < value || (candidate_step *= backtrack; continue)
            tap = candidate_tap
            snap = _bilevel_lower_solve!(state, tap; differentiate=true)
            accepted = true
            break
        end
        accepted || break
    end
    # A rejected trial leaves the live lower model at that trial parameter.
    # Re-solve the accepted tap before extracting the final network result.
    snap = _bilevel_lower_solve!(state, tap; differentiate=true)
    final_value, _ = _bilevel_utility_value_gradient(tap, snap.exported,
        snap.voltages, Dict{String,Float64}(); voltage_reference=Float64(voltage_reference),
        voltage_band=Float64(voltage_band), voltage_weight=Float64(voltage_weight),
        tap_penalty=Float64(tap_penalty))
    lower_result = BMOPFTools.extract_result(state.context)
    BilevelPVResult(snap.status, lower_level, tap, snap.exported, snap.voltages, final_value,
        history, lower_result, BMOPFTools.opf_differentiability_report(state.context))
end

"""
    solve_single_level_pv_tap(net; transformer_id, pv_ids, monitored_buses, ...)

Solve the centralized comparison: utility tap and all PV operating points are
chosen in one BMOPFTools/JuMP OPF. The objective is the same voltage stress and
tap movement used by the bilevel upper level, with an explicit normalized PV
export reward. This is a benchmark for the coordination gap, not the consumer
hierarchy.
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
                                   export_weight::Real=0.25,
                                   export_scale_W::Real=10_000.0,
                                   optimizer=Ipopt.Optimizer,
                                   verbose::Bool=false)
    tid = String(transformer_id); pids = String.(pv_ids); buses = String.(monitored_buses)
    _bilevel_validation(net, tid, pids, buses)
    refs, bands = _bilevel_voltage_data(net, buses, Float64(voltage_reference),
                                        Float64(voltage_band))
    export_scale_W > 0 || throw(ArgumentError("export_scale_W must be > 0"))
    local_net, _, _ = _bilevel_copy_with_tap_bounds(net, tid, tap_bounds,
                                                    Float64(tap_initial))
    model = JuMP.Model(optimizer)
    verbose || JuMP.set_silent(model)
    ctx = BMOPFTools.build_opf_model(local_net; model=model, per_unit=true,
        add_objective=false, softplus=:builtin)
    p_handles = _bilevel_pv_handles(ctx, pids)
    vm_handles = _bilevel_voltage_handles(ctx, buses)
    tap_multiplier = BMOPFTools.opf_object(ctx,
        BMOPFTools.opf_transformer_tap_key(tid))
    bases = _opf_bases(ctx)
    voltage_scales = Dict(bus => bases === nothing ? 1.0 : bases.v_base[bus]
                          for bus in buses)
    power_scale = bases === nothing ? 1.0 : bases.s_base
    voltage_expr(bus) = voltage_scales[bus] * vm_handles[bus]
    stress = sum(((voltage_expr(bus) - refs[bus]) / bands[bus])^4 for bus in buses)
    export_expr = power_scale * sum(values(p_handles)) / Float64(export_scale_W)
    JuMP.@objective(model, Min,
        Float64(voltage_weight) * stress +
        Float64(tap_penalty) * (tap_multiplier - 1.0)^2 -
        Float64(export_weight) * export_expr)
    BMOPFTools.enforce_kcl!(ctx)
    JuMP.optimize!(model)
    status = _bilevel_status(model)
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
