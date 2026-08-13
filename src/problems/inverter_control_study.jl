# Single-snapshot, network-scale phase-aware inverter-control studies.

"""
    ControlledInverterFleetSpec(devices, requests)

Configuration for replacing selected native BMOPFTools IBRs with
[`ControlledDevice`](@ref)s in one network snapshot. `devices` and `requests`
are dictionaries keyed by the native dataset IBR identifier. The two key sets
must be identical and each key must equal the wrapped device's `device_id`.

Unselected IBRs retain their native BMOPFTools formulation. Selected IBRs are
owned by PowerOptLab through BMOPFTools' per-component `OpfDeviceBuilder`
contract, so the advanced plant replaces rather than supplements the native
device.
"""
struct ControlledInverterFleetSpec
    devices::Dict{String,ControlledDevice}
    requests::Dict{String,InverterControlRequest}
    function ControlledInverterFleetSpec(
            devices::Dict{String,ControlledDevice},
            requests::Dict{String,InverterControlRequest})
        isempty(devices) && throw(ArgumentError(
            "a controlled-inverter fleet needs at least one device"))
        device_ids = Set(keys(devices))
        request_ids = Set(keys(requests))
        device_ids == request_ids || throw(ArgumentError(
            "controlled-inverter device and request identifiers must match; " *
            "devices=$(sort!(collect(device_ids))), " *
            "requests=$(sort!(collect(request_ids)))"))
        for id in sort!(collect(device_ids))
            device_id(devices[id]) == id || throw(ArgumentError(
                "fleet key '$id' does not match wrapped device id " *
                "'$(device_id(devices[id]))'"))
            devices[id].device isa AdvancedInverter || throw(ArgumentError(
                "fleet device '$id' must wrap an AdvancedInverter"))
        end
        new(copy(devices), copy(requests))
    end
end

function ControlledInverterFleetSpec(devices::AbstractDict,
                                     requests::AbstractDict)
    typed_devices = Dict{String,ControlledDevice}()
    for (id, controlled) in devices
        controlled isa ControlledDevice || throw(ArgumentError(
            "fleet device '$(id)' must be a ControlledDevice"))
        canonical_id = String(id)
        haskey(typed_devices, canonical_id) && throw(ArgumentError(
            "duplicate controlled-inverter identifier '$canonical_id' after " *
            "string conversion"))
        typed_devices[canonical_id] = controlled
    end
    typed_requests = Dict{String,InverterControlRequest}()
    for (id, request) in requests
        request isa InverterControlRequest || throw(ArgumentError(
            "fleet request '$(id)' must be an InverterControlRequest"))
        canonical_id = String(id)
        haskey(typed_requests, canonical_id) && throw(ArgumentError(
            "duplicate controlled-inverter request identifier " *
            "'$canonical_id' after string conversion"))
        typed_requests[canonical_id] = request
    end
    ControlledInverterFleetSpec(typed_devices, typed_requests)
end

"""Result of [`solve_controlled_inverter_fleet`](@ref)."""
struct ControlledInverterFleetResult <: AbstractSolveResult
    termination_status::String
    devices::Dict{String,ControlledInverterResult}
    spec::ControlledInverterFleetSpec
    network::Dict{String,Any}
    build_manifest::BMOPFTools.OpfBuildManifest
    solve::SolveStatus
end

solve_status(result::ControlledInverterFleetResult) = result.solve
solve_diagnostics(result::ControlledInverterFleetResult) = (
    controlled_device_count=length(result.devices),
    total_p_poc=sum(device.plant.p_poc for device in values(result.devices)),
    total_q_poc=sum(device.plant.q_poc for device in values(result.devices)),
    total_converter_loss=sum(
        device.plant.p_loss + device.plant.p_cap_loss
        for device in values(result.devices)),
    maximum_converter_current=maximum(
        maximum(device.plant.i_mag) for device in values(result.devices)),
    maximum_grid_current=maximum(
        maximum(device.plant.i_grid_mag) for device in values(result.devices)),
    maximum_exact_smooth_current_residual=maximum(
        device.exact_smooth_current_residual
        for device in values(result.devices)),
)

function _native_ibr_phase_count(data::AbstractDict,
                                 controlled::ControlledDevice)
    topology = uppercase(String(get(data, "topology", "FOUR_LEG")))
    terminals = String.(get(data, "terminal_map", String[]))
    inverter = inverter_spec(controlled)
    if topology == "THREE_LEG"
        terminals == inverter.phase_terminals || throw(ArgumentError(
            "native THREE_LEG IBR '$(inverter.id)' terminal_map $terminals " *
            "must equal the advanced inverter phase_terminals " *
            "$(inverter.phase_terminals) in the same sequence order"))
        inverter.neutral === nothing || throw(ArgumentError(
            "native THREE_LEG IBR '$(inverter.id)' is a line-to-line device; " *
            "its advanced replacement must set neutral=nothing"))
        !_lcl_active(inverter) || throw(ArgumentError(
            "native THREE_LEG IBR '$(inverter.id)' cannot be replaced by an " *
            "advanced inverter with a phase-to-ground LCL midpoint branch"))
        inverter.b_filter_shunt == 0.0 || throw(ArgumentError(
            "native THREE_LEG IBR '$(inverter.id)' cannot be replaced by an " *
            "advanced inverter with a phase-to-ground POC shunt"))
        return length(terminals)
    elseif topology == "FOUR_LEG"
        inverter.neutral === nothing && throw(ArgumentError(
            "native FOUR_LEG IBR '$(inverter.id)' needs an explicit advanced " *
            "inverter neutral terminal"))
        expected = [inverter.phase_terminals; inverter.neutral]
        terminals == expected || throw(ArgumentError(
            "native FOUR_LEG IBR '$(inverter.id)' terminal_map $terminals " *
            "must equal phase_terminals followed by neutral $expected"))
        return length(inverter.phase_terminals)
    end
    throw(ArgumentError(
        "native IBR '$(inverter.id)' has topology '$topology'; the initial " *
        "fleet replacement supports THREE_LEG and FOUR_LEG dataset records"))
end

function _validate_controlled_fleet(net::Dict{String,Any},
                                    spec::ControlledInverterFleetSpec)
    native_ibrs = get(net, "ibr", Dict())
    native_ibrs isa AbstractDict || throw(ArgumentError(
        "network field 'ibr' must be a dictionary"))
    phase_counts = Dict{String,Int}()
    for id in sort!(collect(keys(spec.devices)))
        haskey(native_ibrs, id) || throw(ArgumentError(
            "controlled inverter '$id' is not present in network['ibr']; " *
            "fleet control replaces an existing dataset IBR"))
        native = native_ibrs[id]
        native isa AbstractDict || throw(ArgumentError(
            "network IBR '$id' must be a dictionary"))
        if haskey(native, "dc_bus") || get(native, "dc_link_coupled", false) == true
            throw(ArgumentError(
                "controlled replacement of IBR '$id' is not supported while " *
                "its native record has DC-network coupling"))
        end
        controlled = spec.devices[id]
        inverter = inverter_spec(controlled)
        String(get(native, "bus", "")) == inverter.bus || throw(ArgumentError(
            "native IBR '$id' is connected to bus '$(get(native, "bus", ""))', " *
            "but its advanced replacement is connected to '$(inverter.bus)'"))
        validate_device(controlled, (net,); periods=1)
        phase_counts[id] = _native_ibr_phase_count(native, controlled)
        phase_counts[id] == 3 || throw(ArgumentError(
            "controlled fleet IBR '$id' must expose exactly three phase ports"))
    end
    phase_counts
end

# BMOPFTools declares native IBR current variables before device ownership is
# resolved. A custom owner must make the unused native placeholders harmless;
# fixing them also gives the study harness a direct no-double-stamping invariant.
# Match BMOPFTools' declared-variable arity from its public, resolved neutral
# labels. This can exceed the replacement's three physical phases when a source
# dataset uses an unrecognised neutral spelling; every unused native variable
# must still be fixed because native result/cost code can reference it.
function _native_ibr_declared_current_count(ctx, id::String)
    native = _opf_network(ctx)["ibr"][id]
    topology = uppercase(String(get(native, "topology", "FOUR_LEG")))
    terminals = String.(get(native, "terminal_map", String[]))
    # SINGLE_PHASE is included to keep this arity helper faithful to the public
    # BMOPFTools declaration contract when used independently. Fleet validation
    # currently admits only THREE_LEG and FOUR_LEG records, so the branch is not
    # reached by solve_controlled_inverter_fleet.
    topology == "SINGLE_PHASE" && return 1
    topology == "THREE_LEG" && return length(terminals)
    topology == "FOUR_LEG" || throw(ArgumentError(
        "unsupported native IBR topology '$topology' for '$id'"))
    neutral_labels = Set(String.(BMOPFTools.opf_neutral_labels(ctx)))
    # BMOPFTools' internal _phase_positions excludes its resolved neutral
    # position. Counting every non-neutral terminal is equivalent here because
    # _native_ibr_phase_count has already rejected maps with anything other than
    # exactly three ordered phases followed by one replacement neutral.
    count(terminal -> !(terminal in neutral_labels), terminals)
end

function _disable_native_ibr_port!(ctx, id::String)
    phase_count = _native_ibr_declared_current_count(ctx, id)
    for phase in 1:phase_count, component in (:real, :imag)
        current = _opf_ibr_current(ctx, id, phase; component=component)
        JuMP.fix(current, 0.0; force=true)
        JuMP.set_start_value(current, 0.0)
    end
    nothing
end

function _controlled_inverter_result(
        controlled::ControlledDevice{<:AdvancedInverter},
        request::InverterControlRequest,
        handles::_ControlledHandles,
        status::SolveStatus,
        outcome::SolveOutcome,
        bus::AbstractDict)
    extracted = extract_device(controlled, handles, status)
    grid_current = _extract_grid_current(handles.plant, status)
    exact_control = status.publishable ? evaluate_exact(
        controlled.controller,
        InverterControlMeasurement(collect(extracted.control.phase_voltage)),
        request,
        InverterControlRatings(
            controlled.device, controlled.controller.current_target)) :
        extracted.control
    if status.publishable
        exact_control = _apply_plant_capability_exact(
            exact_control, extracted.converter_terminal.phase_voltage,
            extracted.converter_terminal.phase_current, grid_current,
            controlled.device, controlled.controller.current_target)
    end
    residual = status.publishable ? maximum(abs,
        exact_control.phase_current .- extracted.control.phase_current) : NaN
    ControlledInverterResult(
        string(outcome.termination_status), extracted.plant, extracted.control,
        extracted.converter_terminal, grid_current, exact_control, residual,
        bus, status)
end

function _network_without_replaced_ibrs(result::AbstractDict,
                                        selected_ids)
    published = Dict{String,Any}(result)
    native = Dict{String,Any}(
        String(id) => value for (id, value) in get(result, "ibr", Dict()))
    for id in selected_ids
        delete!(native, id)
    end
    published["ibr"] = native
    published
end

"""
    solve_controlled_inverter_fleet(net, spec; kwargs...)
        -> ControlledInverterFleetResult

Solve one network snapshot containing any number of locally controlled advanced
inverters. Every selected native `network["ibr"][id]` is replaced through
BMOPFTools' typed per-component ownership seam. Unselected IBRs and all other
network devices keep their native formulations.

The controllers are fixed algebraic computation graphs based only on their
local POC voltage phasors. This is therefore a simultaneous controlled power
flow, not an optimization-based controller. `selection_objective=:loss`
minimizes residual converter allocation freedom; `:zero` is provided for
objective-invariance checks. `:network_cost` uses BMOPFTools' native generation
cost to select dispatchable *unselected* devices; it should only be used when
that economic objective is part of the intended snapshot. The smooth nonlinear
formulation currently requires `per_unit=true` and `pwm_strategy=:NONE`.

`result.devices[id]` is authoritative for replaced IBRs. To prevent accidental
comparison with unused native placeholder currents, selected identifiers are
removed from `result.network["ibr"]`; that dictionary contains only unselected
native IBR results. `result.build_manifest` records ownership explicitly.
"""
function solve_controlled_inverter_fleet(
        net::Dict{String,Any},
        spec::ControlledInverterFleetSpec;
        per_unit::Bool=true,
        s_base::Real=1e6,
        optimizer=Ipopt.Optimizer,
        verbose::Bool=false,
        selection_objective::Symbol=:loss,
        solver_options=())
    selection_objective in (:loss, :zero, :network_cost) || throw(ArgumentError(
        "selection_objective must be :loss, :zero, or :network_cost"))
    per_unit || throw(ArgumentError(
        "solve_controlled_inverter_fleet requires per_unit=true; raw-SI " *
        "scaling is not supported for the nonlinear controller formulation"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError(
        "s_base must be finite and > 0"))
    _validate_controlled_fleet(net, spec)

    handles = Dict{String,_ControlledHandles}()
    builder = BMOPFTools.OpfDeviceBuilder(:PowerOptLab, function (ctx, ids)
        for id in ids
            _disable_native_ibr_port!(ctx, id)
            handles[id] = stamp_device!(
                ctx, spec.devices[id]; request=spec.requests[id])
        end
        ctx
    end)
    build_spec = BMOPFTools.OpfBuildSpec(component_builders=Dict(
        (:ibr, id) => builder for id in keys(spec.devices)))
    hook! = ctx -> begin
        if selection_objective === :loss
            JuMP.@objective(_opf_model(ctx), Min,
                sum(h.plant.p_loss + h.plant.p_cap_loss
                    for h in values(handles)))
        elseif selection_objective === :network_cost
            JuMP.@objective(_opf_model(ctx), Min,
                BMOPFTools.generation_cost(ctx))
        else
            JuMP.@objective(_opf_model(ctx), Min, 0.0)
        end
    end
    ctx = build_opf_model(
        net; per_unit=true, s_base=Float64(s_base), add_objective=false,
        build_spec=build_spec, model_hook! = hook!, optimizer=optimizer,
        verbose=verbose)
    _set_solver_options!(_opf_model(ctx), solver_options)
    enforce_kcl!(ctx)
    JuMP.optimize!(_opf_model(ctx))

    outcome = _solve_outcome(_opf_model(ctx))
    status = SolveStatus(outcome)
    extracted_network = _extract_result(ctx, outcome)
    # Canonicalize once. Every device and the network result intentionally share
    # this read-only snapshot instead of allocating one top-level bus dictionary
    # per controlled inverter.
    bus = Dict{String,Any}(
        String(id) => value for (id, value) in extracted_network["bus"])
    devices = Dict{String,ControlledInverterResult}()
    for id in sort!(collect(keys(spec.devices)))
        devices[id] = _controlled_inverter_result(
            spec.devices[id], spec.requests[id], handles[id], status, outcome,
            bus)
    end
    network = _network_without_replaced_ibrs(
        extracted_network, keys(spec.devices))
    network["bus"] = bus
    ControlledInverterFleetResult(
        string(outcome.termination_status), devices, deepcopy(spec), network,
        BMOPFTools.opf_build_manifest(ctx), status)
end

"""
    controlled_inverter_rows(result)

Return deterministic, one-row-per-device scalar records suitable for
`DataFrame` construction or CSV/Arrow output. All electrical quantities are SI.
The ordinary VUF column is `|U₂|/|U₁|`; the controller's regularized ratio is
reported separately as `regularized_voltage_unbalance`.
"""
function controlled_inverter_rows(result::ControlledInverterFleetResult)
    rows = nothing
    for id in sort!(collect(keys(result.devices)))
        device = result.devices[id]
        plant = device.plant
        control = device.control
        controlled = result.spec.devices[id]
        inverter = inverter_spec(controlled)
        request = result.spec.requests[id]
        converter_total_current = maximum(hypot.(
            plant.i_mag, inverter.pwm_ac_converter_reserve))
        grid_total_current = maximum(hypot.(
            plant.i_grid_mag, inverter.pwm_ac_grid_reserve))
        u1 = abs(control.voltage_sequence[2])
        u2 = abs(control.voltage_sequence[3])
        vuf = isfinite(u1) && u1 > 0 ? u2/u1 : NaN
        stored_energy = inverter.c_dc === nothing || inverter.v_dc === nothing ?
            NaN : 0.5*inverter.c_dc*inverter.v_dc^2
        row = (
            device_id=id,
            bus=inverter.bus,
            positive_sequence_policy=
                string(nameof(typeof(controlled.controller.positive))),
            unbalance_policy=
                string(nameof(typeof(controlled.controller.unbalance))),
            current_target=
                string(nameof(typeof(controlled.controller.current_target))),
            termination_status=device.termination_status,
            publishable=device.solve.publishable,
            p_available_W=request.p_available,
            p_rated_W=request.p_rated,
            q_scale_var=request.q_scale,
            p_poc_W=plant.p_poc,
            q_poc_var=plant.q_poc,
            p_converter_W=plant.p_conv,
            q_converter_var=plant.q_conv,
            converter_loss_W=plant.p_loss + plant.p_cap_loss,
            voltage_min_V=control.voltage_min,
            voltage_max_V=control.voltage_max,
            positive_sequence_voltage_V=u1,
            negative_sequence_voltage_V=u2,
            voltage_unbalance_factor=vuf,
            regularized_voltage_unbalance=control.eta,
            requested_p_W=control.p_request,
            requested_q_var=control.q_request,
            command_curtailment_W=request.p_available-control.p_request,
            poc_active_power_shortfall_W=request.p_available-plant.p_poc,
            power_scale=control.power_scale,
            current_scale=control.current_scale,
            apparent_power_rating_VA=inverter.s_max,
            converter_current_rating_A=
                inverter.i_max === nothing ? NaN : inverter.i_max,
            grid_current_rating_A=
                inverter.i_grid_max === nothing ? NaN : inverter.i_grid_max,
            maximum_converter_current_A=maximum(plant.i_mag),
            maximum_grid_current_A=maximum(plant.i_grid_mag),
            maximum_converter_total_current_A=converter_total_current,
            maximum_grid_total_current_A=grid_total_current,
            converter_current_utilization=
                inverter.i_max === nothing ? NaN :
                converter_total_current/inverter.i_max,
            grid_current_utilization=
                inverter.i_grid_max === nothing ? NaN :
                grid_total_current/inverter.i_grid_max,
            positive_sequence_current_A=plant.i_positive,
            negative_sequence_current_A=plant.i_negative,
            ripple_power_VA=plant.ripple,
            dc_ripple_voltage_V=plant.dv2,
            dc_capacitance_F=
                inverter.c_dc === nothing ? NaN : inverter.c_dc,
            dc_nominal_voltage_V=
                inverter.v_dc === nothing ? NaN : inverter.v_dc,
            dc_stored_energy_J=stored_energy,
            capacitor_current_A=plant.i_cap,
            capacitor_thermal_current_A=plant.i_cap_thermal,
            capacitor_current_utilization=
                inverter.i_cap_max === nothing ? NaN :
                plant.i_cap_thermal/inverter.i_cap_max,
            dc_ripple_utilization=
                inverter.dv2_max === nothing ? NaN : plant.dv2/inverter.dv2_max,
            exact_smooth_current_residual_A=
                device.exact_smooth_current_residual,
        )
        if rows === nothing
            rows = [row]
        else
            push!(rows, row)
        end
    end
    rows === nothing ? NamedTuple[] : rows
end

"""
    controlled_inverter_phase_rows(result)

Return deterministic long-form records with one row for each controlled device
and phase. Complex voltage and current components, plus magnitudes, are retained
so phase-current ratings and current-location effects can be studied without
parsing nested solver results.
"""
function controlled_inverter_phase_rows(
        result::ControlledInverterFleetResult)
    rows = nothing
    for id in sort!(collect(keys(result.devices)))
        device = result.devices[id]
        for phase in 1:3
            voltage = device.control.phase_voltage[phase]
            converter_current = device.converter_terminal.phase_current[phase]
            grid_current = device.grid_phase_current[phase]
            row = (
                device_id=id,
                phase=phase,
                voltage_real_V=real(voltage),
                voltage_imag_V=imag(voltage),
                voltage_magnitude_V=abs(voltage),
                converter_current_real_A=real(converter_current),
                converter_current_imag_A=imag(converter_current),
                converter_current_magnitude_A=abs(converter_current),
                grid_current_real_A=real(grid_current),
                grid_current_imag_A=imag(grid_current),
                grid_current_magnitude_A=abs(grid_current),
            )
            if rows === nothing
                rows = [row]
            else
                push!(rows, row)
            end
        end
    end
    rows === nothing ? NamedTuple[] : rows
end
