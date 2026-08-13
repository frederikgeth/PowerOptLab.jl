# Counterfactual hardware grids and first-pass sizing diagnostics.

"""
    InverterHardwareSweepPoint(; id, converter_current_scale=nothing,
        grid_current_scale=nothing, dc_capacitance_scale=nothing,
        capacitor_current_scale=nothing)

One counterfactual hardware point for an outer inverter-control study. Scales
are applied to each device's installed rating, preserving heterogeneous fleet
nameplates. `nothing` leaves the corresponding rating unchanged; a numerical
scale requires that the base rating is explicitly present. Consequently a point
constructed with only `id` is an identity transformation, including for legal
plants whose optional ratings are absent.

Apparent-power rating, DC voltage, ripple allowance, controller parameters, and
all other plant data remain fixed. This isolates the four named hardware axes.
"""
struct InverterHardwareSweepPoint
    id::String
    converter_current_scale::Union{Nothing,Float64}
    grid_current_scale::Union{Nothing,Float64}
    dc_capacitance_scale::Union{Nothing,Float64}
    capacitor_current_scale::Union{Nothing,Float64}
end

function _positive_sizing_scale(value, name::String; optional::Bool=false)
    optional && value === nothing && return nothing
    scale = Float64(value)
    isfinite(scale) && scale > 0 || throw(ArgumentError(
        "$name must be finite and > 0"))
    scale
end

function InverterHardwareSweepPoint(;
        id,
        converter_current_scale::Union{Nothing,Real}=nothing,
        grid_current_scale::Union{Nothing,Real}=nothing,
        dc_capacitance_scale::Union{Nothing,Real}=nothing,
        capacitor_current_scale::Union{Nothing,Real}=nothing)
    point_id = String(id)
    isempty(strip(point_id)) && throw(ArgumentError(
        "hardware sweep point id must be non-empty"))
    InverterHardwareSweepPoint(
        point_id,
        _positive_sizing_scale(
            converter_current_scale, "converter_current_scale";
            optional=true),
        _positive_sizing_scale(
            grid_current_scale, "grid_current_scale"; optional=true),
        _positive_sizing_scale(
            dc_capacitance_scale, "dc_capacitance_scale"; optional=true),
        _positive_sizing_scale(
            capacitor_current_scale, "capacitor_current_scale";
            optional=true),
    )
end

function _scaled_required_rating(value, scale, field::String, id::String)
    scale === nothing && return value
    value === nothing && throw(ArgumentError(
        "inverter '$id' has no explicit $field to scale"))
    scaled = value * scale
    isfinite(scaled) && scaled > 0 || throw(ArgumentError(
        "inverter '$id' $field=$value scaled by $scale produces invalid " *
        "absolute rating $scaled; the result must be finite and > 0"))
    scaled
end

function _replace_advanced_inverter(
        inverter::AdvancedInverter,
        point::InverterHardwareSweepPoint)
    names = fieldnames(AdvancedInverter)
    values = Tuple(getfield(inverter, name) for name in names)
    parameters = NamedTuple{names}(values)
    replacements = (
        i_max=_scaled_required_rating(
            inverter.i_max, point.converter_current_scale, "i_max",
            inverter.id),
        i_grid_max=_scaled_required_rating(
            inverter.i_grid_max, point.grid_current_scale, "i_grid_max",
            inverter.id),
        c_dc=_scaled_required_rating(
            inverter.c_dc, point.dc_capacitance_scale, "c_dc", inverter.id),
        i_cap_max=_scaled_required_rating(
            inverter.i_cap_max, point.capacitor_current_scale, "i_cap_max",
            inverter.id),
    )
    AdvancedInverter(; merge(parameters, replacements)...)
end

"""
    resize_controlled_inverter_fleet(spec, point)

Return a new fleet with every advanced inverter scaled by `point`. Controllers
and operating-point requests are retained; the input fleet is not mutated.
Optional grid-side and capacitor-current scales require explicit base ratings,
preventing an absent physical constraint from being silently introduced.
"""
function resize_controlled_inverter_fleet(
        spec::ControlledInverterFleetSpec,
        point::InverterHardwareSweepPoint)
    devices = Dict{String,ControlledDevice}()
    for id in sort!(collect(keys(spec.devices)))
        controlled = spec.devices[id]
        inverter = _replace_advanced_inverter(inverter_spec(controlled), point)
        devices[id] = ControlledDevice(inverter, controlled.controller)
    end
    ControlledInverterFleetSpec(devices, spec.requests)
end

const _HARDWARE_METADATA_KEYS = (
    "base_scenario_id",
    "hardware_point_id",
    "converter_current_scale",
    "grid_current_scale",
    "dc_capacitance_scale",
    "capacitor_current_scale",
)

"""
    validate_inverter_control_campaign(cases; require_shared_network=true)

Validate the fixed-hardware matched-design contract used for control-law
comparisons. Every scenario must contain the same variant set. Within a
scenario, variants must have identical duration, weight, metadata, controlled
device identifiers, physical `AdvancedInverter` specifications, and requests.
Only the controller may differ.

By default, variants must also share the exact network object, an inexpensive
and strong guard against accidental snapshot drift. Set
`require_shared_network=false` only when equivalent networks were deliberately
materialized separately and validated by the dataset adapter.
"""
function validate_inverter_control_campaign(
        cases::AbstractVector;
        require_shared_network::Bool=true)
    _validate_study_cases(cases)
    grouped = Dict{String,Vector{InverterControlStudyCase}}()
    for case in cases
        push!(get!(grouped, case.scenario_id, InverterControlStudyCase[]), case)
    end
    scenarios = sort!(collect(keys(grouped)))
    expected_variants = Set(
        case.variant_id for case in grouped[first(scenarios)])
    for scenario in scenarios
        scenario_cases = grouped[scenario]
        variants = Set(case.variant_id for case in scenario_cases)
        variants == expected_variants || throw(ArgumentError(
            "scenario '$scenario' has variants $(sort!(collect(variants))); " *
            "expected $(sort!(collect(expected_variants)))"))
        reference = first(scenario_cases)
        reference_ids = Set(keys(reference.fleet.devices))
        for case in scenario_cases
            case.duration_h == reference.duration_h || throw(ArgumentError(
                "scenario '$scenario' has variant-dependent duration_h"))
            case.weight == reference.weight || throw(ArgumentError(
                "scenario '$scenario' has variant-dependent weight"))
            isequal(case.metadata, reference.metadata) || throw(ArgumentError(
                "scenario '$scenario' has variant-dependent metadata"))
            if require_shared_network && case.network !== reference.network
                throw(ArgumentError(
                    "scenario '$scenario' variants must share one network " *
                    "object; pass require_shared_network=false only after " *
                    "external network-equivalence validation"))
            end
            Set(keys(case.fleet.devices)) == reference_ids ||
                throw(ArgumentError(
                    "scenario '$scenario' has variant-dependent fleet ids"))
            for id in reference_ids
                _same_struct_fields(
                    inverter_spec(case.fleet.devices[id]),
                    inverter_spec(reference.fleet.devices[id])) ||
                    throw(ArgumentError(
                        "scenario '$scenario' device '$id' has " *
                        "variant-dependent physical hardware"))
                _same_struct_fields(
                    case.fleet.requests[id], reference.fleet.requests[id]) ||
                    throw(ArgumentError(
                        "scenario '$scenario' device '$id' has " *
                        "variant-dependent control request"))
            end
        end
    end
    nothing
end

"""
    expand_inverter_hardware_cases(cases, points)

Expand matched inverter-control cases over a counterfactual hardware grid.
Each hardware point is appended to `scenario_id`, while `variant_id` remains the
controller identifier. Consequently [`inverter_control_paired_rows`](@ref)
continues to compare controllers at identical hardware.

The original scenario identifier and effective hardware scales are added to
metadata (`1.0` for an unchanged axis).
Expanded cases retain the original network by reference. Hardware points are
counterfactual alternatives: aggregate summaries must group by
`"hardware_point_id"` and must never sum their energy as if they were successive
time intervals.
"""
function expand_inverter_hardware_cases(
        cases::AbstractVector,
        points::AbstractVector)
    validate_inverter_control_campaign(cases)
    isempty(points) && throw(ArgumentError(
        "a hardware sweep needs at least one point"))
    all(point -> point isa InverterHardwareSweepPoint, points) ||
        throw(ArgumentError(
            "every hardware-grid entry must be an InverterHardwareSweepPoint"))
    point_ids = [point.id for point in points]
    length(unique(point_ids)) == length(point_ids) || throw(ArgumentError(
        "hardware sweep point identifiers must be unique"))
    for case in cases, key in _HARDWARE_METADATA_KEYS
        haskey(case.metadata, key) && throw(ArgumentError(
            "study metadata key '$key' is reserved by the hardware sweep"))
    end

    ordered_cases = sort!(collect(cases);
        by=case -> (case.scenario_id, case.variant_id))
    ordered_points = sort!(collect(points); by=point -> point.id)
    expanded = InverterControlStudyCase[]
    for case in ordered_cases, point in ordered_points
        metadata = merge(case.metadata, Dict{String,Any}(
            "base_scenario_id" => case.scenario_id,
            "hardware_point_id" => point.id,
            "converter_current_scale" => something(
                point.converter_current_scale, 1.0),
            "grid_current_scale" => something(point.grid_current_scale, 1.0),
            "dc_capacitance_scale" => something(
                point.dc_capacitance_scale, 1.0),
            "capacitor_current_scale" => something(
                point.capacitor_current_scale, 1.0),
        ))
        push!(expanded, InverterControlStudyCase(
            scenario_id="$(case.scenario_id)::hardware[$(point.id)]",
            variant_id=case.variant_id,
            network=case.network,
            fleet=resize_controlled_inverter_fleet(case.fleet, point),
            duration_h=case.duration_h,
            weight=case.weight,
            metadata=metadata,
        ))
    end
    validate_inverter_control_campaign(expanded)
    expanded
end

function _sizing_assessment(publishable::Bool, installed, required)
    publishable && isfinite(installed) && isfinite(required) || return missing
    installed >= required
end

function _sizing_utilization(publishable::Bool, installed, required)
    publishable && isfinite(installed) && installed > 0 && isfinite(required) ||
        return NaN
    required / installed
end

function _sizing_binding(publishable::Bool, limit, achieved, tolerance::Float64)
    publishable && isfinite(limit) && limit > 0 && isfinite(achieved) ||
        return missing
    achieved >= (1 - tolerance) * limit
end

function _sizing_curtailment_active(
        publishable::Bool, available, requested, tolerance::Float64)
    publishable && isfinite(available) && isfinite(requested) || return missing
    available - requested > tolerance * max(abs(available), 1.0)
end

function _sizing_scale_active(publishable::Bool, scale, tolerance::Float64)
    publishable && isfinite(scale) || return missing
    scale < 1 - tolerance
end

"""
    inverter_control_hardware_requirement_rows(study;
        allowed_dc_ripple_fraction, current_margin=1,
        binding_tolerance=1e-4)

Return one first-pass hardware-requirement row per extracted inverter. For a
publishable point, converter/grid total RMS currents and the capacitor's
thermally weighted equivalent current are multiplied by `current_margin`.
Converter and grid requirements combine the fundamental phase currents with any
declared manual carrier-current reserves in RMS quadrature. The monolithic
three-leg DC-link quantities are

```math
I_{2\\omega,rms}=|\\widetilde S|/(\\sqrt{2}V_{dc}),\\qquad
C_{2\\omega,req}=|\\widetilde S|/(2\\omega V_{dc}\\Delta V_{allow}),\\qquad
E_{2\\omega,req}=\\tfrac12 C_{2\\omega,req}V_{dc}^2.
```

`allowed_dc_ripple_fraction` defines
`ΔV_allow = allowed_dc_ripple_fraction*Vdc` and must lie in `(0, 1)`. It is
the zero-to-peak sinusoidal amplitude, so `0.02` means ±2% or 4% peak-to-peak.
Non-publishable numerical requirements are NaN and compliance fields are
`missing`.

Binding diagnostics cover the installed converter/grid/capacitor/ripple
ratings and the converter apparent-power circle. Separate controller flags
identify command curtailment and activation of its power/current allocation
scales. A `false` value for one rating says only that particular rating is not
active; it never proves that the operating point is free of other limits.
Sampled modulation rails, internal-voltage bounds, sequence/neutral limits, and
other optional plant constraints are not exhaustively classified here.

This diagnostic is a sizing result only when evaluated at a common hardware
point that does not alter dispatch or activate a limiter. Otherwise use it to
seed the explicit outer hardware grid. It assumes the modeled three-leg link
supplies the 2ω ripple; it is not a split-link, hold-up, lifetime, or EMT sizing
calculation.
"""
function inverter_control_hardware_requirement_rows(
        study::InverterControlStudyResult;
        allowed_dc_ripple_fraction::Real,
        current_margin::Real=1.0,
        binding_tolerance::Real=1e-4)
    fraction = Float64(allowed_dc_ripple_fraction)
    isfinite(fraction) && 0 < fraction < 1 || throw(ArgumentError(
        "allowed_dc_ripple_fraction must be finite and in (0, 1)"))
    margin = Float64(current_margin)
    isfinite(margin) && margin >= 1 || throw(ArgumentError(
        "current_margin must be finite and >= 1"))
    tolerance = Float64(binding_tolerance)
    isfinite(tolerance) && 0 <= tolerance < 1 || throw(ArgumentError(
        "binding_tolerance must be finite and in [0, 1)"))

    rows = nothing
    for case_result in study.cases
        case_result.result === nothing && continue
        publishable = _case_publishable(case_result)
        for id in sort!(collect(keys(case_result.result.devices)))
            device = case_result.result.devices[id]
            inverter = inverter_spec(case_result.result.spec.devices[id])
            plant = device.plant
            vdc = something(inverter.v_dc, NaN)
            installed_capacitance = something(inverter.c_dc, NaN)
            installed_converter_current = something(inverter.i_max, NaN)
            installed_grid_current = something(inverter.i_grid_max, NaN)
            installed_capacitor_current = something(inverter.i_cap_max, NaN)
            enforced_ripple_limit = something(inverter.dv2_max, NaN)
            allowed_ripple = fraction * vdc
            achieved_converter_current = publishable ? maximum(hypot.(
                plant.i_mag, inverter.pwm_ac_converter_reserve)) : NaN
            achieved_grid_current = publishable ? maximum(hypot.(
                plant.i_grid_mag, inverter.pwm_ac_grid_reserve)) : NaN
            achieved_capacitor_current = publishable ? plant.i_cap_thermal : NaN
            achieved_ripple = publishable ? plant.dv2 : NaN
            achieved_apparent_power = publishable ?
                hypot(plant.p_conv, plant.q_conv) : NaN
            requested_power = publishable ?
                device.control.p_request : NaN
            converter_requirement = margin * achieved_converter_current
            grid_requirement = margin * achieved_grid_current
            capacitor_current_requirement = publishable ?
                margin * achieved_capacitor_current : NaN
            dc_2omega_current = publishable ?
                plant.ripple / (sqrt(2) * vdc) : NaN
            capacitance_requirement = publishable ?
                plant.ripple /
                (2 * (2pi*inverter.f) * vdc * allowed_ripple) : NaN
            installed_stored_energy =
                isfinite(installed_capacitance) && isfinite(vdc) ?
                0.5 * installed_capacitance * vdc^2 : NaN
            stored_energy_requirement = publishable ?
                0.5 * capacitance_requirement * vdc^2 : NaN
            row = (
                scenario_id=case_result.case.scenario_id,
                variant_id=case_result.case.variant_id,
                device_id=id,
                base_scenario_id=get(
                    case_result.case.metadata, "base_scenario_id",
                    case_result.case.scenario_id),
                hardware_point_id=get(
                    case_result.case.metadata, "hardware_point_id", nothing),
                duration_h=case_result.case.duration_h,
                weight=case_result.case.weight,
                publishable=publishable,
                current_margin=margin,
                binding_tolerance=tolerance,
                allowed_dc_ripple_fraction=fraction,
                allowed_dc_ripple_voltage_V=allowed_ripple,
                enforced_dc_ripple_limit_V=enforced_ripple_limit,
                achieved_dc_ripple_voltage_V=achieved_ripple,
                achieved_converter_apparent_power_VA=
                    achieved_apparent_power,
                installed_converter_apparent_power_rating_VA=inverter.s_max,
                converter_current_requirement_A=converter_requirement,
                grid_current_requirement_A=grid_requirement,
                capacitor_thermal_current_requirement_A=
                    capacitor_current_requirement,
                dc_2omega_capacitor_current_A=dc_2omega_current,
                dc_capacitance_2omega_requirement_F=
                    capacitance_requirement,
                dc_stored_energy_2omega_requirement_J=
                    stored_energy_requirement,
                installed_converter_current_rating_A=
                    installed_converter_current,
                installed_grid_current_rating_A=installed_grid_current,
                installed_capacitor_current_rating_A=
                    installed_capacitor_current,
                installed_dc_capacitance_F=installed_capacitance,
                installed_dc_stored_energy_J=installed_stored_energy,
                converter_current_binding=_sizing_binding(
                    publishable, installed_converter_current,
                    achieved_converter_current, tolerance),
                grid_current_binding=_sizing_binding(
                    publishable, installed_grid_current,
                    achieved_grid_current, tolerance),
                capacitor_current_binding=_sizing_binding(
                    publishable, installed_capacitor_current,
                    achieved_capacitor_current, tolerance),
                dc_ripple_binding=_sizing_binding(
                    publishable, enforced_ripple_limit, achieved_ripple,
                    tolerance),
                apparent_power_binding=_sizing_binding(
                    publishable, inverter.s_max, achieved_apparent_power,
                    tolerance),
                command_curtailment_active=_sizing_curtailment_active(
                    publishable,
                    case_result.result.spec.requests[id].p_available,
                    requested_power, tolerance),
                controller_power_limiter_active=_sizing_scale_active(
                    publishable, device.control.power_scale, tolerance),
                controller_current_limiter_active=_sizing_scale_active(
                    publishable, device.control.current_scale, tolerance),
                converter_current_requirement_utilization=_sizing_utilization(
                    publishable, installed_converter_current,
                    converter_requirement),
                grid_current_requirement_utilization=_sizing_utilization(
                    publishable, installed_grid_current, grid_requirement),
                capacitor_current_requirement_utilization=_sizing_utilization(
                    publishable, installed_capacitor_current,
                    capacitor_current_requirement),
                dc_capacitance_2omega_requirement_utilization=
                    _sizing_utilization(
                        publishable, installed_capacitance,
                        capacitance_requirement),
                converter_current_compliant=_sizing_assessment(
                    publishable, installed_converter_current,
                    converter_requirement),
                grid_current_compliant=_sizing_assessment(
                    publishable, installed_grid_current, grid_requirement),
                capacitor_current_compliant=_sizing_assessment(
                    publishable, installed_capacitor_current,
                    capacitor_current_requirement),
                dc_capacitance_2omega_compliant=_sizing_assessment(
                    publishable, installed_capacitance,
                    capacitance_requirement),
                metadata=case_result.case.metadata,
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
