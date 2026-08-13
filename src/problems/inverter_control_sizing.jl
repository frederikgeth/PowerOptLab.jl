# Counterfactual hardware grids and first-pass sizing diagnostics.

"""
    InverterHardwareSweepPoint(; id, converter_current_scale=1,
        grid_current_scale=nothing, dc_capacitance_scale=1,
        capacitor_current_scale=nothing)

One counterfactual hardware point for an outer inverter-control study. Scales
are applied to each device's installed rating, preserving heterogeneous fleet
nameplates. `nothing` leaves the optional grid-current or capacitor-current
rating absent; a numerical scale requires that the corresponding base rating
is explicitly present.

Apparent-power rating, DC voltage, ripple allowance, controller parameters, and
all other plant data remain fixed. This isolates the four named hardware axes.
"""
struct InverterHardwareSweepPoint
    id::String
    converter_current_scale::Float64
    grid_current_scale::Union{Nothing,Float64}
    dc_capacitance_scale::Float64
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
        converter_current_scale::Real=1.0,
        grid_current_scale::Union{Nothing,Real}=nothing,
        dc_capacitance_scale::Real=1.0,
        capacitor_current_scale::Union{Nothing,Real}=nothing)
    point_id = String(id)
    isempty(strip(point_id)) && throw(ArgumentError(
        "hardware sweep point id must be non-empty"))
    InverterHardwareSweepPoint(
        point_id,
        _positive_sizing_scale(
            converter_current_scale, "converter_current_scale"),
        _positive_sizing_scale(
            grid_current_scale, "grid_current_scale"; optional=true),
        _positive_sizing_scale(
            dc_capacitance_scale, "dc_capacitance_scale"),
        _positive_sizing_scale(
            capacitor_current_scale, "capacitor_current_scale";
            optional=true),
    )
end

function _scaled_required_rating(value, scale, field::String, id::String)
    scale === nothing && return value
    value === nothing && throw(ArgumentError(
        "inverter '$id' has no explicit $field to scale"))
    value * scale
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

function _same_struct_fields(left, right)
    typeof(left) === typeof(right) || return false
    all(name -> isequal(getfield(left, name), getfield(right, name)),
        fieldnames(typeof(left)))
end

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
            case.metadata == reference.metadata || throw(ArgumentError(
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

The original scenario identifier and hardware scales are added to metadata.
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
            "converter_current_scale" => point.converter_current_scale,
            "grid_current_scale" => point.grid_current_scale,
            "dc_capacitance_scale" => point.dc_capacitance_scale,
            "capacitor_current_scale" => point.capacitor_current_scale,
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

"""
    inverter_control_hardware_requirement_rows(study;
        allowed_dc_ripple_fraction, current_margin=1)

Return one first-pass hardware-requirement row per extracted inverter. For a
publishable point, converter/grid/capacitor current requirements are the
achieved RMS quantities multiplied by `current_margin`. Converter and grid
requirements combine the fundamental phase currents with any declared manual
carrier-current reserves in RMS quadrature. The monolithic three-leg DC-link
requirement is

```math
C_{2\\omega,req}=|\\widetilde S|/(2\\omega V_{dc}\\Delta V_{allow}).
```

`allowed_dc_ripple_fraction` defines
`ΔV_allow = allowed_dc_ripple_fraction*Vdc` and must lie in `(0, 1)`.
Non-publishable numerical requirements are NaN and compliance fields are
`missing`.

This diagnostic is a sizing result only when evaluated at a common hardware
point that does not alter dispatch or activate a limiter. Otherwise use it to
seed the explicit outer hardware grid. It assumes the modeled three-leg link
supplies the 2ω ripple; it is not a split-link, hold-up, lifetime, or EMT sizing
calculation.
"""
function inverter_control_hardware_requirement_rows(
        study::InverterControlStudyResult;
        allowed_dc_ripple_fraction::Real,
        current_margin::Real=1.0)
    fraction = Float64(allowed_dc_ripple_fraction)
    isfinite(fraction) && 0 < fraction < 1 || throw(ArgumentError(
        "allowed_dc_ripple_fraction must be finite and in (0, 1)"))
    margin = Float64(current_margin)
    isfinite(margin) && margin >= 1 || throw(ArgumentError(
        "current_margin must be finite and >= 1"))

    rows = NamedTuple[]
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
            allowed_ripple = fraction * vdc
            converter_requirement = publishable ?
                margin * maximum(hypot.(
                    plant.i_mag, inverter.pwm_ac_converter_reserve)) : NaN
            grid_requirement = publishable ?
                margin * maximum(hypot.(
                    plant.i_grid_mag, inverter.pwm_ac_grid_reserve)) : NaN
            capacitor_current_requirement = publishable ?
                margin * plant.i_cap_thermal : NaN
            capacitance_requirement = publishable ?
                plant.ripple /
                (2 * (2pi*inverter.f) * vdc * allowed_ripple) : NaN
            push!(rows, (
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
                allowed_dc_ripple_fraction=fraction,
                allowed_dc_ripple_voltage_V=allowed_ripple,
                converter_current_requirement_A=converter_requirement,
                grid_current_requirement_A=grid_requirement,
                capacitor_thermal_current_requirement_A=
                    capacitor_current_requirement,
                dc_capacitance_2omega_requirement_F=
                    capacitance_requirement,
                installed_converter_current_rating_A=
                    installed_converter_current,
                installed_grid_current_rating_A=installed_grid_current,
                installed_capacitor_current_rating_A=
                    installed_capacitor_current,
                installed_dc_capacitance_F=installed_capacitance,
                converter_current_utilization=_sizing_utilization(
                    publishable, installed_converter_current,
                    converter_requirement),
                grid_current_utilization=_sizing_utilization(
                    publishable, installed_grid_current, grid_requirement),
                capacitor_current_utilization=_sizing_utilization(
                    publishable, installed_capacitor_current,
                    capacitor_current_requirement),
                dc_capacitance_2omega_utilization=_sizing_utilization(
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
            ))
        end
    end
    rows
end
