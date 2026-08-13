# Deterministic outer orchestration for matched inverter-control experiments.

"""
    InverterControlStudyCase(; scenario_id, variant_id, network, fleet,
        duration_h=1, weight=1, metadata=Dict())

One independently solved point in a matched inverter-control experiment.
`scenario_id` identifies the common feeder/load/PV realization reused across
control or hardware `variant_id`s. `duration_h` converts powers to energy;
`weight` is a positive scenario multiplicity/probability weight. Arbitrary
dataset keys belong in `metadata` and are copied defensively.

The network is deliberately retained by reference: large matched studies may
reuse one prepared snapshot across variants. Treat it as immutable while a
study is running.
"""
struct InverterControlStudyCase
    scenario_id::String
    variant_id::String
    network::Dict{String,Any}
    fleet::ControlledInverterFleetSpec
    duration_h::Float64
    weight::Float64
    metadata::Dict{String,Any}
end

function InverterControlStudyCase(;
        scenario_id,
        variant_id,
        network::Dict{String,Any},
        fleet::ControlledInverterFleetSpec,
        duration_h::Real=1.0,
        weight::Real=1.0,
        metadata=Dict())
    sid = String(scenario_id)
    vid = String(variant_id)
    isempty(strip(sid)) && throw(ArgumentError("scenario_id must be non-empty"))
    isempty(strip(vid)) && throw(ArgumentError("variant_id must be non-empty"))
    duration = Float64(duration_h)
    scenario_weight = Float64(weight)
    isfinite(duration) && duration > 0 || throw(ArgumentError(
        "duration_h must be finite and > 0"))
    isfinite(scenario_weight) && scenario_weight > 0 || throw(ArgumentError(
        "weight must be finite and > 0"))
    metadata isa AbstractDict || throw(ArgumentError(
        "metadata must be a dictionary"))
    canonical_metadata = Dict{String,Any}()
    for (key, value) in metadata
        skey = String(key)
        haskey(canonical_metadata, skey) && throw(ArgumentError(
            "duplicate metadata key '$skey' after string conversion"))
        canonical_metadata[skey] = deepcopy(value)
    end
    InverterControlStudyCase(
        sid, vid, network, fleet, duration, scenario_weight,
        canonical_metadata)
end

"""Outcome of one [`InverterControlStudyCase`](@ref)."""
struct InverterControlStudyCaseResult
    case::InverterControlStudyCase
    result::Union{Nothing,ControlledInverterFleetResult}
    error_class::Union{Nothing,Symbol}
    error_type::Union{Nothing,String}
    error_message::Union{Nothing,String}
    elapsed_seconds::Float64
end

_case_publishable(case_result::InverterControlStudyCaseResult) =
    case_result.result !== nothing && solve_status(case_result.result).publishable

_case_termination(case_result::InverterControlStudyCaseResult) =
    case_result.result === nothing ? "ERROR" :
    solve_status(case_result.result).termination_status

"""Result of [`run_inverter_control_study`](@ref)."""
struct InverterControlStudyResult <: AbstractSolveResult
    cases::Vector{InverterControlStudyCaseResult}
    solve::SolveStatus
    settings::NamedTuple
end

solve_status(result::InverterControlStudyResult) = result.solve
solve_diagnostics(result::InverterControlStudyResult) = (
    case_count=length(result.cases),
    publishable_case_count=count(_case_publishable, result.cases),
    unpublished_case_count=count(
        case -> case.result !== nothing && !_case_publishable(case),
        result.cases),
    error_case_count=count(case -> case.result === nothing, result.cases),
    validation_error_case_count=count(
        case -> case.error_class === :validation, result.cases),
    unexpected_error_case_count=count(
        case -> case.error_class === :unexpected, result.cases),
    elapsed_seconds=sum(case.elapsed_seconds for case in result.cases),
)

function _study_status(cases::Vector{InverterControlStudyCaseResult})
    complete = all(case -> case.result !== nothing, cases)
    statuses = [solve_status(case.result) for case in cases
                if case.result !== nothing]
    all_publishable = complete && all(status -> status.publishable, statuses)
    all_feasible = complete && all(status -> status.feasible, statuses)
    all_optimal = complete && all(status -> status.optimal, statuses)
    all_primal = complete && all(status -> status.has_primal, statuses)
    termination = all_publishable ? "ALL_PUBLISHABLE" :
        "MIXED[" * join(sort!(unique!(_case_termination.(cases))), ",") * "]"
    primal = all_primal ? "ALL_HAVE_PRIMAL" : "INCOMPLETE_OR_MIXED_PRIMAL"
    SolveStatus(termination, primal, all_primal, all_feasible, all_optimal,
                all_publishable)
end

function _validate_study_cases(cases)
    isempty(cases) && throw(ArgumentError(
        "an inverter-control study needs at least one case"))
    all(case -> case isa InverterControlStudyCase, cases) || throw(ArgumentError(
        "every study entry must be an InverterControlStudyCase"))
    keys = [(case.scenario_id, case.variant_id) for case in cases]
    length(unique(keys)) == length(keys) || throw(ArgumentError(
        "each (scenario_id, variant_id) pair must be unique"))
    nothing
end

function _same_struct_fields(left, right)
    typeof(left) === typeof(right) || return false
    all(name -> isequal(getfield(left, name), getfield(right, name)),
        fieldnames(typeof(left)))
end

function _study_solver_options(solver_options)
    options = solver_options isa NamedTuple ? pairs(solver_options) :
        solver_options
    canonical = Pair{String,Any}[]
    for (name, value) in options
        push!(canonical, string(name) => deepcopy(value))
    end
    canonical
end

"""
    run_inverter_control_study(cases; kwargs...)
        -> InverterControlStudyResult

Solve independent matched cases in deterministic `(scenario_id, variant_id)`
order. Each case invokes [`solve_controlled_inverter_fleet`](@ref); snapshots
are not coupled into one NLP. `continue_on_error=true` (the default, equivalent
to `:validation`) records case-specific `ArgumentError`s from dataset and model
validation and continues. `continue_on_error=false` aborts on every exception.
The explicit unattended-campaign mode `continue_on_error=:all` also records
non-validation exceptions as `error_class=:unexpected`; fatal process-level
exceptions are always rethrown. Use `:all` deliberately: unexpected errors are
software/data failures, not physical infeasibility, and must be reported
separately.

This implementation is deliberately serial and reproducible. Large studies
should partition the case vector across processes outside this function and
concatenate the resulting rows. All fleet-solver keyword arguments are shared
across cases.
"""
function run_inverter_control_study(
        cases::AbstractVector;
        continue_on_error::Union{Bool,Symbol}=true,
        per_unit::Bool=true,
        s_base::Real=1e6,
        optimizer=Ipopt.Optimizer,
        verbose::Bool=false,
        selection_objective::Symbol=:loss,
        solver_options=())
    _validate_study_cases(cases)
    selection_objective in (:loss, :zero, :network_cost) || throw(ArgumentError(
        "selection_objective must be :loss, :zero, or :network_cost"))
    per_unit || throw(ArgumentError(
        "run_inverter_control_study requires per_unit=true"))
    isfinite(s_base) && s_base > 0 || throw(ArgumentError(
        "s_base must be finite and > 0"))
    error_mode = continue_on_error === true ? :validation :
        continue_on_error === false ? :none : continue_on_error
    error_mode in (:none, :validation, :all) || throw(ArgumentError(
        "continue_on_error must be false, true/:validation, or :all"))
    ordered = sort!(collect(cases);
        by=case -> (case.scenario_id, case.variant_id))
    canonical_options = _study_solver_options(solver_options)
    outcomes = InverterControlStudyCaseResult[]
    for case in ordered
        started = time_ns()
        try
            result = solve_controlled_inverter_fleet(
                case.network, case.fleet; per_unit=per_unit, s_base=s_base,
                optimizer=optimizer, verbose=verbose,
                selection_objective=selection_objective,
                solver_options=canonical_options)
            elapsed = (time_ns() - started) / 1e9
            push!(outcomes, InverterControlStudyCaseResult(
                case, result, nothing, nothing, nothing, elapsed))
        catch err
            err isa Union{InterruptException,OutOfMemoryError,StackOverflowError} &&
                rethrow()
            error_class = err isa ArgumentError ? :validation : :unexpected
            (error_mode === :all ||
             (error_mode === :validation && error_class === :validation)) ||
                rethrow()
            elapsed = (time_ns() - started) / 1e9
            push!(outcomes, InverterControlStudyCaseResult(
                case, nothing, error_class, string(typeof(err)),
                sprint(showerror, err), elapsed))
        end
    end
    settings = (
        per_unit=per_unit,
        s_base=Float64(s_base),
        optimizer=string(optimizer),
        verbose=verbose,
        selection_objective=selection_objective,
        continue_on_error=error_mode,
        solver_options=Dict{String,Any}(canonical_options),
    )
    InverterControlStudyResult(outcomes, _study_status(outcomes), settings)
end

"""Return one deterministic status/provenance row per study case."""
function inverter_control_study_case_rows(result::InverterControlStudyResult)
    [(
        scenario_id=case_result.case.scenario_id,
        variant_id=case_result.case.variant_id,
        duration_h=case_result.case.duration_h,
        weight=case_result.case.weight,
        termination_status=_case_termination(case_result),
        publishable=_case_publishable(case_result),
        error_class=case_result.error_class,
        error_type=case_result.error_type,
        error_message=case_result.error_message,
        elapsed_seconds=case_result.elapsed_seconds,
        controlled_device_count=length(case_result.case.fleet.devices),
        metadata=case_result.case.metadata,
    ) for case_result in result.cases]
end

"""
    inverter_control_study_device_rows(result)

Return one row per extracted controlled inverter, prefixed by case identity,
duration, weight, and metadata. Cases that throw before producing a fleet result
remain visible in [`inverter_control_study_case_rows`](@ref), but have no device
row. Non-publishable fleet solves retain device rows with masked numerical data.
"""
function inverter_control_study_device_rows(
        result::InverterControlStudyResult)
    [merge((
        scenario_id=case_result.case.scenario_id,
        variant_id=case_result.case.variant_id,
        duration_h=case_result.case.duration_h,
        weight=case_result.case.weight,
    ), device_row, (metadata=case_result.case.metadata,))
     for case_result in result.cases if case_result.result !== nothing
     for device_row in controlled_inverter_rows(case_result.result)]
end

"""Return one case-prefixed long-form row per controlled inverter phase."""
function inverter_control_study_phase_rows(
        result::InverterControlStudyResult)
    [merge((
        scenario_id=case_result.case.scenario_id,
        variant_id=case_result.case.variant_id,
        duration_h=case_result.case.duration_h,
        weight=case_result.case.weight,
    ), phase_row, (metadata=case_result.case.metadata,))
     for case_result in result.cases if case_result.result !== nothing
     for phase_row in controlled_inverter_phase_rows(case_result.result)]
end

function _weighted_quantile(values, weights, probability::Float64)
    isempty(values) && return NaN
    order = sortperm(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    threshold = probability * sum(sorted_weights)
    cumulative = 0.0
    for (value, weight) in zip(sorted_values, sorted_weights)
        cumulative += weight
        cumulative >= threshold && return Float64(value)
    end
    Float64(last(sorted_values))
end

function _study_group_value(case::InverterControlStudyCase, key::String)
    key == "scenario_id" && return case.scenario_id
    key == "variant_id" && return case.variant_id
    key == "duration_h" && return case.duration_h
    key == "weight" && return case.weight
    haskey(case.metadata, key) || throw(ArgumentError(
        "study case ('$(case.scenario_id)', '$(case.variant_id)') lacks " *
        "grouping metadata '$key'"))
    case.metadata[key]
end

_study_sort_value(value::Real) = (0, Float64(value))
_study_sort_value(value::AbstractString) = (1, String(value))
_study_sort_value(value::Symbol) = (2, String(value))
_study_sort_value(value::Nothing) = (3, "")
_study_sort_value(value::Missing) = (4, "")
_study_sort_value(value) = (5, repr(value))
_study_group_sort_key(group::Tuple) = Tuple(_study_sort_value.(group))

"""
    inverter_control_study_summary_rows(result; group_by=String[])

Return scalar cohort summaries as dictionaries. Every cohort is always grouped
by `variant_id`; `group_by` may add standard case fields or metadata keys such
as `"penetration"` and `"feeder"`. If every case carries
`"hardware_point_id"` metadata, that counterfactual axis is added
automatically; partial hardware tagging is rejected. Failure fractions use all
cases. Electrical and energy metrics use publishable device rows only.

Energy totals weight power by `duration_h*weight`. Current/VUF/capacitor tails
are duration-and-weight empirical inverse-CDF quantiles at 50, 95, and 99 %.
Every tail metric reports its own `*_finite_points` denominator; this is
essential for utilization metrics whose optional rating may be absent.
"""
function inverter_control_study_summary_rows(
        result::InverterControlStudyResult;
        group_by=String[])
    requested = group_by isa AbstractString ? [String(group_by)] :
        String.(collect(group_by))
    hardware_tagged = [haskey(case.case.metadata, "hardware_point_id")
                       for case in result.cases]
    if any(hardware_tagged)
        all(hardware_tagged) || throw(ArgumentError(
            "hardware_point_id metadata must be present on every case or none"))
        "hardware_point_id" in requested || push!(
            requested, "hardware_point_id")
    end
    group_keys = unique(["variant_id"; requested])
    groups = Dict{Tuple,Vector{InverterControlStudyCaseResult}}()
    for case_result in result.cases
        group = Tuple(
            _study_group_value(case_result.case, key) for key in group_keys)
        push!(get!(groups, group, InverterControlStudyCaseResult[]), case_result)
    end

    summaries = Dict{String,Any}[]
    for group in sort!(collect(Base.keys(groups)); by=_study_group_sort_key)
        case_results = groups[group]
        row = Dict{String,Any}(group_keys[index] => group[index]
                               for index in eachindex(group_keys))
        case_count = length(case_results)
        publishable = count(_case_publishable, case_results)
        errors = count(case -> case.result === nothing, case_results)
        unpublished = case_count - publishable - errors
        total_weight = sum(
            case.case.duration_h * case.case.weight for case in case_results)
        publishable_weight = sum((
            case.case.duration_h * case.case.weight for case in case_results
            if _case_publishable(case)); init=0.0)
        row["case_count"] = case_count
        row["publishable_case_count"] = publishable
        row["unpublished_case_count"] = unpublished
        row["error_case_count"] = errors
        row["validation_error_case_count"] = count(
            case -> case.error_class === :validation, case_results)
        row["unexpected_error_case_count"] = count(
            case -> case.error_class === :unexpected, case_results)
        row["publishable_fraction"] = publishable / case_count
        row["weighted_publishable_fraction"] = publishable_weight / total_weight

        device_rows = nothing
        exposures = Float64[]
        for case_result in case_results
            _case_publishable(case_result) || continue
            local_rows = controlled_inverter_rows(case_result.result)
            if device_rows === nothing
                device_rows = copy(local_rows)
            else
                append!(device_rows, local_rows)
            end
            append!(exposures, fill(
                case_result.case.duration_h * case_result.case.weight,
                length(local_rows)))
        end
        device_rows === nothing && (device_rows = NamedTuple[])
        row["publishable_device_points"] = length(device_rows)
        row["weighted_available_energy_kWh"] = sum((
            device_rows[index].p_available_W * exposures[index]
            for index in eachindex(device_rows)); init=0.0) / 1000
        row["weighted_command_curtailment_kWh"] = sum((
            device_rows[index].command_curtailment_W * exposures[index]
            for index in eachindex(device_rows)); init=0.0) / 1000
        row["weighted_poc_shortfall_kWh"] = sum((
            device_rows[index].poc_active_power_shortfall_W * exposures[index]
            for index in eachindex(device_rows)); init=0.0) / 1000
        row["weighted_converter_loss_kWh"] = sum((
            device_rows[index].converter_loss_W * exposures[index]
            for index in eachindex(device_rows)); init=0.0) / 1000

        metrics = (
            ("converter_current_A", :maximum_converter_current_A),
            ("grid_current_A", :maximum_grid_current_A),
            ("converter_total_current_A",
             :maximum_converter_total_current_A),
            ("grid_total_current_A", :maximum_grid_total_current_A),
            ("converter_current_utilization",
             :converter_current_utilization),
            ("grid_current_utilization", :grid_current_utilization),
            ("voltage_unbalance_factor", :voltage_unbalance_factor),
            ("capacitor_thermal_current_A", :capacitor_thermal_current_A),
            ("capacitor_current_utilization",
             :capacitor_current_utilization),
            ("dc_ripple_voltage_V", :dc_ripple_voltage_V),
            ("dc_ripple_utilization", :dc_ripple_utilization),
        )
        for (label, field) in metrics
            values = Float64[getproperty(device_row, field)
                             for device_row in device_rows]
            finite = findall(isfinite, values)
            finite_values = values[finite]
            finite_weights = exposures[finite]
            row["$(label)_finite_points"] = length(finite_values)
            for (suffix, probability) in (("p50", 0.50), ("p95", 0.95),
                                          ("p99", 0.99))
                row["$(label)_$suffix"] = _weighted_quantile(
                    finite_values, finite_weights, probability)
            end
            row["$(label)_max"] = isempty(finite_values) ? NaN :
                maximum(finite_values)
        end
        push!(summaries, row)
    end
    summaries
end

"""
    inverter_control_paired_rows(result, baseline_variant;
        require_shared_network=true)

Return matched device-level differences `variant - baseline` for every
scenario and non-baseline variant in the study. Missing, errored, or
non-publishable pairs are retained with `publishable_pair=false` and NaN
differences. A pair is matched only when exposure metadata, physical inverter
specifications, requests, device identifiers, and (by default) the network
object agree. This prevents both attrition and hardware/network confounding.
Set `require_shared_network=false` only after an external adapter has proved
separately materialized networks equivalent.
"""
function inverter_control_paired_rows(
        result::InverterControlStudyResult,
        baseline_variant::AbstractString;
        require_shared_network::Bool=true)
    baseline = String(baseline_variant)
    lookup = Dict(
        (case.case.scenario_id, case.case.variant_id) => case
        for case in result.cases)
    scenarios = sort!(unique(case.case.scenario_id for case in result.cases))
    variants = sort!(setdiff(
        unique(case.case.variant_id for case in result.cases), [baseline]))
    isempty(variants) && throw(ArgumentError(
        "paired comparison needs at least one non-baseline variant"))
    all(haskey(lookup, (scenario, baseline)) for scenario in scenarios) ||
        throw(ArgumentError(
            "baseline variant '$baseline' is missing from one or more scenarios"))

    rows = nothing
    for scenario in scenarios, variant in variants
        base_case = lookup[(scenario, baseline)]
        variant_case = get(lookup, (scenario, variant), nothing)
        matched_case_definition = variant_case !== nothing &&
            variant_case.case.duration_h == base_case.case.duration_h &&
            variant_case.case.weight == base_case.case.weight &&
            isequal(variant_case.case.metadata, base_case.case.metadata) &&
            (!require_shared_network ||
             variant_case.case.network === base_case.case.network) &&
            Set(keys(variant_case.case.fleet.devices)) ==
                Set(keys(base_case.case.fleet.devices)) &&
            all(id -> _same_struct_fields(
                    inverter_spec(variant_case.case.fleet.devices[id]),
                    inverter_spec(base_case.case.fleet.devices[id])) &&
                _same_struct_fields(
                    variant_case.case.fleet.requests[id],
                    base_case.case.fleet.requests[id]),
                keys(base_case.case.fleet.devices))
        base_rows = _case_publishable(base_case) ? Dict(
            row.device_id => row
            for row in controlled_inverter_rows(base_case.result)) :
            Dict{String,NamedTuple}()
        variant_rows = variant_case !== nothing &&
                       _case_publishable(variant_case) ? Dict(
            row.device_id => row
            for row in controlled_inverter_rows(variant_case.result)) :
            Dict{String,NamedTuple}()
        device_ids = sort!(collect(union(
            Set(keys(base_case.case.fleet.devices)),
            variant_case === nothing ? Set{String}() :
                Set(keys(variant_case.case.fleet.devices)))))
        for id in device_ids
            pair_publishable = _case_publishable(base_case) &&
                variant_case !== nothing && _case_publishable(variant_case) &&
                matched_case_definition && haskey(base_rows, id) &&
                haskey(variant_rows, id)
            if pair_publishable
                base_row = base_rows[id]
                variant_row = variant_rows[id]
                baseline_value = field -> getproperty(base_row, field)
                variant_value = field -> getproperty(variant_row, field)
                delta = field -> getproperty(variant_row, field) -
                                   getproperty(base_row, field)
            else
                baseline_value = _ -> NaN
                variant_value = _ -> NaN
                delta = _ -> NaN
            end
            row = (
                scenario_id=scenario,
                device_id=id,
                baseline_variant=baseline,
                variant_id=variant,
                duration_h=base_case.case.duration_h,
                weight=base_case.case.weight,
                metadata=base_case.case.metadata,
                matched_case_definition=matched_case_definition,
                publishable_pair=pair_publishable,
                baseline_command_curtailment_W=
                    baseline_value(:command_curtailment_W),
                variant_command_curtailment_W=
                    variant_value(:command_curtailment_W),
                delta_command_curtailment_W=delta(:command_curtailment_W),
                baseline_poc_active_power_shortfall_W=
                    baseline_value(:poc_active_power_shortfall_W),
                variant_poc_active_power_shortfall_W=
                    variant_value(:poc_active_power_shortfall_W),
                delta_poc_active_power_shortfall_W=
                    delta(:poc_active_power_shortfall_W),
                baseline_converter_loss_W=baseline_value(:converter_loss_W),
                variant_converter_loss_W=variant_value(:converter_loss_W),
                delta_converter_loss_W=delta(:converter_loss_W),
                baseline_converter_current_A=
                    baseline_value(:maximum_converter_current_A),
                variant_converter_current_A=
                    variant_value(:maximum_converter_current_A),
                delta_converter_current_A=
                    delta(:maximum_converter_current_A),
                baseline_grid_current_A=
                    baseline_value(:maximum_grid_current_A),
                variant_grid_current_A=
                    variant_value(:maximum_grid_current_A),
                delta_grid_current_A=delta(:maximum_grid_current_A),
                baseline_voltage_unbalance_factor=
                    baseline_value(:voltage_unbalance_factor),
                variant_voltage_unbalance_factor=
                    variant_value(:voltage_unbalance_factor),
                delta_voltage_unbalance_factor=
                    delta(:voltage_unbalance_factor),
                baseline_capacitor_thermal_current_A=
                    baseline_value(:capacitor_thermal_current_A),
                variant_capacitor_thermal_current_A=
                    variant_value(:capacitor_thermal_current_A),
                delta_capacitor_thermal_current_A=
                    delta(:capacitor_thermal_current_A),
                baseline_dc_ripple_voltage_V=
                    baseline_value(:dc_ripple_voltage_V),
                variant_dc_ripple_voltage_V=
                    variant_value(:dc_ripple_voltage_V),
                delta_dc_ripple_voltage_V=delta(:dc_ripple_voltage_V),
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

function _paired_group_value(row, key::String)
    key == "scenario_id" && return row.scenario_id
    key == "variant_id" && return row.variant_id
    key == "duration_h" && return row.duration_h
    key == "weight" && return row.weight
    haskey(row.metadata, key) || throw(ArgumentError(
        "paired row for scenario '$(row.scenario_id)' lacks " *
        "grouping metadata '$key'"))
    row.metadata[key]
end

"""
    inverter_control_paired_summary_rows(result, baseline_variant;
        group_by=String[], require_shared_network=true,
        delta_relative_tolerance=1e-6)

Summarize matched device-level `variant - baseline` rows without hiding
attrition. Each cohort reports all candidate, definition-matched, publishable,
and dropped device pairs. Delta means, p05/p50/p95 quantiles, and direction
fractions use `duration_h*weight` exposure over finite publishable pairs only.
A delta is indistinguishable when its magnitude is at most
`delta_relative_tolerance*max(abs(baseline), abs(variant), 1)` in the metric's
reported unit; the unit floor gives zero-valued comparisons a numerical band.
Negative, positive, and indistinguishable fractions are all reported so ties
and solver-scale noise are not conflated with one-sided effects. Per-metric
finite-pair counts expose structurally unavailable quantities.

`group_by` accepts the same standard case fields as
[`inverter_control_study_summary_rows`](@ref), or baseline metadata keys.
Unmatched pairs are assigned to the baseline case's cohort because their row
retains baseline exposure/provenance.
"""
function inverter_control_paired_summary_rows(
        result::InverterControlStudyResult,
        baseline_variant::AbstractString;
        group_by=String[],
        require_shared_network::Bool=true,
        delta_relative_tolerance::Real=1e-6)
    delta_rtol = Float64(delta_relative_tolerance)
    isfinite(delta_rtol) && delta_rtol >= 0 || throw(ArgumentError(
        "delta_relative_tolerance must be finite and >= 0"))
    requested = group_by isa AbstractString ? [String(group_by)] :
        String.(collect(group_by))
    paired = inverter_control_paired_rows(
        result, baseline_variant;
        require_shared_network=require_shared_network)
    hardware_tagged = [haskey(row.metadata, "hardware_point_id") for row in paired]
    if any(hardware_tagged)
        all(hardware_tagged) || throw(ArgumentError(
            "hardware_point_id metadata must be present on every pair or none"))
        "hardware_point_id" in requested || push!(
            requested, "hardware_point_id")
    end
    group_keys = unique(["variant_id"; requested])
    PairRow = eltype(paired)
    groups = Dict{Tuple,Vector{PairRow}}()
    for row in paired
        group = Tuple(_paired_group_value(row, key) for key in group_keys)
        push!(get!(groups, group, PairRow[]), row)
    end

    metrics = (
        ("command_curtailment_W", :delta_command_curtailment_W,
         :baseline_command_curtailment_W, :variant_command_curtailment_W),
        ("poc_active_power_shortfall_W",
         :delta_poc_active_power_shortfall_W,
         :baseline_poc_active_power_shortfall_W,
         :variant_poc_active_power_shortfall_W),
        ("converter_loss_W", :delta_converter_loss_W,
         :baseline_converter_loss_W, :variant_converter_loss_W),
        ("converter_current_A", :delta_converter_current_A,
         :baseline_converter_current_A, :variant_converter_current_A),
        ("grid_current_A", :delta_grid_current_A,
         :baseline_grid_current_A, :variant_grid_current_A),
        ("voltage_unbalance_factor", :delta_voltage_unbalance_factor,
         :baseline_voltage_unbalance_factor,
         :variant_voltage_unbalance_factor),
        ("capacitor_thermal_current_A",
         :delta_capacitor_thermal_current_A,
         :baseline_capacitor_thermal_current_A,
         :variant_capacitor_thermal_current_A),
        ("dc_ripple_voltage_V", :delta_dc_ripple_voltage_V,
         :baseline_dc_ripple_voltage_V, :variant_dc_ripple_voltage_V),
    )
    summaries = Dict{String,Any}[]
    for group in sort!(collect(keys(groups)); by=_study_group_sort_key)
        cohort = groups[group]
        row = Dict{String,Any}(group_keys[index] => group[index]
                               for index in eachindex(group_keys))
        publishable = filter(item -> item.publishable_pair, cohort)
        exposures = [item.duration_h * item.weight for item in publishable]
        row["candidate_pair_count"] = length(cohort)
        row["matched_definition_pair_count"] = count(
            item -> item.matched_case_definition, cohort)
        row["publishable_pair_count"] = length(publishable)
        row["dropped_pair_count"] = length(cohort) - length(publishable)
        row["publishable_pair_fraction"] = length(publishable) / length(cohort)
        for (label, delta_field, baseline_field, variant_field) in metrics
            values = Float64[getproperty(item, delta_field)
                             for item in publishable]
            baselines = Float64[getproperty(item, baseline_field)
                                for item in publishable]
            variants = Float64[getproperty(item, variant_field)
                               for item in publishable]
            finite = findall(index -> isfinite(values[index]) &&
                isfinite(baselines[index]) && isfinite(variants[index]),
                eachindex(values))
            finite_values = values[finite]
            finite_baselines = baselines[finite]
            finite_variants = variants[finite]
            finite_exposures = exposures[finite]
            finite_exposure = sum(finite_exposures; init=0.0)
            tolerances = delta_rtol .* max.(
                abs.(finite_baselines), abs.(finite_variants), 1.0)
            row["$(label)_finite_pair_count"] = length(finite_values)
            row["mean_delta_$label"] = isempty(finite_values) ? NaN :
                sum(finite_values .* finite_exposures) / finite_exposure
            for (suffix, probability) in (("p05", 0.05), ("p50", 0.50),
                                          ("p95", 0.95))
                row["delta_$(label)_$suffix"] = _weighted_quantile(
                    finite_values, finite_exposures, probability)
            end
            for (direction, predicate) in (
                    ("negative", (value, tolerance) -> value < -tolerance),
                    ("positive", (value, tolerance) -> value > tolerance),
                    ("indistinguishable",
                     (value, tolerance) -> abs(value) <= tolerance))
                row["$(direction)_delta_fraction_$label"] =
                    isempty(finite_values) ? NaN :
                    sum(finite_exposures[index]
                        for index in eachindex(finite_values)
                        if predicate(finite_values[index], tolerances[index]);
                        init=0.0) / finite_exposure
            end
        end
        push!(summaries, row)
    end
    summaries
end
