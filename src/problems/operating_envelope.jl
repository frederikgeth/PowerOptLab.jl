# Dynamic operating envelopes (DOEs).
#
# The implementation deliberately distinguishes a feasible allocation at the
# simultaneous bound (`security=:bound_point`) from an envelope whose corners
# have all been represented in the AC model (`security=:corners`).  The latter
# is still a local nonlinear-OPF result, not a proof that a non-convex feasible
# region has no interior holes; this limitation is carried into the result
# diagnostics instead of being hidden behind the word "robust".

"""
    ConnectionPoint(; id, bus, export_max=0.0, import_max=0.0, kwargs...)

An active connection point whose positive operating-envelope capacity is
calculated by [`solve_operating_envelope`](@ref).

- `export_max` / `import_max` are the connection's active-power nameplate limits
  (W, both non-negative). Select which one is used with `direction` on the solve.
- `ibr_id=nothing` retains the lightweight legacy port: aggregate unity-PF power
  is stamped directly at `bus` using `phase_terminals` and `neutral`.
- `ibr_id="pv1"` binds the envelope to an existing BMOPFTools `ibr`. This is the
  recommended representation for PV and batteries: the IBR keeps its prescribed
  Volt-VAr/Volt-Watt or fixed-power-factor control law, apparent/current limits,
  phase topology, and DC coupling while the DOE controls active power only.
- `requested` is an optional requested/forecast capacity (W), used by
  `FairnessPolicy(normalization=:request)`.
- `normalization` is an optional custom fairness reference (W), used by
  `FairnessPolicy(normalization=:custom)`.
"""
Base.@kwdef struct ConnectionPoint
    id::String
    bus::String
    phase_terminals::Vector{String} = ["1"]
    neutral::Union{String,Nothing} = "n"
    export_max::Float64 = 0.0
    import_max::Float64 = 0.0
    ibr_id::Union{String,Nothing} = nothing
    requested::Union{Float64,Nothing} = nothing
    normalization::Union{Float64,Nothing} = nothing
end

"""
    FairnessPolicy(; kind=:equal, normalization=:none, weights=Dict(),
                     alpha=1.0, epsilon=1e-6)

Parameterized policy for allocating active-power envelope capacity.

`normalization` defines the reference in `xᵢ = capacityᵢ/referenceᵢ`:

- `:none` — absolute watts (the reference is 1 W);
- `:capacity` — the connection's export/import nameplate for this direction;
- `:request` — `ConnectionPoint.requested`;
- `:custom` — `ConnectionPoint.normalization`.

Supported `kind` values:

- `:equal` — require equal normalized allocations and maximize their level;
- `:max_total` — maximize the weighted sum of normalized allocations;
- `:proportional` — weighted proportional fairness, `Σwᵢ log(xᵢ+ε)`;
- `:alpha` — weighted alpha fairness (`alpha=0` is weighted sum and `alpha=1`
  is proportional fairness);
- `:max_min` — maximize the minimum normalized allocation, then maximize total
  allocation while retaining that locally optimal minimum;
- `:equal_curtailment` — require equal normalized curtailment from nameplate and
  minimize it.

`weights` is keyed by connection-point id and defaults to one. All weights must
be finite and strictly positive. `epsilon` regularizes logarithmic/negative-power
utilities at zero; it does not create physical capacity.
"""
Base.@kwdef struct FairnessPolicy
    kind::Symbol = :equal
    normalization::Symbol = :none
    weights::Dict{String,Float64} = Dict{String,Float64}()
    alpha::Float64 = 1.0
    epsilon::Float64 = 1e-6
end

const _DOE_CONTROL_STAGES = (:issue, :scenario, :local_law, :context)
const _DOE_FREE_CONTROL_STAGES = (:issue, :scenario, :context)
const _DOE_UNCLASSIFIED_ACTIONS = (:error, :context)
const _DOE_NATIVE_CLASSIFICATIONS = (:free, :fixed_data, :local_law, :envelope)

"""
    DOEControlRegistration(; component, id, quantity, handle, kwargs...)

Register a control exposed by a BMOPFTools public `OpfModelKey` so that a DOE
policy can audit and link it without inspecting JuMP variable names.

`native_classification` is `:free`, `:fixed_data`, `:local_law`, or `:envelope`.
For a local law, provide `automatic_law` and immutable controller provenance in
`metadata`. `model_to_canonical_scale` converts the registered scalar model
value to `canonical_unit`; it must be identical for every context represented
by this registration. Native IBR powers and transformer taps are registered by
the solver automatically.
"""
struct DOEControlRegistration
    component::Symbol
    id::String
    quantity::Symbol
    position::Int
    handle::BMOPFTools.OpfModelKey
    native_classification::Symbol
    automatic_law::Union{Nothing,Symbol}
    canonical_unit::Symbol
    model_to_canonical_scale::Float64
    metadata::Dict{String,Any}

    function DOEControlRegistration(; component::Symbol,
                                    id::AbstractString,
                                    quantity::Symbol,
                                    position::Int=1,
                                    handle::BMOPFTools.OpfModelKey,
                                    native_classification::Symbol=:free,
                                    automatic_law::Union{Nothing,Symbol}=nothing,
                                    canonical_unit::Symbol=:dimensionless,
                                    model_to_canonical_scale::Real=1.0,
                                    metadata=Dict{String,Any}())
        isempty(id) && throw(ArgumentError(
            "DOE control-registration id must be non-empty"))
        position >= 1 || throw(ArgumentError(
            "DOE control-registration position must be positive"))
        handle.kind == :variable || throw(ArgumentError(
            "DOE control-registration handle must identify a variable"))
        native_classification in _DOE_NATIVE_CLASSIFICATIONS ||
            throw(ArgumentError(
                "native_classification must be one of $(_DOE_NATIVE_CLASSIFICATIONS)"))
        native_classification == :local_law && automatic_law === nothing &&
            throw(ArgumentError(
                "a :local_law registration must name automatic_law"))
        scale = Float64(model_to_canonical_scale)
        isfinite(scale) && !iszero(scale) || throw(ArgumentError(
            "model_to_canonical_scale must be finite and nonzero"))
        new(component, String(id), quantity, position, handle,
            native_classification, automatic_law, canonical_unit, scale,
            Dict{String,Any}(string(key) => value for (key, value) in pairs(metadata)))
    end
end

"""
    DOEControlRule(; component, id, quantity, stage)

Override the information stage of one DOE control family. Rules currently use
the stable device-level identity `(component, id, quantity)` and apply to every
phase or regulator position belonging to that device.

Supported native identities are transformer `:tap` and IBR `:active_power` /
`:reactive_power`. A stage is one of:

- `:issue` — one value is shared by every scenario and utilisation context;
- `:scenario` — one value is shared by the utilisation contexts within each
  scenario;
- `:local_law` — the value is determined by a prescribed automatic controller;
- `:context` — independent pointwise recourse (perfect information).

The solver rejects a rule that attempts to replace a prescribed local law or
the envelope's own active-power allocation with a different information stage.
"""
struct DOEControlRule
    component::Symbol
    id::String
    quantity::Symbol
    stage::Symbol

    function DOEControlRule(; component::Symbol, id::AbstractString,
                            quantity::Symbol, stage::Symbol)
        isempty(id) && throw(ArgumentError("DOE control-rule id must be non-empty"))
        stage in _DOE_CONTROL_STAGES || throw(ArgumentError(
            "DOE control-rule stage must be one of $(_DOE_CONTROL_STAGES), got :$stage"))
        new(component, String(id), quantity, stage)
    end
end

"""
    DOEControlPolicy(; name=:custom, default_stage=:context,
                       rules=DOEControlRule[], registrations=DOEControlRegistration[],
                       on_unclassified=:error)

Information structure for controllable assets in a DOE formulation.
`default_stage` applies to discovered free controls that have no explicit rule.
`on_unclassified=:error` fails closed when the model contains a free native
control that the current implementation cannot link safely;
`on_unclassified=:context` retains it as pointwise recourse and reports it.

Use [`PerfectRecourse`](@ref) to reproduce anticipative formulations and
[`IssuePlusLocalLaws`](@ref) for an operational policy in which free setpoints
are issued before uncertainty while prescribed automatic laws still respond
locally.
"""
struct DOEControlPolicy
    name::Symbol
    default_stage::Symbol
    rules::Vector{DOEControlRule}
    registrations::Vector{DOEControlRegistration}
    on_unclassified::Symbol

    function DOEControlPolicy(; name::Symbol=:custom,
                              default_stage::Symbol=:context,
                              rules=DOEControlRule[],
                              registrations=DOEControlRegistration[],
                              on_unclassified::Symbol=:error)
        default_stage in _DOE_FREE_CONTROL_STAGES || throw(ArgumentError(
            "DOE default control stage must be one of $(_DOE_FREE_CONTROL_STAGES), got :$default_stage"))
        on_unclassified in _DOE_UNCLASSIFIED_ACTIONS || throw(ArgumentError(
            "on_unclassified must be one of $(_DOE_UNCLASSIFIED_ACTIONS), got :$on_unclassified"))
        typed_rules = DOEControlRule[rule for rule in rules]
        identities = [(rule.component, rule.id, rule.quantity) for rule in typed_rules]
        length(unique(identities)) == length(identities) || throw(ArgumentError(
            "DOE control policy contains duplicate component/id/quantity rules"))
        typed_registrations = DOEControlRegistration[
            registration for registration in registrations]
        registration_identities = [(registration.component, registration.id,
                                    registration.quantity, registration.position)
                                   for registration in typed_registrations]
        length(unique(registration_identities)) == length(registration_identities) ||
            throw(ArgumentError(
                "DOE control policy contains duplicate control registrations"))
        new(name, default_stage, typed_rules, typed_registrations,
            on_unclassified)
    end
end

"""Policy preset retaining independent controls in every scenario/utilisation context."""
PerfectRecourse(; rules=DOEControlRule[], registrations=DOEControlRegistration[]) = DOEControlPolicy(
    name=:perfect_recourse, default_stage=:context, rules=rules,
    registrations=registrations,
    on_unclassified=:context)

"""Policy preset fixing all discovered free setpoints at DOE issue time."""
IssueFixedControls(; rules=DOEControlRule[], registrations=DOEControlRegistration[]) = DOEControlPolicy(
    name=:issue_fixed_controls, default_stage=:issue, rules=rules,
    registrations=registrations,
    on_unclassified=:error)

"""Operational preset: issue-time free setpoints plus prescribed local control laws."""
IssuePlusLocalLaws(; rules=DOEControlRule[], registrations=DOEControlRegistration[]) = DOEControlPolicy(
    name=:issue_plus_local_laws, default_stage=:issue, rules=rules,
    registrations=registrations,
    on_unclassified=:error)

const _DOE_SCENARIO_ROLES =
    (:train, :calibration, :validation, :test, :stress, :unspecified)

"""
    DOEScenario(; id, network, role=:unspecified, weight=nothing, kwargs...)

One network realization with explicit uncertainty provenance. `role` is one of
`:train`, `:calibration`, `:validation`, `:test`, `:stress`, or `:unspecified`.
`weight` is an optional positive relative probability weight; a scenario set
must either weight every scenario in an interval or none of them.

`source`, `generation_method`, `seed`, `timestamp`, and `metadata` describe how
the realization was obtained. They do not themselves establish calibration or
independence.
"""
struct DOEScenario
    id::String
    network::Dict{String,Any}
    role::Symbol
    weight::Union{Nothing,Float64}
    source::String
    generation_method::Symbol
    seed::Union{Nothing,Int}
    timestamp::Union{Nothing,DateTime}
    metadata::Dict{String,Any}

    function DOEScenario(; id::AbstractString,
                         network::Dict{String,Any},
                         role::Symbol=:unspecified,
                         weight::Union{Nothing,Real}=nothing,
                         source::AbstractString="unspecified",
                         generation_method::Symbol=:unspecified,
                         seed::Union{Nothing,Integer}=nothing,
                         timestamp::Union{Nothing,DateTime}=nothing,
                         metadata=Dict{String,Any}())
        isempty(id) && throw(ArgumentError("DOE scenario id must be non-empty"))
        role in _DOE_SCENARIO_ROLES || throw(ArgumentError(
            "DOE scenario role must be one of $(_DOE_SCENARIO_ROLES)"))
        weight_ = weight === nothing ? nothing : Float64(weight)
        weight_ === nothing || (isfinite(weight_) && weight_ > 0) ||
            throw(ArgumentError("DOE scenario weight must be finite and positive"))
        isempty(source) && throw(ArgumentError(
            "DOE scenario source must be non-empty"))
        new(String(id), network, role, weight_, String(source),
            generation_method, seed === nothing ? nothing : Int(seed),
            timestamp,
            Dict{String,Any}(string(key) => value
                             for (key, value) in pairs(metadata)))
    end
end

"""
    DOEScenarioSet(intervals; dataset_id, metadata=Dict())

Typed scenario ensemble. Pass a vector of [`DOEScenario`](@ref) for one interval
or a vector of non-empty scenario vectors for multiple intervals. Scenario IDs
must be unique within each interval. Optional weights must be complete within an
interval; they are treated as relative weights and normalized only when a
coverage statistic is calculated.
"""
struct DOEScenarioSet
    intervals::Vector{Vector{DOEScenario}}
    dataset_id::String
    metadata::Dict{String,Any}

    function DOEScenarioSet(intervals;
                            dataset_id::AbstractString,
                            metadata=Dict{String,Any}())
        isempty(dataset_id) && throw(ArgumentError(
            "DOE scenario-set dataset_id must be non-empty"))
        groups = if intervals isa AbstractVector &&
                    all(item -> item isa DOEScenario, intervals)
            [DOEScenario[item for item in intervals]]
        elseif intervals isa AbstractVector &&
               all(group -> group isa AbstractVector && !isempty(group) &&
                            all(item -> item isa DOEScenario, group), intervals)
            [DOEScenario[item for item in group] for group in intervals]
        else
            throw(ArgumentError(
                "DOEScenarioSet requires scenarios or non-empty scenario vectors"))
        end
        isempty(groups) && throw(ArgumentError(
            "DOE scenario set needs at least one interval"))
        all(!isempty, groups) || throw(ArgumentError(
            "every DOE scenario-set interval must be non-empty"))
        for (interval, group) in enumerate(groups)
            ids = [scenario.id for scenario in group]
            length(unique(ids)) == length(ids) || throw(ArgumentError(
                "DOE scenario IDs must be unique within interval $interval"))
            weighted = [scenario.weight !== nothing for scenario in group]
            all(weighted) || all(!, weighted) || throw(ArgumentError(
                "DOE scenario weights must be complete within interval $interval"))
        end
        new(groups, String(dataset_id),
            Dict{String,Any}(string(key) => value
                             for (key, value) in pairs(metadata)))
    end
end

"""Return a typed subset containing scenarios whose role is in `roles`."""
function select_doe_scenarios(scenarios::DOEScenarioSet;
                              roles=(:test,))
    role_values = roles isa Symbol ? (roles,) : roles
    selected_roles = Set(Symbol(role) for role in role_values)
    isempty(selected_roles) && throw(ArgumentError(
        "scenario-role selection must be non-empty"))
    issubset(selected_roles, Set(_DOE_SCENARIO_ROLES)) || throw(ArgumentError(
        "unknown DOE scenario role in $(collect(selected_roles))"))
    groups = Vector{Vector{DOEScenario}}()
    for (interval, group) in enumerate(scenarios.intervals)
        selected = [scenario for scenario in group
                    if scenario.role in selected_roles]
        isempty(selected) && throw(ArgumentError(
            "no selected DOE scenarios in interval $interval"))
        push!(groups, selected)
    end
    return DOEScenarioSet(groups;
        dataset_id=scenarios.dataset_id,
        metadata=merge(copy(scenarios.metadata),
            Dict{String,Any}("selected_roles" => sort!(collect(selected_roles)))))
end

"""
    DOEScenarioTimeSplit

Chronological calibration/test split returned by
[`split_doe_scenarios_by_time`](@ref). Scenarios between `calibration_end` and
`test_start` are excluded as a leakage-control gap and listed explicitly.
"""
struct DOEScenarioTimeSplit
    calibration::DOEScenarioSet
    test::DOEScenarioSet
    calibration_end::DateTime
    test_start::DateTime
    excluded_scenario_ids::Vector{String}
    diagnostics::Dict{String,Any}
end

function _doe_scenario_with_role(scenario::DOEScenario, role::Symbol,
                                 split_name::AbstractString)
    metadata = copy(scenario.metadata)
    metadata["original_role"] = scenario.role
    metadata["assigned_by_split"] = String(split_name)
    return DOEScenario(
        id=scenario.id,
        network=scenario.network,
        role=role,
        weight=scenario.weight,
        source=scenario.source,
        generation_method=scenario.generation_method,
        seed=scenario.seed,
        timestamp=scenario.timestamp,
        metadata=metadata)
end

"""
    split_doe_scenarios_by_time(scenarios;
        calibration_end, test_start=calibration_end, split_name="time_block")

Create a chronological calibration/test split for a one-interval typed scenario
ensemble. Timestamps before `calibration_end` are assigned `:calibration`;
timestamps on or after `test_start` are assigned `:test`; the intervening gap
is excluded. Every scenario must have a timestamp and each retained block must
be non-empty.

The one-interval restriction is deliberate in this first slice: it avoids
silently confusing forecast intervals with statistical sample time. Build one
split per forecast interval when their historical samples differ.
"""
function split_doe_scenarios_by_time(
        scenarios::DOEScenarioSet;
        calibration_end::DateTime,
        test_start::DateTime=calibration_end,
        split_name::AbstractString="time_block")
    length(scenarios.intervals) == 1 || throw(ArgumentError(
        "time-block DOE splitting currently requires one scenario interval"))
    test_start >= calibration_end || throw(ArgumentError(
        "test_start must not precede calibration_end"))
    isempty(split_name) && throw(ArgumentError("split_name must be non-empty"))
    group = only(scenarios.intervals)
    all(scenario -> scenario.timestamp !== nothing, group) ||
        throw(ArgumentError(
            "every scenario needs a timestamp for a chronological split"))
    calibration = DOEScenario[]
    test = DOEScenario[]
    excluded = String[]
    for scenario in group
        timestamp = something(scenario.timestamp)
        if timestamp < calibration_end
            push!(calibration,
                _doe_scenario_with_role(scenario, :calibration, split_name))
        elseif timestamp >= test_start
            push!(test, _doe_scenario_with_role(scenario, :test, split_name))
        else
            push!(excluded, scenario.id)
        end
    end
    isempty(calibration) && throw(ArgumentError(
        "chronological split produced no calibration scenarios"))
    isempty(test) && throw(ArgumentError(
        "chronological split produced no test scenarios"))
    common_metadata = merge(copy(scenarios.metadata), Dict{String,Any}(
        "parent_dataset_id" => scenarios.dataset_id,
        "split_name" => String(split_name),
        "calibration_end" => calibration_end,
        "test_start" => test_start,
        "excluded_scenario_ids" => copy(excluded)))
    calibration_set = DOEScenarioSet(calibration;
        dataset_id=scenarios.dataset_id * "-calibration",
        metadata=merge(copy(common_metadata),
            Dict{String,Any}("assigned_role" => :calibration)))
    test_set = DOEScenarioSet(test;
        dataset_id=scenarios.dataset_id * "-test",
        metadata=merge(copy(common_metadata),
            Dict{String,Any}("assigned_role" => :test)))
    diagnostics = Dict{String,Any}(
        "method" => :chronological_holdout,
        "calibration_count" => length(calibration),
        "test_count" => length(test),
        "gap_excluded_count" => length(excluded),
        "temporal_overlap" => false,
        "group_or_site_leakage_assessed" => false)
    return DOEScenarioTimeSplit(
        calibration_set, test_set, calibration_end, test_start,
        excluded, diagnostics)
end

"""
    DOEStudySpec

Immutable top-level provenance record for a reproducible DOE experiment. Use
the network-aware constructor `DOEStudySpec(nets, connection_points; kwargs...)`
and serialize [`doe_study_manifest`](@ref) beside research outputs.
"""
struct DOEStudySpec
    study_id::String
    network_hashes::Vector{Vector{String}}
    connection_points::Vector{Dict{String,Any}}
    direction::Symbol
    security::Symbol
    utilizations
    control_policy::Dict{String,Any}
    fairness::Dict{String,Any}
    solver::String
    solver_options::Dict{String,Any}
    seeds::Dict{String,Int}
    software_versions::Dict{String,String}
    metadata::Dict{String,Any}
    scenario_provenance::Vector{Vector{Dict{String,Any}}}
    scenario_set_metadata::Dict{String,Any}
end

# Source compatibility for the original positional provenance record.
DOEStudySpec(study_id, network_hashes, connection_points, direction, security,
             utilizations, control_policy, fairness, solver, solver_options,
             seeds, software_versions, metadata) =
    DOEStudySpec(study_id, network_hashes, connection_points, direction,
                 security, utilizations, control_policy, fairness, solver,
                 solver_options, seeds, software_versions, metadata,
                 Vector{Vector{Dict{String,Any}}}(), Dict{String,Any}())

function _doe_canonical(value)
    value === nothing && return "null"
    value isa Bool && return value ? "true" : "false"
    value isa Number && return repr(value)
    value isa AbstractString && return repr(String(value))
    value isa Symbol && return ":" * String(value)
    value isa DateTime && return repr(string(value))
    if value isa AbstractDict || value isa NamedTuple
        entries = [(string(key), item) for (key, item) in pairs(value)]
        sort!(entries; by=first)
        return "{" * join((repr(key) * ":" * _doe_canonical(item)
                            for (key, item) in entries), ",") * "}"
    elseif value isa Tuple || value isa AbstractVector
        return "[" * join((_doe_canonical(item) for item in value), ",") * "]"
    end
    return repr(value)
end

_doe_sha256(value) = bytes2hex(SHA.sha256(codeunits(_doe_canonical(value))))

function _doe_policy_manifest(policy::DOEControlPolicy)
    Dict{String,Any}(
        "name" => policy.name,
        "default_stage" => policy.default_stage,
        "on_unclassified" => policy.on_unclassified,
        "rules" => [Dict{String,Any}(
            "component" => rule.component, "id" => rule.id,
            "quantity" => rule.quantity, "stage" => rule.stage)
            for rule in policy.rules],
        "registrations" => [Dict{String,Any}(
            "component" => registration.component,
            "id" => registration.id,
            "quantity" => registration.quantity,
            "position" => registration.position,
            "handle" => Dict("kind" => registration.handle.kind,
                             "family" => registration.handle.family,
                             "index" => registration.handle.index),
            "native_classification" => registration.native_classification,
            "automatic_law" => registration.automatic_law,
            "canonical_unit" => registration.canonical_unit,
            "model_to_canonical_scale" =>
                registration.model_to_canonical_scale,
            "metadata" => copy(registration.metadata))
            for registration in policy.registrations])
end

function _doe_fairness_manifest(policy::FairnessPolicy)
    Dict{String,Any}(
        "kind" => policy.kind,
        "normalization" => policy.normalization,
        "weights" => copy(policy.weights),
        "alpha" => policy.alpha,
        "epsilon" => policy.epsilon)
end

"""
    DOEStudySpec(nets, connection_points; kwargs...)

Create a deterministic study manifest with SHA-256 hashes of every interval /
scenario network and the full connection-point declaration. The manifest also
records typed scenario IDs, roles, weights and construction provenance when
`nets` is a [`DOEScenarioSet`](@ref), plus coverage, control and fairness
policies, solver options, named random seeds, extension metadata, and package
versions. Functions such as custom hooks cannot be serialized; identify their
committed implementation in `metadata`.
"""
function DOEStudySpec(nets, connection_points::AbstractVector{ConnectionPoint};
                      direction::Symbol=:export,
                      security::Symbol=:bound_point,
                      utilizations=nothing,
                      control_policy::DOEControlPolicy=PerfectRecourse(),
                      fairness=:equal,
                      optimizer=Ipopt.Optimizer,
                      solver_options=NamedTuple(),
                      seeds=Dict{String,Int}(),
                      metadata=Dict{String,Any}())
    groups = _scenario_groups(nets)
    scenario_provenance = _scenario_provenance_groups(nets, groups)
    cps = collect(connection_points)
    fairness_policy = _as_policy(fairness)
    network_hashes = [[_doe_sha256(net) for net in group] for group in groups]
    scenario_set_metadata = nets isa DOEScenarioSet ? Dict{String,Any}(
        "dataset_id" => nets.dataset_id,
        "metadata" => copy(nets.metadata),
        "typed" => true) : Dict{String,Any}(
        "dataset_id" => nothing,
        "metadata" => Dict{String,Any}(),
        "typed" => false)
    cp_records = [Dict{String,Any}(
        "id" => cp.id, "bus" => cp.bus,
        "phase_terminals" => copy(cp.phase_terminals),
        "neutral" => cp.neutral, "export_max" => cp.export_max,
        "import_max" => cp.import_max, "ibr_id" => cp.ibr_id,
        "requested" => cp.requested, "normalization" => cp.normalization)
        for cp in cps]
    options = Dict{String,Any}(string(key) => value
                              for (key, value) in pairs(solver_options))
    seed_records = Dict{String,Int}(string(key) => Int(value)
                                    for (key, value) in pairs(seeds))
    metadata_records = Dict{String,Any}(string(key) => value
                                        for (key, value) in pairs(metadata))
    versions = Dict{String,String}(
        "julia" => string(VERSION),
        "PowerOptLab" => string(Base.pkgversion(@__MODULE__)),
        "BMOPFTools" => string(Base.pkgversion(BMOPFTools)),
        "JuMP" => string(Base.pkgversion(JuMP)),
        "Ipopt" => string(Base.pkgversion(Ipopt)))
    manifest = Dict{String,Any}(
        "network_hashes" => network_hashes,
        "scenario_provenance" => scenario_provenance,
        "scenario_set_metadata" => scenario_set_metadata,
        "connection_points" => cp_records,
        "direction" => direction,
        "security" => security,
        "utilizations" => utilizations,
        "control_policy" => _doe_policy_manifest(control_policy),
        "fairness" => _doe_fairness_manifest(fairness_policy),
        "solver" => string(optimizer),
        "solver_options" => options,
        "seeds" => seed_records,
        "software_versions" => versions,
        "metadata" => metadata_records)
    study_id = _doe_sha256(manifest)
    return DOEStudySpec(study_id, network_hashes, cp_records,
        direction, security, deepcopy(utilizations),
        manifest["control_policy"], manifest["fairness"],
        manifest["solver"], options, seed_records, versions, metadata_records,
        scenario_provenance, scenario_set_metadata)
end

"""Return a serialization-friendly dictionary for a [`DOEStudySpec`](@ref)."""
function doe_study_manifest(spec::DOEStudySpec)
    Dict{String,Any}(
        "study_id" => spec.study_id,
        "network_hashes" => deepcopy(spec.network_hashes),
        "scenario_provenance" => deepcopy(spec.scenario_provenance),
        "scenario_set_metadata" => deepcopy(spec.scenario_set_metadata),
        "connection_points" => deepcopy(spec.connection_points),
        "direction" => spec.direction,
        "security" => spec.security,
        "utilizations" => deepcopy(spec.utilizations),
        "control_policy" => deepcopy(spec.control_policy),
        "fairness" => deepcopy(spec.fairness),
        "solver" => spec.solver,
        "solver_options" => deepcopy(spec.solver_options),
        "seeds" => copy(spec.seeds),
        "software_versions" => copy(spec.software_versions),
        "metadata" => deepcopy(spec.metadata))
end

"""
    OperatingEnvelopeResult

Result of [`solve_operating_envelope`](@ref). Capacities are positive SI watts
for the selected `direction`.

`snapshots[t]` is the first scenario at the all-active upper corner. When an
interval has no feasible primal solution it contains only status metadata and
all capacities for that interval are `NaN`; infeasible solver iterates are never
published as envelopes.

`total_export` is retained as a backward-compatible alias of `total_capacity`.
For a new import study use `total_capacity` and inspect `direction`.

`fairness_metrics` reports allocation and curtailment metrics for the published
capacity at each interval. `schedule` records the issue/validity metadata and
whether a non-optimised fallback was published.
"""
struct OperatingEnvelopeResult <: AbstractSolveResult
    termination_status::Vector{String}
    envelope::Dict{String,Vector{Float64}}
    total_export::Vector{Float64}
    snapshots::Vector{Dict{String,Any}}
    direction::Symbol
    total_capacity::Vector{Float64}
    diagnostics::Vector{Dict{String,Any}}
    fairness_metrics::Vector{Dict{String,Any}}
    schedule::Vector{Dict{String,Any}}
end

function solve_status(result::OperatingEnvelopeResult)
    published = !isempty(result.total_capacity) && all(isfinite, result.total_capacity)
    has_primal = !isempty(result.diagnostics) && all(diag ->
        get(diag, "primal_status", "NO_SOLUTION") != "NO_SOLUTION",
        result.diagnostics)
    feasible = !isempty(result.diagnostics) ?
        all(diag -> get(diag, "feasible", false), result.diagnostics) : published
    optimal = feasible && all(status -> status in ("OPTIMAL", "LOCALLY_SOLVED"),
                              result.termination_status)
    primal = feasible ? "FEASIBLE_POINT" : published ? "FALLBACK_POINT" : "NO_SOLUTION"
    SolveStatus(_termination_summary(result.termination_status), primal,
                has_primal, feasible, optimal, published)
end

solve_diagnostics(result::OperatingEnvelopeResult) =
    (interval_count=length(result.total_capacity),
     feasible_count=count(diag -> get(diag, "feasible", false), result.diagnostics),
     published_count=count(isfinite, result.total_capacity),
     direction=result.direction)

# Source compatibility for callers that constructed the original four-field
# result directly. New solves always populate the richer fields above.
OperatingEnvelopeResult(status, envelope, total, snapshots) =
    OperatingEnvelopeResult(status, envelope, copy(total), snapshots, :export,
                            total, Dict{String,Any}[], Dict{String,Any}[],
                            Dict{String,Any}[])

"""
    OperatingEnvelopeContextResult

Evidence for one forecast scenario and participant-utilisation point in an
operating-envelope solve or verification. `feasible=nothing` means that a
multi-context model failed before feasibility of this individual context could
be determined. The nested diagnostics distinguish joint-model evidence from an
independent fixed-control replay.
"""
struct OperatingEnvelopeContextResult
    scenario::Int
    utilization_index::Int
    utilization::Vector{Float64}
    termination_status::String
    primal_status::String
    feasible::Union{Bool,Nothing}
    snapshot::Dict{String,Any}
    diagnostics::Dict{String,Any}
    control_values::Vector{Dict{String,Any}}
end

"""
    OperatingEnvelopeVerification

Result of [`verify_operating_envelope`](@ref). Each interval reports whether a
fixed, already-issued allocation was locally feasible at every requested
utilisation point and forecast/model scenario. This is a verification result,
not a new allocation.
"""
struct OperatingEnvelopeVerification <: AbstractSolveResult
    termination_status::Vector{String}
    feasible::Vector{Bool}
    snapshots::Vector{Dict{String,Any}}
    diagnostics::Vector{Dict{String,Any}}
    context_results::Vector{Vector{OperatingEnvelopeContextResult}}
end

# Source compatibility for callers constructing the original four-field result.
OperatingEnvelopeVerification(status, feasible, snapshots, diagnostics) =
    OperatingEnvelopeVerification(status, feasible, snapshots, diagnostics,
        [OperatingEnvelopeContextResult[] for _ in eachindex(status)])


function solve_status(result::OperatingEnvelopeVerification)
    has_primal = !isempty(result.diagnostics) && all(diag ->
        get(diag, "primal_status", "NO_SOLUTION") != "NO_SOLUTION",
        result.diagnostics)
    feasible = !isempty(result.feasible) && all(result.feasible)
    optimal = feasible && all(status -> status in ("OPTIMAL", "LOCALLY_SOLVED"),
                              result.termination_status)
    SolveStatus(_termination_summary(result.termination_status),
                feasible ? "FEASIBLE_POINT" : "NO_SOLUTION",
                has_primal, feasible, optimal, feasible)
end

solve_diagnostics(result::OperatingEnvelopeVerification) =
    (interval_count=length(result.feasible),
     feasible_count=count(identity, result.feasible), verification=true)

"""
    OperatingEnvelopeSearchResult

Result of [`search_operating_envelope_utilizations`](@ref). `outcome` is
`:search_stable`, `:candidate_counterexample`, or `:inconclusive`. Search-stable
means only that every generated point passed the configured local AC
verification; `global_certificate` remains false.
"""
struct OperatingEnvelopeSearchResult <: AbstractSolveResult
    outcome::Symbol
    utilization_points::Vector{Vector{Float64}}
    verification::OperatingEnvelopeVerification
    candidate_contexts::Vector{OperatingEnvelopeContextResult}
    diagnostics::Dict{String,Any}
end

solve_status(result::OperatingEnvelopeSearchResult) =
    solve_status(result.verification)
solve_diagnostics(result::OperatingEnvelopeSearchResult) =
    (outcome=result.outcome,
     tested_point_count=length(result.utilization_points),
     global_certificate=false)

"""
    DOEAdversarialSearchResult

Result of [`search_operating_envelope_adversarial`](@ref). The search starts
from a deterministic coverage set and performs coordinate refinement around
the utilization points with the smallest normalized network-constraint
headroom. It retains every joint verification round and the final score for
each tested point.

This is a black-box falsification heuristic. `outcome=:search_stable` means no
counterexample was found within the recorded budget; it is not a continuous-set
or global-optimality certificate.
"""
struct DOEAdversarialSearchResult <: AbstractSolveResult
    outcome::Symbol
    utilization_points::Vector{Vector{Float64}}
    point_scores::Vector{Float64}
    verifications::Vector{OperatingEnvelopeVerification}
    candidate_contexts::Vector{OperatingEnvelopeContextResult}
    worst_interval::Union{Nothing,Int}
    worst_context::Union{Nothing,OperatingEnvelopeContextResult}
    diagnostics::Dict{String,Any}
end

solve_status(result::DOEAdversarialSearchResult) =
    isempty(result.verifications) ?
        SolveStatus("NOT_RUN", "NO_SOLUTION", false, false, false, false) :
        solve_status(last(result.verifications))
solve_diagnostics(result::DOEAdversarialSearchResult) =
    (outcome=result.outcome,
     tested_point_count=length(result.utilization_points),
     rounds=length(result.verifications),
     global_certificate=false)

"""
    DOECounterexampleConfirmationResult

Multistart replay evidence returned by
[`confirm_operating_envelope_counterexample`](@ref). `outcome` is
`:repeated_candidate`, `:not_reproduced`, or `:inconclusive`. Repeated local
infeasibility is stronger numerical evidence than one failed solve, but remains
distinct from a globally certified physical violation.
"""
struct DOECounterexampleConfirmationResult <: AbstractSolveResult
    outcome::Symbol
    utilization::Vector{Float64}
    verifications::Vector{OperatingEnvelopeVerification}
    start_scales::Vector{Float64}
    diagnostics::Dict{String,Any}
end

function solve_status(result::DOECounterexampleConfirmationResult)
    feasible_index = findfirst(
        verification -> all(verification.feasible), result.verifications)
    return solve_status(feasible_index === nothing ?
        last(result.verifications) : result.verifications[feasible_index])
end
solve_diagnostics(result::DOECounterexampleConfirmationResult) =
    (outcome=result.outcome,
     run_count=length(result.verifications),
     global_certificate=false)

"""
    AdversarialSearchStableOperatingEnvelopeResult

Evidence from [`solve_adversarial_search_stable_operating_envelope`](@ref),
including every allocation and adaptive falsification search. `envelope` is the
last allocation and `utilization_points` is the final accumulated finite set.
"""
struct AdversarialSearchStableOperatingEnvelopeResult <: AbstractSolveResult
    outcome::Symbol
    envelope::OperatingEnvelopeResult
    allocations::Vector{OperatingEnvelopeResult}
    searches::Vector{DOEAdversarialSearchResult}
    utilization_points::Vector{Vector{Float64}}
    rounds::Int
    diagnostics::Dict{String,Any}
end

solve_status(result::AdversarialSearchStableOperatingEnvelopeResult) =
    solve_status(result.envelope)
solve_diagnostics(result::AdversarialSearchStableOperatingEnvelopeResult) =
    (outcome=result.outcome,
     rounds=result.rounds,
     tested_point_count=length(result.utilization_points),
     global_certificate=false)

"""
    DOECoverageResult

Held-out scenario evaluation returned by
[`evaluate_operating_envelope_coverage`](@ref). `scenario_rows` contains one
safe/unsafe/unresolved record per selected interval/scenario, while `metrics`
reports empirical context and scenario rates. A confidence bound is present
only when the caller explicitly declares the selected scenarios i.i.d.
"""
struct DOECoverageResult <: AbstractSolveResult
    outcome::Symbol
    selected_roles::Vector{Symbol}
    verification::OperatingEnvelopeVerification
    scenario_rows::Vector{NamedTuple}
    metrics::Dict{String,Any}
    diagnostics::Dict{String,Any}
end

solve_status(result::DOECoverageResult) = solve_status(result.verification)
solve_diagnostics(result::DOECoverageResult) =
    (outcome=result.outcome,
     scenario_count=length(result.scenario_rows),
     global_certificate=false)

"""
    DOECoverageCurveResult

Capacity-scaling sensitivity returned by
[`evaluate_operating_envelope_coverage_curve`](@ref). Each element of
`coverages` is a full held-out evaluation at the corresponding scale. The
curve is finite numerical evidence, not a continuous security certificate.
"""
struct DOECoverageCurveResult <: AbstractSolveResult
    scales::Vector{Float64}
    coverages::Vector{DOECoverageResult}
    rows::Vector{NamedTuple}
    diagnostics::Dict{String,Any}
end

solve_status(result::DOECoverageCurveResult) =
    solve_status(last(result.coverages))
solve_diagnostics(result::DOECoverageCurveResult) =
    (scale_count=length(result.scales),
     critical_scale=get(result.diagnostics, "first_candidate_scale", missing),
     global_certificate=false)

"""
    DOECoverageShiftResult

Descriptive comparison returned by [`compare_doe_coverage_shift`](@ref).
`metric_deltas` are shifted minus reference. This record can reveal a change in
observed performance, but does not constitute a statistical distribution-shift
test.
"""
struct DOECoverageShiftResult <: AbstractSolveResult
    outcome::Symbol
    reference::DOECoverageResult
    shifted::DOECoverageResult
    metric_deltas::Dict{String,Any}
    diagnostics::Dict{String,Any}
end

solve_status(result::DOECoverageShiftResult) = solve_status(result.shifted)
solve_diagnostics(result::DOECoverageShiftResult) =
    (outcome=result.outcome,
     distribution_shift_detected=false,
     global_certificate=false)

"""
    SearchStableOperatingEnvelopeResult

Result of [`solve_search_stable_operating_envelope`](@ref). It retains the
final allocation, every utilization-screening round, and the accumulated test
set. `outcome=:search_stable` is a finite-search result rather than a global
robust certificate.
"""
struct SearchStableOperatingEnvelopeResult <: AbstractSolveResult
    outcome::Symbol
    envelope::OperatingEnvelopeResult
    searches::Vector{OperatingEnvelopeSearchResult}
    utilization_points::Vector{Vector{Float64}}
    rounds::Int
    diagnostics::Dict{String,Any}
end

solve_status(result::SearchStableOperatingEnvelopeResult) =
    solve_status(result.envelope)
solve_diagnostics(result::SearchStableOperatingEnvelopeResult) =
    (outcome=result.outcome, rounds=result.rounds,
     tested_point_count=length(result.utilization_points),
     global_certificate=false)

"""
    OperatingEnvelopeMultistartResult

Collection returned by [`solve_operating_envelope_multistart`](@ref). `selected`
is the best publishable run under the declared selection rule; all starts and
capacity spreads remain available for branch-sensitivity reporting.
"""
struct OperatingEnvelopeMultistartResult <: AbstractSolveResult
    selected::OperatingEnvelopeResult
    runs::Vector{OperatingEnvelopeResult}
    start_scales::Vector{Float64}
    selected_index::Int
    diagnostics::Dict{String,Any}
end

solve_status(result::OperatingEnvelopeMultistartResult) =
    solve_status(result.selected)
solve_diagnostics(result::OperatingEnvelopeMultistartResult) =
    (run_count=length(result.runs), selected_index=result.selected_index,
     global_certificate=false,
     maximum_capacity_spread=get(result.diagnostics,
                                 "maximum_capacity_spread_W", NaN))

const _FAIRNESS_KINDS =
    (:equal, :max_total, :proportional, :alpha, :max_min, :equal_curtailment)
const _NORMALIZATIONS = (:none, :capacity, :request, :custom)
const _SECURITY_MODES = (:bound_point, :corners)
const _DIRECTIONS = (:export, :import)

_power_base(ctx) = begin
    bases = _opf_bases(ctx)
    bases === nothing ? 1.0 : Float64(bases.s_base)
end
_capacity_limit(cp::ConnectionPoint, direction::Symbol) =
    direction == :export ? cp.export_max : cp.import_max

function _as_policy(fairness)
    fairness isa FairnessPolicy && return fairness
    fairness isa Symbol || throw(ArgumentError(
        "fairness must be a Symbol or FairnessPolicy, got $(typeof(fairness))"))
    fairness == :sum          && return FairnessPolicy(kind=:max_total)
    fairness == :proportional && return FairnessPolicy(kind=:proportional)
    fairness == :equal        && return FairnessPolicy(kind=:equal)
    fairness in _FAIRNESS_KINDS && return FairnessPolicy(kind=fairness)
    throw(ArgumentError("unknown fairness policy :$fairness"))
end

function _scenario_groups(nets)
    nets isa DOEScenarioSet && return [
        [scenario.network for scenario in group] for group in nets.intervals]
    nets isa Dict{String,Any} && return [[nets]]
    nets isa AbstractVector || throw(ArgumentError(
        "nets must be a network Dict, interval/scenario vectors, or DOEScenarioSet"))
    isempty(nets) && throw(ArgumentError("need at least one interval"))
    if all(n -> n isa Dict{String,Any}, nets)
        return [[n] for n in nets]
    elseif all(g -> g isa AbstractVector && !isempty(g) &&
                    all(n -> n isa Dict{String,Any}, g), nets)
        return [collect(g) for g in nets]
    end
    throw(ArgumentError(
        "nets must contain only network Dicts or only non-empty vectors of network Dicts"))
end

function _scenario_provenance_groups(nets, groups=_scenario_groups(nets))
    if nets isa DOEScenarioSet
        return [[Dict{String,Any}(
            "id" => scenario.id,
            "role" => scenario.role,
            "weight" => scenario.weight,
            "source" => scenario.source,
            "generation_method" => scenario.generation_method,
            "seed" => scenario.seed,
            "timestamp" => scenario.timestamp,
            "metadata" => copy(scenario.metadata),
            "network_hash" => _doe_sha256(scenario.network))
            for scenario in group] for group in nets.intervals]
    end
    return [[Dict{String,Any}(
        "id" => "interval_$(interval)_scenario_$(scenario)",
        "role" => :unspecified,
        "weight" => nothing,
        "source" => "untyped_network_input",
        "generation_method" => :unspecified,
        "seed" => nothing,
        "timestamp" => nothing,
        "metadata" => Dict{String,Any}(),
        "network_hash" => _doe_sha256(net))
        for (scenario, net) in enumerate(group)]
        for (interval, group) in enumerate(groups)]
end

function _validate_policy(policy::FairnessPolicy, cps, direction)
    policy.kind in _FAIRNESS_KINDS || throw(ArgumentError(
        "fairness kind must be one of $(_FAIRNESS_KINDS), got :$(policy.kind)"))
    policy.normalization in _NORMALIZATIONS || throw(ArgumentError(
        "normalization must be one of $(_NORMALIZATIONS), got :$(policy.normalization)"))
    isfinite(policy.alpha) || throw(ArgumentError("fairness alpha must be finite"))
    isfinite(policy.epsilon) && policy.epsilon > 0 || throw(ArgumentError(
        "fairness epsilon must be finite and > 0"))
    unknown = setdiff(Set(keys(policy.weights)), Set(cp.id for cp in cps))
    isempty(unknown) || throw(ArgumentError("fairness weights contain unknown ids: $(collect(unknown))"))
    for cp in cps
        w = get(policy.weights, cp.id, 1.0)
        isfinite(w) && w > 0 || throw(ArgumentError(
            "fairness weight for '$(cp.id)' must be finite and > 0"))
        _fairness_reference(cp, policy.normalization, direction)
    end
end

function _fairness_reference(cp::ConnectionPoint, normalization::Symbol,
                             direction::Symbol)
    ref = if normalization == :none
        1.0
    elseif normalization == :capacity
        _capacity_limit(cp, direction)
    elseif normalization == :request
        cp.requested === nothing && throw(ArgumentError(
            "connection '$(cp.id)' needs requested for normalization=:request"))
        cp.requested
    else
        cp.normalization === nothing && throw(ArgumentError(
            "connection '$(cp.id)' needs normalization for normalization=:custom"))
        cp.normalization
    end
    isfinite(ref) && ref > 0 || throw(ArgumentError(
        "fairness reference for '$(cp.id)' must be finite and > 0, got $ref"))
    return Float64(ref)
end

function _validate_connection_points(groups, cps, policy, direction, security,
                                     max_exact_corners)
    direction in _DIRECTIONS || throw(ArgumentError(
        "direction must be one of $(_DIRECTIONS), got :$direction"))
    security in _SECURITY_MODES || throw(ArgumentError(
        "security must be one of $(_SECURITY_MODES), got :$security"))
    isempty(cps) && throw(ArgumentError("need at least one connection point"))
    ids = [cp.id for cp in cps]
    allunique(ids) || throw(ArgumentError("connection-point ids must be unique: $ids"))
    for cp in cps
        isempty(cp.id) && throw(ArgumentError("connection-point id must not be empty"))
        isempty(cp.bus) && throw(ArgumentError("connection '$(cp.id)' bus must not be empty"))
        isempty(cp.phase_terminals) && throw(ArgumentError(
            "connection '$(cp.id)' needs at least one phase terminal"))
        allunique(cp.phase_terminals) || throw(ArgumentError(
            "connection '$(cp.id)' phase terminals must be unique"))
        cp.neutral in cp.phase_terminals && throw(ArgumentError(
            "connection '$(cp.id)' neutral cannot also be a phase terminal"))
        for (name, value) in (("export_max", cp.export_max), ("import_max", cp.import_max))
            isfinite(value) && value >= 0 || throw(ArgumentError(
                "connection '$(cp.id)' $name must be finite and >= 0, got $value"))
        end
        for (name, value) in (("requested", cp.requested),
                              ("normalization", cp.normalization))
            value === nothing || (isfinite(value) && value > 0) || throw(ArgumentError(
                "connection '$(cp.id)' $name must be finite and > 0"))
        end
    end
    all(_capacity_limit(cp, direction) == 0 for cp in cps) && throw(ArgumentError(
        "all connection points have zero $(direction) capacity"))
    max_exact_corners >= 1 || throw(ArgumentError("max_exact_corners must be >= 1"))
    security == :corners && length(cps) > max_exact_corners && throw(ArgumentError(
        "security=:corners needs 2^N AC contexts; got N=$(length(cps)) > " *
        "max_exact_corners=$max_exact_corners"))
    _validate_policy(policy, cps, direction)

    for (t, group) in enumerate(groups), (s, net) in enumerate(group)
        buses = get(net, "bus", Dict())
        for cp in cps
            haskey(buses, cp.bus) || throw(ArgumentError(
                "interval $t scenario $s: connection '$(cp.id)' bus '$(cp.bus)' not found"))
            terminals = Set(String.(get(buses[cp.bus], "terminal_names", String[])))
            if cp.ibr_id === nothing
                for term in cp.phase_terminals
                    term in terminals || throw(ArgumentError(
                        "interval $t scenario $s: terminal '$term' for '$(cp.id)' not found at bus '$(cp.bus)'"))
                end
                cp.neutral === nothing || cp.neutral in terminals || throw(ArgumentError(
                    "interval $t scenario $s: neutral '$(cp.neutral)' for '$(cp.id)' not found"))
            else
                invs = get(net, "ibr", Dict())
                haskey(invs, cp.ibr_id) || throw(ArgumentError(
                    "interval $t scenario $s: IBR '$(cp.ibr_id)' for '$(cp.id)' not found"))
                inv = invs[cp.ibr_id]
                get(inv, "bus", nothing) == cp.bus || throw(ArgumentError(
                    "interval $t scenario $s: IBR '$(cp.ibr_id)' is not at bus '$(cp.bus)'"))
                topo = uppercase(String(get(inv, "topology", "FOUR_LEG")))
                topo in ("SINGLE_PHASE", "FOUR_LEG") || throw(ArgumentError(
                    "connection-bound IBR '$(cp.ibr_id)' topology '$topo' is not supported; " *
                    "use SINGLE_PHASE or FOUR_LEG for prescribed Q-V control"))
            end
        end
    end
end

# Legacy connection port. It is intentionally aggregate unity-PF and exists for
# backward compatibility and simple teaching examples. Real PV/battery studies
# should bind `ConnectionPoint.ibr_id` to the engine's prescribed-control model.
function _stamp_legacy_port!(ctx, cp::ConnectionPoint)
    m = _opf_model(ctx)
    bus = cp.bus
    P = zero(JuMP.QuadExpr)
    Q = zero(JuMP.QuadExpr)
    for ph in cp.phase_terminals
        cr = JuMP.@variable(m, base_name="cr_doe_$(cp.id)_$(ph)")
        ci = JuMP.@variable(m, base_name="ci_doe_$(cp.id)_$(ph)")
        dvr, dvi = _dv(ctx, bus, ph, cp.neutral)
        P += JuMP.@expression(m, dvr*cr + dvi*ci)
        Q += JuMP.@expression(m, dvi*cr - dvr*ci)
        BMOPFTools.add_terminal_injection!(ctx, bus, ph, cr, ci)
        if cp.neutral !== nothing
            BMOPFTools.add_terminal_injection!(ctx, bus, cp.neutral, -cr, -ci)
        end
    end
    JuMP.@constraint(m, Q == 0.0)
    return P
end

# Recover the active-power expression already stamped by a BMOPFTools IBR. This
# adds no reactive decision: the engine's own constant-PF / Volt-VAr equality is
# retained unchanged.
function _ibr_active_power(ctx, cp::ConnectionPoint)
    m = _opf_model(ctx)
    inv = _opf_network(ctx)["ibr"][cp.ibr_id]
    bus = String(inv["bus"])
    tm = String.(inv["terminal_map"])
    topo = uppercase(String(get(inv, "topology", "FOUR_LEG")))
    vr, vi = _opf_voltage_maps(ctx)

    if topo == "SINGLE_PHASE"
        ph, ref = tm[1], tm[2]
        dvr = JuMP.@expression(m, vr[(bus,ph)] - vr[(bus,ref)])
        dvi = JuMP.@expression(m, vi[(bus,ph)] - vi[(bus,ref)])
        cr = _opf_ibr_current(ctx, cp.ibr_id, 1)
        ci = _opf_ibr_current(ctx, cp.ibr_id, 1; component=:imag)
        return JuMP.@expression(m,
            dvr*cr + dvi*ci)
    end

    neutral = cp.neutral !== nothing && cp.neutral in tm ? cp.neutral : tm[end]
    phases = [term for term in tm if term != neutral]
    terms = JuMP.QuadExpr[]
    for (idx, ph) in enumerate(phases)
        dvr = JuMP.@expression(m, vr[(bus,ph)] - vr[(bus,neutral)])
        dvi = JuMP.@expression(m, vi[(bus,ph)] - vi[(bus,neutral)])
        cr = _opf_ibr_current(ctx, cp.ibr_id, idx)
        ci = _opf_ibr_current(ctx, cp.ibr_id, idx; component=:imag)
        push!(terms, JuMP.@expression(m,
            dvr*cr + dvi*ci))
    end
    return sum(terms)
end

_connection_active_power(ctx, cp) =
    cp.ibr_id === nothing ? _stamp_legacy_port!(ctx, cp) : _ibr_active_power(ctx, cp)

function _dispatch_patterns(n::Int, security::Symbol)
    security == :bound_point && return [ones(Float64, n)]
    return [[Float64((mask >> (i-1)) & 1) for i in 1:n]
            for mask in 0:(Int(2)^n - 1)]
end

function _set_fairness_objective!(model, cap, cps, policy, direction, power_base;
                                  temporal_history=nothing, temporal_dt_h=1.0)
    refs = Dict(cp.id => _fairness_reference(cp, policy.normalization, direction) /
                              power_base for cp in cps)
    x = Dict(cp.id => JuMP.@expression(model, cap[cp.id] / refs[cp.id]) for cp in cps)
    weights = Dict(cp.id => get(policy.weights, cp.id, 1.0) for cp in cps)
    limits = Dict(cp.id => _capacity_limit(cp, direction) / power_base for cp in cps)

    if temporal_history !== nothing
        level = JuMP.@variable(model, base_name="doe_cumulative_fairness", lower_bound=0.0)
        for cp in cps
            prior = get(temporal_history, cp.id, 0.0)
            JuMP.@constraint(model, prior + temporal_dt_h * x[cp.id] / weights[cp.id] >= level)
        end
        JuMP.@objective(model, Max, level)
        return (kind=:cumulative_max_min, level=level, x=x, weights=weights)
    elseif policy.kind == :equal
        level = JuMP.@variable(model, base_name="doe_equal_level", lower_bound=0.0)
        for cp in cps
            JuMP.@constraint(model, x[cp.id] == level)
        end
        JuMP.@objective(model, Max, level)
        return (kind=:single_stage, level=nothing)
    elseif policy.kind == :max_total
        JuMP.@objective(model, Max, sum(weights[cp.id] * x[cp.id] for cp in cps))
        return (kind=:single_stage, level=nothing)
    elseif policy.kind in (:proportional, :alpha)
        α = policy.kind == :proportional ? 1.0 : policy.alpha
        ε = policy.epsilon
        if isapprox(α, 1.0; atol=1e-12, rtol=0.0)
            JuMP.@objective(model, Max,
                sum(weights[cp.id] * log(x[cp.id] + ε) for cp in cps))
        elseif isapprox(α, 0.0; atol=1e-12, rtol=0.0)
            JuMP.@objective(model, Max,
                sum(weights[cp.id] * x[cp.id] for cp in cps))
        else
            JuMP.@objective(model, Max,
                sum(weights[cp.id] * (x[cp.id] + ε)^(1.0-α) / (1.0-α)
                    for cp in cps))
        end
        return (kind=:single_stage, level=nothing)
    elseif policy.kind == :max_min
        level = JuMP.@variable(model, base_name="doe_max_min_level", lower_bound=0.0)
        for cp in cps
            JuMP.@constraint(model, x[cp.id] >= level)
        end
        JuMP.@objective(model, Max, level)
        return (kind=:max_min, level=level, x=x, weights=weights)
    else
        level = JuMP.@variable(model, base_name="doe_curtailment", lower_bound=0.0)
        for cp in cps
            JuMP.@constraint(model, (limits[cp.id] - cap[cp.id]) / refs[cp.id] == level)
        end
        JuMP.@objective(model, Min, level)
        return (kind=:single_stage, level=nothing)
    end
end

_has_primal(model) = _publishable(_solve_outcome(model))

function _result_margins(result, net)
    # Store (normalized margin, physical margin, location). Normalization by the
    # corresponding declared limit makes different instances of one constraint
    # family comparable while preserving physical-unit margins for reporting.
    best = Dict{String,Tuple{Float64,Float64,String}}()
    consider!(kind, margin, label, scale) = begin
        isfinite(margin) && isfinite(scale) && scale > 0 || return
        normalized = margin / scale
        if !haskey(best, kind) || normalized < best[kind][1]
            best[kind] = (normalized, margin, label)
        end
    end

    # Phase-to-ground voltage bounds. More specialized BMOPFTools voltage
    # constraints remain enforced even when not reducible to one scalar margin.
    for (bus_id, bus) in get(net, "bus", Dict())
        rb = get(get(result, "bus", Dict()), bus_id, nothing)
        rb isa Dict || continue
        terminals = String.(get(bus, "terminal_names", String[]))
        grounded = Set(String.(get(bus, "perfectly_grounded_terminals", String[])))
        neutral_term = findfirst(t -> lowercase(t) in ("n", "neutral"), terminals)
        phases = [t for (idx, t) in enumerate(terminals)
                  if !(t in grounded) && idx != neutral_term]
        for (field, sense) in (("v_min", :lower), ("v_max", :upper))
            limits = get(bus, field, nothing)
            limits isa AbstractVector || continue
            for (idx, limit) in enumerate(limits)
                idx <= length(phases) || break
                term = phases[idx]
                haskey(rb, term) || continue
                vm = Float64(rb[term]["vm"])
                margin = sense == :lower ? vm - Float64(limit) : Float64(limit) - vm
                consider!("voltage", margin, "bus:$bus_id:$term:$field",
                          max(abs(Float64(limit)), eps(Float64)))
            end
        end

        vneg_max = get(bus, "vneg_max", nothing)
        if vneg_max isa Number && length(phases) == 3
            vn = if neutral_term === nothing || !haskey(rb, terminals[neutral_term])
                0.0 + 0.0im
            else
                rn = rb[terminals[neutral_term]]
                Float64(rn["vr"]) + im*Float64(rn["vi"])
            end
            V = ComplexF64[]
            for term in phases
                rt = rb[term]
                push!(V, Float64(rt["vr"]) + im*Float64(rt["vi"]) - vn)
            end
            a = cis(2pi/3)
            V2 = (V[1] + a^2*V[2] + a*V[3]) / 3
            consider!("negative_sequence", Float64(vneg_max) - abs(V2),
                      "bus:$bus_id:vneg_max",
                      max(abs(Float64(vneg_max)), eps(Float64)))
        end
    end

    # Per-conductor line ampacity inherited from the referenced linecode.
    for (line_id, line) in get(net, "line", Dict())
        rl = get(get(result, "line", Dict()), line_id, nothing)
        rl isa Dict || continue
        lc = get(get(net, "linecode", Dict()), get(line, "linecode", ""), Dict())
        raw = get(line, "i_max", get(lc, "i_max", nothing))
        raw === nothing && continue
        limits = raw isa Number ? fill(Float64(raw), length(get(line, "terminal_map_from", []))) :
                                  Float64.(raw)
        terms = String.(get(line, "terminal_map_from", String[]))
        for (idx, limit) in enumerate(limits)
            idx <= length(terms) || break
            term = terms[idx]
            haskey(rl, term) || continue
            current = max(Float64(get(rl[term], "cm_fr", 0.0)),
                          Float64(get(rl[term], "cm_to", 0.0)))
            consider!("thermal", limit - current, "line:$line_id:$term:i_max",
                      max(abs(limit), eps(Float64)))
        end
    end
    return best
end

function _merge_margins(results_and_nets)
    worst = Dict{String,Tuple{Float64,Float64,String}}()
    for (result, net) in results_and_nets
        for (kind, item) in _result_margins(result, net)
            if !haskey(worst, kind) || item[1] < worst[kind][1]
                worst[kind] = item
            end
        end
    end
    tolerances = Dict("voltage"=>0.05, "thermal"=>0.01,
                      "negative_sequence"=>0.01)
    binding = [item[3] for (kind, item) in worst
               if item[2] <= get(tolerances, kind, 0.0)]
    sort!(binding)
    return Dict{String,Any}(
        "minimum_margins" => Dict(kind=>item[2] for (kind, item) in worst),
        "minimum_normalized_margins" =>
            Dict(kind=>item[1] for (kind, item) in worst),
        "minimum_margin_locations" => Dict(kind=>item[3] for (kind, item) in worst),
        "binding_constraints" => binding)
end

function _optimize_fairness!(model, stage; max_min_tolerance)
    JuMP.optimize!(model)
    stage.kind in (:max_min, :cumulative_max_min) || return
    _has_primal(model) || return
    best = JuMP.value(stage.level)
    JuMP.@constraint(model, stage.level >= best - max_min_tolerance)
    JuMP.@objective(model, Max,
        sum(stage.weights[id] * stage.x[id] for id in keys(stage.x)))
    JuMP.optimize!(model)
end

function _fairness_metrics(alloc, cps, policy, direction; cumulative=nothing)
    normalized = Dict{String,Float64}()
    curtailment = Dict{String,Float64}()
    for cp in cps
        value = Float64(alloc[cp.id])
        ref = _fairness_reference(cp, policy.normalization, direction)
        normalized[cp.id] = value / ref
        limit = _capacity_limit(cp, direction)
        curtailment[cp.id] = limit > 0 ? 1 - value / limit : 0.0
    end
    values_ = collect(values(normalized))
    denominator = length(values_) * sum(x^2 for x in values_)
    # Roundoff can push the mathematically bounded index marginally above one
    # (for example 1.0000000000000002 for equal allocations).
    jain = denominator > 0 ? clamp(sum(values_)^2 / denominator, 0.0, 1.0) : 1.0
    out = Dict{String,Any}(
        "total_capacity_W" => sum(values(alloc)),
        "normalized_allocations" => normalized,
        "curtailment_fraction" => curtailment,
        "jain_index" => jain,
        "min_normalized" => minimum(values_),
        "max_normalized" => maximum(values_),
        "mean_normalized" => sum(values_) / length(values_))
    cumulative === nothing || (out["cumulative_normalized"] = copy(cumulative))
    return out
end

function _validate_temporal_fairness(mode, history, dt_h, cps)
    mode in (:none, :cumulative_max_min) || throw(ArgumentError(
        "temporal_fairness must be :none or :cumulative_max_min"))
    isfinite(dt_h) && dt_h > 0 || throw(ArgumentError("temporal_dt_h must be finite and > 0"))
    unknown = setdiff(Set(keys(history)), Set(cp.id for cp in cps))
    isempty(unknown) || throw(ArgumentError("fairness_history contains unknown ids: $(collect(unknown))"))
    for (id, value) in history
        isfinite(value) && value >= 0 || throw(ArgumentError(
            "fairness_history for '$id' must be finite and >= 0"))
    end
end

function _capacity_trajectory(capacities, cps, T, direction)
    source = capacities isa OperatingEnvelopeResult ? capacities.envelope : capacities
    source isa AbstractDict || throw(ArgumentError(
        "capacities must be an OperatingEnvelopeResult or a dictionary keyed by connection-point id"))
    result = Vector{Dict{String,Float64}}(undef, T)
    for t in 1:T
        item = Dict{String,Float64}()
        for cp in cps
            haskey(source, cp.id) || throw(ArgumentError("capacities missing id '$(cp.id)'"))
            raw = source[cp.id]
            value = raw isa AbstractVector ? (length(raw) == T || throw(ArgumentError(
                "capacity vector for '$(cp.id)' must have $T entries")); raw[t]) : raw
            value isa Number && isfinite(value) && value >= 0 || throw(ArgumentError(
                "capacity for '$(cp.id)' at interval $t must be finite and >= 0"))
            value <= _capacity_limit(cp, direction) + 1e-8 || throw(ArgumentError(
                "capacity for '$(cp.id)' exceeds its declared $(direction) limit"))
            item[cp.id] = Float64(value)
        end
        result[t] = item
    end
    return result
end

function _doe_issued_values_from_result(capacities, interval,
                                        policy::DOEControlPolicy,
                                        stages)
    capacities isa OperatingEnvelopeResult ||
        return Dict{_DOEPolicyValueKey,Float64}(), :capacity_values_only
    interval <= length(capacities.diagnostics) ||
        return Dict{_DOEPolicyValueKey,Float64}(), :missing_interval_diagnostics
    diagnostics = capacities.diagnostics[interval]
    get(diagnostics, "control_policy_signature", nothing) ==
        _doe_sha256(_doe_policy_manifest(policy)) ||
        return Dict{_DOEPolicyValueKey,Float64}(), :policy_mismatch
    records = get(diagnostics, "issued_control_values", nothing)
    records isa AbstractVector ||
        return Dict{_DOEPolicyValueKey,Float64}(), :not_recorded
    values_ = Dict{_DOEPolicyValueKey,Float64}()
    for record in records
        record isa AbstractDict || continue
        stage = get(record, "stage", nothing)
        stage in stages || continue
        key = _doe_control_key(
            record["component"], record["id"], record["quantity"],
            record["position"])
        information_index = stage == :issue ? 0 : Int(record["scenario"])
        value = Float64(record["value"])
        isfinite(value) || throw(ArgumentError(
            "recorded issued control value is not finite"))
        values_[(key, information_index)] = value
    end
    source = if !isempty(values_)
        :operating_envelope_result
    elseif isempty(records)
        :no_issued_controls
    else
        :no_controls_at_selected_stages
    end
    return values_, source
end

function _verification_patterns(utilizations, n)
    utilizations == :bound_point && return _dispatch_patterns(n, :bound_point)
    utilizations == :corners && return _dispatch_patterns(n, :corners)
    utilizations isa AbstractVector || throw(ArgumentError(
        "utilizations must be :bound_point, :corners, or a vector of utilization vectors"))
    patterns = Vector{Vector{Float64}}()
    for pattern in utilizations
        pattern isa AbstractVector && length(pattern) == n || throw(ArgumentError(
            "each utilization vector must contain one value per connection point"))
        values_ = Float64.(pattern)
        all(x -> isfinite(x) && 0 <= x <= 1, values_) || throw(ArgumentError(
            "utilization values must be finite and lie in [0, 1]"))
        push!(patterns, values_)
    end
    isempty(patterns) && throw(ArgumentError("need at least one utilization point"))
    return patterns
end

const _DOEControlKey = NamedTuple{
    (:component, :id, :quantity, :position),
    Tuple{Symbol,String,Symbol,Int}}
const _DOEPolicyValueKey = Tuple{_DOEControlKey,Int}

_doe_control_key(component, id, quantity, position) =
    _DOEControlKey((component, String(id), quantity, Int(position)))
_doe_rule_identity(key::_DOEControlKey) = (key.component, key.id, key.quantity)

function _doe_rule_map(policy::DOEControlPolicy)
    Dict((rule.component, rule.id, rule.quantity) => rule for rule in policy.rules)
end

function _doe_numeric_at(value, position::Int)
    value isa Real && return Float64(value)
    value isa AbstractVector && position <= length(value) &&
        value[position] isa Real && return Float64(value[position])
    return nothing
end

function _doe_fixed_bounds(device, quantity::Symbol, position::Int)
    prefix = quantity == :active_power ? "p" : "q"
    lower = _doe_numeric_at(get(device, "$(prefix)_min", nothing), position)
    upper = _doe_numeric_at(get(device, "$(prefix)_max", nothing), position)
    lower === nothing && return false
    upper === nothing && return false
    return isapprox(lower, upper; rtol=1e-12,
                    atol=1e-12 * max(1.0, abs(lower), abs(upper)))
end

function _doe_native_control_handles(record, cps, registrations)
    handles = Dict{_DOEControlKey,Any}()
    declarations = Dict{_DOEControlKey,DOEControlRegistration}()
    constraint_families = Dict{Any,Set{Symbol}}()
    object_keys = BMOPFTools.opf_object_keys(record.ctx)
    for object_key in object_keys
        if object_key.kind == :constraint && object_key.index !== nothing
            push!(get!(constraint_families, object_key.index, Set{Symbol}()),
                  object_key.family)
        elseif object_key.kind == :variable && object_key.family in (:p_ibr, :q_ibr)
            index = object_key.index
            index isa Tuple && length(index) >= 2 || continue
            quantity = object_key.family == :p_ibr ? :active_power : :reactive_power
            key = _doe_control_key(:ibr, index[1], quantity, index[2])
            handles[key] = BMOPFTools.opf_object(record.ctx, object_key)
        elseif object_key.kind == :variable && object_key.family == :tap
            index = object_key.index
            if index isa Tuple
                length(index) >= 2 || continue
                key = _doe_control_key(:transformer, index[1], :tap, index[2])
            else
                key = _doe_control_key(:transformer, index, :tap, 1)
            end
            handles[key] = BMOPFTools.opf_object(record.ctx, object_key)
        end
    end

    available_keys = Set(object_keys)
    for registration in registrations
        registration.handle in available_keys || continue
        key = _doe_control_key(registration.component, registration.id,
                               registration.quantity, registration.position)
        haskey(handles, key) && throw(ArgumentError(
            "custom DOE control registration collides with native control $(key)"))
        handles[key] = BMOPFTools.opf_object(record.ctx, registration.handle)
        declarations[key] = registration
    end

    envelope_ibrs = Set(cp.ibr_id for cp in cps if cp.ibr_id !== nothing)
    classifications = Dict{_DOEControlKey,Tuple{Symbol,Union{Nothing,Symbol}}}()
    for key in keys(handles)
        if haskey(declarations, key)
            registration = declarations[key]
            classifications[key] = (registration.native_classification,
                                    registration.automatic_law)
            continue
        end
        if key.component == :transformer
            classifications[key] = (:free, nothing)
            continue
        end
        device = get(get(record.net, "ibr", Dict()), key.id, Dict())
        index = (key.id, key.position)
        families = get(constraint_families, index, Set{Symbol}())
        if key.quantity == :active_power && key.id in envelope_ibrs
            classifications[key] = (:envelope, nothing)
        elseif key.quantity == :reactive_power && :ibr_power_factor in families
            classifications[key] = (:local_law, :power_factor)
        elseif key.quantity == :reactive_power && :ibr_q_volt_var in families
            classifications[key] = (:local_law, :volt_var)
        elseif _doe_fixed_bounds(device, key.quantity, key.position)
            classifications[key] = (:fixed_data, nothing)
        else
            classifications[key] = (:free, nothing)
        end
    end
    return handles, classifications, declarations
end


function _doe_control_scale(record, key, declarations)
    registration = get(declarations, key, nothing)
    registration !== nothing && return registration.model_to_canonical_scale
    return key.quantity in (:active_power, :reactive_power) ?
        _power_base(record.ctx) : 1.0
end

function _doe_control_unit(key, declarations)
    registration = get(declarations, key, nothing)
    registration !== nothing && return registration.canonical_unit
    key.quantity == :active_power && return :W
    key.quantity == :reactive_power && return :var
    return :effective_ratio_coefficient
end

function _doe_control_metadata(key, declarations)
    registration = get(declarations, key, nothing)
    registration === nothing ? Dict{String,Any}() : copy(registration.metadata)
end

function _doe_policy_stage(policy, rules, used_rules, key, native_states)
    identity = _doe_rule_identity(key)
    rule = get(rules, identity, nothing)
    rule !== nothing && push!(used_rules, identity)
    distinct_states = unique(first.(native_states))

    if all(==(:fixed_data), distinct_states)
        rule === nothing || throw(ArgumentError(
            "control rule for $(identity) targets a value already fixed by network data"))
        return :fixed_data, :network_data
    elseif all(==(:envelope), distinct_states)
        rule === nothing || throw(ArgumentError(
            "control rule for $(identity) targets active power already governed by the DOE allocation"))
        return :envelope, :formulation
    elseif all(==(:local_law), distinct_states)
        if rule !== nothing && rule.stage != :local_law
            throw(ArgumentError(
                "control rule for $(identity) cannot replace a prescribed local law with :$(rule.stage); change the network control profile instead"))
        end
        return :local_law, rule === nothing ? :network_control_profile : :explicit_rule
    end

    stage = rule === nothing ? policy.default_stage : rule.stage
    stage == :local_law && throw(ArgumentError(
        "control rule for $(identity) requests :local_law, but no prescribed local-law equality was discovered"))
    if length(distinct_states) > 1 && stage != :context
        throw(ArgumentError(
            "control $(identity) changes native classification across scenarios ($(distinct_states)); only :context is currently safe"))
    end
    return stage, rule === nothing ? :policy_default : :explicit_rule
end

function _doe_link_groups(records, present_indices, stage::Symbol)
    stage == :issue && return [present_indices]
    stage == :scenario || return Vector{Vector{Int}}()
    groups = Vector{Vector{Int}}()
    for scenario in sort!(unique(records[index].scenario for index in present_indices))
        push!(groups, [index for index in present_indices
                       if records[index].scenario == scenario])
    end
    return groups
end

function _doe_generator_recourses(records)
    recourses = Dict{_DOEControlKey,Set{Int}}()
    for (record_index, record) in enumerate(records)
        for (id, device) in get(record.net, "generator", Dict())
            for quantity in (:active_power, :reactive_power)
                prefix = quantity == :active_power ? "p" : "q"
                lower = get(device, "$(prefix)_min", nothing)
                upper = get(device, "$(prefix)_max", nothing)
                positions = max(lower isa AbstractVector ? length(lower) : lower isa Real ? 1 : 0,
                                upper isa AbstractVector ? length(upper) : upper isa Real ? 1 : 0)
                for position in 1:positions
                    _doe_fixed_bounds(device, quantity, position) && continue
                    key = _doe_control_key(:generator, id, quantity, position)
                    push!(get!(recourses, key, Set{Int}()), record_index)
                end
            end
        end
    end
    return recourses
end

function _apply_doe_control_policy!(model, records, cps,
                                    policy::DOEControlPolicy;
                                    fixed_control_values=Dict{_DOEControlKey,Float64}(),
                                    fixed_policy_control_values=
                                        Dict{_DOEPolicyValueKey,Float64}())
    rules = _doe_rule_map(policy)
    used_rules = Set{Tuple{Symbol,String,Symbol}}()
    handles_by_record = Vector{Dict{_DOEControlKey,Any}}(undef, length(records))
    declarations_by_record = Vector{Dict{_DOEControlKey,DOEControlRegistration}}(
        undef, length(records))
    classes_by_record = Vector{Dict{_DOEControlKey,Tuple{Symbol,Union{Nothing,Symbol}}}}(
        undef, length(records))
    all_keys = Set{_DOEControlKey}()
    for (index, record) in enumerate(records)
        handles, classifications, declarations = _doe_native_control_handles(
            record, cps, policy.registrations)
        handles_by_record[index] = handles
        classes_by_record[index] = classifications
        declarations_by_record[index] = declarations
        union!(all_keys, keys(handles))
    end
    discovered_registrations = Set(
        (registration.component, registration.id, registration.quantity,
         registration.position)
        for declarations in declarations_by_record
        for registration in values(declarations))
    declared_registrations = Set(
        (registration.component, registration.id, registration.quantity,
         registration.position) for registration in policy.registrations)
    missing_registrations = setdiff(
        declared_registrations, discovered_registrations)
    isempty(missing_registrations) || throw(ArgumentError(
        "DOE control registrations did not resolve to a model variable: $(collect(missing_registrations))"))

    audit = Dict{String,Any}[]
    link_count = 0
    fix_count = 0
    used_policy_values = Set{_DOEPolicyValueKey}()
    for key in sort!(collect(all_keys); by=key ->
                     (string(key.component), key.id, string(key.quantity), key.position))
        present = [index for index in eachindex(records)
                   if haskey(handles_by_record[index], key)]
        native_states = [classes_by_record[index][key] for index in present]
        stage, source = _doe_policy_stage(
            policy, rules, used_rules, key, native_states)
        local_laws = unique(last.(native_states))
        filter!(!isnothing, local_laws)
        links = 0
        equality_groups = Vector{Vector{Dict{String,Int}}}()
        if stage in (:issue, :scenario)
            if stage == :issue && length(present) != length(records)
                throw(ArgumentError(
                    "issue-time control $(_doe_rule_identity(key)) is absent from one or more DOE contexts"))
            end
            for group in _doe_link_groups(records, present, stage)
                isempty(group) && continue
                push!(equality_groups, [Dict("scenario" => records[index].scenario,
                                             "pattern" => records[index].pattern)
                                        for index in group])
                reference_index = first(group)
                information_index = stage == :issue ? 0 :
                    records[reference_index].scenario
                # Fixing every member to the recorded issued value already
                # enforces non-anticipativity. Adding the usual linking tree as
                # well creates redundant equalities and can make Ipopt's KKT
                # system singular at an otherwise valid replay point.
                haskey(fixed_policy_control_values,
                       (key, information_index)) && continue
                reference = handles_by_record[reference_index][key]
                reference_scale = _doe_control_scale(
                    records[reference_index], key,
                    declarations_by_record[reference_index])
                for index in Iterators.drop(group, 1)
                    scale = _doe_control_scale(
                        records[index], key, declarations_by_record[index])
                    JuMP.@constraint(model,
                        scale * handles_by_record[index][key] ==
                        reference_scale * reference)
                    links += 1
                end
            end
        end
        if haskey(fixed_control_values, key)
            for index in present
                JuMP.@constraint(model,
                    handles_by_record[index][key] == fixed_control_values[key])
                fix_count += 1
            end
        end
        if stage in (:issue, :scenario)
            for index in present
                information_index = stage == :issue ? 0 : records[index].scenario
                value_key = (key, information_index)
                haskey(fixed_policy_control_values, value_key) || continue
                scale = _doe_control_scale(
                    records[index], key, declarations_by_record[index])
                JuMP.@constraint(model,
                    scale * handles_by_record[index][key] ==
                    fixed_policy_control_values[value_key])
                push!(used_policy_values, value_key)
                fix_count += 1
            end
        end
        link_count += links
        push!(audit, Dict{String,Any}(
            "component" => key.component,
            "id" => key.id,
            "quantity" => key.quantity,
            "position" => key.position,
            "native_classification" => length(unique(first.(native_states))) == 1 ?
                first(first(native_states)) : :mixed,
            "automatic_laws" => Symbol[law for law in local_laws],
            "canonical_unit" => _doe_control_unit(
                key, declarations_by_record[first(present)]),
            "metadata" => _doe_control_metadata(
                key, declarations_by_record[first(present)]),
            "stage" => stage,
            "stage_source" => source,
            "linkage_supported" => true,
            "contexts_present" => length(present),
            "context_count" => length(records),
            "equality_groups" => equality_groups,
            "link_constraints" => links))
    end

    unknown_fixed = setdiff(Set(keys(fixed_control_values)), all_keys)
    isempty(unknown_fixed) || throw(ArgumentError(
        "fixed replay controls were not discovered in this context: $(collect(unknown_fixed))"))
    unused_policy_values = setdiff(
        Set(keys(fixed_policy_control_values)), used_policy_values)
    isempty(unused_policy_values) || throw(ArgumentError(
        "issued policy-control values did not match the verification policy/contexts: $(collect(unused_policy_values))"))

    # Generator P/Q is currently represented by current variables whose power
    # depends on voltage. Linking those currents would impose the wrong policy,
    # so only perfect/context recourse is accepted until stable P/Q handles exist.
    for (key, present_set) in sort!(collect(_doe_generator_recourses(records));
                                    by=pair -> (pair.first.id,
                                                string(pair.first.quantity),
                                                pair.first.position))
        identity = _doe_rule_identity(key)
        rule = get(rules, identity, nothing)
        rule !== nothing && push!(used_rules, identity)
        stage = if rule !== nothing
            rule.stage
        elseif policy.on_unclassified == :context
            :context
        else
            throw(ArgumentError(
                "free generator control $(identity) has no stable P/Q linkage handle; use PerfectRecourse(), an explicit :context rule, or remove/fix the control"))
        end
        stage == :context || throw(ArgumentError(
            "generator control $(identity) cannot yet be linked at :$stage because the engine exposes current, not voltage-independent P/Q, handles"))
        push!(audit, Dict{String,Any}(
            "component" => key.component,
            "id" => key.id,
            "quantity" => key.quantity,
            "position" => key.position,
            "native_classification" => :free,
            "automatic_laws" => Symbol[],
            "canonical_unit" => key.quantity == :active_power ? :W : :var,
            "metadata" => Dict{String,Any}(),
            "stage" => :context,
            "stage_source" => rule === nothing ? :unclassified_fallback : :explicit_rule,
            "linkage_supported" => false,
            "contexts_present" => length(present_set),
            "context_count" => length(records),
            "equality_groups" => Vector{Vector{Dict{String,Int}}}(),
            "link_constraints" => 0))
    end

    unused = setdiff(Set(keys(rules)), used_rules)
    isempty(unused) || throw(ArgumentError(
        "DOE control rules did not match a discovered control: $(sort!(collect(unused)))"))
    adaptive = [item for item in audit if item["stage"] == :context &&
                item["native_classification"] in (:free, :mixed)]
    shared = [item for item in audit if item["stage"] in (:issue, :scenario)]
    diagnostics = Dict{String,Any}(
        "control_policy" => policy.name,
        "control_policy_signature" =>
            _doe_sha256(_doe_policy_manifest(policy)),
        "control_default_stage" => policy.default_stage,
        "control_audit" => audit,
        "control_link_constraints" => link_count,
        "control_fix_constraints" => fix_count,
        "shared_control_count" => length(shared),
        "adaptive_control_count" => length(adaptive),
        "control_nonanticipativity" => !isempty(shared),
        "nonanticipativity_enforced" =>
            link_count > 0 || !isempty(used_policy_values),
        "perfect_recourse_controls_present" => !isempty(adaptive),
        "ideal_recourse_used" => !isempty(adaptive),
        "all_discovered_free_controls_classified" => true,
        "custom_control_registration_count" => length(policy.registrations),
        "control_discovery_scope" =>
            :native_controls_plus_explicit_registered_extension_controls)
    return (diagnostics=diagnostics,
            handles_by_record=handles_by_record,
            declarations_by_record=declarations_by_record,
            audit=audit)
end

function _doe_key_from_audit(item)
    _doe_control_key(item["component"], item["id"], item["quantity"],
                     item["position"])
end

function _doe_context_control_values(record_index, records, handles_by_record,
                                     declarations_by_record, audit)
    values_ = Dict{String,Any}[]
    handles = handles_by_record[record_index]
    declarations = declarations_by_record[record_index]
    for item in audit
        key = _doe_key_from_audit(item)
        entry = Dict{String,Any}(
            "component" => key.component,
            "id" => key.id,
            "quantity" => key.quantity,
            "position" => key.position,
            "stage" => item["stage"])
        if haskey(handles, key)
            model_value = JuMP.value(handles[key])
            entry["value"] = model_value * _doe_control_scale(
                records[record_index], key, declarations)
            entry["unit"] = _doe_control_unit(key, declarations)
            entry["model_value"] = model_value
            entry["available"] = true
        else
            entry["available"] = false
            entry["reason"] = :no_stable_value_handle
        end
        push!(values_, entry)
    end
    return values_
end

function _doe_policy_control_values(records, handles_by_record,
                                    declarations_by_record, audit)
    values_ = Dict{_DOEPolicyValueKey,Float64}()
    records_ = Dict{String,Any}[]
    for item in audit
        stage = item["stage"]
        stage in (:issue, :scenario) || continue
        key = _doe_key_from_audit(item)
        information_indices = stage == :issue ? [0] :
            sort!(unique(record.scenario for record in records))
        for information_index in information_indices
            record_index = findfirst(eachindex(records)) do index
                haskey(handles_by_record[index], key) &&
                    (stage == :issue ||
                     records[index].scenario == information_index)
            end
            record_index === nothing && continue
            handle = handles_by_record[record_index][key]
            canonical_value = JuMP.value(handle) * _doe_control_scale(
                records[record_index], key,
                declarations_by_record[record_index])
            value_key = (key, information_index)
            values_[value_key] = canonical_value
            push!(records_, Dict{String,Any}(
                "component" => key.component,
                "id" => key.id,
                "quantity" => key.quantity,
                "position" => key.position,
                "stage" => stage,
                "scenario" => stage == :issue ? nothing : information_index,
                "value" => canonical_value,
                "unit" => _doe_control_unit(
                    key, declarations_by_record[record_index])))
        end
    end
    return values_, records_
end

function _doe_context_replay_values(record_index, handles_by_record, audit)
    handles = handles_by_record[record_index]
    values_ = Dict{_DOEControlKey,Float64}()
    unsupported = Dict{String,Any}[]
    for item in audit
        item["native_classification"] in (:free, :mixed) || continue
        key = _doe_key_from_audit(item)
        if haskey(handles, key)
            values_[key] = JuMP.value(handles[key])
        else
            push!(unsupported, item)
        end
    end
    return values_, unsupported
end

function _solve_interval_group(group, cps, policy;
                               direction, security, per_unit, s_base,
                               optimizer, verbose, solver_options,
                               volt_var_watt_eps, max_min_tolerance,
                               control_policy, control_policy_source,
                               context_hook!,
                               start_hook!,
                               temporal_history=nothing,
                               temporal_dt_h=1.0,
                               fixed_capacity=nothing,
                               patterns_override=nothing,
                               fixed_control_values=Dict{_DOEControlKey,Float64}(),
                               fixed_policy_control_values=
                                   Dict{_DOEPolicyValueKey,Float64}())
    model = JuMP.Model(optimizer)
    pb = per_unit ? s_base : 1.0
    cap = Dict{String,Any}()
    for cp in cps
        upper = _capacity_limit(cp, direction) / pb
        if fixed_capacity === nothing
            cap[cp.id] = JuMP.@variable(model, base_name="doe_capacity_$(cp.id)",
                lower_bound=0.0, upper_bound=upper)
        else
            value = fixed_capacity[cp.id] / pb
            cap[cp.id] = JuMP.@variable(model, base_name="doe_capacity_$(cp.id)",
                lower_bound=value, upper_bound=value)
        end
    end

    patterns = patterns_override === nothing ? _dispatch_patterns(length(cps), security) : patterns_override
    sign = direction == :export ? 1.0 : -1.0
    specifications = [(net=net, scenario=scenario_index,
                       pattern=pattern_index, fractions=fractions)
                      for (scenario_index, net) in enumerate(group)
                      for (pattern_index, fractions) in enumerate(patterns)]
    hook_factory = context_index -> begin
        fractions = specifications[context_index].fractions
        ctx -> begin
            for (i, cp) in enumerate(cps)
                p = _connection_active_power(ctx, cp)
                JuMP.@constraint(_opf_model(ctx),
                    p == sign * fractions[i] * cap[cp.id])
            end
            context_hook! === nothing || context_hook!(ctx)
        end
    end
    multi = build_multi_context([spec.net for spec in specifications]; model,
        hook_factory, per_unit, s_base, optimizer, verbose, solver_options,
        context_options=(volt_var_watt_eps=volt_var_watt_eps, verbose=verbose))
    records = [(ctx=multi.contexts[index], net=spec.net,
                scenario=spec.scenario, pattern=spec.pattern,
                fractions=spec.fractions)
               for (index, spec) in enumerate(specifications)]
    control_payload = _apply_doe_control_policy!(
        model, records, cps, control_policy;
        fixed_control_values=fixed_control_values,
        fixed_policy_control_values=fixed_policy_control_values)
    control_diagnostics = control_payload.diagnostics
    foreach(r -> enforce_kcl!(r.ctx), records)
    start_hook! === nothing || start_hook!(multi.contexts)
    solve_time_seconds = @elapsed begin
        if fixed_capacity === nothing
            stage = _set_fairness_objective!(model, cap, cps, policy, direction, pb;
                temporal_history=temporal_history, temporal_dt_h=temporal_dt_h)
            _optimize_fairness!(model, stage; max_min_tolerance=max_min_tolerance)
        else
            JuMP.@objective(model, Min, 0.0)
            JuMP.optimize!(model)
        end
    end

    outcome = _solve_outcome(model)
    status = string(outcome.termination_status)
    primal = string(outcome.primal_status)
    feasible = _publishable(outcome)
    primal_residual_diagnostics = Dict{String,Any}(
        "primal_residual_available" => false,
        "maximum_primal_constraint_violation" => NaN,
        "constraint_violation_count_above_1e_6" => 0)
    if outcome.has_primal
        try
            report = JuMP.primal_feasibility_report(model; atol=0.0)
            violations = Float64.(collect(values(report)))
            primal_residual_diagnostics["primal_residual_available"] = true
            primal_residual_diagnostics["maximum_primal_constraint_violation"] =
                maximum(violations; init=0.0)
            primal_residual_diagnostics["constraint_violation_count_above_1e_6"] =
                count(>(1e-6), violations)
        catch error_
            primal_residual_diagnostics["primal_residual_error"] =
                sprint(showerror, error_)
        end
    end
    total = NaN
    alloc = Dict(cp.id => NaN for cp in cps)
    snapshot = Dict{String,Any}("termination_status"=>status,
                                "primal_status"=>primal)
    objective = NaN
    context_results = OperatingEnvelopeContextResult[]
    replay_values = Vector{Tuple{Int,Int,Dict{_DOEControlKey,Float64},Vector{Dict{String,Any}}}}()
    policy_control_values = Dict{_DOEPolicyValueKey,Float64}()
    policy_control_records = Dict{String,Any}[]
    margin_diagnostics = Dict{String,Any}(
        "minimum_margins"=>Dict{String,Float64}(),
        "minimum_normalized_margins"=>Dict{String,Float64}(),
        "minimum_margin_locations"=>Dict{String,String}(),
        "binding_constraints"=>String[])
    if feasible
        policy_control_values, policy_control_records = _doe_policy_control_values(
            records, control_payload.handles_by_record,
            control_payload.declarations_by_record, control_payload.audit)
        for cp in cps
            raw = JuMP.value(cap[cp.id]) * pb
            alloc[cp.id] = clamp(raw, 0.0, _capacity_limit(cp, direction))
        end
        total = sum(values(alloc))
        objective = JuMP.objective_value(model)
        representative_index = findfirst(r -> r.scenario == 1 && all(isone, r.fractions), records)
        representative = records[something(representative_index, 1)]
        snapshot = extract_result(representative.ctx)
        checked = Tuple{Dict{String,Any},Dict{String,Any}}[]
        for (record_index, record) in enumerate(records)
            result = record.scenario == 1 && all(isone, record.fractions) ?
                     snapshot : extract_result(record.ctx)
            push!(checked, (result, record.net))
            local_margins = _merge_margins([(result, record.net)])
            control_values = _doe_context_control_values(
                record_index, records, control_payload.handles_by_record,
                control_payload.declarations_by_record,
                control_payload.audit)
            values_, unsupported = _doe_context_replay_values(
                record_index, control_payload.handles_by_record,
                control_payload.audit)
            push!(replay_values,
                (record.scenario, record.pattern, values_, unsupported))
            evidence = Dict{String,Any}(
                "evidence_scope" => :joint_policy_model,
                "margin_evaluation" =>
                    :independent_snapshot_quantity_recomputation,
                "joint_solve_time_seconds" => solve_time_seconds,
                "independent_replay" => nothing)
            merge!(evidence, local_margins)
            merge!(evidence, primal_residual_diagnostics)
            push!(context_results, OperatingEnvelopeContextResult(
                record.scenario, record.pattern, copy(record.fractions),
                status, primal, true, result, evidence, control_values))
        end
        margin_diagnostics = _merge_margins(checked)
    else
        individually_resolved = length(records) == 1
        for record in records
            push!(context_results, OperatingEnvelopeContextResult(
                record.scenario, record.pattern, copy(record.fractions),
                status, primal, individually_resolved ? false : nothing,
                Dict{String,Any}(),
                Dict{String,Any}(
                    "evidence_scope" => individually_resolved ?
                        :single_context_model : :unresolved_joint_model,
                    "joint_solve_time_seconds" => solve_time_seconds,
                    "independent_replay" => nothing),
                Dict{String,Any}[]))
        end
    end

    security_scope = if security == :corners
        :all_box_corners
    elseif length(patterns) == 1 && all(isone, only(patterns))
        :simultaneous_upper_bound_only
    else
        :explicit_utilization_points
    end
    diag = Dict{String,Any}(
        "feasible" => feasible,
        "primal_status" => primal,
        "objective" => objective,
        "solve_time_seconds" => solve_time_seconds,
        "direction" => direction,
        "security" => security,
        "security_scope" => security_scope,
        "guarantee" => :local_ac_feasibility_at_tested_dispatches,
        "global_certificate" => false,
        "solver_class" => :local_nonlinear,
        "uncertainty_semantics" => length(group) == 1 ?
            :single_declared_snapshot : :finite_scenario_set,
        "allocation_coupling" => :shared_across_scenarios_and_dispatch_points,
        "control_recourse" => control_policy.name,
        "prescribed_ibr_controls" => :retained,
        "control_policy_source" => control_policy_source,
        "issued_policy_controls_fixed" =>
            !isempty(fixed_policy_control_values),
        "issued_control_values" => policy_control_records,
        "scenario_count" => length(group),
        "dispatch_points_per_scenario" => length(patterns),
        "fairness_kind" => policy.kind,
        "normalization" => policy.normalization,
        "verification" => fixed_capacity !== nothing,
        "temporal_fairness" => temporal_history === nothing ? :none : :cumulative_max_min)
    merge!(diag, control_diagnostics)
    merge!(diag, margin_diagnostics)
    merge!(diag, primal_residual_diagnostics)
    return (status=status, alloc=alloc, total=total,
            snapshot=snapshot, diagnostics=diag,
            context_results=context_results,
            replay_values=replay_values,
            policy_control_values=policy_control_values)
end

"""
    solve_operating_envelope(nets, connection_points; kwargs...)
        -> OperatingEnvelopeResult

Calculate active-power operating-envelope capacity for each interval.

`nets` accepts three shapes:

- one network `Dict` — one interval, one forecast scenario;
- `Vector{Dict}` — several intervals, one scenario each (backward compatible);
- `Vector{Vector{Dict}}` — several intervals, each containing one or more
  forecast/model scenarios. One capacity allocation is shared by every scenario
  in an interval;
- [`DOEScenarioSet`](@ref) — the same interval/scenario structure with typed
  IDs, train/calibration/validation/test roles, optional weights, timestamps,
  generation method, seeds, and provenance metadata.

Important keywords:

- `direction=:export` or `:import`; returned capacities are positive magnitudes;
- `fairness=:equal`, the legacy symbols `:sum` / `:proportional`, or a
  [`FairnessPolicy`](@ref);
- `security=:bound_point` enforces only simultaneous full utilisation;
- `security=:corners` embeds all `2^N` zero/full-utilisation corners for every
  scenario. It is deliberately capped by `max_exact_corners=10` and is reported
  as local AC feasibility at tested points, not a global robust certificate;
- `control_policy=nothing` retains the historical pointwise-perfect-recourse
  formulation and reports it as a legacy default. Pass [`PerfectRecourse`](@ref)
  explicitly for published-work replication, or [`IssuePlusLocalLaws`](@ref)
  to share free setpoints across all scenario/utilisation contexts while
  retaining prescribed automatic IBR laws;
- `utilizations=nothing` may instead be an explicit collection of participant
  utilisation vectors in `[0,1]^N`. This supports reproducible adaptive or
  counterexample-guided studies and is reported as an explicit test set;
- `context_hook! = nothing` may register research-extension variables through
  BMOPFTools public `OpfModelKey`s before the control policy is applied. Pair
  those handles with [`DOEControlRegistration`](@ref);
- `start_hook! = nothing` receives the completed OPF contexts immediately
  before optimization and can set reproducible JuMP start values;
- `volt_var_watt_eps=2e-3` controls the engine's smooth approximation of
  mandatory IBR Volt-VAr/Volt-Watt curve corners.

Loads retain their known P/Q from each network snapshot. Connection-bound IBRs
retain their prescribed Q-V law. Other network devices, including BMOPFTools
STATCOM IBRs, remain available to the OPF; comparing otherwise identical nets
with and without a STATCOM therefore quantifies its impact on active-power DOEs.
Network voltage, phase-to-neutral, negative-sequence (`vneg_max`), branch thermal,
neutral-current, and device limits declared by BMOPFTools remain in force.

Diagnostics label the result as a local nonlinear solve, never as a global
certificate. They also record whether uncertainty is represented by one
declared snapshot or a finite scenario set, the selected control policy, and a
per-control audit of shared, automatic, fixed, and context-adaptive variables.
Connection-bound prescribed IBR controls remain enforced.
"""
function solve_operating_envelope(nets,
                                  connection_points::AbstractVector{ConnectionPoint};
                                  fairness=:equal,
                                  direction::Symbol=:export,
                                  security::Symbol=:bound_point,
                                  utilizations=nothing,
                                  per_unit::Bool=true,
                                  s_base::Float64=1e6,
                                  optimizer=Ipopt.Optimizer,
                                  verbose::Bool=false,
                                  solver_options=(),
                                  volt_var_watt_eps::Float64=2e-3,
                                  max_exact_corners::Int=10,
                                  max_min_tolerance::Float64=1e-7,
                                  control_policy::Union{Nothing,DOEControlPolicy}=nothing,
                                  context_hook! = nothing,
                                  start_hook! = nothing,
                                  temporal_fairness::Symbol=:none,
                                  fairness_history::AbstractDict=Dict{String,Float64}(),
                                  temporal_dt_h::Float64=1.0,
                                  issued_at::Union{Nothing,DateTime}=nothing,
                                  interval_seconds::Float64=300.0,
                                  validity_seconds::Float64=interval_seconds,
                                  fallback::Symbol=:missing)
    groups = _scenario_groups(nets)
    scenario_provenance = _scenario_provenance_groups(nets, groups)
    cps = collect(connection_points)
    policy = _as_policy(fairness)
    resolved_control_policy = control_policy === nothing ?
        PerfectRecourse() : control_policy
    control_policy_source = control_policy === nothing ? :legacy_default : :explicit
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be finite and > 0"))
    isfinite(volt_var_watt_eps) && volt_var_watt_eps > 0 || throw(ArgumentError(
        "volt_var_watt_eps must be finite and > 0"))
    isfinite(max_min_tolerance) && max_min_tolerance >= 0 || throw(ArgumentError(
        "max_min_tolerance must be finite and >= 0"))
    _validate_connection_points(groups, cps, policy, direction, security,
                                max_exact_corners)
    utilization_patterns = if utilizations === nothing
        nothing
    else
        security == :bound_point || throw(ArgumentError(
            "explicit utilizations cannot be combined with security=:corners"))
        _verification_patterns(utilizations, length(cps))
    end
    _validate_temporal_fairness(temporal_fairness, fairness_history, temporal_dt_h, cps)
    isfinite(interval_seconds) && interval_seconds > 0 || throw(ArgumentError(
        "interval_seconds must be finite and > 0"))
    isfinite(validity_seconds) && validity_seconds > 0 || throw(ArgumentError(
        "validity_seconds must be finite and > 0"))
    fallback in (:missing, :zero, :last_feasible) || throw(ArgumentError(
        "fallback must be :missing, :zero, or :last_feasible"))

    T = length(groups)
    ids = [cp.id for cp in cps]
    envelope = Dict(id => fill(NaN, T) for id in ids)
    total = fill(NaN, T)
    statuses = Vector{String}(undef, T)
    snapshots = Vector{Dict{String,Any}}(undef, T)
    diagnostics = Vector{Dict{String,Any}}(undef, T)
    metrics = Vector{Dict{String,Any}}(undef, T)
    schedule = Vector{Dict{String,Any}}(undef, T)
    history = Dict{String,Float64}(cp.id => Float64(get(fairness_history, cp.id, 0.0)) for cp in cps)
    last_feasible = nothing

    for t in 1:T
        solved = _solve_interval_group(groups[t], cps, policy;
            direction=direction, security=security, per_unit=per_unit,
            s_base=s_base, optimizer=optimizer, verbose=verbose,
            solver_options=solver_options, volt_var_watt_eps=volt_var_watt_eps,
            max_min_tolerance=max_min_tolerance,
            control_policy=resolved_control_policy,
            control_policy_source=control_policy_source,
            context_hook! = context_hook!,
            start_hook! = start_hook!,
            patterns_override=utilization_patterns,
            temporal_history=temporal_fairness == :cumulative_max_min ? history : nothing,
            temporal_dt_h=temporal_dt_h)
        solved.diagnostics["scenario_provenance"] =
            deepcopy(scenario_provenance[t])
        solved.diagnostics["scenario_dataset_id"] =
            nets isa DOEScenarioSet ? nets.dataset_id : nothing
        nets isa DOEScenarioSet &&
            (solved.diagnostics["uncertainty_semantics"] =
                :typed_finite_scenario_set)
        statuses[t] = solved.status
        snapshots[t] = solved.snapshot
        diagnostics[t] = solved.diagnostics
        total[t] = solved.total
        publication_source = :optimized
        published = solved.alloc
        if !solved.diagnostics["feasible"]
            if fallback == :zero
                published = Dict(id => 0.0 for id in ids)
                total[t] = 0.0
                publication_source = :zero_fallback
            elseif fallback == :last_feasible && last_feasible !== nothing
                published = copy(last_feasible)
                total[t] = sum(values(published))
                publication_source = :last_feasible_fallback
            else
                publication_source = :missing
            end
            diagnostics[t]["fallback_network_safe"] = false
        else
            last_feasible = copy(published)
            if temporal_fairness == :cumulative_max_min
                for cp in cps
                    history[cp.id] += temporal_dt_h * published[cp.id] /
                        _fairness_reference(cp, policy.normalization, direction)
                end
            end
        end
        for id in ids
            envelope[id][t] = published[id]
        end
        metrics[t] = if publication_source == :missing
            Dict{String,Any}("available"=>false, "publication_source"=>publication_source)
        else
            outcome = _fairness_metrics(published, cps, policy, direction;
                cumulative=temporal_fairness == :cumulative_max_min ? history : nothing)
            outcome["available"] = true
            outcome["publication_source"] = publication_source
            outcome
        end
        valid_from = issued_at === nothing ? nothing :
            issued_at + Millisecond(round(Int, 1000 * interval_seconds * (t - 1)))
        valid_until = valid_from === nothing ? nothing :
            valid_from + Millisecond(round(Int, 1000 * validity_seconds))
        schedule[t] = Dict{String,Any}("interval_index"=>t, "issued_at"=>issued_at,
            "valid_from"=>valid_from, "valid_until"=>valid_until,
            "publication_source"=>publication_source)
    end

    OperatingEnvelopeResult(statuses, envelope, copy(total), snapshots,
                            direction, total, diagnostics, metrics, schedule)
end

function _doe_scale_registered_starts!(contexts, scale::Float64)
    changed = 0
    seen = Set{Any}()
    for ctx in contexts
        for key in BMOPFTools.opf_object_keys(ctx; kind=:variable)
            variable = BMOPFTools.opf_object(ctx, key)
            variable isa JuMP.VariableRef || continue
            index = JuMP.index(variable)
            index in seen && continue
            push!(seen, index)
            start = JuMP.start_value(variable)
            start === nothing && continue
            target = Float64(start) * scale
            JuMP.has_lower_bound(variable) &&
                (target = max(target, JuMP.lower_bound(variable)))
            JuMP.has_upper_bound(variable) &&
                (target = min(target, JuMP.upper_bound(variable)))
            JuMP.is_fixed(variable) && (target = JuMP.fix_value(variable))
            JuMP.set_start_value(variable, target)
            changed += 1
        end
    end
    return changed
end

"""
    solve_operating_envelope_multistart(nets, connection_points;
        start_scales=(1.0, 0.9, 1.1), base_start_hook! = nothing, kwargs...)

Run the same DOE formulation from several deterministic perturbations of the
native registered-variable starting point. The selected result maximizes total
published capacity among runs with feasible primal solutions, while diagnostics
retain every termination status and the per-interval capacity spread.

This is a branch-sensitivity diagnostic, not a proof of global optimality.
`base_start_hook!`, when supplied, runs before each scale perturbation. Remaining
keywords are forwarded to [`solve_operating_envelope`](@ref).
"""
function solve_operating_envelope_multistart(
        nets, connection_points::AbstractVector{ConnectionPoint};
        start_scales=(1.0, 0.9, 1.1),
        base_start_hook! = nothing,
        kwargs...)
    scales = Float64.(collect(start_scales))
    isempty(scales) && throw(ArgumentError("start_scales must be non-empty"))
    all(scale -> isfinite(scale) && scale > 0, scales) ||
        throw(ArgumentError("start_scales must be finite and positive"))
    runs = OperatingEnvelopeResult[]
    changed_counts = Int[]
    for scale in scales
        changed = Ref(0)
        start_hook! = contexts -> begin
            base_start_hook! === nothing || base_start_hook!(contexts)
            changed[] = _doe_scale_registered_starts!(contexts, scale)
        end
        push!(runs, solve_operating_envelope(
            nets, connection_points; start_hook! = start_hook!, kwargs...))
        push!(changed_counts, changed[])
    end
    accepted = [index for index in eachindex(runs)
                if all(get(diag, "feasible", false)
                       for diag in runs[index].diagnostics)]
    selected_index = isempty(accepted) ? 1 : accepted[argmax(
        [sum(runs[index].total_capacity) for index in accepted])]
    interval_count = length(runs[selected_index].total_capacity)
    spreads = fill(NaN, interval_count)
    if !isempty(accepted)
        for interval in 1:interval_count
            values_ = [runs[index].total_capacity[interval] for index in accepted]
            spreads[interval] = maximum(values_) - minimum(values_)
        end
    end
    diagnostics = Dict{String,Any}(
        "selection" => :maximum_total_capacity_among_feasible_primals,
        "accepted_run_indices" => accepted,
        "start_changed_variable_counts" => changed_counts,
        "termination_statuses" => [run.termination_status for run in runs],
        "maximum_primal_constraint_violation" =>
            [[get(diag, "maximum_primal_constraint_violation", NaN)
              for diag in run.diagnostics] for run in runs],
        "capacity_spread_W" => spreads,
        "maximum_capacity_spread_W" =>
            isempty(accepted) ? NaN : maximum(spreads),
        "global_certificate" => false,
        "claim" => :deterministic_multistart_branch_sensitivity)
    return OperatingEnvelopeMultistartResult(
        runs[selected_index], runs, scales, selected_index, diagnostics)
end

"""
    compare_operating_envelope_policies(nets, connection_points, policies; kwargs...)

Solve the same DOE study under several fairness policies. `policies` is a
dictionary or vector of `label => fairness` pairs. The returned dictionary maps
each label to an [`OperatingEnvelopeResult`](@ref), whose `fairness_metrics`
make the capacity/fairness trade-off directly comparable.
"""
function compare_operating_envelope_policies(nets, connection_points, policies; kwargs...)
    entries = policies isa AbstractDict ? collect(pairs(policies)) : collect(policies)
    isempty(entries) && throw(ArgumentError("need at least one fairness policy"))
    result = Dict{String,OperatingEnvelopeResult}()
    for entry in entries
        entry isa Pair || throw(ArgumentError("policies must contain label => policy pairs"))
        label = string(first(entry))
        haskey(result, label) && throw(ArgumentError("duplicate policy label '$label'"))
        result[label] = solve_operating_envelope(nets, connection_points;
            fairness=last(entry), kwargs...)
    end
    return result
end

function _doe_snapshot_voltage_difference(left, right)
    buses_left = get(left, "bus", Dict())
    buses_right = get(right, "bus", Dict())
    difference = 0.0
    found = false
    for bus in intersect(Set(keys(buses_left)), Set(keys(buses_right)))
        left_bus = buses_left[bus]
        right_bus = buses_right[bus]
        for terminal in intersect(Set(keys(left_bus)), Set(keys(right_bus)))
            left_terminal = left_bus[terminal]
            right_terminal = right_bus[terminal]
            all(haskey(left_terminal, field) && haskey(right_terminal, field)
                for field in ("vr", "vi")) || continue
            delta = hypot(Float64(left_terminal["vr"]) - Float64(right_terminal["vr"]),
                          Float64(left_terminal["vi"]) - Float64(right_terminal["vi"]))
            difference = max(difference, delta)
            found = true
        end
    end
    return found ? difference : NaN
end

function _doe_relabel_context(result::OperatingEnvelopeContextResult,
                              scenario, pattern, fractions)
    OperatingEnvelopeContextResult(
        scenario, pattern, copy(fractions), result.termination_status,
        result.primal_status, result.feasible, result.snapshot,
        result.diagnostics, result.control_values)
end

"""
    verify_operating_envelope(nets, connection_points, capacities; kwargs...)
        -> OperatingEnvelopeVerification

Check an already-issued capacity trajectory without optimising it. The active
power capacities are fixed and the normal network physics, including prescribed
Q-V IBR controls and any STATCOM model present in the network, is solved at
every requested scenario and utilisation point. `utilizations` is
`:bound_point`, `:corners`, or explicit vectors in `[0, 1]^N`.

The same local-solve and finite-test-set semantics as
[`solve_operating_envelope`](@ref) apply. `control_policy` also applies during
verification, so an issued envelope can be checked under the same information
structure used to quantify it. Custom vectors are reported as
`security_scope=:explicit_utilization_points`. When `capacities` is an
`OperatingEnvelopeResult` produced under the same policy,
`replay_issued_controls=true` fixes its recorded `:issue` and `:scenario`
control values during verification. A capacity dictionary has no such issuance
record and is labelled accordingly.
"""
function verify_operating_envelope(nets,
                                   connection_points::AbstractVector{ConnectionPoint},
                                   capacities;
                                   direction::Symbol=:export,
                                   utilizations=:bound_point,
                                   per_unit::Bool=true,
                                   s_base::Float64=1e6,
                                   optimizer=Ipopt.Optimizer,
                                   verbose::Bool=false,
                                   solver_options=(),
                                   volt_var_watt_eps::Float64=2e-3,
                                   control_policy::Union{Nothing,DOEControlPolicy}=nothing,
                                   context_hook! = nothing,
                                   start_hook! = nothing,
                                   independent_replay::Bool=true,
                                   diagnose_infeasible_contexts::Bool=true,
                                   replay_issued_controls::Bool=true,
                                   replay_control_stages=(:issue, :scenario),
                                   max_exact_corners::Int=10)
    groups = _scenario_groups(nets)
    scenario_provenance = _scenario_provenance_groups(nets, groups)
    cps = collect(connection_points)
    policy = FairnessPolicy(kind=:max_total)
    resolved_control_policy = control_policy === nothing ?
        PerfectRecourse() : control_policy
    control_policy_source = control_policy === nothing ? :legacy_default : :explicit
    isfinite(s_base) && s_base > 0 || throw(ArgumentError("s_base must be finite and > 0"))
    isfinite(volt_var_watt_eps) && volt_var_watt_eps > 0 || throw(ArgumentError(
        "volt_var_watt_eps must be finite and > 0"))
    security = utilizations == :corners ? :corners : :bound_point
    replay_stage_values = replay_control_stages isa Symbol ?
        (replay_control_stages,) : replay_control_stages
    replay_stages = Set(Symbol(stage) for stage in replay_stage_values)
    isempty(replay_stages) && throw(ArgumentError(
        "replay_control_stages must be non-empty"))
    issubset(replay_stages, Set((:issue, :scenario))) || throw(ArgumentError(
        "replay_control_stages may contain only :issue and :scenario"))
    _validate_connection_points(groups, cps, policy, direction, security, max_exact_corners)
    patterns = _verification_patterns(utilizations, length(cps))
    trajectory = _capacity_trajectory(capacities, cps, length(groups), direction)
    statuses = String[]
    feasible = Bool[]
    snapshots = Dict{String,Any}[]
    diagnostics = Dict{String,Any}[]
    context_results = Vector{OperatingEnvelopeContextResult}[]
    for (t, group) in enumerate(groups)
        issued_values, issued_source = replay_issued_controls ?
            _doe_issued_values_from_result(
                capacities, t, resolved_control_policy, replay_stages) :
            (Dict{_DOEPolicyValueKey,Float64}(), :disabled)
        solved = _solve_interval_group(group, cps, policy;
            direction=direction, security=security, per_unit=per_unit, s_base=s_base,
            optimizer=optimizer, verbose=verbose, solver_options=solver_options,
            volt_var_watt_eps=volt_var_watt_eps, max_min_tolerance=1e-7,
            control_policy=resolved_control_policy,
            control_policy_source=control_policy_source,
            context_hook! = context_hook!,
            start_hook! = start_hook!,
            fixed_capacity=trajectory[t], patterns_override=patterns,
            fixed_policy_control_values=issued_values)
        solved.diagnostics["issued_control_replay_source"] = issued_source
        solved.diagnostics["issued_control_replay_count"] = length(issued_values)
        solved.diagnostics["replay_control_stages"] =
            sort!(collect(replay_stages))
        solved.diagnostics["scenario_provenance"] =
            deepcopy(scenario_provenance[t])
        solved.diagnostics["scenario_dataset_id"] =
            nets isa DOEScenarioSet ? nets.dataset_id : nothing
        nets isa DOEScenarioSet &&
            (solved.diagnostics["uncertainty_semantics"] =
                :typed_finite_scenario_set)
        push!(statuses, solved.status)
        push!(snapshots, solved.snapshot)
        interval_contexts = OperatingEnvelopeContextResult[]
        if solved.diagnostics["feasible"]
            replay_all_feasible = true
            if independent_replay
                for (context_index, (scenario, net, pattern_index, fractions)) in enumerate(
                        (scenario, net, pattern_index, fractions)
                        for (scenario, net) in enumerate(group)
                        for (pattern_index, fractions) in enumerate(patterns))
                    stored_scenario, stored_pattern, fixed_values, unsupported =
                        solved.replay_values[context_index]
                    stored_scenario == scenario && stored_pattern == pattern_index ||
                        error("internal DOE context ordering mismatch")
                    replay = _solve_interval_group([net], cps, policy;
                        direction=direction, security=:bound_point,
                        per_unit=per_unit, s_base=s_base,
                        optimizer=optimizer, verbose=verbose,
                        solver_options=solver_options,
                        volt_var_watt_eps=volt_var_watt_eps,
                        max_min_tolerance=1e-7,
                        control_policy=resolved_control_policy,
                        control_policy_source=control_policy_source,
                        context_hook! = context_hook!,
                        start_hook! = start_hook!,
                        fixed_capacity=trajectory[t],
                        patterns_override=[fractions],
                        fixed_control_values=fixed_values)
                    replay_context = only(replay.context_results)
                    replay_feasible = replay_context.feasible === true
                    replay_all_feasible &= replay_feasible
                    joint_context = solved.context_results[context_index]
                    joint_context.diagnostics["independent_replay"] = Dict{String,Any}(
                        "feasible" => replay_feasible,
                        "termination_status" => replay_context.termination_status,
                        "primal_status" => replay_context.primal_status,
                        "control_replay_complete" => isempty(unsupported),
                        "unreplayed_controls" => unsupported,
                        "maximum_voltage_difference_V" => replay_feasible ?
                            _doe_snapshot_voltage_difference(
                                joint_context.snapshot, replay_context.snapshot) : NaN,
                        "snapshot" => replay_context.snapshot,
                        "diagnostics" => replay_context.diagnostics)
                    push!(interval_contexts, joint_context)
                end
            else
                append!(interval_contexts, solved.context_results)
            end
            solved.diagnostics["joint_policy_feasible"] = true
            solved.diagnostics["independent_replay_all_feasible"] =
                independent_replay ? replay_all_feasible : nothing
            solved.diagnostics["verification_evidence"] =
                independent_replay ? :joint_model_plus_fixed_control_replay :
                :joint_model_only
            push!(feasible, replay_all_feasible)
        elseif diagnose_infeasible_contexts
            any_individually_infeasible = false
            all_individually_feasible = true
            for (scenario, net) in enumerate(group)
                for (pattern_index, fractions) in enumerate(patterns)
                    individual = _solve_interval_group([net], cps, policy;
                        direction=direction, security=:bound_point,
                        per_unit=per_unit, s_base=s_base,
                        optimizer=optimizer, verbose=verbose,
                        solver_options=solver_options,
                        volt_var_watt_eps=volt_var_watt_eps,
                        max_min_tolerance=1e-7,
                        control_policy=resolved_control_policy,
                        control_policy_source=control_policy_source,
                        context_hook! = context_hook!,
                        start_hook! = start_hook!,
                        fixed_capacity=trajectory[t],
                        patterns_override=[fractions],
                        fixed_policy_control_values=issued_values)
                    context = _doe_relabel_context(
                        only(individual.context_results), scenario,
                        pattern_index, fractions)
                    context.diagnostics["evidence_scope"] =
                        :individual_context_diagnostic
                    push!(interval_contexts, context)
                    individually_feasible = context.feasible === true
                    all_individually_feasible &= individually_feasible
                    any_individually_infeasible |= context.feasible === false
                end
            end
            solved.diagnostics["joint_policy_feasible"] = false
            solved.diagnostics["independent_contexts_all_feasible"] =
                all_individually_feasible
            solved.diagnostics["infeasibility_interpretation"] =
                any_individually_infeasible ? :at_least_one_context_infeasible :
                all_individually_feasible ?
                    :shared_control_incompatibility_or_joint_nlp_failure :
                    :individual_context_diagnostics_inconclusive
            solved.diagnostics["verification_evidence"] =
                :joint_model_plus_individual_context_diagnostics
            push!(feasible, false)
        else
            append!(interval_contexts, solved.context_results)
            solved.diagnostics["joint_policy_feasible"] = false
            solved.diagnostics["verification_evidence"] = :joint_model_only
            push!(feasible, false)
        end
        solved.diagnostics["context_result_count"] = length(interval_contexts)
        push!(context_results, interval_contexts)
        push!(diagnostics, solved.diagnostics)
    end
    return OperatingEnvelopeVerification(
        statuses, feasible, snapshots, diagnostics, context_results)
end

"""
    evaluate_operating_envelope_coverage(
        scenarios, connection_points, capacities;
        roles=(:test,), utilizations=:bound_point,
        iid_assumption=false, confidence=0.95, kwargs...)

Evaluate an issued envelope on a role-selected subset of a typed
[`DOEScenarioSet`](@ref). A scenario is a candidate violation when any tested
utilization context is individually infeasible, unresolved when numerical or
joint-policy evidence cannot classify every context, and passed otherwise.

Empirical context/scenario rates are always returned. A one-sided Hoeffding
upper bound is returned only when `iid_assumption=true`; that flag is an explicit
scientific assertion by the caller, not something inferred from metadata.
Issue-time controls from an `OperatingEnvelopeResult` are replayed, while
scenario-stage controls may adapt to each newly observed held-out scenario.
"""
function evaluate_operating_envelope_coverage(
        scenarios::DOEScenarioSet,
        connection_points::AbstractVector{ConnectionPoint}, capacities;
        roles=(:test,),
        utilizations=:bound_point,
        iid_assumption::Bool=false,
        confidence::Float64=0.95,
        kwargs...)
    0 < confidence < 1 || throw(ArgumentError(
        "confidence must lie strictly between zero and one"))
    selected = select_doe_scenarios(scenarios; roles=roles)
    selected_roles = sort!(unique(
        scenario.role for group in selected.intervals for scenario in group))
    verification = verify_operating_envelope(
        selected, connection_points, capacities;
        utilizations=utilizations,
        replay_control_stages=(:issue,),
        kwargs...)

    rows = NamedTuple[]
    candidate_scenarios = 0
    unresolved_scenarios = 0
    passed_scenarios = 0
    candidate_contexts = 0
    unresolved_contexts = 0
    passed_contexts = 0
    for (interval, group) in enumerate(selected.intervals)
        for (scenario_index, scenario) in enumerate(group)
            contexts = [context for context in
                        verification.context_results[interval]
                        if context.scenario == scenario_index]
            failed_count = count(context -> context.feasible === false, contexts)
            unresolved_count = count(
                context -> context.feasible === nothing, contexts)
            passed_count = count(context -> context.feasible === true, contexts)
            candidate_contexts += failed_count
            unresolved_contexts += unresolved_count
            passed_contexts += passed_count
            status = failed_count > 0 ? :candidate_violation :
                unresolved_count > 0 ? :unresolved : :passed
            candidate_scenarios += status == :candidate_violation
            unresolved_scenarios += status == :unresolved
            passed_scenarios += status == :passed
            normalized_margins = Float64[]
            for context in contexts
                margins = get(context.diagnostics,
                    "minimum_normalized_margins", Dict())
                append!(normalized_margins,
                    Float64(value) for value in values(margins)
                    if value isa Real && isfinite(value))
            end
            push!(rows, (
                dataset_id=selected.dataset_id,
                interval=interval,
                scenario_index=scenario_index,
                scenario_id=scenario.id,
                role=scenario.role,
                weight=scenario.weight,
                timestamp=scenario.timestamp,
                source=scenario.source,
                generation_method=scenario.generation_method,
                utilization_count=length(contexts),
                status=status,
                candidate_context_count=failed_count,
                unresolved_context_count=unresolved_count,
                minimum_normalized_margin=isempty(normalized_margins) ?
                    NaN : minimum(normalized_margins),
            ))
        end
    end

    scenario_count = length(rows)
    context_count = candidate_contexts + unresolved_contexts + passed_contexts
    resolved_scenarios = candidate_scenarios + passed_scenarios
    resolved_contexts = candidate_contexts + passed_contexts
    scenario_frequency = resolved_scenarios == 0 ? NaN :
        candidate_scenarios / resolved_scenarios
    context_frequency = resolved_contexts == 0 ? NaN :
        candidate_contexts / resolved_contexts
    conservative_scenario_frequency = scenario_count == 0 ? NaN :
        (candidate_scenarios + unresolved_scenarios) / scenario_count

    weights_available = all(row.weight !== nothing for row in rows)
    weighted_frequency = missing
    weighted_conservative_frequency = missing
    if weights_available
        interval_estimates = Float64[]
        interval_conservative = Float64[]
        for interval in eachindex(selected.intervals)
            interval_rows = [row for row in rows if row.interval == interval]
            total_weight = sum(Float64(row.weight) for row in interval_rows)
            resolved_weight = sum(Float64(row.weight) for row in interval_rows
                                  if row.status != :unresolved)
            resolved_weight > 0 && push!(interval_estimates,
                sum(Float64(row.weight) for row in interval_rows
                    if row.status == :candidate_violation) / resolved_weight)
            push!(interval_conservative,
                sum(Float64(row.weight) for row in interval_rows
                    if row.status != :passed) / total_weight)
        end
        weighted_frequency = isempty(interval_estimates) ? missing :
            sum(interval_estimates) / length(interval_estimates)
        weighted_conservative_frequency =
            sum(interval_conservative) / length(interval_conservative)
    end

    hoeffding_upper = missing
    if iid_assumption
        radius = sqrt(log(1 / (1 - confidence)) / (2 * scenario_count))
        hoeffding_upper = min(1.0, conservative_scenario_frequency + radius)
    end
    outcome = candidate_scenarios > 0 ? :candidate_violations_observed :
        unresolved_scenarios > 0 ? :inconclusive :
        :no_candidate_violations_observed
    metrics = Dict{String,Any}(
        "scenario_count" => scenario_count,
        "context_count" => context_count,
        "passed_scenario_count" => passed_scenarios,
        "candidate_scenario_count" => candidate_scenarios,
        "unresolved_scenario_count" => unresolved_scenarios,
        "passed_context_count" => passed_contexts,
        "candidate_context_count" => candidate_contexts,
        "unresolved_context_count" => unresolved_contexts,
        "empirical_candidate_scenario_frequency" => scenario_frequency,
        "empirical_candidate_context_frequency" => context_frequency,
        "conservative_scenario_frequency" =>
            conservative_scenario_frequency,
        "weighted_candidate_scenario_frequency" => weighted_frequency,
        "weighted_conservative_scenario_frequency" =>
            weighted_conservative_frequency,
        "confidence" => confidence,
        "one_sided_hoeffding_upper_bound" => hoeffding_upper)
    diagnostics = Dict{String,Any}(
        "dataset_id" => selected.dataset_id,
        "selected_roles" => selected_roles,
        "iid_assumption" => iid_assumption,
        "statistical_bound" => iid_assumption ?
            :one_sided_hoeffding_candidate_frequency : :none,
        "unresolved_treatment" =>
            :counted_as_candidate_in_conservative_rate_and_bound,
        "scenario_weights" => weights_available ? :declared_relative_weights :
            :not_available,
        "issue_controls" => :replayed_when_available,
        "scenario_controls" => :adaptive_to_selected_scenario,
        "distribution_shift_assessed" => false,
        "global_certificate" => false,
        "claim" => iid_assumption ?
            :empirical_coverage_with_declared_iid_concentration_bound :
            :empirical_coverage_only)
    return DOECoverageResult(
        outcome, selected_roles, verification, rows, metrics, diagnostics)
end

function _doe_scaled_capacities(result::OperatingEnvelopeResult,
                                scale::Float64)
    envelope = Dict{String,Vector{Float64}}(
        id => scale .* values for (id, values) in result.envelope)
    return OperatingEnvelopeResult(
        copy(result.termination_status), envelope,
        scale .* result.total_export, deepcopy(result.snapshots),
        result.direction, scale .* result.total_capacity,
        deepcopy(result.diagnostics), deepcopy(result.fairness_metrics),
        deepcopy(result.schedule))
end

function _doe_scaled_capacities(capacities::AbstractDict, scale::Float64)
    scaled = Dict{String,Any}()
    for (id, value) in pairs(capacities)
        scaled[string(id)] = value isa AbstractVector ?
            scale .* Float64.(value) : scale * Float64(value)
    end
    return scaled
end

"""
    evaluate_operating_envelope_coverage_curve(
        scenarios, connection_points, capacities;
        scales=(0.5, 0.75, 1.0), kwargs...)

Evaluate the same selected held-out scenarios over an ascending sequence of
capacity multipliers. When `capacities` is an [`OperatingEnvelopeResult`](@ref),
its recorded issue-time controls are retained at every scale; they are not
silently re-optimized after test scenarios are observed. Other coverage
keywords, including `roles`, `utilizations`, and `iid_assumption`, are forwarded
to [`evaluate_operating_envelope_coverage`](@ref).

The first scale with a candidate violation is a finite-grid diagnostic only.
Candidate-count reversals are reported rather than hidden: they can reflect
physical non-monotonicity, discrete/control effects, or local-solver behavior.
"""
function evaluate_operating_envelope_coverage_curve(
        scenarios::DOEScenarioSet,
        connection_points::AbstractVector{ConnectionPoint}, capacities;
        scales=(0.5, 0.75, 1.0),
        kwargs...)
    capacities isa OperatingEnvelopeResult || capacities isa AbstractDict ||
        throw(ArgumentError(
            "capacities must be an OperatingEnvelopeResult or dictionary"))
    scale_values = sort!(unique(Float64.(collect(scales))))
    isempty(scale_values) && throw(ArgumentError("scales must be non-empty"))
    all(scale -> isfinite(scale) && scale >= 0, scale_values) ||
        throw(ArgumentError("scales must be finite and non-negative"))

    coverages = DOECoverageResult[]
    rows = NamedTuple[]
    for scale in scale_values
        coverage = evaluate_operating_envelope_coverage(
            scenarios, connection_points,
            _doe_scaled_capacities(capacities, scale); kwargs...)
        push!(coverages, coverage)
        push!(rows, (
            scale=scale,
            total_capacity_W=capacities isa OperatingEnvelopeResult ?
                scale .* capacities.total_capacity : missing,
            outcome=coverage.outcome,
            scenario_count=coverage.metrics["scenario_count"],
            candidate_scenario_count=
                coverage.metrics["candidate_scenario_count"],
            unresolved_scenario_count=
                coverage.metrics["unresolved_scenario_count"],
            empirical_candidate_scenario_frequency=
                coverage.metrics["empirical_candidate_scenario_frequency"],
            empirical_candidate_context_frequency=
                coverage.metrics["empirical_candidate_context_frequency"],
            conservative_scenario_frequency=
                coverage.metrics["conservative_scenario_frequency"],
            weighted_candidate_scenario_frequency=
                coverage.metrics["weighted_candidate_scenario_frequency"],
            weighted_conservative_scenario_frequency=
                coverage.metrics["weighted_conservative_scenario_frequency"],
            one_sided_hoeffding_upper_bound=
                coverage.metrics["one_sided_hoeffding_upper_bound"],
        ))
    end

    first_candidate = findfirst(
        row -> row.candidate_scenario_count > 0, rows)
    candidate_count_reversals = NamedTuple[]
    for index in 2:length(rows)
        rows[index].candidate_scenario_count >=
            rows[index - 1].candidate_scenario_count && continue
        push!(candidate_count_reversals, (
            lower_scale=rows[index - 1].scale,
            higher_scale=rows[index].scale,
            lower_candidate_count=rows[index - 1].candidate_scenario_count,
            higher_candidate_count=rows[index].candidate_scenario_count,
        ))
    end
    diagnostics = Dict{String,Any}(
        "first_candidate_scale" => first_candidate === nothing ? missing :
            rows[first_candidate].scale,
        "scale_order" => :ascending,
        "capacity_scaling" => :uniform,
        "issue_control_treatment" => capacities isa OperatingEnvelopeResult ?
            :retained_from_issued_result : :not_available,
        "control_reoptimization" => :none,
        "candidate_count_reversals" => candidate_count_reversals,
        "candidate_count_monotonic_nondecreasing" =>
            isempty(candidate_count_reversals),
        "continuous_threshold_estimated" => false,
        "global_certificate" => false,
        "claim" => :finite_capacity_sensitivity_curve)
    return DOECoverageCurveResult(
        scale_values, coverages, rows, diagnostics)
end

function _doe_coverage_metric_delta(reference::DOECoverageResult,
                                    shifted::DOECoverageResult,
                                    key::String)
    left = get(reference.metrics, key, missing)
    right = get(shifted.metrics, key, missing)
    if left isa Real && right isa Real && isfinite(left) && isfinite(right)
        return Float64(right - left)
    end
    return missing
end

"""
    compare_doe_coverage_shift(reference, shifted;
        reference_label="reference", shifted_label="shifted", tolerance=0.0)

Compare two already evaluated held-out ensembles descriptively. Deltas are
`shifted - reference`; the outcome uses the conservative scenario-frequency
delta so unresolved cases are not treated as safe. This function deliberately
does not infer that the scenario distributions differ. Establishing covariate,
concept, or temporal distribution shift requires a separately declared test
and assumptions appropriate to the data-generating process.
"""
function compare_doe_coverage_shift(
        reference::DOECoverageResult,
        shifted::DOECoverageResult;
        reference_label::AbstractString="reference",
        shifted_label::AbstractString="shifted",
        tolerance::Float64=0.0)
    isempty(reference_label) && throw(ArgumentError(
        "reference_label must be non-empty"))
    isempty(shifted_label) && throw(ArgumentError(
        "shifted_label must be non-empty"))
    isfinite(tolerance) && tolerance >= 0 || throw(ArgumentError(
        "tolerance must be finite and non-negative"))
    metric_keys = (
        "empirical_candidate_scenario_frequency",
        "empirical_candidate_context_frequency",
        "conservative_scenario_frequency",
        "weighted_candidate_scenario_frequency",
        "weighted_conservative_scenario_frequency",
        "one_sided_hoeffding_upper_bound",
    )
    deltas = Dict{String,Any}(key =>
        _doe_coverage_metric_delta(reference, shifted, key)
        for key in metric_keys)
    primary_delta = deltas["conservative_scenario_frequency"]
    outcome = primary_delta === missing ? :inconclusive :
        primary_delta > tolerance ? :higher_candidate_frequency :
        primary_delta < -tolerance ? :lower_candidate_frequency :
        :no_observed_change
    diagnostics = Dict{String,Any}(
        "reference_label" => String(reference_label),
        "shifted_label" => String(shifted_label),
        "reference_dataset_id" =>
            get(reference.diagnostics, "dataset_id", "unknown"),
        "shifted_dataset_id" =>
            get(shifted.diagnostics, "dataset_id", "unknown"),
        "reference_scenario_count" =>
            get(reference.metrics, "scenario_count", missing),
        "shifted_scenario_count" =>
            get(shifted.metrics, "scenario_count", missing),
        "comparison_metric" => :conservative_scenario_frequency,
        "tolerance" => tolerance,
        "performance_shift_observed" => outcome in
            (:higher_candidate_frequency, :lower_candidate_frequency),
        "distribution_shift_test" => :not_performed,
        "distribution_shift_detected" => false,
        "global_certificate" => false,
        "claim" => :descriptive_held_out_performance_comparison)
    return DOECoverageShiftResult(
        outcome, reference, shifted, deltas, diagnostics)
end

function _doe_first_primes(count::Int)
    primes = Int[]
    candidate = 2
    while length(primes) < count
        isprime = all(candidate % divisor != 0 for divisor in 2:isqrt(candidate))
        isprime && push!(primes, candidate)
        candidate += 1
    end
    return primes
end

function _doe_radical_inverse(index::Int, base::Int)
    value = 0.0
    factor = 1.0 / base
    remaining = index
    while remaining > 0
        digit = remaining % base
        value += digit * factor
        remaining ÷= base
        factor /= base
    end
    return value
end

function _doe_search_points(participant_count::Int, samples::Int,
                            sequence_offset::Int;
                            include_zero, include_bound, include_corners,
                            max_exact_corners)
    participant_count >= 1 || throw(ArgumentError(
        "utilization search needs at least one participant"))
    samples >= 0 || throw(ArgumentError("samples must be non-negative"))
    sequence_offset >= 0 || throw(ArgumentError(
        "sequence_offset must be non-negative"))
    include_corners && participant_count > max_exact_corners &&
        throw(ArgumentError(
            "include_corners=true requires participant count <= max_exact_corners"))
    points = Vector{Vector{Float64}}()
    include_zero && push!(points, zeros(participant_count))
    include_bound && push!(points, ones(participant_count))
    include_corners && append!(points,
        _dispatch_patterns(participant_count, :corners))
    bases = _doe_first_primes(participant_count)
    for index in (sequence_offset + 1):(sequence_offset + samples)
        push!(points, [_doe_radical_inverse(index, base) for base in bases])
    end
    unique_points = Vector{Vector{Float64}}()
    seen = Set{Tuple}()
    for point in points
        identity = Tuple(point)
        identity in seen && continue
        push!(seen, identity)
        push!(unique_points, point)
    end
    return unique_points
end

"""
    search_operating_envelope_utilizations(nets, connection_points, capacities;
        samples=16, sequence_offset=0, include_zero=true,
        include_bound=true, include_corners=false, kwargs...)

Screen an issued operating envelope at a deterministic Halton low-discrepancy
set of continuous participant-utilisation vectors. The generated set is solved
jointly, so the selected [`DOEControlPolicy`](@ref) and its non-anticipativity
constraints apply across the search points.

The result is deliberately labelled `:search_stable`,
`:candidate_counterexample`, or `:inconclusive`; this finite local-NLP search is
not a robust-feasibility certificate. `sequence_offset` makes disjoint,
reproducible batches possible. Remaining keywords are forwarded to
[`verify_operating_envelope`](@ref).
"""
function search_operating_envelope_utilizations(
        nets, connection_points::AbstractVector{ConnectionPoint}, capacities;
        samples::Int=16,
        sequence_offset::Int=0,
        include_zero::Bool=true,
        include_bound::Bool=true,
        include_corners::Bool=false,
        max_exact_corners::Int=10,
        independent_replay::Bool=false,
        kwargs...)
    points = _doe_search_points(length(connection_points), samples,
        sequence_offset; include_zero, include_bound, include_corners,
        max_exact_corners)
    isempty(points) && throw(ArgumentError(
        "utilization search generated no points"))
    verification = verify_operating_envelope(
        nets, connection_points, capacities;
        utilizations=points,
        independent_replay=independent_replay,
        max_exact_corners=max_exact_corners,
        kwargs...)
    candidates = OperatingEnvelopeContextResult[]
    unresolved = false
    for interval_contexts in verification.context_results
        for context in interval_contexts
            context.feasible === false && push!(candidates, context)
            context.feasible === nothing && (unresolved = true)
        end
    end
    outcome = if all(verification.feasible)
        :search_stable
    elseif !isempty(candidates)
        :candidate_counterexample
    else
        :inconclusive
    end
    diagnostics = Dict{String,Any}(
        "method" => :halton_low_discrepancy_screening,
        "requested_samples" => samples,
        "sequence_offset" => sequence_offset,
        "tested_point_count" => length(points),
        "include_zero" => include_zero,
        "include_bound" => include_bound,
        "include_corners" => include_corners,
        "candidate_context_count" => length(candidates),
        "unresolved_contexts_present" => unresolved,
        "global_certificate" => false,
        "claim" => outcome == :search_stable ?
            :local_ac_feasibility_at_generated_points :
            outcome == :candidate_counterexample ?
                :candidate_violation_requires_confirmation :
                :numerically_or_policy_inconclusive)
    return OperatingEnvelopeSearchResult(
        outcome, points, verification, candidates, diagnostics)
end

function _doe_context_adversarial_score(context::OperatingEnvelopeContextResult)
    context.feasible === false && return Inf
    context.feasible === nothing && return NaN
    margins = get(context.diagnostics, "minimum_normalized_margins", nothing)
    margins isa AbstractDict && !isempty(margins) || return NaN
    values_ = Float64[value for value in values(margins) if value isa Real &&
                       isfinite(value)]
    isempty(values_) && return NaN
    # Higher is worse: zero is a binding declared limit and positive values
    # represent negative normalized headroom within solver tolerances.
    return -minimum(values_)
end

function _doe_adversarial_point_scores(points, verification)
    point_index = Dict(Tuple(point) => index for (index, point) in enumerate(points))
    scores = fill(NaN, length(points))
    worst_contexts = Vector{Union{Nothing,OperatingEnvelopeContextResult}}(
        nothing, length(points))
    worst_intervals = Vector{Union{Nothing,Int}}(nothing, length(points))
    for (interval, contexts) in enumerate(verification.context_results)
        for context in contexts
            index = get(point_index, Tuple(context.utilization), nothing)
            index === nothing && continue
            score = _doe_context_adversarial_score(context)
            if !isnan(score) && (isnan(scores[index]) || score > scores[index])
                scores[index] = score
                worst_contexts[index] = context
                worst_intervals[index] = interval
            end
        end
    end
    return scores, worst_contexts, worst_intervals
end

function _doe_coordinate_refinement(points, scores, restarts, step)
    ranked = [index for index in eachindex(points) if isfinite(scores[index])]
    sort!(ranked; by=index -> (-scores[index], Tuple(points[index])))
    resize!(ranked, min(restarts, length(ranked)))
    seen = Set(Tuple(point) for point in points)
    additions = Vector{Vector{Float64}}()
    for index in ranked
        for coordinate in eachindex(points[index]), direction in (-1.0, 1.0)
            candidate = copy(points[index])
            candidate[coordinate] = clamp(
                candidate[coordinate] + direction * step, 0.0, 1.0)
            key = Tuple(candidate)
            key in seen && continue
            push!(seen, key)
            push!(additions, candidate)
        end
    end
    return additions, ranked
end

"""
    search_operating_envelope_adversarial(
        nets, connection_points, capacities;
        seed_samples=8, refinement_rounds=3, restarts=3,
        initial_step=0.25, step_decay=0.5, kwargs...)

Search for unsafe partial utilization using deterministic Halton seeds followed
by coordinate refinement around points with the least normalized voltage,
thermal, or negative-sequence headroom. Every round jointly verifies the full
accumulated point set, preserving the selected control-recourse policy and its
non-anticipativity constraints.

The score is the negative of the smallest declared normalized constraint
margin, maximized over intervals and scenarios for each utilization point.
An individually infeasible context receives an infinite candidate score;
unresolved contexts remain unscored and make the outcome inconclusive. Because
each AC verification and the outer point search are local and finite, this
function is a falsification heuristic, not the continuous
violation-maximizing oracle or a robust certificate.
"""
function search_operating_envelope_adversarial(
        nets, connection_points::AbstractVector{ConnectionPoint}, capacities;
        seed_samples::Int=8,
        sequence_offset::Int=0,
        include_zero::Bool=true,
        include_bound::Bool=true,
        include_corners::Bool=false,
        max_exact_corners::Int=10,
        refinement_rounds::Int=3,
        restarts::Int=3,
        initial_step::Float64=0.25,
        step_decay::Float64=0.5,
        minimum_step::Float64=1e-3,
        independent_replay::Bool=false,
        kwargs...)
    seed_samples >= 0 || throw(ArgumentError(
        "seed_samples must be non-negative"))
    refinement_rounds >= 0 || throw(ArgumentError(
        "refinement_rounds must be non-negative"))
    restarts >= 1 || throw(ArgumentError("restarts must be positive"))
    isfinite(initial_step) && initial_step > 0 || throw(ArgumentError(
        "initial_step must be finite and positive"))
    isfinite(step_decay) && 0 < step_decay <= 1 || throw(ArgumentError(
        "step_decay must be finite and lie in (0, 1]"))
    isfinite(minimum_step) && minimum_step > 0 || throw(ArgumentError(
        "minimum_step must be finite and positive"))

    points = _doe_search_points(length(connection_points), seed_samples,
        sequence_offset; include_zero, include_bound, include_corners,
        max_exact_corners)
    isempty(points) && throw(ArgumentError(
        "adversarial search generated no seed points"))
    verifications = OperatingEnvelopeVerification[]
    candidates = OperatingEnvelopeContextResult[]
    candidate_locations = Dict{String,Any}[]
    scores = fill(NaN, length(points))
    worst_contexts = Vector{Union{Nothing,OperatingEnvelopeContextResult}}(
        nothing, length(points))
    worst_intervals = Vector{Union{Nothing,Int}}(nothing, length(points))
    refinement_sources = Vector{Vector{Int}}()
    step = initial_step
    outcome = :search_stable

    for round in 0:refinement_rounds
        verification = verify_operating_envelope(
            nets, connection_points, capacities;
            utilizations=points,
            independent_replay=independent_replay,
            max_exact_corners=max_exact_corners,
            kwargs...)
        push!(verifications, verification)
        scores, worst_contexts, worst_intervals = _doe_adversarial_point_scores(
            points, verification)
        empty!(candidates)
        empty!(candidate_locations)
        for (interval, contexts) in enumerate(verification.context_results),
            context in contexts
            if context.feasible === false
                push!(candidates, context)
                push!(candidate_locations, Dict{String,Any}(
                    "interval" => interval,
                    "scenario" => context.scenario,
                    "utilization_index" => context.utilization_index,
                    "utilization" => copy(context.utilization),
                    "termination_status" => context.termination_status,
                    "primal_status" => context.primal_status))
            end
        end
        if !isempty(candidates)
            outcome = :candidate_counterexample
            break
        elseif !all(verification.feasible) || all(isnan, scores)
            outcome = :inconclusive
            break
        elseif round == refinement_rounds
            outcome = :search_stable
            break
        end

        additions, ranked = _doe_coordinate_refinement(
            points, scores, restarts, max(step, minimum_step))
        push!(refinement_sources, ranked)
        if isempty(additions)
            outcome = :search_stable
            break
        end
        append!(points, additions)
        append!(scores, fill(NaN, length(additions)))
        append!(worst_contexts,
                fill(nothing, length(additions)))
        append!(worst_intervals, fill(nothing, length(additions)))
        step *= step_decay
    end

    scored = [index for index in eachindex(scores) if !isnan(scores[index])]
    worst_index = isempty(scored) ? nothing :
        scored[argmax(scores[scored])]
    worst_context = worst_index === nothing ? nothing :
        worst_contexts[worst_index]
    worst_interval = worst_index === nothing ? nothing :
        worst_intervals[worst_index]
    diagnostics = Dict{String,Any}(
        "method" => :adaptive_coordinate_normalized_margin_search,
        "seed_method" => :halton_low_discrepancy,
        "seed_samples" => seed_samples,
        "sequence_offset" => sequence_offset,
        "requested_refinement_rounds" => refinement_rounds,
        "completed_verification_rounds" => length(verifications),
        "restarts" => restarts,
        "initial_step" => initial_step,
        "step_decay" => step_decay,
        "minimum_step" => minimum_step,
        "refinement_source_indices" => refinement_sources,
        "tested_point_count" => length(points),
        "round_tested_point_counts" => [
            get(verification.diagnostics[1], "dispatch_points_per_scenario", 0)
            for verification in verifications],
        "candidate_context_count" => length(candidates),
        "candidate_locations" => candidate_locations,
        "worst_point_index" => worst_index,
        "worst_interval" => worst_interval,
        "worst_score" => worst_index === nothing ? NaN : scores[worst_index],
        "score_definition" =>
            :negative_minimum_declared_normalized_constraint_margin,
        "global_certificate" => false,
        "claim" => outcome == :candidate_counterexample ?
            :candidate_violation_requires_confirmation :
            outcome == :search_stable ?
                :no_counterexample_found_within_recorded_search_budget :
                :numerically_or_policy_inconclusive)
    return DOEAdversarialSearchResult(
        outcome, points, scores, verifications, copy(candidates),
        worst_interval, worst_context, diagnostics)
end

"""
    confirm_operating_envelope_counterexample(
        nets, connection_points, capacities, utilization;
        start_scales=(1.0, 0.9, 1.1), kwargs...)

Repeat one candidate utilization point from deterministic perturbations of the
registered OPF starting values. If `capacities` is an
[`OperatingEnvelopeResult`](@ref), recorded issue-time/scenario controls are
replayed by default through [`verify_operating_envelope`](@ref).

`:repeated_candidate` means every start produced at least one individually
infeasible context; `:not_reproduced` means at least one start verified every
interval/scenario; all other combinations are `:inconclusive`. None is a global
infeasibility certificate.
"""
function confirm_operating_envelope_counterexample(
        nets, connection_points::AbstractVector{ConnectionPoint}, capacities,
        utilization;
        start_scales=(1.0, 0.9, 1.1),
        base_start_hook! = nothing,
        independent_replay::Bool=false,
        kwargs...)
    point = only(_verification_patterns(
        [utilization], length(connection_points)))
    scales = Float64.(collect(start_scales))
    isempty(scales) && throw(ArgumentError("start_scales must be non-empty"))
    all(scale -> isfinite(scale) && scale > 0, scales) ||
        throw(ArgumentError("start_scales must be finite and positive"))
    verifications = OperatingEnvelopeVerification[]
    changed_counts = Int[]
    run_candidates = Vector{Vector{Dict{String,Any}}}()
    for scale in scales
        changed = Ref(0)
        start_hook! = contexts -> begin
            base_start_hook! === nothing || base_start_hook!(contexts)
            changed[] = _doe_scale_registered_starts!(contexts, scale)
        end
        verification = verify_operating_envelope(
            nets, connection_points, capacities;
            utilizations=[point],
            start_hook! = start_hook!,
            independent_replay=independent_replay,
            kwargs...)
        push!(verifications, verification)
        push!(changed_counts, changed[])
        locations = Dict{String,Any}[]
        for (interval, contexts) in enumerate(verification.context_results),
            context in contexts
            context.feasible === false || continue
            push!(locations, Dict{String,Any}(
                "interval" => interval,
                "scenario" => context.scenario,
                "termination_status" => context.termination_status,
                "primal_status" => context.primal_status))
        end
        push!(run_candidates, locations)
    end
    reproduced = [!isempty(locations) for locations in run_candidates]
    fully_feasible = [all(verification.feasible)
                      for verification in verifications]
    outcome = any(fully_feasible) ? :not_reproduced :
        all(reproduced) ? :repeated_candidate : :inconclusive
    diagnostics = Dict{String,Any}(
        "method" => :deterministic_multistart_candidate_replay,
        "start_changed_variable_counts" => changed_counts,
        "termination_statuses" =>
            [verification.termination_status for verification in verifications],
        "fully_feasible_runs" => fully_feasible,
        "candidate_reproduced_runs" => reproduced,
        "candidate_locations_by_run" => run_candidates,
        "global_certificate" => false,
        "claim" => outcome == :repeated_candidate ?
            :candidate_repeated_across_recorded_starts :
            outcome == :not_reproduced ?
                :candidate_not_reproduced_from_every_start :
                :candidate_confirmation_inconclusive)
    return DOECounterexampleConfirmationResult(
        outcome, point, verifications, scales, diagnostics)
end

function _doe_append_unique_points!(points, additions)
    seen = Set(Tuple(point) for point in points)
    added = 0
    for point in additions
        values_ = Float64.(point)
        key = Tuple(values_)
        key in seen && continue
        push!(seen, key)
        push!(points, values_)
        added += 1
    end
    return added
end

"""
    solve_adversarial_search_stable_operating_envelope(
        nets, connection_points; initial_utilizations=:bound_point,
        max_rounds=3, solve_keywords=NamedTuple(),
        search_keywords=NamedTuple())

Allocate on an accumulated finite utilization set, run
[`search_operating_envelope_adversarial`](@ref) against the actual recorded
issue/scenario controls, add all discovered points, and repeat. Every allocation
and search is retained. A final reallocation is performed when the search budget
is exhausted so the returned envelope represents every accumulated point,
although that final allocation has not itself passed another adaptive screen.

Outcomes are `:search_stable`, `:budget_exhausted`, `:inconclusive`, or
`:allocation_failed`. They remain finite local-NLP evidence and never a global
certificate.
"""
function solve_adversarial_search_stable_operating_envelope(
        nets, connection_points::AbstractVector{ConnectionPoint};
        initial_utilizations=:bound_point,
        max_rounds::Int=3,
        control_policy::Union{Nothing,DOEControlPolicy}=nothing,
        context_hook! = nothing,
        solve_keywords=NamedTuple(),
        search_keywords=NamedTuple())
    max_rounds >= 1 || throw(ArgumentError("max_rounds must be positive"))
    solve_keywords isa NamedTuple || throw(ArgumentError(
        "solve_keywords must be a NamedTuple"))
    search_keywords isa NamedTuple || throw(ArgumentError(
        "search_keywords must be a NamedTuple"))
    points = _verification_patterns(
        initial_utilizations, length(connection_points))
    allocations = OperatingEnvelopeResult[]
    searches = DOEAdversarialSearchResult[]
    added_by_round = Int[]

    for round in 1:max_rounds
        allocation = solve_operating_envelope(
            nets, connection_points;
            utilizations=points,
            control_policy=control_policy,
            context_hook! = context_hook!,
            solve_keywords...)
        push!(allocations, allocation)
        if any(!isfinite, allocation.total_capacity)
            return AdversarialSearchStableOperatingEnvelopeResult(
                :allocation_failed, allocation, allocations, searches, points,
                length(searches), Dict{String,Any}(
                    "failure_round" => round,
                    "global_certificate" => false,
                    "claim" => :no_publishable_allocation))
        end
        search = search_operating_envelope_adversarial(
            nets, connection_points, allocation;
            control_policy=control_policy,
            context_hook! = context_hook!,
            search_keywords...)
        push!(searches, search)
        added = _doe_append_unique_points!(
            points, search.utilization_points)
        push!(added_by_round, added)
        if search.outcome == :search_stable
            return AdversarialSearchStableOperatingEnvelopeResult(
                :search_stable, allocation, allocations, searches, points,
                round, Dict{String,Any}(
                    "points_added_by_round" => added_by_round,
                    "final_allocation_screened" => true,
                    "global_certificate" => false,
                    "claim" =>
                        :no_counterexample_found_within_recorded_search_budget))
        elseif search.outcome == :inconclusive && added == 0
            return AdversarialSearchStableOperatingEnvelopeResult(
                :inconclusive, allocation, allocations, searches, points,
                round, Dict{String,Any}(
                    "points_added_by_round" => added_by_round,
                    "final_allocation_screened" => false,
                    "global_certificate" => false,
                    "claim" => :numerically_or_policy_inconclusive))
        end
    end

    final_allocation = solve_operating_envelope(
        nets, connection_points;
        utilizations=points,
        control_policy=control_policy,
        context_hook! = context_hook!,
        solve_keywords...)
    push!(allocations, final_allocation)
    outcome = any(!isfinite, final_allocation.total_capacity) ?
        :allocation_failed : :budget_exhausted
    return AdversarialSearchStableOperatingEnvelopeResult(
        outcome, final_allocation, allocations, searches, points,
        length(searches), Dict{String,Any}(
            "points_added_by_round" => added_by_round,
            "final_allocation_screened" => false,
            "global_certificate" => false,
            "claim" => outcome == :allocation_failed ?
                :no_publishable_final_allocation :
                :finite_search_budget_exhausted))
end

"""
    solve_search_stable_operating_envelope(nets, connection_points;
        initial_utilizations=:bound_point, samples_per_round=16,
        max_rounds=3, solve_keywords=NamedTuple(),
        verification_keywords=NamedTuple(), kwargs...)

Counterexample-guided research wrapper for a box DOE. Each round allocates
capacity on the accumulated utilization set, screens that allocation on an
increasing deterministic low-discrepancy set, and adds the full screened set to
the next allocation when a candidate violation or joint-policy conflict is
found.

The result is `:search_stable`, `:budget_exhausted`, or `:allocation_failed`.
No outcome is a global certificate: this procedure improves finite tested-point
coverage while retaining local nonlinear-solver limitations. Put fairness and
other allocation-only options in `solve_keywords`; put replay/solver options
used only for screening in `verification_keywords`.
"""
function solve_search_stable_operating_envelope(
        nets, connection_points::AbstractVector{ConnectionPoint};
        initial_utilizations=:bound_point,
        samples_per_round::Int=16,
        max_rounds::Int=3,
        sequence_offset::Int=0,
        control_policy::Union{Nothing,DOEControlPolicy}=nothing,
        context_hook! = nothing,
        solve_keywords=NamedTuple(),
        verification_keywords=NamedTuple())
    samples_per_round >= 1 || throw(ArgumentError(
        "samples_per_round must be positive"))
    max_rounds >= 1 || throw(ArgumentError("max_rounds must be positive"))
    sequence_offset >= 0 || throw(ArgumentError(
        "sequence_offset must be non-negative"))
    solve_keywords isa NamedTuple || throw(ArgumentError(
        "solve_keywords must be a NamedTuple"))
    verification_keywords isa NamedTuple || throw(ArgumentError(
        "verification_keywords must be a NamedTuple"))
    points = _verification_patterns(
        initial_utilizations, length(connection_points))
    searches = OperatingEnvelopeSearchResult[]
    allocation = solve_operating_envelope(
        nets, connection_points;
        utilizations=points,
        control_policy=control_policy,
        context_hook! = context_hook!,
        solve_keywords...)
    for round in 1:max_rounds
        if any(!isfinite, allocation.total_capacity)
            return SearchStableOperatingEnvelopeResult(
                :allocation_failed, allocation, searches, points,
                length(searches), Dict{String,Any}(
                    "global_certificate" => false,
                    "failure_round" => round,
                    "claim" => :no_publishable_allocation))
        end
        search = search_operating_envelope_utilizations(
            nets, connection_points, allocation;
            samples=round * samples_per_round,
            sequence_offset=sequence_offset,
            control_policy=control_policy,
            context_hook! = context_hook!,
            independent_replay=false,
            verification_keywords...)
        push!(searches, search)
        points = search.utilization_points
        if search.outcome == :search_stable
            return SearchStableOperatingEnvelopeResult(
                :search_stable, allocation, searches, points, round,
                Dict{String,Any}(
                    "global_certificate" => false,
                    "claim" => :local_ac_feasibility_at_generated_points,
                    "samples_per_round" => samples_per_round,
                    "sequence_offset" => sequence_offset))
        end
        allocation = solve_operating_envelope(
            nets, connection_points;
            utilizations=points,
            control_policy=control_policy,
            context_hook! = context_hook!,
            solve_keywords...)
    end
    return SearchStableOperatingEnvelopeResult(
        :budget_exhausted, allocation, searches, points,
        length(searches), Dict{String,Any}(
            "global_certificate" => false,
            "claim" => :finite_search_budget_exhausted,
            "samples_per_round" => samples_per_round,
            "sequence_offset" => sequence_offset))
end

"""
    doe_benchmark_rows(spec, result; method_label="DOE")

Return stable, tidy interval rows for a DOE allocation. The rows retain the
study identity, formulation labels, solver evidence, capacity, and claim
boundary without embedding full network snapshots.
"""
function doe_benchmark_rows(spec::DOEStudySpec,
                            result::OperatingEnvelopeResult;
                            method_label::AbstractString="DOE")
    interval_count = length(result.total_capacity)
    length(spec.network_hashes) == interval_count || throw(ArgumentError(
        "DOEStudySpec interval count does not match the result"))
    return [(
        study_id=spec.study_id,
        method=String(method_label),
        interval=interval,
        termination_status=result.termination_status[interval],
        primal_status=string(get(result.diagnostics[interval],
                                 "primal_status", "UNKNOWN")),
        feasible=Bool(get(result.diagnostics[interval], "feasible", false)),
        published=isfinite(result.total_capacity[interval]),
        direction=result.direction,
        total_capacity_W=result.total_capacity[interval],
        fairness_kind=get(result.diagnostics[interval], "fairness_kind", missing),
        control_policy=get(result.diagnostics[interval], "control_policy", missing),
        control_policy_signature=get(
            result.diagnostics[interval], "control_policy_signature", missing),
        issued_control_count=length(get(
            result.diagnostics[interval], "issued_control_values", Any[])),
        scenario_dataset_id=get(
            result.diagnostics[interval], "scenario_dataset_id", missing),
        security_scope=get(result.diagnostics[interval], "security_scope", missing),
        scenario_count=Int(get(result.diagnostics[interval], "scenario_count", 0)),
        dispatch_point_count=Int(get(result.diagnostics[interval],
                                     "dispatch_points_per_scenario", 0)),
        solve_time_seconds=Float64(get(result.diagnostics[interval],
                                       "solve_time_seconds", NaN)),
        maximum_primal_constraint_violation=Float64(get(
            result.diagnostics[interval],
            "maximum_primal_constraint_violation", NaN)),
        global_certificate=Bool(get(result.diagnostics[interval],
                                    "global_certificate", false)),
    ) for interval in 1:interval_count]
end

"""
    doe_context_benchmark_rows(spec, verification; method_label="verification")

Return one tidy row per interval/scenario/utilization context, including replay
status and minimum-margin summaries.
"""
function doe_context_benchmark_rows(
        spec::DOEStudySpec, verification::OperatingEnvelopeVerification;
        method_label::AbstractString="verification")
    length(spec.network_hashes) == length(verification.context_results) ||
        throw(ArgumentError(
            "DOEStudySpec interval count does not match the verification"))
    rows = NamedTuple[]
    for (interval, contexts) in enumerate(verification.context_results)
        for context in contexts
            replay = get(context.diagnostics, "independent_replay", nothing)
            margins = get(context.diagnostics, "minimum_margins",
                          Dict{String,Float64}())
            normalized_margins = get(
                context.diagnostics, "minimum_normalized_margins",
                Dict{String,Float64}())
            push!(rows, (
                study_id=spec.study_id,
                method=String(method_label),
                interval=interval,
                scenario=context.scenario,
                utilization_index=context.utilization_index,
                utilization=copy(context.utilization),
                termination_status=context.termination_status,
                primal_status=context.primal_status,
                feasible=context.feasible,
                replay_feasible=replay isa AbstractDict ?
                    get(replay, "feasible", missing) : missing,
                control_replay_complete=replay isa AbstractDict ?
                    get(replay, "control_replay_complete", missing) : missing,
                replay_voltage_difference_V=replay isa AbstractDict ?
                    get(replay, "maximum_voltage_difference_V", NaN) : NaN,
                minimum_margins=copy(margins),
                minimum_normalized_margins=copy(normalized_margins),
                issued_control_replay_source=get(
                    verification.diagnostics[interval],
                    "issued_control_replay_source", missing),
                issued_control_replay_count=get(
                    verification.diagnostics[interval],
                    "issued_control_replay_count", 0),
                global_certificate=false,
            ))
        end
    end
    return rows
end
