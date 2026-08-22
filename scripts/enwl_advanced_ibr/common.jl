module ENWLAdvancedIBR

using PowerOptLab
using BMOPFTools: parse_bmopf

export default_data_root, parse_options, snapshot_paths, snapshot_id,
       feeder_id, inventory_rows, build_study_cases, write_csv,
       print_inventory, print_campaign_plan

const DEFAULT_VARIANTS = ["baseline", "worst_phase", "sequence_droop"]
const PHASES = ["a", "b", "c"]

"Default to a sibling BMOPFDraftData checkout, without depending on the CWD."
default_data_root() = normpath(joinpath(
    @__DIR__, "..", "..", "..", "BMOPFDraftData", "benchmarks",
    "ENWLsnapshots"))

function parse_options(args; defaults=Dict{String,String}())
    options = copy(defaults)
    for argument in args
        startswith(argument, "--") || throw(ArgumentError(
            "expected --name=value, got '$argument'"))
        argument == "--help" && (options["help"] = "true"; continue)
        fields = split(argument[3:end], "="; limit=2)
        length(fields) == 2 || throw(ArgumentError(
            "expected --name=value, got '$argument'"))
        isempty(fields[1]) && throw(ArgumentError("option name cannot be empty"))
        options[fields[1]] = fields[2]
    end
    options
end

function snapshot_paths(root::AbstractString; feeder::Union{Nothing,String}=nothing,
                        snapshots::Union{Nothing,Vector{String}}=nothing)
    isdir(root) || throw(ArgumentError("ENWL snapshot root does not exist: $root"))
    feeders = feeder === nothing ?
        sort!(filter(isdir, readdir(root; join=true))) : [joinpath(root, feeder)]
    all(isdir, feeders) || throw(ArgumentError(
        "requested feeder directory does not exist under $root"))
    wanted = snapshots === nothing ? nothing : Set(snapshots)
    matched = Set{String}()
    paths = String[]
    for directory in feeders
        for path in sort!(readdir(directory; join=true))
            endswith(path, ".bmopf.json") || continue
            sid = snapshot_id(path)
            wanted === nothing && (push!(paths, path); continue)
            suffix = replace(sid, feeder_id(path) * "_" => "")
            keys_here = sid == suffix ? (sid,) : (sid, suffix)
            hit = false
            for key in keys_here
                key in wanted || continue
                push!(matched, key)
                hit = true
            end
            hit && push!(paths, path)
        end
    end
    if wanted !== nothing
        # Dropping a typo silently would yield a short campaign that still
        # looks complete once the CSVs are written.
        unmatched = sort!(collect(setdiff(wanted, matched)))
        isempty(unmatched) || throw(ArgumentError(
            "no snapshot matches: $(join(unmatched, ", "))"))
    end
    isempty(paths) && throw(ArgumentError("no matching .bmopf.json snapshots found"))
    paths
end

snapshot_id(path::AbstractString) = replace(basename(path), r"\.bmopf\.json$" => "")
feeder_id(path::AbstractString) = basename(dirname(path))

_number(value) = Float64(value)
_scalar(value) = value isa AbstractVector ? _number(only(value)) : _number(value)
_total(value) = value isa AbstractVector ? sum(_number, value) : _number(value)

# The schema types `p_avail` as a scalar and `p_max` as a per-phase array,
# while ENWL writes both as scalars, so total over whatever shape is present.
function _available_power(record, fallback::Float64)
    haskey(record, "p_avail") && return _total(record["p_avail"])
    haskey(record, "p_max") && return _total(record["p_max"])
    fallback
end

function _phase(record)
    terminals = String.(record["terminal_map"])
    phases = filter(terminal -> terminal in PHASES, terminals)
    length(phases) == 1 ? only(phases) : "other"
end

function _snapshot_inventory(path)
    net = parse_bmopf(path)
    ibrs = get(net, "ibr", Dict())
    phase_counts = Dict(phase => 0 for phase in ("a", "b", "c", "other"))
    topology_counts = Dict{String,Int}()
    total_rating = 0.0
    total_available = 0.0
    for record in values(ibrs)
        topology = uppercase(String(get(record, "topology", "FOUR_LEG")))
        topology_counts[topology] = get(topology_counts, topology, 0) + 1
        phase_counts[_phase(record)] += 1
        total_rating += sum(_number.(record["s_max"]))
        total_available += _available_power(record, 0.0)
    end
    (
        feeder=feeder_id(path), snapshot=snapshot_id(path), path=abspath(path),
        buses=length(get(net, "bus", Dict())),
        lines=length(get(net, "line", Dict())),
        loads=length(get(net, "load", Dict())),
        ibrs=length(ibrs), single_phase_ibrs=get(topology_counts, "SINGLE_PHASE", 0),
        phase_a_ibrs=phase_counts["a"], phase_b_ibrs=phase_counts["b"],
        phase_c_ibrs=phase_counts["c"], other_phase_ibrs=phase_counts["other"],
        total_ibr_rating_kVA=total_rating / 1e3,
        total_p_available_kW=total_available / 1e3,
    )
end

inventory_rows(paths) = [_snapshot_inventory(path) for path in paths]

function print_inventory(rows)
    println("feeder\tsnapshots\tbuses\tlines\tloads\tIBRs\t1ph IBRs\tPV kVA\tmax available kW\tphases a/b/c")
    for feeder in sort!(unique(row.feeder for row in rows))
        group = filter(row -> row.feeder == feeder, rows)
        first_row = first(group)
        @assert all(row.buses == first_row.buses && row.lines == first_row.lines &&
                    row.loads == first_row.loads && row.ibrs == first_row.ibrs
                    for row in group)
        println(join((
            feeder, length(group), first_row.buses, first_row.lines,
            first_row.loads, first_row.ibrs, first_row.single_phase_ibrs,
            round(first_row.total_ibr_rating_kVA; digits=2),
            round(maximum(row.total_p_available_kW for row in group); digits=2),
            "$(first_row.phase_a_ibrs)/$(first_row.phase_b_ibrs)/$(first_row.phase_c_ibrs)",
        ), '\t'))
    end
end

function _selected_ids(net, count::Int)
    ids = sort!(collect(String.(keys(get(net, "ibr", Dict()))));
                by=id -> (String(net["ibr"][id]["bus"]), id))
    isempty(ids) && throw(ArgumentError("snapshot contains no IBR records"))
    requested = count == 0 ? length(ids) : count
    1 <= requested <= length(ids) || throw(ArgumentError(
        "devices must be 0 (all) or lie in 1:$(length(ids)); got $count"))
    requested == length(ids) && return ids
    # Evenly span the feeder rather than selecting a lexicographic prefix.
    indices = requested == 1 ? [cld(length(ids), 2)] :
        unique(round.(Int, range(1, length(ids); length=requested)))
    length(indices) == requested || error("internal device selection was not unique")
    ids[indices]
end

const PER_PHASE_REFERENCES = ("PG_PER_PHASE", "PN_PER_PHASE")

function _profile_field(section, name::String, allowed)
    haskey(section, name) || throw(ArgumentError(
        "ENWL control profile must declare '$name'"))
    value = uppercase(String(section[name]))
    value in allowed || throw(ArgumentError(
        "unsupported $name '$value'; this study supports $(join(allowed, ", "))"))
    value
end

function _as4777_laws(net)
    profiles = get(net, "control_profile", Dict())
    length(profiles) == 1 || throw(ArgumentError(
        "expected exactly one ENWL control profile, found $(length(profiles))"))
    profile = only(values(profiles))
    vv = profile["volt_var"]
    vw = profile["volt_watt"]

    # The ordinate units decide what the curve fractions multiply. Checking
    # them keeps the retrofitted plants on the same law as the native PVs
    # solved alongside them: q_scale carries VAR_MAX and p_rated carries S_MAX.
    _profile_field(vv, "q_unit", ("VA_FRACTION",))
    _profile_field(vv, "q_ref", ("VAR_MAX",))
    _profile_field(vw, "p_unit", ("VA_FRACTION",))
    volt_watt_basis =
        _profile_field(vw, "p_ref", ("S_MAX", "P_AVAILABLE")) == "S_MAX" ?
        :rated : :available

    # Both curves must monitor the same per-phase magnitude; the caller decides
    # what to do about the phase-to-neutral/phase-to-ground gap.
    reference = _profile_field(vv, "voltage_reference", PER_PHASE_REFERENCES)
    vw_reference = _profile_field(vw, "voltage_reference", PER_PHASE_REFERENCES)
    reference == vw_reference || throw(ArgumentError(
        "volt-var and volt-watt must monitor the same voltage reference; " *
        "got '$reference' and '$vw_reference'"))

    vv_breakpoints = Float64.(vv["breakpoints"])
    q_limits = Float64.(vv["q_limits"])
    length(vv_breakpoints) == 4 && length(q_limits) == 2 || throw(ArgumentError(
        "expected the four-point AS/NZS volt-var profile used by ENWL"))
    q_low, q_high = extrema(q_limits)
    volt_var = PiecewiseLinearLaw(
        vv_breakpoints, [q_high, 0.0, 0.0, q_low]; smoothing_epsilon=0.05)
    p_limits = Float64.(vw["p_limits"])
    length(p_limits) == 2 || throw(ArgumentError(
        "expected the two-point AS/NZS volt-watt profile used by ENWL"))
    volt_watt = PiecewiseLinearLaw(
        Float64.(vw["breakpoints"]), reverse(sort(p_limits));
        smoothing_epsilon=0.05)
    volt_watt, volt_var, volt_watt_basis, reference
end

function _controller(variant::String, volt_watt, volt_var, volt_watt_basis::Symbol)
    positive = if variant == "baseline"
        AverageVoltageVoltVarWatt(
            volt_watt=volt_watt, volt_var=volt_var,
            volt_watt_basis=volt_watt_basis)
    elseif variant == "worst_phase"
        WorstPhaseVoltVarWatt(
            volt_watt=volt_watt, volt_var=volt_var,
            extrema_epsilon=0.5, conflict_epsilon=0.05,
            volt_watt_basis=volt_watt_basis)
    elseif variant == "sequence_droop"
        PositiveSequenceVoltVarWatt(
            volt_watt=volt_watt, volt_var=volt_var,
            volt_watt_basis=volt_watt_basis)
    else
        throw(ArgumentError("unknown variant '$variant'; choose from " *
                            join(DEFAULT_VARIANTS, ", ")))
    end
    unbalance = if variant == "sequence_droop"
        gain = PiecewiseLinearLaw(
            [0.0, 0.005, 0.02, 0.05], [0.0, 0.0, 0.08, 0.16];
            smoothing_epsilon=2.5e-4)
        NegativeSequenceAdmittanceDroop(gain; ripple_blend=0.25)
    else
        NoUnbalanceControl()
    end
    SequenceController(positive; unbalance=unbalance)
end

function _retrofit!(net, id::String; nominal_voltage::Float64=230.0)
    record = net["ibr"][id]
    uppercase(String(get(record, "topology", ""))) == "SINGLE_PHASE" ||
        throw(ArgumentError("ENWL retrofit expects SINGLE_PHASE IBR '$id'"))
    original_terminals = String.(record["terminal_map"])
    original_phase = _phase(record)
    original_phase != "other" || throw(ArgumentError(
        "IBR '$id' does not have one a/b/c phase terminal"))
    bus = String(record["bus"])
    bus_terminals = Set(String.(net["bus"][bus]["terminal_names"]))
    issubset(Set([PHASES; "n"]), bus_terminals) || throw(ArgumentError(
        "IBR '$id' bus '$bus' is not an a-b-c-n bus"))

    rating = sum(_number.(record["s_max"]))
    available = _available_power(record, rating)
    q_scale = max(abs(_scalar(get(record, "q_min", -rating))),
                  abs(_scalar(get(record, "q_max", rating))))
    per_phase_rating = rating / 3

    # The native record is only a typed placeholder after ownership transfers
    # to PowerOptLab, but keeping it internally coherent makes the adapted
    # network inspectable and validates the intended three-wire connection.
    record["topology"] = "THREE_LEG"
    record["terminal_map"] = copy(PHASES)
    record["s_max"] = fill(per_phase_rating, 3)
    record["p_min"] = zeros(3)
    record["p_max"] = fill(per_phase_rating, 3)
    record["q_min"] = fill(-per_phase_rating, 3)
    record["q_max"] = fill(per_phase_rating, 3)

    current_rating = 1.20rating / (3nominal_voltage)
    inverter = AdvancedInverter(
        id=id, bus=bus, phase_terminals=copy(PHASES), neutral=nothing,
        topology=:THREE_LEG, s_max=rating, i_max=current_rating,
        v_dc=700.0, c_dc=1.1e-3 * rating / 20e3,
        m_max=0.96, r_filter=0.05, x_filter=0.15)
    request = InverterControlRequest(
        p_available=available, p_rated=rating, q_scale=q_scale)
    provenance = Dict{String,Any}(
        "original_topology" => "SINGLE_PHASE",
        "original_terminal_map" => original_terminals,
        "original_phase" => original_phase,
        "synthetic_topology" => "THREE_LEG",
        "rating_VA" => rating,
        "available_W" => available,
    )
    inverter, request, provenance
end

function _build_snapshot_cases(path::String, variants, device_count::Int)
    net = parse_bmopf(path)
    selected = _selected_ids(net, device_count)
    volt_watt, volt_var, volt_watt_basis, reference = _as4777_laws(net)
    if reference == "PN_PER_PHASE"
        # A three-wire retrofit has no neutral to sense against, so the curve
        # input carries the local neutral displacement on these feeders.
        message = "profile monitors phase-to-neutral voltage, but the " *
                  "retrofitted three-wire plants sense phase-to-ground"
        @warn message reference maxlog=1
    end
    inverters = Dict{String,AdvancedInverter}()
    requests = Dict{String,InverterControlRequest}()
    provenance = Dict{String,Any}()
    for id in selected
        inverter, request, source = _retrofit!(net, id)
        inverters[id] = inverter
        requests[id] = request
        provenance[id] = source
    end
    metadata = Dict{String,Any}(
        "feeder" => feeder_id(path),
        "source_path" => abspath(path),
        "retrofit" => "single_phase_to_three_leg_total_rating_preserved",
        "controlled_device_count" => length(selected),
        "selected_device_ids" => selected,
        "device_provenance" => provenance,
        "profile_voltage_reference" => reference,
        "sensed_voltage_reference" => "PG_PER_PHASE",
        "volt_watt_basis" => String(volt_watt_basis),
    )
    [begin
        devices = Dict(id => ControlledDevice(
            inverter, _controller(variant, volt_watt, volt_var, volt_watt_basis))
            for (id, inverter) in inverters)
        InverterControlStudyCase(
            scenario_id=snapshot_id(path), variant_id=variant,
            network=net, fleet=ControlledInverterFleetSpec(devices, requests),
            duration_h=0.5, metadata=metadata)
     end for variant in variants]
end

function build_study_cases(paths; variants=DEFAULT_VARIANTS, device_count=3)
    isempty(variants) && throw(ArgumentError("at least one variant is required"))
    reduce(vcat, (_build_snapshot_cases(path, variants, device_count)
                  for path in paths); init=InverterControlStudyCase[])
end

function print_campaign_plan(cases)
    println("ENWL advanced-IBR campaign")
    println("  cases: ", length(cases))
    println("  snapshots: ", length(unique(case.scenario_id for case in cases)))
    println("  variants: ", join(sort!(unique(case.variant_id for case in cases)), ", "))
    println("  controlled devices/case: ",
            join(sort!(unique(length(case.fleet.devices) for case in cases)), ", "))
    println("  retrofit: single-phase PV -> synthetic balanced THREE_LEG plant")
    println("  invariant: original total VA rating and available W are preserved")
end

_csv_names(row::NamedTuple) = collect(propertynames(row))
_csv_names(row::AbstractDict) = sort!(collect(keys(row)); by=string)
_csv_get(row::NamedTuple, name) = getproperty(row, name)
_csv_get(row::AbstractDict, name) = row[name]

function _csv_field(value)
    value === nothing && return ""
    value === missing && return ""
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    text
end

function write_csv(path::AbstractString, rows)
    isempty(rows) && return nothing
    mkpath(dirname(path))
    names = _csv_names(first(rows))
    open(path, "w") do io
        println(io, join(_csv_field.(string.(names)), ','))
        for row in rows
            println(io, join((_csv_field(_csv_get(row, name)) for name in names), ','))
        end
    end
    path
end

end # module
