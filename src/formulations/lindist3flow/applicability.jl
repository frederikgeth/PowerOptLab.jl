const _L3F_UNSUPPORTED_FAMILIES = (
    "switch", "capacitor", "ibr",
    "dc_bus", "dc_line", "dc_load", "dc_source", "dc_converter",
)

_l3f_error!(out, code, component, id, message; evidence=Dict{String,Any}()) =
    push!(out, L3FFinding(code, :error, component,
        isnothing(id) ? nothing : String(id), message, Dict{String,Any}(evidence)))

_l3f_warning!(out, code, component, id, message; evidence=Dict{String,Any}()) =
    push!(out, L3FFinding(code, :warning, component,
        isnothing(id) ? nothing : String(id), message, Dict{String,Any}(evidence)))

function _l3f_input(input)
    if input isa AbstractString
        return BMOPFTools.parse_bmopf(String(input); from_string=true)
    elseif input isa AbstractDict
        return deepcopy(Dict{String,Any}(string(k) => v for (k, v) in input))
    end
    throw(ArgumentError("L3F input must be a BMOPF dictionary or JSON string"))
end

function _l3f_shunt_matrix(data, n::Int)
    value(prefix, i, j) = if haskey(data, "$(prefix)$(i)_$(j)")
        Float64(data["$(prefix)$(i)_$(j)"])
    elseif haskey(data, "$(prefix)$(j)_$(i)")
        Float64(data["$(prefix)$(j)_$(i)"])
    else
        0.0
    end
    [value("G_", i, j) + im * value("B_", i, j) for i in 1:n, j in 1:n]
end

function _l3f_connection_incidence(configuration::String, nterminal::Int, nchannel::Int)
    cfg = uppercase(configuration)
    if cfg in ("WYE", "SINGLE_PHASE")
        nterminal == nchannel || throw(DimensionMismatch(
            "$cfg requires one retained terminal per power channel"))
        return Matrix{Float64}(I, nterminal, nterminal)
    elseif cfg == "DELTA"
        if nterminal == 2 && nchannel == 1
            return reshape([1.0, -1.0], 1, 2)
        elseif nterminal == 3 && nchannel == 3
            return [1.0 -1.0 0.0; 0.0 1.0 -1.0; -1.0 0.0 1.0]
        end
        throw(DimensionMismatch("DELTA requires two terminals/one channel or three terminals/three channels"))
    end
    throw(ArgumentError("unsupported connection configuration '$configuration'"))
end

function _l3f_validate_shunts!(findings, net)
    buses = get(net, "bus", Dict())
    for (id, data) in get(net, "shunt", Dict())
        busid = String(get(data, "bus", "")); bus = get(buses, busid, nothing)
        tm = string.(get(data, "terminal_map", String[]))
        valid = bus isa AbstractDict && !isempty(tm) && allunique(tm) &&
            all(in(string.(get(bus, "terminal_names", String[]))), tm)
        valid || _l3f_error!(findings, "E.L3F.TERMINAL_MAP_INVALID", :shunt, id,
            "fixed shunt terminal_map must contain declared retained bus terminals")
        try
            Y = _l3f_shunt_matrix(data, length(tm))
            all(isfinite, real.(Y)) && all(isfinite, imag.(Y)) ||
                throw(ArgumentError("shunt admittance must be finite"))
        catch err
            _l3f_error!(findings, "E.L3F.SHUNT_INVALID", :shunt, id, sprint(showerror, err))
        end
        haskey(data, "time_series") && _l3f_error!(findings,
            "E.L3F.TIME_SERIES_UNSUPPORTED", :shunt, id,
            "LinDist3Flow accepts fixed shunt admittance only")
    end
end

function _l3f_transformers(net)
    out = Tuple{String,String,AbstractDict}[]
    for (subtype, table) in get(net, "transformer", Dict())
        table isa AbstractDict || continue
        for (id, data) in table
            data isa AbstractDict || continue
            push!(out, (String(subtype), String(id), data))
        end
    end
    out
end

function _l3f_transformer_neff(subtype::String, data)
    if subtype == "single_phase"
        Float64(data["v_nom_from"]) / Float64(data["v_nom_to"]) *
            Float64(get(data, "tap", 1.0))
    elseif subtype == "single_phase_autotransformer"
        tap = Float64(get(data, "tap_ratio", 1.0))
        uppercase(String(get(data, "regulator_type", "B"))) == "A" ? inv(tap) : tap
    else
        throw(ArgumentError("unsupported transformer subtype '$subtype'"))
    end
end

function _l3f_validate_transformers!(findings, net)
    buses = get(net, "bus", Dict())
    for (subtype, id, data) in _l3f_transformers(net)
        subtype in ("single_phase", "single_phase_autotransformer") || begin
            _l3f_error!(findings, "E.L3F.TRANSFORMER_UNSUPPORTED", :transformer, id,
                "only fixed ideal single_phase and single_phase_autotransformer devices are supported")
            continue
        end
        bf, bt = String(get(data, "bus_from", "")), String(get(data, "bus_to", ""))
        mf, mt = string.(get(data, "terminal_map_from", String[])),
                 string.(get(data, "terminal_map_to", String[]))
        valid = haskey(buses, bf) && haskey(buses, bt) && !isempty(mf) &&
            length(mf) == length(mt) && allunique(mf) && allunique(mt) &&
            all(in(string.(get(buses[bf], "terminal_names", String[]))), mf) &&
            all(in(string.(get(buses[bt], "terminal_names", String[]))), mt)
        valid || _l3f_error!(findings, "E.L3F.TERMINAL_MAP_INVALID", :transformer, id,
            "transformer terminal maps must be aligned retained conductors on declared buses")
        length(mf) == 1 || _l3f_error!(findings, "E.L3F.DEVICE_ARITY", :transformer, id,
            "single-phase transformer/regulator support requires one retained power channel")
        lo_key, hi_key = subtype == "single_phase" ? ("tap_min", "tap_max") :
                                                    ("tap_ratio_min", "tap_ratio_max")
        if haskey(data, lo_key) || haskey(data, hi_key)
            lo, hi = Float64(get(data, lo_key, NaN)), Float64(get(data, hi_key, NaN))
            (!isfinite(lo) || !isfinite(hi) || lo != hi) && _l3f_error!(findings,
                "E.L3F.ADJUSTABLE_TAP_UNSUPPORTED", :transformer, id,
                "LinDist3Flow accepts fixed regulator settings only; adjustable tap intervals are excluded")
        end
        try
            neff = _l3f_transformer_neff(subtype, data)
            isfinite(neff) && neff > 0 || throw(ArgumentError("effective ratio must be positive and finite"))
        catch err
            _l3f_error!(findings, "E.L3F.TRANSFORMER_RATIO_INVALID", :transformer, id,
                        sprint(showerror, err))
        end
        for key in ("r_series_from", "x_series_from", "r_series_to", "x_series_to",
                    "g_no_load", "b_no_load")
            abs(Float64(get(data, key, 0.0))) <= 1e-12 || _l3f_error!(findings,
                "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED", :transformer, id,
                "fixed-ratio support is ideal; nonzero '$key' must be represented separately")
        end
        for key in ("s_rating", "i_max_from", "i_max_to")
            haskey(data, key) || continue
            value = data[key]
            value isa Real && isfinite(value) && value > 0 || _l3f_error!(findings,
                "E.L3F.LIMIT_INVALID", :transformer, id,
                "$key must be a positive finite scalar")
        end
    end
end

function _l3f_neutral_labels(net)
    labels = Set{String}()
    tc = get(net, "terminal_conventions", nothing)
    if tc isa AbstractDict && haskey(tc, "neutral")
        foreach(x -> push!(labels, string(x)), tc["neutral"])
    end
    for bus in values(get(net, "bus", Dict()))
        bus isa AbstractDict || continue
        nt = get(bus, "neutral_terminal", nothing)
        nt === nothing || push!(labels, string(nt))
    end
    isempty(labels) && push!(labels, "n")
    labels
end

function _l3f_has_explicit_neutral(net)
    labels = _l3f_neutral_labels(net)
    any(bus isa AbstractDict && any(t -> string(t) in labels || lowercase(string(t)) == "n",
                                   get(bus, "terminal_names", String[]))
        for bus in values(get(net, "bus", Dict())))
end

function _l3f_has_matrix_key(data, prefixes)
    any(key -> any(prefix -> startswith(String(key), prefix), prefixes), keys(data))
end

function _l3f_pattern_matrix(data::AbstractDict, real_prefix::String,
                             imag_prefix::String, n::Int, label::String)
    n > 0 || throw(ArgumentError("$label has no conductors"))
    pattern = r"^(R_series_|X_series_)(\d+)_(\d+)$"
    seen = false
    for key in keys(data)
        m = match(pattern, String(key))
        m === nothing && continue
        seen = true
        max(parse(Int, m.captures[2]), parse(Int, m.captures[3])) <= n ||
            throw(DimensionMismatch("$label matrix key '$key' exceeds terminal-map arity $n"))
    end
    seen || throw(ArgumentError("$label has no series impedance matrix"))
    value(prefix, i, j) = if haskey(data, "$(prefix)$(i)_$(j)")
        Float64(data["$(prefix)$(i)_$(j)"])
    elseif haskey(data, "$(prefix)$(j)_$(i)")
        Float64(data["$(prefix)$(j)_$(i)"])
    else
        0.0
    end
    Z = Matrix{ComplexF64}(undef, n, n)
    for i in 1:n, j in 1:n
        Z[i, j] = value(real_prefix, i, j) + im * value(imag_prefix, i, j)
    end
    all(isfinite, real.(Z)) && all(isfinite, imag.(Z)) ||
        throw(ArgumentError("$label series matrix contains non-finite values"))
    Z
end

function _l3f_series_matrix(net, line, n::Int, id::String)
    inline = _l3f_has_matrix_key(line, ("R_series_", "X_series_"))
    if inline
        haskey(line, "linecode") && throw(ArgumentError(
            "line '$id' defines both inline impedance and a linecode"))
        return _l3f_pattern_matrix(line, "R_series_", "X_series_", n, "line '$id'")
    end
    lcid = get(line, "linecode", nothing)
    lcid isa AbstractString || throw(ArgumentError("line '$id' has no impedance source"))
    lc = get(get(net, "linecode", Dict()), String(lcid), nothing)
    lc isa AbstractDict || throw(ArgumentError("line '$id' references unknown linecode '$lcid'"))
    length = Float64(get(line, "length", 1.0))
    isfinite(length) && length > 0 || throw(ArgumentError("line '$id' length must be positive"))
    length .* _l3f_pattern_matrix(lc, "R_series_", "X_series_", n,
                                  "linecode '$lcid'")
end

function _l3f_bus_bound(bus, key::String, k::Int, n::Int)
    value = get(bus, key, nothing)
    value === nothing && return nothing
    if value isa Real
        return Float64(value)
    elseif value isa AbstractVector
        length(value) == n || throw(DimensionMismatch(
            "bus $key length $(length(value)) does not match $n terminals"))
        return Float64(value[k])
    end
    throw(ArgumentError("bus $key must be a scalar or vector"))
end

function _l3f_component_time_series!(findings, net)
    for family in ("bus", "line", "load", "generator", "voltage_source")
        for (id, component) in get(net, family, Dict())
            component isa AbstractDict || continue
            haskey(component, "time_series") || continue
            _l3f_error!(findings, "E.L3F.TIME_SERIES_UNSUPPORTED", Symbol(family), id,
                "L3F version 0.1 accepts one resolved snapshot, not component time-series references")
        end
    end
end

function _l3f_validate_components!(findings, net)
    buses = get(net, "bus", Dict())
    isempty(buses) && _l3f_error!(findings, "E.L3F.EMPTY_NETWORK", :network, nothing,
                                  "the network has no buses")

    for family in _L3F_UNSUPPORTED_FAMILIES
        table = get(net, family, nothing)
        table isa AbstractDict && !isempty(table) || continue
        for id in keys(table)
            code = startswith(family, "dc_") ? "E.L3F.DC_SUBSYSTEM_UNSUPPORTED" :
                   family == "transformer" ? "E.L3F.TRANSFORMER_UNSUPPORTED" :
                   family == "switch" ? "E.L3F.SWITCH_UNSUPPORTED" :
                   "E.L3F.COMPONENT_UNSUPPORTED"
            _l3f_error!(findings, code, Symbol(family), id,
                "component family '$family' is outside the minimal L3F implementation slice")
        end
    end
    for key in ("control_profile", "time_series")
        table = get(net, key, nothing)
        table isa AbstractDict && !isempty(table) || continue
        _l3f_error!(findings, "E.L3F.CONTROL_PROFILE_UNSUPPORTED", :network, nothing,
            "top-level '$key' data must be resolved before building an L3F snapshot")
    end
    _l3f_component_time_series!(findings, net)

    for (bid, bus) in buses
        bus isa AbstractDict || begin
            _l3f_error!(findings, "E.L3F.BUS_INVALID", :bus, bid,
                        "bus must be an object"); continue
        end
        terminals = string.(get(bus, "terminal_names", String[]))
        isempty(terminals) && _l3f_error!(findings, "E.L3F.TERMINAL_MAP_INVALID", :bus, bid,
                                          "bus has no terminal_names")
        allunique(terminals) || _l3f_error!(findings, "E.L3F.TERMINAL_MAP_INVALID", :bus, bid,
                                            "bus terminal_names are not unique")
        for key in ("vpp_min", "vpp_max", "vpos_min", "vpos_max", "vneg_max",
                    "vzero_max", "vm_unbalance_max")
            haskey(bus, key) || continue
            _l3f_error!(findings, "E.L3F.LIMIT_UNSUPPORTED", :bus, bid,
                "bus limit '$key' is not assessed by the minimal L3F implementation")
        end
        try
            for k in eachindex(terminals)
                lo = _l3f_bus_bound(bus, "v_min", k, length(terminals))
                hi = _l3f_bus_bound(bus, "v_max", k, length(terminals))
                lo === nothing || (isfinite(lo) && lo >= 0) || throw(ArgumentError("invalid v_min"))
                hi === nothing || (isfinite(hi) && hi >= 0) || throw(ArgumentError("invalid v_max"))
                lo !== nothing && hi !== nothing && lo > hi &&
                    throw(ArgumentError("v_min exceeds v_max at terminal $(terminals[k])"))
            end
        catch err
            _l3f_error!(findings, "E.L3F.VOLTAGE_BOUND_INVALID", :bus, bid,
                        sprint(showerror, err))
        end
    end

    for (sid, source) in get(net, "voltage_source", Dict())
        busid = String(get(source, "bus", ""))
        bus = get(buses, busid, nothing)
        bus isa AbstractDict || begin
            _l3f_error!(findings, "E.L3F.BUS_UNKNOWN", :voltage_source, sid,
                        "source references unknown bus '$busid'"); continue
        end
        cfg = uppercase(String(get(source, "configuration", "WYE")))
        cfg in ("WYE", "SINGLE_PHASE") || _l3f_error!(findings,
            "E.L3F.CONNECTION_UNSUPPORTED", :voltage_source, sid,
            "only grounded-wye/single-phase voltage sources are supported")
        tm = string.(get(source, "terminal_map", String[]))
        bt = string.(get(bus, "terminal_names", String[]))
        !isempty(tm) && allunique(tm) && all(in(bt), tm) || _l3f_error!(findings,
            "E.L3F.TERMINAL_MAP_INVALID", :voltage_source, sid,
            "source terminal_map must contain distinct declared bus terminals")
        Set(tm) == Set(bt) || _l3f_error!(findings, "E.L3F.REFERENCE_MISSING",
            :voltage_source, sid, "the source must reference every terminal in its root bus")
        vm = Float64.(get(source, "v_magnitude", Float64[]))
        va = Float64.(get(source, "v_angle", Float64[]))
        length(vm) == length(tm) == length(va) || _l3f_error!(findings,
            "E.L3F.REFERENCE_MISSING", :voltage_source, sid,
            "v_magnitude/v_angle must cover the source terminal_map")
        length(vm) == length(tm) && any(x -> !isfinite(x) || x <= 0, vm) &&
            _l3f_error!(findings, "E.L3F.REFERENCE_ZERO_WINDING", :voltage_source, sid,
                        "source reference magnitudes must be finite and nonzero")
        for field in ("i_max", "s_max")
            haskey(source, field) || continue
            value = source[field]
            value isa AbstractVector && length(value) == length(tm) &&
                all(x -> x isa Real && isfinite(x) && x > 0, value) ||
                _l3f_error!(findings, "E.L3F.LIMIT_INVALID", :voltage_source, sid,
                            "$field must contain one positive finite rating per terminal")
        end
    end

    for (lid, load) in get(net, "load", Dict())
        busid = String(get(load, "bus", ""))
        bus = get(buses, busid, nothing)
        bus isa AbstractDict || begin
            _l3f_error!(findings, "E.L3F.BUS_UNKNOWN", :load, lid,
                        "load references unknown bus '$busid'"); continue
        end
        cfg = uppercase(String(get(load, "configuration", "WYE")))
        cfg in ("WYE", "SINGLE_PHASE", "DELTA") || _l3f_error!(findings,
            "E.L3F.CONNECTION_UNSUPPORTED", :load, lid,
            "only grounded-wye, single-phase, and delta constant-power loads are supported")
        lowercase(String(get(load, "model", "constant_power"))) == "constant_power" ||
            _l3f_error!(findings, "E.L3F.LOAD_MODEL_UNSUPPORTED", :load, lid,
                        "only model=constant_power is supported in the minimal slice")
        tm = string.(get(load, "terminal_map", String[]))
        bt = string.(get(bus, "terminal_names", String[]))
        !isempty(tm) && allunique(tm) && all(in(bt), tm) || _l3f_error!(findings,
            "E.L3F.TERMINAL_MAP_INVALID", :load, lid,
            "load terminal_map must contain distinct declared phase terminals")
        p = get(load, "p_nom", Any[]); q = get(load, "q_nom", Any[])
        length(p) == length(q) || _l3f_error!(findings,
            "E.L3F.DEVICE_ARITY", :load, lid,
            "p_nom and q_nom must have equal channel counts")
        try
            _l3f_connection_incidence(cfg, length(tm), length(p))
            values = vcat(Float64.(p), Float64.(q))
            all(isfinite, values) || throw(ArgumentError("nominal powers must be finite"))
        catch err
            _l3f_error!(findings, "E.L3F.DEVICE_DATA_INVALID", :load, lid,
                        sprint(showerror, err))
        end
    end

    for (gid, gen) in get(net, "generator", Dict())
        busid = String(get(gen, "bus", ""))
        bus = get(buses, busid, nothing)
        bus isa AbstractDict || begin
            _l3f_error!(findings, "E.L3F.BUS_UNKNOWN", :generator, gid,
                        "generator references unknown bus '$busid'"); continue
        end
        cfg = uppercase(String(get(gen, "configuration", "WYE")))
        cfg in ("WYE", "SINGLE_PHASE", "DELTA") || _l3f_error!(findings,
            "E.L3F.CONNECTION_UNSUPPORTED", :generator, gid,
            "only grounded-wye, single-phase, and delta generators are supported")
        tm = string.(get(gen, "terminal_map", String[]))
        bt = string.(get(bus, "terminal_names", String[]))
        !isempty(tm) && allunique(tm) && all(in(bt), tm) || _l3f_error!(findings,
            "E.L3F.TERMINAL_MAP_INVALID", :generator, gid,
            "generator terminal_map must contain distinct declared phase terminals")
        nch = length(get(gen, "p_min", Any[]))
        try
            _l3f_connection_incidence(cfg, length(tm), nch)
        catch err
            _l3f_error!(findings, "E.L3F.DEVICE_ARITY", :generator, gid, sprint(showerror, err))
        end
        for field in ("p_min", "p_max", "q_min", "q_max")
            value = get(gen, field, nothing)
            value isa AbstractVector && length(value) == nch ||
                _l3f_error!(findings, "E.L3F.DEVICE_ARITY", :generator, gid,
                            "$field is required and must have one value per physical power channel")
        end
        try
            pmin, pmax = Float64.(gen["p_min"]), Float64.(gen["p_max"])
            qmin, qmax = Float64.(gen["q_min"]), Float64.(gen["q_max"])
            all(isfinite, vcat(pmin, pmax, qmin, qmax)) ||
                throw(ArgumentError("generator bounds must be finite"))
            all(pmin .<= pmax) || throw(ArgumentError("p_min exceeds p_max"))
            all(qmin .<= qmax) || throw(ArgumentError("q_min exceeds q_max"))
        catch err
            _l3f_error!(findings, "E.L3F.DEVICE_DATA_INVALID", :generator, gid,
                        sprint(showerror, err))
        end
        for field in ("i_max", "s_max")
            haskey(gen, field) || continue
            value = gen[field]
            value isa AbstractVector && length(value) == nch &&
                all(x -> x isa Real && isfinite(x) && x > 0, value) ||
                _l3f_error!(findings, "E.L3F.LIMIT_INVALID", :generator, gid,
                            "$field must contain one positive finite rating per physical channel")
        end
    end
end

function _l3f_validate_lines!(findings, net)
    buses = get(net, "bus", Dict())
    linecodes = get(net, "linecode", Dict())
    for (id_raw, line) in get(net, "line", Dict())
        id = String(id_raw)
        bf, bt = String(get(line, "bus_from", "")), String(get(line, "bus_to", ""))
        haskey(buses, bf) && haskey(buses, bt) || begin
            _l3f_error!(findings, "E.L3F.BUS_UNKNOWN", :line, id,
                        "line endpoints must reference declared buses"); continue
        end
        mf = string.(get(line, "terminal_map_from", String[]))
        mt = string.(get(line, "terminal_map_to", String[]))
        valid = !isempty(mf) && length(mf) == length(mt) && allunique(mf) && allunique(mt) &&
            all(in(string.(get(buses[bf], "terminal_names", String[]))), mf) &&
            all(in(string.(get(buses[bt], "terminal_names", String[]))), mt)
        valid || _l3f_error!(findings, "E.L3F.TERMINAL_MAP_INVALID", :line, id,
                             "line terminal maps must be aligned, unique, and declared at each bus")
        if valid
            try
                _l3f_series_matrix(net, line, length(mf), id)
            catch err
                _l3f_error!(findings, "E.L3F.LINE_MATRIX_INVALID", :line, id,
                            sprint(showerror, err))
            end
        end
        sources = Any[line]
        lcid = get(line, "linecode", nothing)
        lcid isa AbstractString && haskey(linecodes, lcid) && push!(sources, linecodes[lcid])
        if any(data -> _l3f_has_matrix_key(data,
                    ("G_from_", "B_from_", "G_to_", "B_to_")), sources)
            _l3f_error!(findings, "E.L3F.LINE_SHUNT_UNSUPPORTED", :line, id,
                        "line shunts require the Phase 3 affine shunt kernel")
        end
        for field in ("i_max", "s_max")
            data = haskey(line, field) ? line :
                   (lcid isa AbstractString && haskey(linecodes, lcid) &&
                    haskey(linecodes[lcid], field) ? linecodes[lcid] : nothing)
            data === nothing && continue
            value = data[field]
            value isa AbstractVector && length(value) == length(mf) &&
                all(x -> x isa Real && isfinite(x) && x > 0, value) ||
                _l3f_error!(findings, "E.L3F.LIMIT_INVALID", :line, id,
                            "$field must contain one positive finite rating per conductor")
        end
    end
end

function _l3f_topology!(findings, net)
    buses = sort!(String.(collect(keys(get(net, "bus", Dict())))))
    adjacency = Dict(bus => Tuple{String,String,Symbol,String}[] for bus in buses)
    pair_ids = Dict{Tuple{String,String},Vector{String}}()
    for (id_raw, line) in get(net, "line", Dict())
        id = String(id_raw)
        a, b = String(get(line, "bus_from", "")), String(get(line, "bus_to", ""))
        haskey(adjacency, a) && haskey(adjacency, b) || continue
        push!(adjacency[a], (b, id, :line, "")); push!(adjacency[b], (a, id, :line, ""))
        pair = a < b ? (a, b) : (b, a)
        push!(get!(pair_ids, pair, String[]), "line/$id")
    end
    for (subtype, id, transformer) in _l3f_transformers(net)
        a, b = String(get(transformer, "bus_from", "")), String(get(transformer, "bus_to", ""))
        haskey(adjacency, a) && haskey(adjacency, b) || continue
        push!(adjacency[a], (b, id, :transformer, subtype))
        push!(adjacency[b], (a, id, :transformer, subtype))
        pair = a < b ? (a, b) : (b, a)
        push!(get!(pair_ids, pair, String[]), "transformer/$subtype/$id")
    end
    for (pair, ids) in pair_ids
        length(ids) <= 1 || _l3f_error!(findings, "E.L3F.TOPOLOGY_NOT_RADIAL",
            :network, nothing, "parallel energized lines connect $(pair[1]) and $(pair[2])";
            evidence=Dict("lines" => ids))
    end

    islands = Vector{Vector{String}}()
    unseen = Set(buses)
    while !isempty(unseen)
        start = minimum(unseen)
        queue = [start]; delete!(unseen, start); component = String[]
        while !isempty(queue)
            node = popfirst!(queue); push!(component, node)
            for (neighbor, _, _, _) in adjacency[node]
                neighbor in unseen || continue
                delete!(unseen, neighbor); push!(queue, neighbor)
            end
        end
        push!(islands, sort!(component))
    end

    source_buses = Dict{String,Vector{String}}()
    for (sid, source) in get(net, "voltage_source", Dict())
        push!(get!(source_buses, String(get(source, "bus", "")), String[]), String(sid))
    end
    roots = String[]
    topology = L3FOrientedLine[]
    lines = get(net, "line", Dict())
    for island in islands
        members = Set(island)
        sources = [(bus, sid) for (bus, ids) in source_buses if bus in members for sid in ids]
        if length(sources) != 1
            code = isempty(sources) ? "E.L3F.SOURCE_MISSING" : "E.L3F.MULTIPLE_SOURCES"
            _l3f_error!(findings, code, :network, nothing,
                "energized island $(join(island, ", ")) has $(length(sources)) active sources")
            continue
        end
        root = first(first(sources)); push!(roots, root)
        edge_ids = Set((family, subtype, id) for bus in island
                       for (_, id, family, subtype) in adjacency[bus])
        if length(edge_ids) != length(island) - 1
            _l3f_error!(findings, "E.L3F.TOPOLOGY_NOT_RADIAL", :network, nothing,
                "island rooted at '$root' has $(length(island)) buses and $(length(edge_ids)) lines")
            continue
        end
        visited = Set([root]); queue = [root]
        while !isempty(queue)
            parent = popfirst!(queue)
            for (child, id, family, subtype) in sort(adjacency[parent]; by=x -> (x[3], x[2]))
                child in visited && continue
                push!(visited, child); push!(queue, child)
                data = family == :line ? lines[id] : net["transformer"][subtype][id]
                original = String(get(data, "bus_from", "")) == parent
                pm = string.(get(data, original ? "terminal_map_from" : "terminal_map_to", String[]))
                cm = string.(get(data, original ? "terminal_map_to" : "terminal_map_from", String[]))
                push!(topology, L3FOrientedLine(id, family, subtype,
                    parent, child, pm, cm, !original))
            end
        end
    end
    topology, roots, islands
end

function _l3f_prepare(input; options::L3FOptions=L3FOptions(), reference=nothing)
    net = _l3f_input(input)
    findings = L3FFinding[]
    reduced = false
    if _l3f_has_explicit_neutral(net)
        if options.kron_reduce
            try
                net = kron_reduce_bmopf(net)
                reduced = true
            catch err
                _l3f_error!(findings, "E.L3F.KRON_REDUCTION_FAILED", :network, nothing,
                            sprint(showerror, err))
            end
        else
            _l3f_error!(findings, "E.L3F.EXPLICIT_NEUTRAL_UNSUPPORTED", :network, nothing,
                "explicit neutrals must be Kron-reduced before L3F model construction")
        end
    end
    if options.require_neutral_provenance &&
       !haskey(get(net, "_meta", Dict()), "kron_reduction") && reference === nothing
        _l3f_error!(findings, "E.L3F.NEUTRAL_REDUCTION_UNDECLARED", :network, nothing,
            "neutral-reduction provenance or an explicit reference is required by options")
    end
    options.reference_policy == :explicit && reference === nothing && _l3f_error!(findings,
        "E.L3F.REFERENCE_MISSING", :network, nothing,
        "reference_policy=:explicit requires a reference argument")
    _l3f_validate_components!(findings, net)
    _l3f_validate_shunts!(findings, net)
    _l3f_validate_transformers!(findings, net)
    _l3f_validate_lines!(findings, net)
    topology, roots, islands = _l3f_topology!(findings, net)
    report = L3FApplicabilityReport(
        any(f -> f.severity == :error, findings) ? :inapplicable : :applicable,
        findings, roots, islands, reduced)
    (network=net, applicability=report, topology=topology)
end

"""
    check_l3f_applicability(net; options=L3FOptions(), reference=nothing)

Validate and classify a BMOPF snapshot without constructing a JuMP model. If
`options.kron_reduce` is true, an explicit neutral is reduced on a deep copy
before the supported-domain checks. The returned report preserves every
finding and never mutates the caller's data.
"""
function check_l3f_applicability(net; options::L3FOptions=L3FOptions(), reference=nothing)
    _l3f_prepare(net; options, reference).applicability
end
