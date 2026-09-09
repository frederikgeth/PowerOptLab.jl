# Yd / Dy banks. Their voltage map is singular, so which winding faces the
# source is part of the model rather than a bookkeeping detail.
const _L3F_DELTA_SUBTYPES = ("wye_delta", "delta_wye")

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
    # Every declared entry must lie inside the terminal-map arity. Silently
    # truncating an out-of-range entry would drop physics the caller declared,
    # which the applicability contract forbids.
    for key in keys(data)
        m = match(r"^([GB])_(\d+)_(\d+)$", String(key))
        m === nothing && continue
        i, j = parse(Int, m.captures[2]), parse(Int, m.captures[3])
        min(i, j) >= 1 && max(i, j) <= n || throw(
            DimensionMismatch("shunt admittance key '$key' falls outside terminal indices 1:$n"))
    end
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
        uppercase(String(get(data, "regulator_type", "B"))) == "A" ? tap : inv(tap)
    elseif subtype == "open_delta_regulator"
        tap = Float64.(get(data, "tap_ratio", [1.0, 1.0]))
        uppercase(String(get(data, "regulator_type", "B"))) == "A" ? tap : inv.(tap)
    elseif subtype in _L3F_DELTA_SUBTYPES
        # Both v_nom are phase-to-neutral equivalents, so this is the nominal
        # from/to ratio; the sqrt(3) of the delta coil lives in the gain, not here.
        Float64(data["v_nom_from"]) / Float64(data["v_nom_to"])
    else
        throw(ArgumentError("unsupported transformer subtype '$subtype'"))
    end
end

"""
    _l3f_delta_transformer_gain(subtype, data) -> (g, delta_is_from)

Coil relation of an ideal Yd/Dy bank as ``D v_{delta} = g\\,v_{wye}``, where `D`
is the delta incidence and `v_wye` is measured against the (grounded, reduced)
star point.

This is BMOPFTools' executable convention. Its `wye_delta` uses
``n_{eff}=\\sqrt3/N`` and its `delta_wye` uses ``n_{eff}=N\\sqrt3``, both with
``N = v^{nom}_{from}/v^{nom}_{to}`` and both `v_nom` given as phase-to-neutral
equivalents. The two collapse to the same statement,
``g=\\sqrt3\\,v^{nom}_{delta}/v^{nom}_{wye}``: the delta coil spans a
line-to-line voltage while its nominal is quoted phase-to-neutral, and the
``\\sqrt3`` is exactly that difference.
"""
function _l3f_delta_transformer_gain(subtype::String, data)
    ratio = Float64(data["v_nom_from"]) / Float64(data["v_nom_to"])
    isfinite(ratio) && ratio > 0 ||
        throw(ArgumentError("v_nom_from / v_nom_to must be positive and finite"))
    subtype == "wye_delta" ? (sqrt(3.0) / ratio, false) : (sqrt(3.0) * ratio, true)
end

"""
    _l3f_transformer_map(subtype, data, parent_is_from) -> Matrix{Float64}

The fixed voltage map ``v_{child} = T v_{parent}`` for one orientation.

Every supported subtype is built directly for the orientation asked for rather
than by inverting a forward map, because a delta winding's map is singular: the
delta incidence has rank 2, so the transform cannot be inverted and the
orientation genuinely matters.

With the delta winding upstream the map is ``T = D/g`` as an ideal connection
map — the
downstream wye voltages are fully determined, and their zero-sequence component
is zero, which is the correct behaviour of an ideal bank with no zero-sequence
impedance. With the wye winding upstream, the delta terminals are determined
only up to a common offset; that case needs a gauge and is handled by
[`_l3f_delta_gauge_map`](@ref).
"""
function _l3f_transformer_map(subtype::String, data, parent_is_from::Bool)
    if subtype in ("single_phase", "single_phase_autotransformer")
        n = _l3f_transformer_neff(subtype, data)
        return reshape([parent_is_from ? inv(n) : n], 1, 1)
    elseif subtype == "open_delta_regulator"
        A = regulator_gain_matrix("OPEN_DELTA", get(data, "tap_ratio", [1.0, 1.0]);
            connection=String(get(data, "connection", "")),
            regulator_type=String(get(data, "regulator_type", "B")))
        return parent_is_from ? inv(A) : A
    elseif subtype in _L3F_DELTA_SUBTYPES
        g, delta_is_from = _l3f_delta_transformer_gain(subtype, data)
        D = _l3f_connection_incidence("DELTA", 3, 3)
        return parent_is_from == delta_is_from ? D ./ g : _l3f_delta_gauge_map(D, g)
    end
    throw(ArgumentError("unsupported transformer subtype '$subtype'"))
end

"""
    _l3f_delta_gauge_map(D, g) -> Matrix{Float64}

Voltage map for a Yd/Dy bank traversed with the **wye** winding upstream.

``D v_{delta} = g\\,v_{wye}`` leaves the delta terminal voltages free in their
common (zero-sequence) component, because a delta winding neither imposes nor
carries one. The pseudo-inverse selects the minimum-norm solution, which is the
one orthogonal to `D`'s null space `span{1}` — that is, the solution with zero
zero-sequence voltage at the delta bus.

That is a modelling assumption, not physics: in a phase-to-ground formulation
the delta bus's ground reference actually comes from capacitive coupling the
formulation has already discarded. It is offered only under
`L3FOptions(unsupported=:approximate)`, with `A.L3F.DELTA_ZERO_SEQUENCE_GAUGE`
recording it.
"""
_l3f_delta_gauge_map(D::AbstractMatrix, g::Real) = g .* pinv(Matrix{Float64}(D))

function _l3f_transformer_forward_map(subtype::String, data)
    _l3f_transformer_map(subtype, data, true)
end

_l3f_transformer_oriented_map(edge::L3FOrientedLine, data) =
    _l3f_transformer_map(edge.subtype, data, !edge.reversed)

"""Whether `edge` presents a delta winding to its parent bus."""
function _l3f_delta_parent(edge::L3FOrientedLine, data)
    edge.subtype in _L3F_DELTA_SUBTYPES || return nothing
    _, delta_is_from = _l3f_delta_transformer_gain(edge.subtype, data)
    (!edge.reversed) == delta_is_from
end

function _l3f_validate_transformers!(findings, net)
    buses = get(net, "bus", Dict())
    for (subtype, id, data) in _l3f_transformers(net)
        subtype in ("single_phase", "single_phase_autotransformer",
                    "open_delta_regulator", _L3F_DELTA_SUBTYPES...) || begin
            _l3f_error!(findings, "E.L3F.TRANSFORMER_UNSUPPORTED", :transformer, id,
                "only fixed ideal single_phase, single_phase_autotransformer, " *
                "wye_delta, delta_wye, and open_delta_regulator devices are supported")
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
        expected_arity = subtype == "open_delta_regulator" ||
                         subtype in _L3F_DELTA_SUBTYPES ? 3 : 1
        length(mf) == expected_arity || _l3f_error!(findings,
            "E.L3F.DEVICE_ARITY", :transformer, id,
            subtype == "open_delta_regulator" ?
                "a Kron-reduced open-delta regulator requires three retained phase terminals" :
            subtype in _L3F_DELTA_SUBTYPES ?
                "a Kron-reduced $subtype bank requires three retained phase terminals " *
                "on each side; the wye star point is the eliminated conductor" :
                "single-phase transformer/regulator support requires one retained power channel")
        lo_key, hi_key = subtype == "single_phase" ? ("tap_min", "tap_max") :
                                                    ("tap_ratio_min", "tap_ratio_max")
        if haskey(data, lo_key) || haskey(data, hi_key)
            raw_lo, raw_hi = get(data, lo_key, NaN), get(data, hi_key, NaN)
            lo = raw_lo isa AbstractVector ? Float64.(raw_lo) : [Float64(raw_lo)]
            hi = raw_hi isa AbstractVector ? Float64.(raw_hi) : [Float64(raw_hi)]
            (!all(isfinite, lo) || !all(isfinite, hi) || length(lo) != length(hi) ||
             any(lo .!= hi)) && _l3f_error!(findings,
                "E.L3F.ADJUSTABLE_TAP_UNSUPPORTED", :transformer, id,
                "LinDist3Flow accepts fixed regulator settings only; adjustable tap intervals are excluded")
        end
        try
            neff = _l3f_transformer_neff(subtype, data)
            neff_values = neff isa AbstractVector ? neff : [neff]
            expected_ratio_count = subtype == "open_delta_regulator" ? 2 : 1
            subtype in _L3F_DELTA_SUBTYPES && _l3f_delta_transformer_gain(subtype, data)
            length(neff_values) == expected_ratio_count || throw(DimensionMismatch(
                "expected $expected_ratio_count effective ratio(s)"))
            all(x -> isfinite(x) && x > 0, neff_values) ||
                throw(ArgumentError("effective ratio must be positive and finite"))
            T = _l3f_transformer_forward_map(subtype, data)
            all(isfinite, T) || throw(ArgumentError("voltage-gain matrix must be finite"))
        catch err
            _l3f_error!(findings, "E.L3F.TRANSFORMER_RATIO_INVALID", :transformer, id,
                        sprint(showerror, err))
        end
        for key in ("r_series_from", "x_series_from", "r_series_to", "x_series_to",
                    "g_no_load", "b_no_load", "r_neutral_from", "x_neutral_from",
                    "r_neutral_to", "x_neutral_to")
            abs(Float64(get(data, key, 0.0))) <= 1e-12 || _l3f_error!(findings,
                "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED", :transformer, id,
                "fixed-ratio support is ideal; nonzero '$key' must be represented separately")
        end
        # BMOPF types `i_max_from`/`i_max_to` as `number[]` (per conductor) and
        # `s_rating` as a scalar nameplate. Accept a bare scalar for the
        # single-conductor subtypes too, since that shape reaches the wild.
        for key in ("s_rating", "i_max_from", "i_max_to")
            haskey(data, key) || continue
            value = data[key]
            positive(x) = x isa Real && isfinite(x) && x > 0
            # Per-conductor, except for an open-delta bank whose two entries
            # are its two units rather than its three conductors.
            expected = subtype == "open_delta_regulator" ? 2 : length(mf)
            valid_limit = if key == "s_rating"
                positive(value)
            else
                (value isa AbstractVector && length(value) == expected &&
                 all(positive, value)) || (expected == 1 && positive(value))
            end
            valid_limit || _l3f_error!(findings,
                "E.L3F.LIMIT_INVALID", :transformer, id,
                key == "s_rating" ? "$key must be a positive finite scalar" :
                subtype == "open_delta_regulator" ?
                    "$key must contain two positive finite winding ratings" :
                expected == 1 ?
                    "$key must be a positive finite scalar or one-element vector" :
                    "$key must contain $expected positive finite per-conductor ratings")
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
    # A perfectly grounded terminal is an explicit return conductor whatever it
    # is named. Matching only on the label set would skip Kron reduction and
    # then surface as a confusing zero-reference-phasor error instead.
    any(bus isa AbstractDict &&
        (any(t -> string(t) in labels || lowercase(string(t)) == "n",
             get(bus, "terminal_names", String[])) ||
         !isempty(get(bus, "perfectly_grounded_terminals", String[])))
        for bus in values(get(net, "bus", Dict())))
end

function _l3f_has_matrix_key(data, prefixes)
    any(key -> any(prefix -> startswith(String(key), prefix), prefixes), keys(data))
end

# BMOPFTools selects an inline series-coefficient source from these canonical
# markers, not from an arbitrary sparse R/X entry.  A malformed inline matrix
# may still be reported by `_l3f_pattern_matrix`, but it must not silently
# displace a valid referenced linecode merely because it contains e.g.
# `R_series_2_2`.
_l3f_has_inline_z(data) = haskey(data, "R_series_1_1") ||
                          haskey(data, "X_series_1_1")

function _l3f_pattern_matrix(data::AbstractDict, real_prefix::String,
                             imag_prefix::String, n::Int, label::String)
    n > 0 || throw(ArgumentError("$label has no conductors"))
    pattern = r"^(R_series_|X_series_)(\d+)_(\d+)$"
    seen = false
    for key in keys(data)
        m = match(pattern, String(key))
        m === nothing && continue
        seen = true
        i, j = parse(Int, m.captures[2]), parse(Int, m.captures[3])
        min(i, j) >= 1 && max(i, j) <= n ||
            throw(DimensionMismatch("$label matrix key '$key' falls outside terminal indices 1:$n"))
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
    inline = _l3f_has_inline_z(line)
    if inline
        # BMOPFTools uses the inline absolute matrix as the active coefficient
        # source when both representations happen to be present.  Keep that
        # deterministic precedence here as well; the linecode is not merged.
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
        # Only phase-to-ground v_min/v_max are stamped by the L3F model.  All
        # other BMOPF voltage-bound flavours must remain visible as errors;
        # dropping them under `unsupported=:approximate` would silently relax
        # an engineering constraint.
        for key in ("vpn_min", "vpn_max", "vn_max",
                    "vpp_min", "vpp_max", "vpos_min", "vpos_max",
                    "vuf_max", "vneg_max", "vzero_max", "vm_unbalance_max",
                    "va_diff_min", "va_diff_max")
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
            "only grounded-wye, single-phase, and delta loads are supported")
        load_model = lowercase(String(get(load, "model", "constant_power")))
        load_model in ("constant_power", "constant_impedance", "zip") ||
            _l3f_error!(findings, "E.L3F.LOAD_MODEL_UNSUPPORTED", :load, lid,
                        "only constant-power, constant-impedance, and pure ZP ZIP loads are supported")
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
            nch = length(p)
            _l3f_connection_incidence(cfg, length(tm), nch)
            values = vcat(Float64.(p), Float64.(q))
            all(isfinite, values) || throw(ArgumentError("nominal powers must be finite"))

            if load_model in ("constant_impedance", "zip")
                vnom_raw = get(load, "v_nom", nothing)
                vnom_raw === nothing && throw(ArgumentError(
                    "model=$load_model requires v_nom"))
                vnom = vnom_raw isa AbstractVector ? Float64.(vnom_raw) : [Float64(vnom_raw)]
                length(vnom) in (1, nch) || throw(ArgumentError(
                    "v_nom must be scalar or have one entry per load channel"))
                all(x -> isfinite(x) && x > 0.0, vnom) || throw(ArgumentError(
                    "v_nom entries must be positive and finite"))
            end

            if load_model == "zip"
                for field in ("alpha_z", "alpha_i", "alpha_p",
                              "beta_z", "beta_i", "beta_p")
                    raw = get(load, field, nothing)
                    raw === nothing && continue
                    coeff = raw isa AbstractVector ? Float64.(raw) : [Float64(raw)]
                    length(coeff) in (1, nch) || throw(ArgumentError(
                        "$field must be scalar or have one entry per load channel"))
                    all(isfinite, coeff) || throw(ArgumentError(
                        "$field entries must be finite"))
                end
                for field in ("alpha_i", "beta_i")
                    raw = get(load, field, nothing)
                    raw === nothing && continue
                    coeff = raw isa AbstractVector ? Float64.(raw) : [Float64(raw)]
                    any(!iszero, coeff) && _l3f_error!(findings,
                        "E.L3F.ZIP_CURRENT_UNSUPPORTED", :load, lid,
                        "$field must be identically zero; constant-current ZIP terms are not affine in squared voltage")
                end
            end
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

_l3f_is_dispatchable(component) = begin
    lo = Float64.(get(component, "p_min", Float64[]))
    hi = Float64.(get(component, "p_max", Float64[]))
    length(lo) == length(hi) && any(hi .> lo)
end

"""
Objective-dependent data checks.

BMOPF requires `cost` on a generator. Defaulting a missing vector to zero
   leaves the objective flat in that unit's dispatch direction.
"""
function _l3f_validate_objective!(findings, net, options::L3FOptions)
    options.objective == :cost || return
    for (family, symbol) in (("generator", :generator), ("voltage_source", :voltage_source))
        for (id, component) in get(net, family, Dict())
            component isa AbstractDict || continue
            cost = get(component, "cost", nothing)
            if cost === nothing
                (family == "voltage_source" || _l3f_is_dispatchable(component)) &&
                    _l3f_warning!(findings, "W.L3F.COST_MISSING", symbol, id,
                        "objective=:cost but no 'cost' vector is declared; this " *
                        "unit is priced at zero and the optimum may be non-unique")
                continue
            end
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
        lcid = get(line, "linecode", nothing)
        if _l3f_has_inline_z(line) && lcid isa AbstractString
            _l3f_error!(findings, "E.L3F.LINE_IMPEDANCE_SOURCE", :line, id,
                "line declares both an inline series-impedance matrix and " *
                "linecode '$lcid'; BMOPF requires exactly one impedance source")
        end
        # After recording the integrity error above, keep deterministic inline
        # precedence for downstream diagnostics. Never merge an inline source
        # with coefficients from the referenced linecode.
        coefficient_source = if _l3f_has_inline_z(line)
            line
        elseif lcid isa AbstractString && haskey(linecodes, lcid)
            linecodes[lcid]
        else
            nothing
        end
        if coefficient_source isa AbstractDict && _l3f_has_matrix_key(coefficient_source,
                    ("G_from_", "B_from_", "G_to_", "B_to_"))
            _l3f_error!(findings, "E.L3F.LINE_SHUNT_UNSUPPORTED", :line, id,
                        "line shunts require the Phase 3 affine shunt kernel")
        end
        for field in ("va_diff_min", "va_diff_max")
            haskey(line, field) || continue
            _l3f_error!(findings, "E.L3F.LIMIT_UNSUPPORTED", :line, id,
                "line bound '$field' is not assessed by the minimal L3F implementation")
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

"""
Two-port devices as `(family, subtype, id, bus_from, map_from, bus_to, map_to)`,
in a deterministic order. Both `line` and every supported `transformer` subtype
carry the same aligned conductor maps, so radiality, parallelism, and
orientation are all decided on the terminal graph rather than the bus graph.
"""
function _l3f_branches(net)
    out = Tuple{Symbol,String,String,String,Vector{String},String,Vector{String}}[]
    buses = get(net, "bus", Dict())
    add(family, subtype, id, data) = begin
        a, b = String(get(data, "bus_from", "")), String(get(data, "bus_to", ""))
        haskey(buses, a) && haskey(buses, b) || return
        ma = string.(get(data, "terminal_map_from", String[]))
        mb = string.(get(data, "terminal_map_to", String[]))
        length(ma) == length(mb) && !isempty(ma) || return
        push!(out, (family, subtype, id, a, ma, b, mb))
    end
    for (id, line) in get(net, "line", Dict())
        line isa AbstractDict && add(:line, "", String(id), line)
    end
    for (subtype, id, transformer) in _l3f_transformers(net)
        add(:transformer, subtype, id, transformer)
    end
    sort!(out; by=x -> (x[1], x[2], x[3]))
end

_l3f_branch_label(family, subtype, id) =
    family == :line ? "line/$id" : "transformer/$subtype/$id"

"""
    _l3f_topology!(findings, net)

Orient every two-port device away from its island's unique voltage source.

Radiality is a **per-conductor** property: the graph whose nodes are
`(bus, terminal)` pairs and whose edges are the aligned conductor pairs of each
device must be a forest. A three-unit single-phase regulator bank on one bus
pair is therefore radial — the units occupy disjoint conductors — even though
the bus graph shows three parallel edges. Islands and source assignment stay at
bus granularity, which is the conservative choice for a single-source model.
"""
function _l3f_topology!(findings, net)
    buses = sort!(String.(collect(keys(get(net, "bus", Dict())))))
    branches = _l3f_branches(net)

    # Terminal graph: nodes are (bus, terminal), edges are aligned conductors.
    terminal_adjacency = Dict{Tuple{String,String},Vector{Tuple{Tuple{String,String},Int}}}()
    terminal_pairs = Dict{Tuple{Tuple{String,String},Tuple{String,String}},Vector{String}}()
    bus_adjacency = Dict(bus => Tuple{String,Int}[] for bus in buses)
    for (e, (family, subtype, id, bf, mf, bt, mt)) in enumerate(branches)
        label = _l3f_branch_label(family, subtype, id)
        push!(bus_adjacency[bf], (bt, e)); push!(bus_adjacency[bt], (bf, e))
        for k in eachindex(mf)
            u, v = (bf, mf[k]), (bt, mt[k])
            push!(get!(terminal_adjacency, u, valtype(terminal_adjacency)()), (v, e))
            push!(get!(terminal_adjacency, v, valtype(terminal_adjacency)()), (u, e))
            key = u <= v ? (u, v) : (v, u)
            push!(get!(terminal_pairs, key, String[]), label)
        end
    end
    for (key, labels) in sort!(collect(terminal_pairs); by=first)
        length(labels) <= 1 || _l3f_error!(findings, "E.L3F.TOPOLOGY_NOT_RADIAL",
            :network, nothing,
            "parallel conductors connect $(key[1][1]).$(key[1][2]) and " *
            "$(key[2][1]).$(key[2][2])";
            evidence=Dict("branches" => sort(labels)))
    end

    # Per-conductor radiality: every terminal component must be a tree.
    terminal_nodes = sort!(collect(keys(terminal_adjacency)))
    unseen_terminals = Set(terminal_nodes)
    while !isempty(unseen_terminals)
        start = minimum(unseen_terminals)
        queue = [start]; delete!(unseen_terminals, start)
        nodes = 0; incidences = 0
        while !isempty(queue)
            node = popfirst!(queue); nodes += 1
            incidences += length(terminal_adjacency[node])
            for (neighbor, _) in terminal_adjacency[node]
                neighbor in unseen_terminals || continue
                delete!(unseen_terminals, neighbor); push!(queue, neighbor)
            end
        end
        edges = incidences ÷ 2
        edges <= nodes - 1 || _l3f_error!(findings, "E.L3F.TOPOLOGY_NOT_RADIAL",
            :network, nothing,
            "conductor component containing $(start[1]).$(start[2]) has $nodes " *
            "terminals and $edges conductors, so it is not radial";
            evidence=Dict("terminals" => nodes, "conductors" => edges))
    end

    islands = Vector{Vector{String}}()
    unseen = Set(buses)
    while !isempty(unseen)
        start = minimum(unseen)
        queue = [start]; delete!(unseen, start); component = String[]
        while !isempty(queue)
            node = popfirst!(queue); push!(component, node)
            for (neighbor, _) in bus_adjacency[node]
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
    any(f -> f.code == "E.L3F.TOPOLOGY_NOT_RADIAL", findings) &&
        return topology, roots, islands
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
        # Breadth-first over conductors: a device is oriented by the first of its
        # conductors that is traversed. In a forest every conductor of a device
        # agrees, because a disagreeing conductor would close a cycle.
        oriented = Dict{Int,Bool}()
        visited = Set((root, terminal)
                      for terminal in string.(get(net["bus"][root], "terminal_names", String[])))
        queue = sort!(collect(visited))
        while !isempty(queue)
            node = popfirst!(queue)
            for (neighbor, e) in sort(get(terminal_adjacency, node, []); by=x -> (x[1], x[2]))
                neighbor in visited && continue
                push!(visited, neighbor); push!(queue, neighbor)
                get!(oriented, e, branches[e][4] == node[1])
            end
        end
        for e in sort!(collect(keys(oriented)))
            family, subtype, id, bf, mf, bt, mt = branches[e]
            original = oriented[e]
            push!(topology, L3FOrientedLine(id, family, subtype,
                original ? bf : bt, original ? bt : bf,
                original ? mf : mt, original ? mt : mf, !original))
        end
    end
    topology, roots, islands
end

"""
Gate Yd/Dy banks on which winding faces the source.

The delta incidence has rank 2, so `D v_delta = g v_wye` determines the wye
voltages from the delta ones but not the reverse: a delta winding neither
imposes nor carries a zero-sequence terminal voltage. With the delta winding
upstream the ideal connection map is directly determined. With the wye winding
upstream the delta bus's common
voltage is genuinely undetermined by the transformer, and only a gauge can close
it — so that orientation is an error unless `unsupported=:approximate` accepts
the zero-zero-sequence assumption.
"""
function _l3f_validate_delta_orientation!(findings, net, topology, options::L3FOptions)
    for edge in topology
        edge.family == :transformer || continue
        edge.subtype in _L3F_DELTA_SUBTYPES || continue
        data = net["transformer"][edge.subtype][edge.id]
        delta_parent = try
            _l3f_delta_parent(edge, data)
        catch
            continue                      # ratio already reported as invalid
        end
        delta_parent === true && continue
        wye_bus, delta_bus = edge.parent, edge.child
        if options.unsupported == :approximate
            _l3f_warning!(findings, "A.L3F.DELTA_ZERO_SEQUENCE_GAUGE", :transformer,
                edge.id,
                "the wye winding faces the source, so the delta terminals at bus " *
                "'$delta_bus' are determined only up to a common offset; the " *
                "zero-sequence voltage there is assumed zero. A different ground " *
                "reference at that bus would give a different answer";
                evidence=Dict("wye_bus" => wye_bus, "delta_bus" => delta_bus))
        else
            _l3f_error!(findings, "E.L3F.DELTA_ORIENTATION_UNSUPPORTED", :transformer,
                edge.id,
                "the wye winding faces the source at bus '$wye_bus', so the delta " *
                "terminals at '$delta_bus' are determined only up to a common " *
                "offset — a delta winding neither imposes nor carries a " *
                "zero-sequence terminal voltage. Put the delta winding upstream, " *
                "or accept the zero-zero-sequence gauge with " *
                "L3FOptions(unsupported=:approximate)";
                evidence=Dict("wye_bus" => wye_bus, "delta_bus" => delta_bus))
        end
    end
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
    # `kron_reduce_bmopf` eliminates the conductors it recognises as neutrals.
    # A perfectly grounded terminal under any other label survives, and the model
    # would then meet it as a zero reference phasor — a true but unhelpful
    # diagnostic. Name the actual obstacle instead.
    for (bid, bus) in get(net, "bus", Dict())
        bus isa AbstractDict || continue
        retained = intersect(string.(get(bus, "perfectly_grounded_terminals", String[])),
                             string.(get(bus, "terminal_names", String[])))
        isempty(retained) && continue
        _l3f_error!(findings, "E.L3F.GROUNDED_TERMINAL_RETAINED", :bus, bid,
            "terminal(s) $(join(sort(retained), ", ")) are perfectly grounded but " *
            "still retained; LinDist3Flow models retained conductors only. Eliminate " *
            "them first — `kron_reduce_bmopf` removes conductors declared as neutrals " *
            "via `terminal_conventions` or `neutral_terminal`";
            evidence=Dict("terminals" => sort(retained)))
    end
    if options.require_neutral_provenance &&
       !haskey(get(net, "_meta", Dict()), "kron_reduction") && reference === nothing
        _l3f_error!(findings, "E.L3F.NEUTRAL_REDUCTION_UNDECLARED", :network, nothing,
            "neutral-reduction provenance or an explicit reference is required by options")
    end
    lowered = _l3f_lower!(findings, net, options)
    options.reference_policy == :explicit && reference === nothing && _l3f_error!(findings,
        "E.L3F.REFERENCE_MISSING", :network, nothing,
        "reference_policy=:explicit requires a reference argument")
    _l3f_validate_components!(findings, net)
    _l3f_validate_shunts!(findings, net)
    _l3f_validate_transformers!(findings, net)
    _l3f_validate_lines!(findings, net)
    _l3f_validate_objective!(findings, net, options)
    topology, roots, islands = _l3f_topology!(findings, net)
    _l3f_validate_delta_orientation!(findings, net, topology, options)
    report = L3FApplicabilityReport(
        any(f -> f.severity == :error, findings) ? :inapplicable : :applicable,
        findings, roots, islands, reduced, lowered)
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
