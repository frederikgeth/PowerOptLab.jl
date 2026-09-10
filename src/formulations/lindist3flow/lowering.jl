# Widening the admissible input without widening the formulation.
#
# The model itself only knows about a small vocabulary: series lines, fixed
# shunts, ideal fixed-ratio transformers, and P/Z loads. A lowering is canonical
# only when it stays within the same LinDist3Flow approximation (fixed phase
# angles and omitted series losses). It is never an exact AC re-representation.
# A smaller set genuinely cannot be represented and can only be projected onto
# it at a cost.
#
# `L3FOptions(unsupported=...)` selects how far to go:
#
#   :reject       nothing is rewritten (the default; what the formulation
#                 promises when it says a component is unsupported)
#   :lower        canonical re-representations only, reported as `L.L3F.*` at
#                 severity :info — the same L3F approximation is retained
#   :approximate  also the lossy projections, reported as `A.L3F.*` at severity
#                 :warning — the solved model is a different problem
#   :permissive   also removes otherwise-valid operational constraints that
#                 have no S-W representation, always with original-value evidence
#
# Every rewrite lands in the applicability report, so a result can always be
# traced back to what was actually solved. Nothing here is silent.

_l3f_info!(out, code, component, id, message; evidence=Dict{String,Any}()) =
    push!(out, L3FFinding(code, :info, component,
        isnothing(id) ? nothing : String(id), message, Dict{String,Any}(evidence)))

"""Synthetic id for an element the lowering pass introduces."""
_l3f_derived(kind::AbstractString, parts...) =
    string("_l3f_", kind, "_", join(parts, "_"))

_l3f_table!(net, family) = get!(net, family, Dict{String,Any}())

function _l3f_scalar_impedance(data, r_key, x_key)
    r = Float64(get(data, r_key, 0.0))
    x = Float64(get(data, x_key, 0.0))
    (r, x, hypot(r, x) > 1e-12)
end

"""
Diagonal inline series-impedance fields for `n` conductors sharing one scalar
impedance, the shape a per-coil transformer leakage takes.
"""
function _l3f_diagonal_series(r::Real, x::Real, n::Int)
    out = Dict{String,Any}()
    for k in 1:n
        out["R_series_$(k)_$(k)"] = Float64(r)
        out["X_series_$(k)_$(k)"] = Float64(x)
    end
    out
end

# ---------------------------------------------------------------------------
# Canonical L3F re-representations
# ---------------------------------------------------------------------------

"""
Closed switches become zero-impedance lines; open switches are removed.

Both preserve the canonical L3F approximation. A closed ideal switch is a
branch with no series drop, which the line kernel already represents, and an open
switch carries no current at all. If
removing an open switch leaves a subnetwork with no source, the island check
reports that in the usual way rather than this pass guessing at intent.
"""
function _l3f_lower_switches!(findings, net)
    switches = get(net, "switch", nothing)
    switches isa AbstractDict && !isempty(switches) || return
    lines = _l3f_table!(net, "line")
    for (sid_raw, switch) in sort!(collect(switches); by=first)
        sid = String(sid_raw)
        switch isa AbstractDict || continue
        if Bool(get(switch, "open_switch", false))
            _l3f_info!(findings, "L.L3F.SWITCH_OPEN_REMOVED", :switch, sid,
                "an open switch carries no current and is removed; any subnetwork " *
                "it alone energized is now a separate island")
            continue
        end
        mf = string.(get(switch, "terminal_map_from", String[]))
        id = _l3f_derived("switch", sid)
        line = Dict{String,Any}(
            "bus_from" => String(get(switch, "bus_from", "")),
            "bus_to" => String(get(switch, "bus_to", "")),
            "terminal_map_from" => mf,
            "terminal_map_to" => string.(get(switch, "terminal_map_to", String[])))
        merge!(line, _l3f_diagonal_series(0.0, 0.0, length(mf)))
        haskey(switch, "i_max") && (line["i_max"] = copy(switch["i_max"]))
        haskey(switch, "s_max") && (line["s_max"] = copy(switch["s_max"]))
        lines[id] = line
        _l3f_info!(findings, "L.L3F.SWITCH_LOWERED", :switch, sid,
            "closed switch represented in the canonical L3F model as the " *
            "zero-impedance line '$id'";
            evidence=Dict("line" => id))
    end
    delete!(net, "switch")
end

"""
Fixed capacitor banks become fixed shunt susceptance.

`B = q_rated / v_nom^2` per coil. For a delta bank the coil admittances are
referred to the terminals by ``D^{T}\\operatorname{diag}(b)D``. BMOPF capacitors
carry no switching state, so this preserves the canonical L3F component model;
it does not make the surrounding LinDist3Flow equations exact AC physics.
"""
function _l3f_lower_capacitors!(findings, net)
    capacitors = get(net, "capacitor", nothing)
    capacitors isa AbstractDict && !isempty(capacitors) || return
    shunts = _l3f_table!(net, "shunt")
    for (cid_raw, capacitor) in sort!(collect(capacitors); by=first)
        cid = String(cid_raw)
        capacitor isa AbstractDict || continue
        try
            tm = string.(capacitor["terminal_map"])
            q = Float64.(capacitor["q_rated"])
            v_nom = Float64(capacitor["v_nom"])
            isfinite(v_nom) && v_nom > 0 ||
                throw(ArgumentError("v_nom must be positive and finite"))
            all(isfinite, q) || throw(ArgumentError("q_rated must be finite"))
            cfg = uppercase(String(get(capacitor, "configuration", "WYE")))
            D = _l3f_connection_incidence(cfg, length(tm), length(q))
            Y = transpose(D) * Diagonal(q ./ v_nom^2) * D
            id = _l3f_derived("capacitor", cid)
            shunt = Dict{String,Any}("bus" => String(capacitor["bus"]),
                                     "terminal_map" => tm)
            for i in eachindex(tm), j in i:length(tm)
                iszero(Y[i, j]) && continue
                shunt["B_$(i)_$(j)"] = Y[i, j]
            end
            shunts[id] = shunt
            _l3f_info!(findings, "L.L3F.CAPACITOR_LOWERED", :capacitor, cid,
                "fixed capacitor represented in the canonical L3F model as the " *
                "shunt '$id' with " *
                "B = q_rated / v_nom^2";
                evidence=Dict("shunt" => id))
        catch err
            _l3f_error!(findings, "E.L3F.CAPACITOR_INVALID", :capacitor, cid,
                        sprint(showerror, err))
        end
    end
    delete!(net, "capacitor")
end

"""
Inline line shunt admittance becomes a fixed shunt at each end.

BMOPF already splits a line's shunt into declared from- and to-side halves, so
moving each half onto its own bus preserves the canonical pi data. Linecode entries are per unit length and are
scaled by the line's `length`, matching the series path.
"""
function _l3f_lower_line_shunts!(findings, net)
    lines = get(net, "line", Dict())
    linecodes = get(net, "linecode", Dict())
    shunts = _l3f_table!(net, "shunt")
    for (lid_raw, line) in sort!(collect(lines); by=first)
        lid = String(lid_raw)
        line isa AbstractDict || continue
        # BMOPFTools selects one complete coefficient source.  Inline series
        # matrices select inline shunts as well; do not merge them with a
        # referenced linecode.  A linecode supplies both kinds of coefficients
        # and is scaled by the line length.
        inline = _l3f_has_inline_z(line)
        lcid = get(line, "linecode", nothing)
        source = if inline
            (line, 1.0)
        elseif lcid isa AbstractString && haskey(linecodes, lcid)
            (linecodes[lcid], Float64(get(line, "length", 1.0)))
        else
            nothing
        end
        for (side, bus_key, map_key) in (("from", "bus_from", "terminal_map_from"),
                                         ("to", "bus_to", "terminal_map_to"))
            tm = string.(get(line, map_key, String[]))
            isempty(tm) && continue
            source === nothing && continue
            data, scale = source
            n = length(tm)
            bad_keys = String[]
            for key in keys(data)
                m = match(Regex("^[GB]_$(side)_(\\d+)_(\\d+)\$"), String(key))
                m === nothing && continue
                i, j = parse(Int, m.captures[1]), parse(Int, m.captures[2])
                min(i, j) >= 1 && max(i, j) <= n ||
                    push!(bad_keys, String(key))
            end
            isempty(bad_keys) || _l3f_error!(findings, "E.L3F.LINE_SHUNT_INVALID",
                :line, lid, "line shunt key(s) $(join(sort(bad_keys), ", ")) " *
                "fall outside terminal indices 1:$n")
            # Match BMOPFTools' full-vs-upper-triangular matrix semantics:
            # when both orientations are present, the direct key wins, rather
            # than summing the symmetric off-diagonal twice.
            matrix(prefix) = [begin
                key = "$(prefix)$(side)_$(i)_$(j)"
                reverse_key = "$(prefix)$(side)_$(j)_$(i)"
                value = haskey(data, key) ? data[key] :
                        (haskey(data, reverse_key) ? data[reverse_key] : 0.0)
                Float64(value) * scale
            end for i in 1:n, j in 1:n]
            G, B = matrix("G_"), matrix("B_")
            symmetric = all(isapprox(G[i, j], G[j, i]; atol=1e-12, rtol=1e-10) &&
                            isapprox(B[i, j], B[j, i]; atol=1e-12, rtol=1e-10)
                            for i in 1:n, j in 1:n)
            symmetric || _l3f_warning!(findings,
                "W.L3F.LINE_SHUNT_ASYMMETRIC", :line, lid,
                "$side-side line shunt admittance is asymmetric; it is preserved " *
                "to match BMOPF's direct-key semantics, but a passive reciprocal " *
                "line admittance is expected to be symmetric";
                evidence=Dict("side" => side))
            # Keep every reconstructed entry.  This preserves a deliberately
            # full (possibly asymmetric) input matrix, while a symmetric full
            # matrix is still represented once per matrix position and is not
            # accidentally doubled by the lowering.
            entries = Dict{Tuple{Int,Int},ComplexF64}((i, j) =>
                ComplexF64(G[i, j], B[i, j]) for i in 1:n, j in 1:n)
            filter!(pair_value -> !iszero(pair_value[2]), entries)
            isempty(entries) && continue
            id = _l3f_derived("lineshunt", lid, side)
            shunt = Dict{String,Any}("bus" => String(line[bus_key]), "terminal_map" => tm)
            for ((i, j), value) in entries
                iszero(real(value)) || (shunt["G_$(i)_$(j)"] = real(value))
                iszero(imag(value)) || (shunt["B_$(i)_$(j)"] = imag(value))
            end
            shunts[id] = shunt
            _l3f_info!(findings, "L.L3F.LINE_SHUNT_LOWERED", :line, lid,
                "$side-side line shunt moved onto bus $(line[bus_key]) as '$id'; " *
                "the declared pi halves are unchanged";
                evidence=Dict("shunt" => id, "side" => side))
        end
        for key in collect(keys(line))
            occursin(r"^[GB]_(from|to)_\d+_\d+$", String(key)) && delete!(line, key)
        end
    end
    # A linecode may be shared, so its shunt entries are cleared only after every
    # line that references it has taken its own copy.
    for (_, linecode) in linecodes
        linecode isa AbstractDict || continue
        for key in collect(keys(linecode))
            occursin(r"^[GB]_(from|to)_\d+_\d+$", String(key)) && delete!(linecode, key)
        end
    end
end

"""
A non-ideal single-phase or center-tap transformer becomes an ideal one plus
explicit series and shunt elements. For a center tap, one primary leakage arm
carries aggregate power and two identical secondary arms carry the individual
leg powers. This is the coupled three-winding star equivalent: the common
primary arm, not two independent transformers, retains the leg coupling.

This is a canonical L3F lowering: it preserves the fixed-angle/lossless
approximation, but it is not an exact AC multi-port. Yd/Dy banks use a
wye-referred short-circuit impedance, and autotransformers use their
from-referred regulating-winding impedance and from-side exciting branch.
"""
function _l3f_lower_transformer_impedance!(findings, net)
    tables = get(net, "transformer", Dict())
    tables isa AbstractDict || return
    buses = get(net, "bus", Dict())
    lines = _l3f_table!(net, "line")
    shunts = _l3f_table!(net, "shunt")
    for subtype in sort!(collect(String.(keys(tables))))
        table = tables[subtype]
        table isa AbstractDict || continue
        for (tid_raw, transformer) in sort!(collect(table); by=first)
            tid = String(tid_raw)
            transformer isa AbstractDict || continue
            loss_fields = ("r_series_from", "x_series_from", "r_series_to",
                "x_series_to", "g_no_load", "b_no_load")
            invalid_loss = false
            for field in loss_fields
                haskey(transformer, field) || continue
                value = transformer[field]
                if !(value isa Real && !(value isa Bool) && isfinite(value)) ||
                   (startswith(field, "r_") || field == "g_no_load") && value < 0
                    _l3f_error!(findings, "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED",
                        :transformer, tid,
                        "$field must be finite numeric data" *
                        ((startswith(field, "r_") || field == "g_no_load") ?
                         " and nonnegative" : ""))
                    invalid_loss = true
                end
            end
            invalid_loss && continue
            if subtype == "single_phase_autotransformer"
                fields = ("r_series_from", "x_series_from",
                          "r_series_to", "x_series_to")
                values = try
                    Float64[get(transformer, key, 0.0) for key in
                        (fields..., "g_no_load", "b_no_load")]
                catch err
                    _l3f_error!(findings, "E.L3F.TRANSFORMER_DATA_INVALID",
                        :transformer, tid, sprint(showerror, err)); continue
                end
                all(isfinite, values) || begin
                    _l3f_error!(findings, "E.L3F.TRANSFORMER_DATA_INVALID",
                        :transformer, tid, "leakage and exciting parameters must be finite"); continue
                end
                has_series = any(abs(values[k]) > 1e-12 for k in 1:4)
                has_shunt = hypot(values[5], values[6]) > 1e-12
                has_series || has_shunt || continue
                for side in ("from", "to")
                    get!(transformer, "_l3f_physical_bus_$side",
                         String(get(transformer, "bus_$side", "")))
                    get!(transformer, "_l3f_physical_terminal_map_$side",
                         string.(get(transformer, "terminal_map_$side", String[])))
                end
                if has_series
                    n = _l3f_transformer_neff(subtype, transformer)
                    z = complex(Float64(get(transformer, "r_series_from", 0.0)),
                                Float64(get(transformer, "x_series_from", 0.0))) +
                        n^2 * complex(Float64(get(transformer, "r_series_to", 0.0)),
                                      Float64(get(transformer, "x_series_to", 0.0)))
                    bus = String(transformer["bus_from"])
                    tm = string.(transformer["terminal_map_from"])
                    internal = _l3f_derived("xfmr", tid, "from")
                    buses[internal] = Dict{String,Any}("terminal_names" => copy(tm))
                    line_id = _l3f_derived("leakage", tid, "from")
                    line = Dict{String,Any}(
                        "bus_from" => bus, "bus_to" => internal,
                        "terminal_map_from" => copy(tm), "terminal_map_to" => copy(tm))
                    merge!(line, _l3f_diagonal_series(real(z), imag(z), length(tm)))
                    lines[line_id] = line
                    transformer["bus_from"] = internal
                    foreach(key -> delete!(transformer, key), fields)
                    _l3f_info!(findings, "L.L3F.TRANSFORMER_LEAKAGE_LOWERED",
                        :transformer, tid,
                        "autotransformer leakage referred to the physical from side " *
                        "and represented by series line '$line_id'";
                        evidence=Dict("line" => line_id, "side" => "from",
                                      "r" => real(z), "x" => imag(z)))
                end
                if has_shunt
                    bus = String(transformer["_l3f_physical_bus_from"])
                    tm = string.(transformer["_l3f_physical_terminal_map_from"])
                    id = _l3f_derived("noload", tid, "from")
                    shunt = Dict{String,Any}("bus" => bus, "terminal_map" => copy(tm))
                    ncoil = length(tm)
                    for k in 1:ncoil
                        shunt["G_$(k)_$(k)"] = Float64(get(transformer, "g_no_load", 0.0)) / ncoil
                        shunt["B_$(k)_$(k)"] = Float64(get(transformer, "b_no_load", 0.0)) / ncoil
                    end
                    shunts[id] = shunt
                    delete!(transformer, "g_no_load"); delete!(transformer, "b_no_load")
                    _l3f_info!(findings, "L.L3F.TRANSFORMER_NO_LOAD_LOWERED",
                        :transformer, tid,
                        "autotransformer exciting admittance represented at its " *
                        "external from terminal as shunt '$id'";
                        evidence=Dict("shunt" => id, "bus" => bus))
                end
                continue
            end
            if subtype in _L3F_DELTA_SUBTYPES
                fields = ("r_series_from", "x_series_from",
                          "r_series_to", "x_series_to")
                has_series = any(abs(Float64(get(transformer, key, 0.0))) > 1e-12
                                 for key in fields)
                has_shunt = hypot(Float64(get(transformer, "g_no_load", 0.0)),
                                  Float64(get(transformer, "b_no_load", 0.0))) > 1e-12
                has_series || has_shunt || continue
                for side in ("from", "to")
                    get!(transformer, "_l3f_physical_bus_$side",
                         String(get(transformer, "bus_$side", "")))
                    get!(transformer, "_l3f_physical_terminal_map_$side",
                         string.(get(transformer, "terminal_map_$side", String[])))
                end
                if has_series
                    wye_side = subtype == "wye_delta" ? "from" : "to"
                    wye_bus_key, wye_map_key = "bus_$wye_side", "terminal_map_$wye_side"
                    bus = String(transformer[wye_bus_key]); tm = string.(transformer[wye_map_key])
                    internal = _l3f_derived("xfmr", tid, wye_side)
                    buses[internal] = Dict{String,Any}("terminal_names" => copy(tm))
                    N0 = Float64(transformer["v_nom_from"]) /
                         Float64(transformer["v_nom_to"])
                    g0 = subtype == "wye_delta" ? sqrt(3.0) / N0 : sqrt(3.0) * N0
                    tap = if haskey(transformer, "tap")
                        Float64(transformer["tap"])
                    elseif haskey(transformer, "tap_min") && haskey(transformer, "tap_max") &&
                           Float64(transformer["tap_min"]) == Float64(transformer["tap_max"])
                        Float64(transformer["tap_min"])
                    else
                        1.0
                    end
                    Zw = subtype == "wye_delta" ?
                        complex(Float64(get(transformer, "r_series_from", 0.0)),
                                Float64(get(transformer, "x_series_from", 0.0))) :
                        complex(Float64(get(transformer, "r_series_to", 0.0)),
                                Float64(get(transformer, "x_series_to", 0.0)))
                    Zd = subtype == "wye_delta" ?
                        complex(Float64(get(transformer, "r_series_to", 0.0)),
                                Float64(get(transformer, "x_series_to", 0.0))) :
                        complex(Float64(get(transformer, "r_series_from", 0.0)),
                                Float64(get(transformer, "x_series_from", 0.0)))
                    zsc = (subtype == "wye_delta" ? tap^2 : 1.0) *
                          (Zw + 3 / g0^2 * Zd)
                    line_id = _l3f_derived("leakage", tid, wye_side)
                    line = Dict{String,Any}(
                        "bus_from" => (wye_side == "from" ? bus : internal),
                        "bus_to" => (wye_side == "from" ? internal : bus),
                        "terminal_map_from" => copy(tm), "terminal_map_to" => copy(tm))
                    merge!(line, _l3f_diagonal_series(real(zsc), imag(zsc), length(tm)))
                    lines[line_id] = line; transformer[wye_bus_key] = internal
                    foreach(key -> delete!(transformer, key), fields)
                    _l3f_info!(findings, "L.L3F.TRANSFORMER_LEAKAGE_LOWERED",
                        :transformer, tid,
                        "connection-aware Yd/Dy leakage referred to the wye winding " *
                        "and represented by series line '$line_id'";
                        evidence=Dict("line" => line_id, "side" => wye_side,
                                      "r" => real(zsc), "x" => imag(zsc)))
                end
                if has_shunt
                    bus = String(transformer["_l3f_physical_bus_to"])
                    tm = string.(transformer["_l3f_physical_terminal_map_to"])
                    n = length(tm); G = Float64(get(transformer, "g_no_load", 0.0)) / n
                    B = Float64(get(transformer, "b_no_load", 0.0)) / n
                    Y = subtype == "wye_delta" ? begin
                        D = _l3f_connection_incidence("DELTA", n, n)
                        transpose(D) * Diagonal(fill(complex(G, B), n)) * D
                    end : Diagonal(fill(complex(G, B), n))
                    id = _l3f_derived("noload", tid)
                    shunt = Dict{String,Any}("bus" => bus, "terminal_map" => copy(tm))
                    for i in 1:n, j in 1:n
                        iszero(real(Y[i,j])) || (shunt["G_$(i)_$(j)"] = real(Y[i,j]))
                        iszero(imag(Y[i,j])) || (shunt["B_$(i)_$(j)"] = imag(Y[i,j]))
                    end
                    shunts[id] = shunt
                    delete!(transformer, "g_no_load"); delete!(transformer, "b_no_load")
                    _l3f_info!(findings, "L.L3F.TRANSFORMER_NO_LOAD_LOWERED",
                        :transformer, tid,
                        "Yd/Dy winding-2 exciting admittance represented at its " *
                        "external terminal bus as shunt '$id'";
                        evidence=Dict("shunt" => id, "bus" => bus))
                end
                continue
            end
            # A center tap is a three-winding star: its one from-side arm is
            # common to both legs, while the two identical to-side star arms
            # are diagonal after the ideal 1->2 core map. Autotransformer,
            # Yd/Dy, and open-delta leakage cannot be decomposed this way.
            subtype in ("single_phase", "center_tap") || continue
            # Leakage lines introduce internal buses, but nameplate limits and
            # the exciting branch are defined at the caller's physical winding
            # endpoints. Preserve those identities before changing either side.
            get!(transformer, "_l3f_physical_bus_from",
                 String(get(transformer, "bus_from", "")))
            get!(transformer, "_l3f_physical_bus_to",
                 String(get(transformer, "bus_to", "")))
            get!(transformer, "_l3f_physical_terminal_map_from",
                 string.(get(transformer, "terminal_map_from", String[])))
            get!(transformer, "_l3f_physical_terminal_map_to",
                 string.(get(transformer, "terminal_map_to", String[])))
            for (side, bus_key, map_key, r_key, x_key) in (
                    ("from", "bus_from", "terminal_map_from", "r_series_from", "x_series_from"),
                    ("to", "bus_to", "terminal_map_to", "r_series_to", "x_series_to"))
                r, x, present = _l3f_scalar_impedance(transformer, r_key, x_key)
                present || continue
                bus = String(get(transformer, bus_key, ""))
                haskey(buses, bus) || continue
                tm = string.(get(transformer, map_key, String[]))
                internal = _l3f_derived("xfmr", tid, side)
                buses[internal] = Dict{String,Any}("terminal_names" => copy(tm))
                line_id = _l3f_derived("leakage", tid, side)
                line = Dict{String,Any}(
                    "bus_from" => bus, "terminal_map_from" => copy(tm),
                    "bus_to" => internal, "terminal_map_to" => copy(tm))
                merge!(line, _l3f_diagonal_series(r, x, length(tm)))
                lines[line_id] = line
                transformer[bus_key] = internal
                delete!(transformer, r_key); delete!(transformer, x_key)
                _l3f_info!(findings, "L.L3F.TRANSFORMER_LEAKAGE_LOWERED",
                    :transformer, tid,
                "$side-winding leakage $(r) + j$(x) ohm represented in the " *
                    "series line '$line_id' through the internal bus '$internal'";
                    evidence=Dict("line" => line_id, "bus" => internal,
                                  "side" => side, "r" => r, "x" => x))
            end
            g, b, present = _l3f_scalar_impedance(transformer, "g_no_load", "b_no_load")
            present || continue
            # BMOPFTools places the exciting branch across winding 2: the
            # ordinary transformer's to-side coil, or center-tap LV leg 1.
            # It is one total admittance and is never duplicated over both legs.
            bus = String(get(transformer, "_l3f_physical_bus_to",
                             get(transformer, "bus_to", "")))
            haskey(buses, bus) || continue
            tm = string.(get(transformer, "_l3f_physical_terminal_map_to",
                             get(transformer, "terminal_map_to", String[])))
            id = _l3f_derived("noload", tid)
            shunt = Dict{String,Any}("bus" => bus, "terminal_map" => copy(tm))
            shunt_positions = subtype == "center_tap" ? (1,) : eachindex(tm)
            for k in shunt_positions
                iszero(g) || (shunt["G_$(k)_$(k)"] = g)
                iszero(b) || (shunt["B_$(k)_$(k)"] = b)
            end
            shunts[id] = shunt
            delete!(transformer, "g_no_load"); delete!(transformer, "b_no_load")
            _l3f_info!(findings, "L.L3F.TRANSFORMER_NO_LOAD_LOWERED", :transformer, tid,
                "no-load admittance $(g) + j$(b) S represented in the " *
                "shunt '$id' across winding 2";
                evidence=Dict("shunt" => id, "g" => g, "b" => b))
        end
    end
end

# ---------------------------------------------------------------------------
# Lossy projections
# ---------------------------------------------------------------------------

"""
Tangent ZP fractions for a voltage exponent.

A term ``(V/V_{nom})^\\gamma`` is matched in value and first derivative at
``V=V_{nom}`` by ``\\alpha_P+\\alpha_Z(V/V_{nom})^2`` with
``\\alpha_Z=\\gamma/2`` and ``\\alpha_P=1-\\gamma/2``. The endpoint classes
``\\gamma=0`` (constant power) and ``\\gamma=2`` (constant impedance) coincide
with the native model classes; for ``\\gamma=1`` it is the familiar half-to-Z,
half-to-P split of a constant-current term. Other exponents are experimental
tangent projections.
"""
_l3f_zp_tangent(gamma::Real) = (gamma / 2, 1 - gamma / 2)

function _l3f_project_load!(findings, load, lid)
    model = lowercase(String(get(load, "model", "constant_power")))
    model in ("constant_current", "zip", "exponential") || return
    nch = length(get(load, "p_nom", Float64[]))
    channels(field, default) = begin
        raw = get(load, field, nothing)
        raw === nothing && return fill(Float64(default), nch)
        values = raw isa AbstractVector ? Float64.(raw) : [Float64(raw)]
        length(values) == 1 ? fill(values[1], nch) : values
    end

    if model == "constant_current"
        z, p = _l3f_zp_tangent(1.0)
        merge!(load, Dict{String,Any}("model" => "zip",
            "alpha_z" => fill(z, nch), "alpha_p" => fill(p, nch),
            "alpha_i" => zeros(nch),
            "beta_z" => fill(z, nch), "beta_p" => fill(p, nch),
            "beta_i" => zeros(nch)))
        _l3f_warning!(findings, "A.L3F.LOAD_LAW_PROJECTED", :load, lid,
            "constant-current law projected onto the ZP tangent at v_nom " *
            "(alpha_Z = beta_Z = 0.5, alpha_P = beta_P = 0.5); the two laws " *
            "agree in value and slope at nominal voltage and diverge away from it";
            evidence=Dict("model" => "constant_current", "alpha_z" => z, "alpha_p" => p))
        return
    end

    if model == "exponential"
        gp, gq = channels("gamma_p", 0.0), channels("gamma_q", 0.0)
        zp = _l3f_zp_tangent.(gp); zq = _l3f_zp_tangent.(gq)
        for key in ("gamma_p", "gamma_q")
            delete!(load, key)
        end
        merge!(load, Dict{String,Any}("model" => "zip",
            "alpha_z" => first.(zp), "alpha_p" => last.(zp), "alpha_i" => zeros(nch),
            "beta_z" => first.(zq), "beta_p" => last.(zq), "beta_i" => zeros(nch)))
        _l3f_warning!(findings, "A.L3F.LOAD_LAW_PROJECTED", :load, lid,
            "exponential law with exponents $(gp) / $(gq) projected onto the ZP " *
            "tangent at v_nom; endpoint classes coincide with native P/Z laws, " *
            "and other exponents are experimental";
            evidence=Dict("model" => "exponential", "gamma_p" => gp, "gamma_q" => gq))
        return
    end

    # ZIP with a non-zero current fraction: reassign that fraction to Z and P.
    # Each P/Q coefficient family defaults independently to constant power when
    # the whole family is absent, matching `_l3f_zip_coefficients`. Resolving
    # those defaults before merging prevents projection of one family from
    # erasing an omitted nominal contribution in the other.
    family(prefix) = begin
        fields = ("$(prefix)_z", "$(prefix)_i", "$(prefix)_p")
        all(field -> !haskey(load, field), fields) &&
            return (zeros(nch), zeros(nch), ones(nch))
        (channels(fields[1], 0.0), channels(fields[2], 0.0),
         channels(fields[3], 0.0))
    end
    az, ai, ap = family("alpha")
    bz, bi, bp = family("beta")
    all(iszero, ai) && all(iszero, bi) && return
    merge!(load, Dict{String,Any}(
        "alpha_z" => az .+ ai ./ 2, "alpha_p" => ap .+ ai ./ 2, "alpha_i" => zeros(nch),
        "beta_z" => bz .+ bi ./ 2, "beta_p" => bp .+ bi ./ 2, "beta_i" => zeros(nch)))
    _l3f_warning!(findings, "A.L3F.LOAD_LAW_PROJECTED", :load, lid,
        "constant-current ZIP fractions $(ai) / $(bi) split evenly onto Z and P, " *
        "the tangent at v_nom; the total fractions are preserved";
        evidence=Dict("model" => "zip", "alpha_i" => ai, "beta_i" => bi))
end

"""
Collapse an adjustable tap interval onto a single fixed setting.

Unlike the canonical lowerings above this removes a decision the caller asked for.
The declared operating tap is kept when it lies inside the interval, otherwise
the midpoint is used, and the discarded interval is recorded so the loss is
visible in the report.
"""
function _l3f_project_taps!(findings, net)
    for (subtype, tid, transformer) in _l3f_transformers(net)
        value_key, lo_key, hi_key = subtype in ("single_phase", "center_tap", _L3F_DELTA_SUBTYPES...,
                                                 "grounded_wye_wye", "delta_delta") ?
            ("tap", "tap_min", "tap_max") :
            ("tap_ratio", "tap_ratio_min", "tap_ratio_max")
        haskey(transformer, lo_key) || haskey(transformer, hi_key) || continue
        vector(x) = x isa AbstractVector ? Float64.(x) : [Float64(x)]
        lo, hi = try
            vector(get(transformer, lo_key, NaN)),
            vector(get(transformer, hi_key, NaN))
        catch err
            _l3f_error!(findings, "E.L3F.TAP_INTERVAL_INVALID", :transformer, tid,
                "tap interval entries must be numeric: $(sprint(showerror, err))")
            continue
        end
        expected = subtype == "open_delta_regulator" ? 2 :
                   subtype == "closed_delta_regulator" ? 3 : 1
        valid = length(lo) == expected && length(hi) == expected &&
                all(isfinite, lo) && all(isfinite, hi) &&
                all(lo .> 0) && all(hi .> 0) && all(lo .<= hi)
        if !valid
            _l3f_error!(findings, "E.L3F.TAP_INTERVAL_INVALID", :transformer, tid,
                "tap interval must contain $expected positive finite lower/upper " *
                "value(s) with minimum no greater than maximum";
                evidence=Dict("tap_min" => lo, "tap_max" => hi))
            continue
        end
        declared = try
            haskey(transformer, value_key) ? vector(transformer[value_key]) : Float64[]
        catch err
            _l3f_error!(findings, "E.L3F.TAP_INTERVAL_INVALID", :transformer, tid,
                "declared tap entries must be numeric: $(sprint(showerror, err))")
            continue
        end
        if all(lo .== hi)
            if !isempty(declared) &&
               (length(declared) != length(lo) || any(declared .!= lo))
                _l3f_error!(findings, "E.L3F.TAP_INTERVAL_INVALID", :transformer, tid,
                    "declared tap $(declared) contradicts the fixed interval $(lo)";
                    evidence=Dict("tap" => declared, "fixed" => lo))
            end
            continue
        end
        fixed = if length(declared) == length(lo) && all(lo .<= declared .<= hi)
            declared
        else
            (lo .+ hi) ./ 2
        end
        transformer[value_key] = length(fixed) == 1 && !(get(transformer, value_key, nothing) isa AbstractVector) ?
            fixed[1] : fixed
        delete!(transformer, lo_key); delete!(transformer, hi_key)
        _l3f_warning!(findings, "A.L3F.ADJUSTABLE_TAP_PROJECTED", :transformer, tid,
            "adjustable tap interval $(lo) to $(hi) collapsed to the fixed setting " *
            "$(fixed); the optimizer no longer selects the tap, so the solution is " *
            "feasible for that setting rather than optimal over the range";
            evidence=Dict("tap_min" => lo, "tap_max" => hi, "fixed" => fixed))
    end
end

function _l3f_permissive_numeric(value)
    values = value isa AbstractVector ? value : [value]
    !isempty(values) && all(x -> x isa Real && isfinite(x), values)
end

"""Drop only well-formed operational data outside the one-shot S-W vocabulary."""
function _l3f_permissive_constraints!(findings, net)
    for key in ("control_profile", "time_series")
        table = get(net, key, nothing)
        table isa AbstractDict && !isempty(table) || continue
        if !_l3f_permissive_wellformed(table)
            _l3f_error!(findings, "E.L3F.CONTROL_PROFILE_UNSUPPORTED", :network,
                nothing, "top-level '$key' contains malformed or non-finite metadata")
            continue
        end
        original = deepcopy(table)
        empty!(table)
        _l3f_warning!(findings, "A.L3F.PERMISSIVE_METADATA_DROPPED", :network,
            nothing, "unused top-level '$key' metadata was dropped";
            evidence=Dict("field" => key, "original_value" => original))
    end
    bus_fields = ("vpn_min", "vpn_max", "vpos_min", "vpos_max",
        "vuf_max", "vneg_max", "vzero_max", "vm_unbalance_max",
        "va_diff_min", "va_diff_max")
    for (bid, bus) in get(net, "bus", Dict())
        bus isa AbstractDict || continue
        for (lo, hi) in (("vpn_min", "vpn_max"), ("vpos_min", "vpos_max"),
                         ("va_diff_min", "va_diff_max"))
            haskey(bus, lo) && haskey(bus, hi) || continue
            lv, hv = bus[lo], bus[hi]
            if _l3f_permissive_numeric(lv) && _l3f_permissive_numeric(hv)
                lvs = lv isa AbstractVector ? Float64.(lv) : [Float64(lv)]
                hvs = hv isa AbstractVector ? Float64.(hv) : [Float64(hv)]
                (length(lvs) in (1, length(hvs)) || length(hvs) == 1) &&
                    any(lvs .> hvs) && _l3f_error!(findings,
                        "E.L3F.VOLTAGE_BOUND_INVALID", :bus, bid,
                        "$lo exceeds $hi before permissive removal")
            end
        end
        for field in bus_fields
            haskey(bus, field) || continue
            value = bus[field]
            if !_l3f_permissive_numeric(value)
                _l3f_error!(findings, "E.L3F.VOLTAGE_BOUND_INVALID", :bus, bid,
                    "$field must contain finite numeric values before it can be dropped")
                continue
            end
            !startswith(field, "va_diff") && any(Float64.(value isa AbstractVector ? value : [value]) .< 0) && begin
                _l3f_error!(findings, "E.L3F.VOLTAGE_BOUND_INVALID", :bus, bid,
                    "$field must be nonnegative before it can be dropped"); continue
            end
            delete!(bus, field)
            _l3f_warning!(findings, "A.L3F.PERMISSIVE_VOLTAGE_LIMIT_DROPPED",
                :bus, bid, "unsupported bus voltage constraint '$field' was dropped";
                evidence=Dict("field" => field, "original_value" => deepcopy(value)))
        end
    end
    for (lid, line) in get(net, "line", Dict())
        line isa AbstractDict || continue
        for field in ("va_diff_min", "va_diff_max")
            haskey(line, field) || continue
            value = line[field]
            if !_l3f_permissive_numeric(value)
                _l3f_error!(findings, "E.L3F.LIMIT_INVALID", :line, lid,
                    "$field must contain finite numeric values before it can be dropped")
                continue
            end
            delete!(line, field)
            _l3f_warning!(findings, "A.L3F.PERMISSIVE_LINE_LIMIT_DROPPED",
                :line, lid, "unsupported line angle constraint '$field' was dropped";
                evidence=Dict("field" => field, "original_value" => deepcopy(value)))
        end
    end

    local_tables = get(net, "transformer", Dict())
    for subtype in _L3F_LOCAL_BANK_SUBTYPES
        for (tid, tx) in get(local_tables, subtype, Dict())
            tx isa AbstractDict || continue
            for field in ("r_series_from", "x_series_from", "r_series_to",
                          "x_series_to", "g_no_load", "b_no_load")
                haskey(tx, field) || continue
                value = tx[field]
                if !(value isa Real && isfinite(value))
                    _l3f_error!(findings, "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED",
                        :transformer, tid,
                        "$field must be finite numeric data before idealization")
                    continue
                end
                delete!(tx, field)
                iszero(Float64(value)) || _l3f_warning!(findings,
                    "A.L3F.PERMISSIVE_TRANSFORMER_IDEALIZED", :transformer, tid,
                    "local bank '$subtype' dropped nonideal field '$field'";
                    evidence=Dict("subtype" => subtype, "field" => field,
                                  "original_value" => deepcopy(value)))
            end
            if haskey(tx, "s_rating") && subtype in ("delta_delta", "closed_delta_regulator")
                value = tx["s_rating"]
                if !(value isa Real && isfinite(value) && value > 0)
                    _l3f_error!(findings, "E.L3F.LIMIT_INVALID", :transformer, tid,
                        "s_rating must be positive finite data before it can be dropped")
                else
                    delete!(tx, "s_rating")
                    _l3f_warning!(findings,
                        "A.L3F.PERMISSIVE_TRANSFORMER_LIMIT_DROPPED", :transformer, tid,
                        "local bank '$subtype' dropped s_rating without a defined coil interpretation";
                        evidence=Dict("subtype" => subtype, "field" => "s_rating",
                                      "original_value" => value))
                end
            end
            if haskey(tx, "no_load_shunt")
                value = tx["no_load_shunt"]
                winding = value isa AbstractDict ? get(value, "winding", 0) : 0
                g = value isa AbstractDict ? get(value, "g", nothing) : nothing
                b = value isa AbstractDict ? get(value, "b", nothing) : nothing
                valid = value isa AbstractDict && winding isa Integer && !(winding isa Bool) &&
                    Int(winding) in (1, 2, 3) && g isa Real && !(g isa Bool) &&
                    isfinite(g) && g >= 0 && b isa Real && !(b isa Bool) && isfinite(b) &&
                    !any(haskey(tx, key) for key in ("g_no_load", "b_no_load"))
                if !valid
                    _l3f_error!(findings, "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED",
                        :transformer, tid, "no_load_shunt must be an object before idealization")
                else
                    delete!(tx, "no_load_shunt")
                    _l3f_warning!(findings, "A.L3F.PERMISSIVE_TRANSFORMER_IDEALIZED",
                        :transformer, tid, "local bank '$subtype' dropped nested no_load_shunt";
                        evidence=Dict("subtype" => subtype, "field" => "no_load_shunt",
                                      "original_value" => deepcopy(value)))
                end
            end
        end
    end

    records = get(get(net, "_meta", Dict()), "explicit_transformer_core_shunts", Dict())
    records isa AbstractDict || return net
    shunts = get(net, "shunt", Dict())
    for (record_id, record) in collect(records)
        record isa AbstractDict || continue
        subtype = String(get(record, "subtype", ""))
        subtype in _L3F_LOCAL_BANK_SUBTYPES || continue
        shunt_id = String(get(record, "shunt_id", ""))
        haskey(shunts, shunt_id) || continue
        original = deepcopy(shunts[shunt_id])
        delete!(shunts, shunt_id); delete!(records, record_id)
        _l3f_warning!(findings, "A.L3F.PERMISSIVE_TRANSFORMER_IDEALIZED",
            :transformer, get(record, "transformer_id", nothing),
            "parser-materialized local-bank core shunt '$shunt_id' was removed";
            evidence=Dict("subtype" => subtype, "shunt_id" => shunt_id,
                          "original_value" => original))
    end
    net
end

function _l3f_permissive_pre_kron_metadata!(findings, net)
    tables = [(family, get(net, family, Dict())) for family in
        ("bus", "line", "load", "generator", "voltage_source", "shunt",
         "switch", "capacitor", "ibr")]
    for (subtype, table) in get(net, "transformer", Dict())
        push!(tables, ("transformer", table))
    end
    for (family, table) in tables, (id, component) in table
        component isa AbstractDict && haskey(component, "time_series") || continue
        value = component["time_series"]
        if !(value isa AbstractDict || value isa AbstractString) ||
           !_l3f_permissive_wellformed(value)
            _l3f_error!(findings, "E.L3F.TIME_SERIES_UNSUPPORTED",
                Symbol(family), id, "time_series metadata must be an object or profile id")
            continue
        end
        delete!(component, "time_series")
        _l3f_warning!(findings, "A.L3F.PERMISSIVE_METADATA_DROPPED",
            Symbol(family), id,
            "time_series metadata was dropped; declared static snapshot values are used";
            evidence=Dict("field" => "time_series", "original_value" => deepcopy(value)))
    end
    net
end

"""Capture neutral-voltage limits before Kron erases their terminal."""
function _l3f_permissive_pre_kron_limits!(findings, net)
    for (bid, bus) in get(net, "bus", Dict())
        bus isa AbstractDict && haskey(bus, "vn_max") || continue
        value = bus["vn_max"]
        if !_l3f_permissive_numeric(value) ||
           any(Float64.(value isa AbstractVector ? value : [value]) .< 0)
            _l3f_error!(findings, "E.L3F.VOLTAGE_BOUND_INVALID", :bus, bid,
                "vn_max must contain nonnegative finite numeric values")
            continue
        end
        delete!(bus, "vn_max")
        _l3f_warning!(findings, "A.L3F.PERMISSIVE_VOLTAGE_LIMIT_DROPPED",
            :bus, bid, "neutral-voltage constraint 'vn_max' was dropped before reduction";
            evidence=Dict("field" => "vn_max", "original_value" => deepcopy(value)))
    end
    net
end

function _l3f_permissive_transformer_residuals!(findings, net)
    supported = ("single_phase", "center_tap", "single_phase_autotransformer",
        "open_delta_regulator", _L3F_DELTA_SUBTYPES..., _L3F_LOCAL_BANK_SUBTYPES...)
    fields = ("r_series_from", "x_series_from", "r_series_to", "x_series_to",
        "g_no_load", "b_no_load", "r_neutral_from", "x_neutral_from",
        "r_neutral_to", "x_neutral_to")
    for subtype in supported, (tid, tx) in get(get(net, "transformer", Dict()), subtype, Dict())
        tx isa AbstractDict || continue
        for field in fields
            haskey(tx, field) || continue
            value = tx[field]
            value isa Real && isfinite(value) || begin
                _l3f_error!(findings, "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED",
                    :transformer, tid, "$field must be finite numeric data before idealization")
                continue
            end
            delete!(tx, field)
            iszero(value) || _l3f_warning!(findings,
                "A.L3F.PERMISSIVE_TRANSFORMER_IDEALIZED", :transformer, tid,
                "supported transformer '$subtype' dropped residual nonideal field '$field'";
                evidence=Dict("subtype" => subtype, "field" => field,
                              "original_value" => deepcopy(value)))
        end
        haskey(tx, "no_load_shunt") || continue
        value = tx["no_load_shunt"]
        winding = value isa AbstractDict ? get(value, "winding", 0) : 0
        g = value isa AbstractDict ? get(value, "g", nothing) : nothing
        b = value isa AbstractDict ? get(value, "b", nothing) : nothing
        maxw = subtype == "center_tap" ? 3 : 2
        valid = value isa AbstractDict && winding isa Integer && !(winding isa Bool) &&
            Int(winding) in 1:maxw && g isa Real && !(g isa Bool) &&
            isfinite(g) && g >= 0 && b isa Real && !(b isa Bool) && isfinite(b) &&
            !any(haskey(tx, key) for key in ("g_no_load", "b_no_load"))
        if valid
            delete!(tx, "no_load_shunt")
            _l3f_warning!(findings, "A.L3F.PERMISSIVE_TRANSFORMER_IDEALIZED",
                :transformer, tid, "supported transformer '$subtype' dropped nested no_load_shunt";
                evidence=Dict("subtype" => subtype, "field" => "no_load_shunt",
                              "original_value" => deepcopy(value)))
        else
            _l3f_error!(findings, "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED",
                :transformer, tid, "no_load_shunt is malformed and cannot be idealized")
        end
    end
end

"""
    _l3f_lower!(findings, net, options) -> Bool

Rewrite `net` in place towards the supported vocabulary, appending one finding
per rewrite. Returns whether anything changed. `net` is already a private copy.
"""
function _l3f_lower!(findings, net, options::L3FOptions)
    options.unsupported == :reject && return false
    before = length(findings)
    _l3f_lower_restricted_controls!(findings, net)
    if options.unsupported == :permissive
        for (_, tid, tx) in _l3f_transformers(net)
            haskey(tx, "no_load_shunt") || continue
            any(haskey(tx, key) for key in ("g_no_load", "b_no_load")) || continue
            _l3f_error!(findings, "E.L3F.TRANSFORMER_NONIDEAL_UNSUPPORTED",
                :transformer, tid,
                "nested no_load_shunt conflicts with scalar g_no_load/b_no_load")
        end
    end
    options.unsupported == :permissive &&
        _l3f_permissive_constraints!(findings, net)
    _l3f_lower_switches!(findings, net)
    _l3f_lower_capacitors!(findings, net)
    _l3f_lower_line_shunts!(findings, net)
    # Leakage referral depends on the chosen fixed tap. Resolve an experimental
    # adjustable interval before synthesizing its series impedance.
    options.unsupported in (:approximate, :permissive) && _l3f_project_taps!(findings, net)
    _l3f_lower_transformer_impedance!(findings, net)
    options.unsupported == :permissive &&
        _l3f_permissive_transformer_residuals!(findings, net)
    if options.unsupported in (:approximate, :permissive)
        for (lid, load) in sort!(collect(get(net, "load", Dict())); by=first)
            load isa AbstractDict && _l3f_project_load!(findings, load, String(lid))
        end
        # In :approximate, unsupported bus/line bounds remain errors. The
        # permissive preprocessor has already removed its explicitly listed
        # operational bounds with structured warnings.
    end
    length(findings) > before
end
