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
A non-ideal single-phase transformer becomes an ideal one plus explicit series
and shunt elements. This is a canonical L3F lowering: it preserves the
single-phase topology and the fixed-angle/lossless approximation, but it is not
an exact AC two-port. Connection-aware Yd/Dy and multi-phase transformer
leakage/no-load lowering is deliberately not attempted; those components remain
unsupported rather than being represented by false diagonal phase lines.
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
            # A diagonal phase-line decomposition is justified only for a
            # single phase-to-ground transformer winding.  In particular,
            # autotransformer, Yd/Dy, and open-delta leakage lives in coupled
            # coil coordinates and cannot be replaced by independent phase
            # lines without losing the connection and shared-winding physics.
            subtype == "single_phase" || continue
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
            # BMOPFTools places the ordinary single-phase exciting branch
            # across winding 2 (the to-side coil), at its total coil
            # admittance.  There is no per-phase splitting for this subtype.
            bus = String(get(transformer, "bus_to", ""))
            haskey(buses, bus) || continue
            tm = string.(get(transformer, "terminal_map_to", String[]))
            id = _l3f_derived("noload", tid)
            shunt = Dict{String,Any}("bus" => bus, "terminal_map" => copy(tm))
            for k in eachindex(tm)
                iszero(g) || (shunt["G_$(k)_$(k)"] = g)
                iszero(b) || (shunt["B_$(k)_$(k)"] = b)
            end
            shunts[id] = shunt
            delete!(transformer, "g_no_load"); delete!(transformer, "b_no_load")
            _l3f_info!(findings, "L.L3F.TRANSFORMER_NO_LOAD_LOWERED", :transformer, tid,
                "no-load admittance $(g) + j$(b) S represented in the " *
                "shunt '$id' across the to-side coil";
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
    ai, bi = channels("alpha_i", 0.0), channels("beta_i", 0.0)
    all(iszero, ai) && all(iszero, bi) && return
    az, ap = channels("alpha_z", 0.0), channels("alpha_p", 0.0)
    bz, bp = channels("beta_z", 0.0), channels("beta_p", 0.0)
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
        value_key, lo_key, hi_key = subtype == "single_phase" ?
            ("tap", "tap_min", "tap_max") :
            ("tap_ratio", "tap_ratio_min", "tap_ratio_max")
        haskey(transformer, lo_key) || haskey(transformer, hi_key) || continue
        vector(x) = x isa AbstractVector ? Float64.(x) : [Float64(x)]
        lo = vector(get(transformer, lo_key, NaN))
        hi = vector(get(transformer, hi_key, NaN))
        length(lo) == length(hi) && all(isfinite, lo) && all(isfinite, hi) || continue
        all(lo .== hi) && continue
        declared = haskey(transformer, value_key) ? vector(transformer[value_key]) : Float64[]
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

"""
Legacy helper retained for source compatibility.

Unimplemented BMOPF voltage limits are now applicability errors under every
policy, so this helper deliberately does not mutate the network. Constraints
must not be silently relaxed by `unsupported=:approximate`.
"""
function _l3f_project_bus_limits!(findings, net)
    nothing
end

"""
    _l3f_lower!(findings, net, options) -> Bool

Rewrite `net` in place towards the supported vocabulary, appending one finding
per rewrite. Returns whether anything changed. `net` is already a private copy.
"""
function _l3f_lower!(findings, net, options::L3FOptions)
    options.unsupported == :reject && return false
    before = length(findings)
    _l3f_lower_switches!(findings, net)
    _l3f_lower_capacitors!(findings, net)
    _l3f_lower_line_shunts!(findings, net)
    _l3f_lower_transformer_impedance!(findings, net)
    if options.unsupported == :approximate
        for (lid, load) in sort!(collect(get(net, "load", Dict())); by=first)
            load isa AbstractDict && _l3f_project_load!(findings, load, String(lid))
        end
        _l3f_project_taps!(findings, net)
        # Unsupported bus/line voltage bounds are applicability errors, never
        # experimental relaxations.  Keep them in the private copy so the
        # validator can report every offending field.
    end
    length(findings) > before
end
