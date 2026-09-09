# Widening the admissible input without widening the formulation.
#
# The model itself only knows about a small vocabulary: series lines, fixed
# shunts, ideal fixed-ratio transformers, and P/Z loads. Much of what BMOPF can
# express sits outside that vocabulary while remaining inside the same
# mathematical class, and can be rewritten into it *exactly*. A smaller set
# genuinely cannot, and can only be projected onto it at a cost.
#
# `L3FOptions(unsupported=...)` selects how far to go:
#
#   :reject       nothing is rewritten (the default; what the formulation
#                 promises when it says a component is unsupported)
#   :lower        exact re-representations only, reported as `L.L3F.*` at
#                 severity :info — the solved model is the same physics
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
# Exact re-representations
# ---------------------------------------------------------------------------

"""
Closed switches become zero-impedance lines; open switches are removed.

Both are exact. A closed ideal switch is a branch with no series drop, which the
line kernel already represents, and an open switch carries no current at all. If
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
        lines[id] = line
        _l3f_info!(findings, "L.L3F.SWITCH_LOWERED", :switch, sid,
            "closed switch represented exactly as the zero-impedance line '$id'";
            evidence=Dict("line" => id))
    end
    delete!(net, "switch")
end

"""
Fixed capacitor banks become fixed shunt susceptance.

`B = q_rated / v_nom^2` per coil. For a delta bank the coil admittances are
referred to the terminals by ``D^{T}\\operatorname{diag}(b)D``, which is exactly
the terminal admittance matrix of the same three coils. BMOPF capacitors carry
no switching state, so this is an exact restatement, not a projection.
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
                "fixed capacitor represented exactly as the shunt '$id' with " *
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
moving each half onto its own bus is a restatement of the same pi model rather
than a lumping approximation. Linecode entries are per unit length and are
scaled by the line's `length`, matching the series path.
"""
function _l3f_lower_line_shunts!(findings, net)
    lines = get(net, "line", Dict())
    linecodes = get(net, "linecode", Dict())
    shunts = _l3f_table!(net, "shunt")
    for (lid_raw, line) in sort!(collect(lines); by=first)
        lid = String(lid_raw)
        line isa AbstractDict || continue
        sources = Tuple{Any,Float64}[(line, 1.0)]
        lcid = get(line, "linecode", nothing)
        if lcid isa AbstractString && haskey(linecodes, lcid)
            push!(sources, (linecodes[lcid], Float64(get(line, "length", 1.0))))
        end
        for (side, bus_key, map_key) in (("from", "bus_from", "terminal_map_from"),
                                         ("to", "bus_to", "terminal_map_to"))
            tm = string.(get(line, map_key, String[]))
            isempty(tm) && continue
            entries = Dict{Tuple{Int,Int},ComplexF64}()
            for (data, scale) in sources, key in keys(data)
                m = match(Regex("^([GB])_$(side)_(\\d+)_(\\d+)\$"), String(key))
                m === nothing && continue
                i, j = parse(Int, m.captures[2]), parse(Int, m.captures[3])
                max(i, j) <= length(tm) || continue
                value = Float64(data[key]) * scale
                pair = i <= j ? (i, j) : (j, i)
                current = get(entries, pair, ComplexF64(0))
                entries[pair] = m.captures[1] == "G" ?
                    current + value : current + im * value
            end
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
A non-ideal transformer becomes an ideal one plus explicit series and shunt
elements.

Winding leakage is a series impedance in the coil's own coordinates and the
no-load admittance is a shunt across the winding-2 coil, so introducing an
internal bus per non-zero winding and stamping those as ordinary line and shunt
elements reproduces the same two-port exactly. Nothing is dropped and no
accuracy is lost — the formulation simply had no vocabulary for a component
carrying its own impedance.
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
                    "$side-winding leakage $(r) + j$(x) ohm represented exactly as " *
                    "the series line '$line_id' through the internal bus '$internal'";
                    evidence=Dict("line" => line_id, "bus" => internal,
                                  "side" => side, "r" => r, "x" => x))
            end
            g, b, present = _l3f_scalar_impedance(transformer, "g_no_load", "b_no_load")
            present || continue
            # The no-load shunt is stamped across the winding-2 (to-side) coil.
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
                "no-load admittance $(g) + j$(b) S represented exactly as the " *
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
``\\alpha_Z=\\gamma/2`` and ``\\alpha_P=1-\\gamma/2``. The construction is exact
for ``\\gamma=0`` (constant power) and ``\\gamma=2`` (constant impedance), and
for ``\\gamma=1`` it is the familiar half-to-Z, half-to-P split of a
constant-current term.
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
            "tangent at v_nom; exact for exponents 0 and 2, approximate otherwise";
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

Unlike the exact rewrites above this removes a decision the caller asked for.
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
Drop bus limits the formulation does not assess.

A phase-to-phase or sequence-voltage limit is representable through the same
cross-voltage closure the model already forms, so this is a gap rather than an
impossibility. Until those constraints exist, `:approximate` drops them and says
so — the solved problem is a relaxation, and its solution may violate them.
"""
function _l3f_project_bus_limits!(findings, net)
    for (bid, bus) in sort!(collect(get(net, "bus", Dict())); by=first)
        bus isa AbstractDict || continue
        dropped = String[]
        for key in ("vpp_min", "vpp_max", "vpos_min", "vpos_max", "vneg_max",
                    "vzero_max", "vm_unbalance_max")
            haskey(bus, key) || continue
            push!(dropped, key); delete!(bus, key)
        end
        isempty(dropped) && continue
        _l3f_warning!(findings, "A.L3F.BUS_LIMIT_DROPPED", :bus, bid,
            "bus limit(s) $(join(dropped, ", ")) dropped; the solved problem is a " *
            "relaxation and its solution may violate them";
            evidence=Dict("limits" => dropped))
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
    _l3f_lower_switches!(findings, net)
    _l3f_lower_capacitors!(findings, net)
    _l3f_lower_line_shunts!(findings, net)
    _l3f_lower_transformer_impedance!(findings, net)
    if options.unsupported == :approximate
        for (lid, load) in sort!(collect(get(net, "load", Dict())); by=first)
            load isa AbstractDict && _l3f_project_load!(findings, load, String(lid))
        end
        _l3f_project_taps!(findings, net)
        _l3f_project_bus_limits!(findings, net)
    end
    length(findings) > before
end
