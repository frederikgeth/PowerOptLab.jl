"""One machine-readable estimator preflight finding (`:error`, `:warning`, or `:info`)."""
struct SEPreflightFinding
    severity::Symbol
    code::Symbol
    element::String
    message::String
end

"""
    SEPreflightReport

Electrical support inventory for one snapshot. `supported` means the input can
be represented by the compiled model, not that it is observable or that a solve
will converge. `elements` records representation and fixed taps; `omitted`
records network injections and operating limits not used by this estimator.
Sources, declared perfect grounds, measurement references, and exact-equation
requests are recorded separately. Compilation retains the report in `s.preflight`.
"""
struct SEPreflightReport
    supported::Bool
    findings::Vector{SEPreflightFinding}
    elements::Vector{NamedTuple}
    omitted::Vector{NamedTuple}
    sources::Vector{NamedTuple}
    grounded_terminals::Vector{Tuple{String,String}}
    measurements::Vector{NamedTuple}
    exact_equations::Vector{NamedTuple}
end

"""Compilation error carrying an inspectable [`SEPreflightReport`](@ref)."""
struct SEUnsupportedNetwork <: Exception
    report::SEPreflightReport
end
function Base.showerror(io::IO, err::SEUnsupportedNetwork)
    print(io, "Unsupported state-estimation network:")
    for f in err.report.findings
        f.severity === :error && print(io, "\n  ", f.code, " [", f.element, "]: ", f.message)
    end
end

# Do not silently accept a primitive exporter warning about a skipped element or
# a singular shunt-only substitute. Retain the diagnostic on the report instead.
function _se_capture_primitive(f)
    buffer = IOBuffer()
    value = with_logger(SimpleLogger(buffer, Logging.Warn)) do
        f()
    end
    value, String(take!(buffer))
end

function _se_preflight(net, measurements; neutral="n", zero_injection=String[], exact_devices=Any[])
    findings = SEPreflightFinding[]; elements = NamedTuple[]; omitted = NamedTuple[]
    sources = NamedTuple[]; refs = NamedTuple[]; exact = NamedTuple[]
    grounds = Tuple{String,String}[]
    issue(severity, code, element, message) = push!(findings,
        SEPreflightFinding(severity, code, element, message))
    buses = get(net, "bus", Dict())
    validnode(bus, terminal) = haskey(buses, bus) && terminal in get(buses[bus], "terminal_names", [])
    function terminals!(element, bus, terminals)
        isempty(terminals) && issue(:error, :empty_terminal_map, element, "Terminal map is empty.")
        length(unique(terminals)) == length(terminals) || issue(:error, :duplicate_terminal, element, "Terminal map repeats a conductor.")
        for t in terminals
            validnode(bus, t) || issue(:error, :absent_terminal, element, "Terminal ($bus, $t) is not declared on a bus.")
        end
    end
    for (id, bus) in buses
        terminals!("bus/$id", id, get(bus, "terminal_names", String[]))
        for terminal in get(bus, "perfectly_grounded_terminals", [])
            terminals!("bus/$id", id, [terminal]); push!(grounds, (String(id), String(terminal)))
        end
        for key in ("v_min", "v_max", "vm_min", "vm_max")
            haskey(bus, key) && push!(omitted, (element="bus/$id", field=key, reason="Operating limit is not an estimator constraint."))
        end
    end
    for (id, source) in get(net, "voltage_source", Dict())
        ts = get(source, "terminal_map", String[]); bus = get(source, "bus", "")
        terminals!("voltage_source/$id", bus, ts)
        vs = get(source, "v_magnitude", Float64[]); angles = get(source, "v_angle", zeros(length(ts)))
        length(vs) == length(ts) && length(angles) == length(ts) &&
            all(v -> isfinite(v) && v >= 0, vs) && all(isfinite, angles) ||
            issue(:error, :invalid_source, "voltage_source/$id", "Supply finite phasors with one magnitude and angle per source conductor.")
        push!(sources, (element="voltage_source/$id", bus=bus, terminals=copy(ts), magnitudes=copy(vs), angles=copy(angles), treatment=:fixed_phasor))
    end
    isempty(sources) && issue(:warning, :no_fixed_source, "network", "No fixed source: supply an initial state and enough reference information in measurements.")
    # This deliberately excludes all ideal transformers, including unity-ratio
    # delta/wye cases: a scalar nameplate ratio is not a conductor identity.
    supported_xfmrs = ("single_phase", "center_tap", "wye_delta", "delta_wye",
                       "single_phase_autotransformer", "open_delta_regulator", "n_winding")
    for (kind, group) in get(net, "transformer", Dict()), (id, x) in group
        label = "transformer/$kind/$id"
        if !(kind in supported_xfmrs)
            issue(:error, :unsupported_transformer, label, "No supported finite-admittance transformer primitive."); continue
        end
        try
            if kind == "n_winding"
                raw = get(x, "windings", Any[])
                windings = BMOPFTools._nw_windings(x)
                length(windings) >= 2 && length(windings) == length(raw) ||
                    issue(:error, :invalid_windings, label, "Supply at least two winding records.")
                phase_counts = Int[]
                for w in windings
                    terminals!(label, w.bus, w.terminal_map)
                    isfinite(w.v_nom) && w.v_nom > 0 ||
                        issue(:error, :invalid_ratio, label, "Every winding nominal voltage must be finite and positive.")
                    w.connection in ("WYE", "DELTA") ||
                        issue(:error, :invalid_connection, label, "Winding configuration must be WYE or DELTA.")
                    phases, ref = BMOPFTools._nw_phase_terminals(w.terminal_map)
                    length(phases) + (ref === nothing ? 0 : 1) == length(w.terminal_map) ||
                        issue(:error, :unsupported_terminal_name, label, "The n-winding exporter recognises only a/b/c/n conductor names.")
                    push!(phase_counts, length(phases))
                    w.connection == "DELTA" && length(phases) < 2 &&
                        issue(:error, :terminal_count, label, "A delta winding requires at least two phases.")
                end
                !isempty(phase_counts) && minimum(phase_counts) > 0 && all(==(first(phase_counts)), phase_counts) ||
                    issue(:error, :terminal_count, label, "All windings must have the same nonzero phase count.")
                any(key -> haskey(x, key), ("tap", "tap_ratio")) &&
                    issue(:error, :unsupported_tap, label, "The n-winding exporter does not apply tap fields; encode known winding ratios in v_nom.")
            else
                for side in ("from", "to")
                    terminals!(label, get(x, "bus_$side", ""), get(x, "terminal_map_$side", String[]))
                end
                if !(kind in ("single_phase_autotransformer", "open_delta_regulator"))
                    all(key -> haskey(x, key) && isfinite(x[key]) && x[key] > 0, ("v_nom_from", "v_nom_to")) ||
                        issue(:error, :invalid_ratio, label, "Both nominal winding voltages must be finite and positive.")
                end
                taps = get(x, kind in ("single_phase_autotransformer", "open_delta_regulator") ? "tap_ratio" : "tap", 1.0)
                tv = taps isa Real ? [taps] : taps
                !isempty(tv) && all(t -> isfinite(t) && t > 0, tv) || issue(:error, :invalid_tap, label, "Fixed tap ratios must be finite and positive.")
                kind in ("single_phase_autotransformer", "open_delta_regulator") &&
                    !(uppercase(get(x, "regulator_type", "B")) in ("A", "B")) &&
                    issue(:error, :invalid_regulator_type, label, "Regulator type must be A or B.")
                if kind == "single_phase_autotransformer"
                    all(side -> length(BMOPFTools._xfmr_winding_pairs(x["terminal_map_$side"])) == 1, ("from", "to")) ||
                        issue(:error, :terminal_count, label, "A single-phase regulator must contain exactly one winding pair on each side.")
                end
                if kind == "open_delta_regulator"
                    all(side -> length(BMOPFTools._phase_positions(x["terminal_map_$side"])) == 3, ("from", "to")) ||
                        issue(:error, :terminal_count, label, "An open-delta regulator requires exactly three phase conductors on each side.")
                    length(tv) == 2 || issue(:error, :invalid_tap, label, "Supply exactly two open-delta tap ratios.")
                elseif !(taps isa Real)
                    issue(:error, :invalid_tap, label, "Supply one scalar fixed tap ratio.")
                end
                z = [Float64(get(x, key, 0.0)) for key in ("r_series_from", "x_series_from", "r_series_to", "x_series_to")]
                all(isfinite, z) || issue(:error, :invalid_impedance, label, "Leakage impedances must be finite.")
                if maximum(abs, z) <= 1e-6 || (kind == "center_tap" && (iszero(complex(z[1],z[2])) || iszero(complex(z[3],z[4]))))
                    issue(:error, :ideal_transformer, label, "Ideal/degenerate coupling is outside the voltage-only support contract. Supply physical finite leakage or use an augmented-state formulation; do not invent a small impedance.")
                    continue
                end
            end
            (primitive, warnings) = _se_capture_primitive() do
                kind == "n_winding" ? BMOPFTools.nwinding_yprim(x) : transformer_yprim(x, kind)
            end
            nodes, Y = primitive
            for (bus, terminal) in nodes; terminals!(label, bus, [terminal]); end
            size(Y) == (length(nodes), length(nodes)) || issue(:error, :primitive_shape, label, "Primitive dimensions do not match its terminals.")
            isempty(warnings) || issue(:error, :primitive_warning, label, strip(warnings))
            isempty(nodes) && issue(:error, :empty_primitive, label, "Transformer was not represented.")
            all(isfinite, Y) || issue(:error, :nonfinite_primitive, label, "Transformer primitive contains nonfinite entries.")
            push!(elements, (element=label, representation=:finite_admittance,
                tap=deepcopy(get(x, "tap_ratio", get(x, "tap", 1.0)))))
        catch err
            issue(:error, :invalid_transformer, label, sprint(showerror, err))
        end
    end
    for kind in ("line", "switch", "shunt", "capacitor"), (id, x) in get(net, kind, Dict())
        label = "$kind/$id"
        try
            if kind in ("line", "switch")
                for side in ("from", "to")
                    terminals!(label, get(x, "bus_$side", ""), get(x, "terminal_map_$side", String[]))
                end
                length(get(x,"terminal_map_from",[])) == length(get(x,"terminal_map_to",[])) ||
                    issue(:error, :terminal_count, label, "Conductor maps must have equal length.")
            else
                terminals!(label, get(x,"bus",""), get(x,"terminal_map",String[]))
            end
            if kind == "line"
                Z, _ = BMOPFTools._line_z_complex(x, get(net,"linecode",Dict()))
                Z === nothing && issue(:error, :missing_line_impedance, label, "No series impedance can be assembled.")
                Z !== nothing && !all(isfinite, Z) && issue(:error, :invalid_impedance, label, "Line impedance must be finite.")
                if Z !== nothing && norm(Z) > 1e-4
                    (primitive, warnings) = _se_capture_primitive(() -> line_yprim(x, get(net,"linecode",Dict())))
                    isempty(warnings) || issue(:error, :primitive_warning, label, strip(warnings))
                    isempty(primitive[1]) && issue(:error,:empty_primitive,label,"Line was not represented.")
                    all(isfinite, primitive[2]) || issue(:error,:nonfinite_primitive,label,"Line primitive is nonfinite.")
                end
            elseif kind in ("shunt", "capacitor")
                primitive = kind == "shunt" ? BMOPFTools._shunt_yprim(x) : BMOPFTools._capacitor_yprim(x)
                nodes, Y = primitive
                isempty(nodes) && issue(:error, :empty_primitive, label, "Element has no admittance data.")
                size(Y) == (length(nodes), length(nodes)) && all(isfinite, Y) ||
                    issue(:error, :invalid_primitive, label, "Admittance must be finite and match the terminal map.")
            elseif kind == "switch"
                get(x,"status","closed") in ("open","closed") || issue(:error,:invalid_switch,label,"Status must be open or closed.")
            end
            push!(elements, (element=label, representation=kind == "switch" ? Symbol(get(x,"status","closed")) : :passive, tap=nothing))
        catch err
            issue(:error, :invalid_element, label, sprint(showerror,err))
        end
    end
    for kind in ("load", "generator", "ibr"), (id, x) in get(net, kind, Dict())
        push!(omitted, (element="$kind/$id", field="injection", reason="Not imported from network data: supply measurements, pseudomeasurements, or exact_devices explicitly."))
    end
    # Controls and operating envelopes are not part of the passive snapshot.
    function inventory_omissions(group, path)
        if group isa AbstractVector
            for (i, record) in enumerate(group)
                inventory_omissions(record, "$path/$i")
            end
            return
        end
        group isa AbstractDict || return
        for (key, value) in group
            field = String(key)
            if endswith(field, "_min") || endswith(field, "_max") || occursin("rating", lowercase(field)) || occursin("control", lowercase(field))
                item = (element=path, field=field, reason="Operating bound or control is not enforced; taps and topology are fixed snapshot data.")
                any(o -> o.element == path && o.field == field, omitted) || push!(omitted, item)
            elseif value isa Union{AbstractDict,AbstractVector}
                inventory_omissions(value, isempty(path) ? field : "$path/$field")
            end
        end
    end
    inventory_omissions(net, "")
    !isempty(omitted) && issue(:warning,:omitted_information,"network","Some injections or operating limits are not used; inspect report.omitted.")
    known = Set(("bus","linecode","line","switch","shunt","capacitor","transformer","voltage_source","load","generator","ibr"))
    for (kind, group) in net
        kind in known && continue
        group isa AbstractDict || continue
        if any(x -> x isa AbstractDict && (haskey(x,"bus") || haskey(x,"bus_from")), values(group))
            issue(:error,:unsupported_element,String(kind),"Unrecognised electrical element group would be omitted.")
        end
    end
    ybus = nothing
    if !any(f -> f.severity === :error, findings)
        try
            ybus, warnings = _se_capture_primitive(() -> ybus_passive(net))
            isempty(warnings) || issue(:error,:primitive_warning,"network",strip(warnings))
            for (node, index) in sort!(collect(ybus.index); by=first)
                if index != 0 && node != ybus.nodes[index]
                    push!(elements, (element="$(node[1])/$(node[2])", representation=:node_alias, tap=nothing))
                end
            end
            all(isfinite, ybus.Y) || issue(:error,:nonfinite_admittance,"network","Passive Ybus is nonfinite.")
            fixed = _source_phasors(net, ybus.index)
            for node in _zero_injection_set(net, zero_injection, neutral; include_neutral=true)
                i = get(ybus.index, node, 0)
                i == 0 && throw(ArgumentError("zero-injection terminal $node is grounded or absent"))
                haskey(fixed, i) && throw(ArgumentError("zero-injection terminal $node is fixed by a source"))
            end
            _compile_measurements(net, measurements, ybus.index, neutral)
            _compile_exact_devices(exact_devices, ybus.index)
        catch err
            issue(:error,:invalid_compilation,"network",sprint(showerror,err))
        end
    end
    for (i,m) in enumerate(measurements)
        if !(m isa Union{Measurement,BranchMeasurement})
            issue(:error, :invalid_measurement, string(i), "Expected Measurement or BranchMeasurement.")
            continue
        end
        push!(refs, m isa BranchMeasurement ? (row=i, kind=m.kind, location="line/$(m.line)/$(m.side)/$(m.terminal)", reference=nothing) :
            (row=i, kind=m.kind, location="$(m.bus)/$(m.terminal)", reference=_resolve_ref(m,neutral)))
    end
    for node in zero_injection; push!(exact,(kind=:zero_injection, specification=deepcopy(node))); end
    for d in exact_devices; push!(exact,(kind=:device, specification=deepcopy(d))); end
    report = SEPreflightReport(!any(f -> f.severity === :error, findings), findings, elements,
        omitted, sources, sort!(grounds), refs, exact)
    report, ybus
end

"""
    state_estimator_preflight(net, measurements=[]; neutral="n", zero_injection=[], exact_devices=[])

Inspect the supported electrical interpretation before compiling. Ideal and
warning-producing transformer primitives are rejected, never regularised or
silently skipped. Fixed taps are data, not estimated controls. This report does
not certify observability, model accuracy, or convergence.
"""
state_estimator_preflight(net::Dict{String,Any}, measurements::AbstractVector=Measurement[]; kwargs...) =
    first(_se_preflight(net, measurements; kwargs...))
