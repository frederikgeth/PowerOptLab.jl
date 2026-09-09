function _l3f_reference_hash(voltage)
    records = ["$(bus)\0$(terminal)\0$(repr(real(value)))\0$(repr(imag(value)))"
               for ((bus, terminal), value) in sort!(collect(voltage); by=first)]
    bytes2hex(SHA.sha256(join(records, "\n")))
end

function _l3f_explicit_reference(reference)
    reference isa L3FReferenceState && return reference
    reference isa AbstractDict || throw(ArgumentError(
        "reference must be L3FReferenceState or a dictionary"))
    voltage = Dict{Tuple{String,String},ComplexF64}()
    if all(key -> key isa Tuple && length(key) == 2, keys(reference))
        for (key, value) in reference
            voltage[(String(key[1]), String(key[2]))] = ComplexF64(value)
        end
    else
        for (bus, terminals) in reference
            terminals isa AbstractDict || throw(ArgumentError(
                "nested reference entries must map terminals to phasors"))
            for (terminal, value) in terminals
                phasor = if value isa AbstractDict
                    ComplexF64(Float64(get(value, "vr", get(value, :vr, 0.0))),
                               Float64(get(value, "vi", get(value, :vi, 0.0))))
                else
                    ComplexF64(value)
                end
                voltage[(String(bus), String(terminal))] = phasor
            end
        end
    end
    L3FReferenceState(voltage, :explicit, _l3f_reference_hash(voltage), :not_run)
end

function _l3f_source_reference(net, topology)
    voltage = Dict{Tuple{String,String},ComplexF64}()
    for (_, source) in get(net, "voltage_source", Dict())
        bus = String(source["bus"])
        tm = string.(source["terminal_map"])
        vm = Float64.(source["v_magnitude"])
        va = Float64.(source["v_angle"])
        for k in eachindex(tm)
            voltage[(bus, tm[k])] = vm[k] * cis(va[k])
        end
    end
    pending = copy(topology)
    while !isempty(pending)
        progressed = false
        for edge in copy(pending)
            all(haskey(voltage, (edge.parent, terminal)) for terminal in edge.parent_map) || continue
            transform = if edge.family == :line
                Matrix{Float64}(I, length(edge.parent_map), length(edge.parent_map))
            else
                data = net["transformer"][edge.subtype][edge.id]
                _l3f_transformer_oriented_map(edge, data)
            end
            for k in eachindex(edge.parent_map)
                key = (edge.child, edge.child_map[k])
                value = sum(transform[k, j] * voltage[(edge.parent, edge.parent_map[j])]
                            for j in eachindex(edge.parent_map))
                if haskey(voltage, key) && !isapprox(voltage[key], value; atol=1e-10, rtol=1e-10)
                    throw(ArgumentError("conflicting propagated reference at $key"))
                end
                voltage[key] = value
            end
            deleteat!(pending, findfirst(==(edge), pending))
            progressed = true
        end
        progressed || throw(ArgumentError("could not propagate source references through the radial topology"))
    end
    L3FReferenceState(voltage, :source_propagated,
                      _l3f_reference_hash(voltage), :not_run)
end

function _l3f_validate_reference(reference::L3FReferenceState, net)
    for (bus, data) in get(net, "bus", Dict())
        for terminal in string.(get(data, "terminal_names", String[]))
            key = (String(bus), terminal)
            haskey(reference.voltage, key) ||
                throw(ArgumentError("reference is missing bus terminal $key"))
            value = reference.voltage[key]
            isfinite(real(value)) && isfinite(imag(value)) && !iszero(value) ||
                throw(ArgumentError("reference at $key must be finite and nonzero"))
        end
    end
    for (sid, source) in get(net, "voltage_source", Dict())
        bus = String(source["bus"])
        tm = string.(source["terminal_map"])
        vm = Float64.(source["v_magnitude"])
        va = Float64.(source["v_angle"])
        for k in eachindex(tm)
            expected = vm[k] * cis(va[k])
            isapprox(reference.voltage[(bus, tm[k])], expected; atol=1e-8, rtol=1e-8) ||
                throw(ArgumentError("reference conflicts with voltage source '$sid' at terminal $(tm[k])"))
        end
    end
    reference
end

function _l3f_reference(net, topology, reference, policy::Symbol=:auto)
    # `:source_propagated` names the flat propagated profile specifically, so it
    # must not be silently overridden by a caller-supplied reference; `:auto`
    # prefers an explicit reference when one is given. `:explicit` requires one,
    # which `_l3f_prepare` has already enforced.
    use_explicit = reference !== nothing && policy != :source_propagated
    result = use_explicit ? _l3f_explicit_reference(reference) :
             _l3f_source_reference(net, topology)
    _l3f_validate_reference(result, net)
end

"""
Restate every voltage-source rating in working coordinates from its SI value.

BMOPFTools' classic per-unit preparation scales a source's `v_magnitude`,
`p_min`/`p_max`/`q_min`/`q_max` and `cost`, but at the pinned revision it leaves
`s_max` and `i_max` in SI. Stamping those raw values against per-unit variables
would silently widen a source nameplate by a factor of `s_base`, so the rating
would never bind in the default coordinates while binding correctly under
`per_unit=false`.

Deriving the working value from the caller's SI network on every build — rather
than patching whatever the working copy happens to hold — keeps this correct and
idempotent if upstream later scales these fields itself.
"""
function _l3f_rescale_source_ratings!(working, physical, bases, options::L3FOptions)
    for (sid, source) in get(working, "voltage_source", Dict())
        source isa AbstractDict || continue
        original = get(get(physical, "voltage_source", Dict()), sid, nothing)
        original isa AbstractDict || continue
        bus = String(get(source, "bus", ""))
        haskey(source, "s_max") && haskey(original, "s_max") &&
            (source["s_max"] = Float64.(original["s_max"]) ./ options.s_base)
        haskey(source, "i_max") && haskey(original, "i_max") &&
            (source["i_max"] = Float64.(original["i_max"]) ./ get(bases.i_base, bus, 1.0))
    end
    working
end

function _l3f_working_coordinates(net, reference::L3FReferenceState,
                                  options::L3FOptions)
    options.per_unit || return deepcopy(net), reference, nothing

    # Reuse BMOPFTools' public coordinate-preparation contract, but none of its
    # nonlinear component builders.  The returned network is a private working
    # copy and the caller's BMOPF dictionary remains in SI.
    context = BMOPFTools.initialize_opf_model(net;
        per_unit=true, s_base=options.s_base, model=JuMP.Model(), kcl_guard=false)
    working = BMOPFTools.opf_network(context)
    bases = BMOPFTools.opf_bases(context)
    _l3f_rescale_source_ratings!(working, net, bases, options)
    voltage = Dict{Tuple{String,String},ComplexF64}(
        key => value / bases.v_base[key[1]] for (key, value) in reference.voltage)
    working_reference = L3FReferenceState(
        voltage, reference.provenance, _l3f_reference_hash(voltage),
        reference.nonlinear_status)
    _l3f_validate_reference(working_reference, working)
    working, working_reference, bases
end

function _l3f_report_error(report::L3FApplicabilityReport, code, message)
    findings = copy(report.findings)
    _l3f_error!(findings, code, :network, nothing, message)
    L3FApplicabilityReport(:inapplicable, findings, report.roots,
                           report.islands, report.kron_reduced, report.lowered)
end

function _l3f_add_bounds!(variable, data, index::Int,
                          lower_key::String, upper_key::String)
    lower = get(data, lower_key, nothing)
    upper = get(data, upper_key, nothing)
    lower_value = lower isa AbstractVector ? (length(lower) >= index ? Float64(lower[index]) : nothing) :
                  lower isa Real ? Float64(lower) : nothing
    upper_value = upper isa AbstractVector ? (length(upper) >= index ? Float64(upper[index]) : nothing) :
                  upper isa Real ? Float64(upper) : nothing
    lower_value === nothing || JuMP.set_lower_bound(variable, lower_value)
    upper_value === nothing || JuMP.set_upper_bound(variable, upper_value)
    variable
end

_l3f_name(parts...) = join(replace.(string.(parts), r"[^A-Za-z0-9_]" => "_"), "__")

function _l3f_register_constraint!(constraints, family::Symbol, key, constraint)
    get!(constraints, family, Dict{Any,Any}())[key] = constraint
    constraint
end

"""Read a scalar or per-conductor rating, returning `nothing` when absent."""
function _l3f_scalar_or_indexed(data, field::String, k::Int)
    haskey(data, field) || return nothing
    raw = data[field]
    raw isa AbstractVector ? Float64(raw[length(raw) == 1 ? 1 : k]) : Float64(raw)
end

function _l3f_rating(component, fallback, field::String, k::Int)
    data = haskey(component, field) ? component : fallback
    data === nothing || !haskey(data, field) ? nothing : Float64(data[field][k])
end

function _l3f_add_power_circle!(model, constraints, family, key, p, q, radius)
    radius === nothing && return
    _l3f_register_constraint!(constraints, family, key,
        @constraint(model, [Float64(radius), p, q] in JuMP.SecondOrderCone()))
end

"""
Stamp ``p^2+q^2 ≤ w I_max^2`` as a native rotated SOC.

`voltage_squared` may be a bus `w` variable or an affine fixed-angle winding-
voltage closure. The cone itself is not outer-linearized.
"""
function _l3f_add_current_cone!(model, constraints, family, key,
                                p, q, voltage_squared, rating)
    rating === nothing && return
    imax = Float64(rating)
    _l3f_register_constraint!(constraints, family, key,
        @constraint(model,
            [voltage_squared, imax^2 / 2, p, q] in JuMP.RotatedSecondOrderCone()))
end

"""Affine fixed-angle approximation of one physical channel's `|D*v|²`."""
function _l3f_winding_voltage_squared(D, k::Int, vbar, bus::String, tm, w)
    winding = winding_voltage_coefficients(view(D, k, :), vbar)
    out = JuMP.AffExpr(winding.constant)
    for terminal in eachindex(tm)
        JuMP.add_to_expression!(out, winding.coefficients[terminal],
                                w[(bus, tm[terminal])])
    end
    out
end

"""Power of a synthetic shunt introduced by canonical component lowering."""
function _l3f_derived_shunt_power(net, kind::String, id::String, side,
                                  bus::String, tm, reference, w, k::Int)
    shunt_id = side === nothing ? _l3f_derived(kind, id) :
                                  _l3f_derived(kind, id, side)
    shunt = get(get(net, "shunt", Dict()), shunt_id, nothing)
    shunt isa AbstractDict || return (JuMP.AffExpr(0.0), JuMP.AffExpr(0.0))
    shunt_tm = string.(get(shunt, "terminal_map", String[]))
    shunt_bus = String(get(shunt, "bus", ""))
    shunt_bus == bus && shunt_tm == tm || throw(ArgumentError(
        "derived shunt '$shunt_id' no longer matches its $kind endpoint"))
    _l3f_shunt_power(shunt, bus, tm, reference, w, k)
end

"""Return `sign*Sseries + Sshunt` as affine real/reactive expressions."""
function _l3f_total_endpoint_power(p, q, sign::Real, ps, qs)
    pt, qt = JuMP.AffExpr(0.0), JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(pt, sign, p)
    JuMP.add_to_expression!(qt, sign, q)
    JuMP.add_to_expression!(pt, ps)
    JuMP.add_to_expression!(qt, qs)
    pt, qt
end

"""Per-coil WYE rating for a Yd/Dy bank in the active coordinates."""
function _l3f_ywye_rating(transformer, subtype::String, n_ph::Int)
    # BMOPFTools stores the scalar nameplate on the from-side base for
    # compatibility, and preserves both side-base values privately after
    # per-unit preparation.  The WYE coil is from-side for Yd and to-side for
    # Dy; use that side's value before taking the total-bank per-coil share.
    key = subtype == "wye_delta" ? "_s_rating_from_pu" : "_s_rating_to_pu"
    value = haskey(transformer, key) ? transformer[key] :
            get(transformer, "s_rating", nothing)
    value === nothing ? nothing : Float64(value) / n_ph
end

function _l3f_shunt_power(data, bus::String, tm, reference, w, phi::Int)
    Y = _l3f_shunt_matrix(data, length(tm))
    p, q = JuMP.AffExpr(0.0), JuMP.AffExpr(0.0)
    for psi in eachindex(tm)
        cross = cross_voltage_coefficients(reference.voltage[(bus, tm[phi])],
                                           reference.voltage[(bus, tm[psi])])
        scale = conj(Y[phi, psi])
        constant = scale * cross.constant
        a = scale * cross.coefficient_phi
        b = scale * cross.coefficient_psi
        JuMP.add_to_expression!(p, real(constant)); JuMP.add_to_expression!(q, imag(constant))
        JuMP.add_to_expression!(p, real(a), w[(bus, tm[phi])])
        JuMP.add_to_expression!(q, imag(a), w[(bus, tm[phi])])
        JuMP.add_to_expression!(p, real(b), w[(bus, tm[psi])])
        JuMP.add_to_expression!(q, imag(b), w[(bus, tm[psi])])
    end
    p, q
end

function _l3f_channel_value(data, field::String, k::Int, default::Float64)
    raw = get(data, field, nothing)
    raw === nothing && return default
    raw isa AbstractVector ? Float64(length(raw) == 1 ? raw[1] : raw[k]) : Float64(raw)
end

function _l3f_zip_coefficients(load, family::Symbol, k::Int)
    fields = family == :p ? ("alpha_z", "alpha_i", "alpha_p") :
                            ("beta_z", "beta_i", "beta_p")
    all(field -> !haskey(load, field), fields) && return (0.0, 0.0, 1.0)
    (_l3f_channel_value(load, fields[1], k, 0.0),
     _l3f_channel_value(load, fields[2], k, 0.0),
     _l3f_channel_value(load, fields[3], k, 0.0))
end

function _l3f_load_channel_power(load, D, vbar, bus::String, tm, w, k::Int)
    model = lowercase(String(get(load, "model", "constant_power")))
    p_nom, q_nom = Float64(load["p_nom"][k]), Float64(load["q_nom"][k])
    model == "constant_power" && return (JuMP.AffExpr(p_nom), JuMP.AffExpr(q_nom))

    v_nom = _l3f_channel_value(load, "v_nom", k, NaN)
    w_winding = _l3f_winding_voltage_squared(D, k, vbar, bus, tm, w)

    if model == "constant_impedance"
        p_coeffs, q_coeffs = (1.0, 0.0, 0.0), (1.0, 0.0, 0.0)
    else
        p_coeffs = _l3f_zip_coefficients(load, :p, k)
        q_coeffs = _l3f_zip_coefficients(load, :q, k)
    end
    p = JuMP.AffExpr(p_nom * p_coeffs[3])
    q = JuMP.AffExpr(q_nom * q_coeffs[3])
    JuMP.add_to_expression!(p, p_nom * p_coeffs[1] / v_nom^2, w_winding)
    JuMP.add_to_expression!(q, q_nom * q_coeffs[1] / v_nom^2, w_winding)
    p, q
end

function _l3f_mapped_terminal_power(H::ConnectionPowerMap, p, q, terminal::Int)
    pt, qt = JuMP.AffExpr(0.0), JuMP.AffExpr(0.0)
    for channel in eachindex(p)
        JuMP.add_to_expression!(pt, H.real_part[terminal, channel], p[channel])
        JuMP.add_to_expression!(pt, -H.imag_part[terminal, channel], q[channel])
        JuMP.add_to_expression!(qt, H.imag_part[terminal, channel], p[channel])
        JuMP.add_to_expression!(qt, H.real_part[terminal, channel], q[channel])
    end
    pt, qt
end

function _l3f_open_delta_limits!(model, constraints, edge, transformer,
                                 child_p, child_q, parent_power, reference, w)
    pairs = _L3F_OPEN_DELTA_PAIRS[uppercase(String(transformer["connection"]))]
    shared = only(intersect(collect(pairs[1]), collect(pairs[2])))
    original_from_is_parent = !edge.reversed
    for (j, pair) in enumerate(pairs)
        other = pair[1] == shared ? pair[2] : pair[1]
        for (side, p, q, bus, tm) in (
            (:parent, parent_power[other][1], parent_power[other][2], edge.parent, edge.parent_map),
            (:child, child_p[other], child_q[other], edge.child, edge.child_map))
            original_side = (side == :parent) == original_from_is_parent ? "from" : "to"
            wterminal = w[(bus, tm[other])]
            vterminal = abs(reference.voltage[(bus, tm[other])])
            winding_voltage = abs(reference.voltage[(bus, tm[pair[1]])] -
                                  reference.voltage[(bus, tm[pair[2]])])
            if haskey(transformer, "s_rating")
                radius = Float64(transformer["s_rating"]) * vterminal / winding_voltage
                _l3f_add_power_circle!(model, constraints,
                    :transformer_winding_apparent_power,
                    (edge.subtype, edge.id, original_side, j), p, q, radius)
            end
            ikey = original_side == "from" ? "i_max_from" : "i_max_to"
            if haskey(transformer, ikey)
                rating = _l3f_scalar_or_indexed(transformer, ikey, j)
                _l3f_add_current_cone!(model, constraints,
                    :transformer_current,
                    (edge.subtype, edge.id, original_side, j),
                    p, q, wterminal, rating)
            end
        end
    end
end

function _l3f_cost_coefficient(component, k::Int, kind::String, id::String)
    cost = get(component, "cost", nothing)
    cost === nothing && return 0.0
    cost isa AbstractVector || throw(ArgumentError("$kind '$id' cost must be a per-phase vector"))
    length(cost) >= k || throw(ArgumentError("$kind '$id' cost does not cover phase $k"))
    Float64(cost[k]) / 1000.0
end

function _l3f_build_model(net, topology, reference, report, optimizer, options;
                          physical_network=net, physical_reference=reference,
                          bases=nothing)
    model = optimizer === nothing ? JuMP.Model() : JuMP.Model(optimizer)
    variables = Dict{Symbol,Any}()
    constraints = Dict{Symbol,Any}()

    w = Dict{Tuple{String,String},JuMP.VariableRef}()
    for (busid_raw, bus) in sort!(collect(get(net, "bus", Dict())); by=first)
        busid = String(busid_raw)
        terminals = string.(bus["terminal_names"])
        for (k, terminal) in enumerate(terminals)
            variable = @variable(model, lower_bound=0.0,
                base_name=_l3f_name("l3f_w", busid, terminal))
            lo = _l3f_bus_bound(bus, "v_min", k, length(terminals))
            hi = _l3f_bus_bound(bus, "v_max", k, length(terminals))
            lo === nothing || JuMP.set_lower_bound(variable, lo^2)
            hi === nothing || JuMP.set_upper_bound(variable, hi^2)
            w[(busid, terminal)] = variable
        end
    end
    variables[:w] = w

    p_line = Dict{Tuple{Symbol,String,Int},JuMP.VariableRef}()
    q_line = Dict{Tuple{Symbol,String,Int},JuMP.VariableRef}()
    for edge in topology, k in eachindex(edge.parent_map)
        p_line[(edge.family, edge.id, k)] = @variable(model,
            base_name=_l3f_name("l3f_p_branch", edge.family, edge.id, k))
        q_line[(edge.family, edge.id, k)] = @variable(model,
            base_name=_l3f_name("l3f_q_branch", edge.family, edge.id, k))
    end
    variables[:p_line] = p_line; variables[:q_line] = q_line

    p_generator = Dict{Tuple{String,Int},JuMP.VariableRef}()
    q_generator = Dict{Tuple{String,Int},JuMP.VariableRef}()
    for (gid_raw, gen) in sort!(collect(get(net, "generator", Dict())); by=first)
        gid = String(gid_raw); tm = string.(gen["terminal_map"])
        nch = length(gen["p_min"])
        D = _l3f_connection_incidence(String(gen["configuration"]), length(tm), nch)
        vbar = ComplexF64[reference.voltage[(String(gen["bus"]), t)] for t in tm]
        for k in 1:nch
            p = @variable(model, base_name=_l3f_name("l3f_pg", gid, k))
            q = @variable(model, base_name=_l3f_name("l3f_qg", gid, k))
            _l3f_add_bounds!(p, gen, k, "p_min", "p_max")
            _l3f_add_bounds!(q, gen, k, "q_min", "q_max")
            p_generator[(gid, k)] = p; q_generator[(gid, k)] = q
            smax = _l3f_rating(gen, nothing, "s_max", k)
            imax = _l3f_rating(gen, nothing, "i_max", k)
            _l3f_add_power_circle!(model, constraints, :generator_apparent_power,
                (gid, k), p, q, smax)
            w_channel = _l3f_winding_voltage_squared(D, k, vbar,
                String(gen["bus"]), tm, w)
            _l3f_add_current_cone!(model, constraints, :generator_current,
                (gid, k), p, q, w_channel, imax)
        end
    end
    variables[:p_generator] = p_generator; variables[:q_generator] = q_generator

    p_source = Dict{Tuple{String,Int},JuMP.VariableRef}()
    q_source = Dict{Tuple{String,Int},JuMP.VariableRef}()
    for (sid_raw, source) in sort!(collect(get(net, "voltage_source", Dict())); by=first)
        sid = String(sid_raw); tm = string.(source["terminal_map"])
        for k in eachindex(tm)
            p = @variable(model, base_name=_l3f_name("l3f_p_source", sid, k))
            q = @variable(model, base_name=_l3f_name("l3f_q_source", sid, k))
            _l3f_add_bounds!(p, source, k, "p_min", "p_max")
            _l3f_add_bounds!(q, source, k, "q_min", "q_max")
            p_source[(sid, k)] = p; q_source[(sid, k)] = q
            key = (String(source["bus"]), tm[k])
            _l3f_register_constraint!(constraints, :source_voltage, (sid, tm[k]),
                @constraint(model, w[key] == abs2(reference.voltage[key])))
            smax = _l3f_rating(source, nothing, "s_max", k)
            imax = _l3f_rating(source, nothing, "i_max", k)
            _l3f_add_power_circle!(model, constraints, :source_apparent_power,
                (sid, k), p, q, smax)
            _l3f_add_current_cone!(model, constraints, :source_current,
                (sid, k), p, q, w[key], imax)
        end
    end
    variables[:p_source] = p_source; variables[:q_source] = q_source

    for edge in topology
        if edge.family == :transformer
            transformer = net["transformer"][edge.subtype][edge.id]
            transform = _l3f_transformer_oriented_map(edge, transformer)
            vbar_parent = ComplexF64[reference.voltage[(edge.parent, terminal)]
                                     for terminal in edge.parent_map]
            H = connection_power_map(transform, vbar_parent)
            child_p = [p_line[(edge.family, edge.id, k)] for k in eachindex(edge.child_map)]
            child_q = [q_line[(edge.family, edge.id, k)] for k in eachindex(edge.child_map)]
            parent_power = [_l3f_mapped_terminal_power(H, child_p, child_q, k)
                            for k in eachindex(edge.parent_map)]
            for phi in eachindex(edge.parent_map)
                winding = winding_voltage_coefficients(view(transform, phi, :), vbar_parent)
                rhs = JuMP.AffExpr(winding.constant)
                for psi in eachindex(edge.parent_map)
                    JuMP.add_to_expression!(rhs, winding.coefficients[psi],
                        w[(edge.parent, edge.parent_map[psi])])
                end
                _l3f_register_constraint!(constraints, :transformer_voltage_ratio,
                    (edge.subtype, edge.id, phi), @constraint(model,
                        w[(edge.child, edge.child_map[phi])] == rhs))
            end
            parent_is_from = !edge.reversed
            if edge.subtype == "open_delta_regulator"
                _l3f_open_delta_limits!(model, constraints, edge, transformer,
                    child_p, child_q, parent_power, reference, w)
            elseif edge.subtype in _L3F_DELTA_SUBTYPES
                # BMOPF's Yd/Dy `s_rating` is the total bank VA.  Its thermal
                # boundary is applied to the three WYE coils, each at the
                # equal per-coil share S_rating / n_ph.  Delta terminal powers
                # are not coil powers and must not receive a second copy of
                # this constraint.
                n_ph = length(edge.parent_map)
                wye_is_from = edge.subtype == "wye_delta"
                wye_is_parent = wye_is_from == parent_is_from
                wye_p, wye_q = wye_is_parent ?
                    ([parent_power[k][1] for k in eachindex(parent_power)],
                     [parent_power[k][2] for k in eachindex(parent_power)]) :
                    (child_p, child_q)
                s_coil = _l3f_ywye_rating(transformer, edge.subtype, n_ph)
                for phi in 1:n_ph
                    _l3f_add_power_circle!(model, constraints,
                        :transformer_apparent_power,
                        (edge.subtype, edge.id, :wye, phi),
                        wye_p[phi], wye_q[phi], s_coil)
                end
                # Current limits remain side-specific.  BMOPFTools defines the
                # delta-side limit on the bushing/terminal current (not the
                # internal coil arm), so use the corresponding terminal-power
                # expression and phase-ground reference voltage here.
                for phi in eachindex(edge.parent_map)
                    parent_p, parent_q = parent_power[phi]
                    child_ep, child_eq = _l3f_total_endpoint_power(
                        child_p[phi], child_q[phi], -1.0,
                        JuMP.AffExpr(0.0), JuMP.AffExpr(0.0))
                    from_p, from_q = parent_is_from ?
                        (parent_p, parent_q) : (child_ep, child_eq)
                    to_p, to_q = parent_is_from ?
                        (child_ep, child_eq) : (parent_p, parent_q)
                    from_bus = String(transformer["bus_from"])
                    to_bus = String(transformer["bus_to"])
                    from_tm = string.(transformer["terminal_map_from"])
                    to_tm = string.(transformer["terminal_map_to"])
                    _l3f_add_current_cone!(model, constraints, :transformer_current,
                        (edge.subtype, edge.id, "from", phi), from_p, from_q,
                        w[(from_bus, from_tm[phi])],
                        _l3f_scalar_or_indexed(transformer, "i_max_from", phi))
                    _l3f_add_current_cone!(model, constraints, :transformer_current,
                        (edge.subtype, edge.id, "to", phi), to_p, to_q,
                        w[(to_bus, to_tm[phi])],
                        _l3f_scalar_or_indexed(transformer, "i_max_to", phi))
                end
            else
                # Reconstruct power entering each physical endpoint. The
                # to-side expression includes a lowered exciting shunt, when
                # present, so its terminal ratings retain BMOPF semantics.
                for phi in eachindex(edge.parent_map)
                    parent_p, parent_q = parent_power[phi]
                    child_ep, child_eq = _l3f_total_endpoint_power(
                        child_p[phi], child_q[phi], -1.0,
                        JuMP.AffExpr(0.0), JuMP.AffExpr(0.0))
                    from_p, from_q = parent_is_from ?
                        (parent_p, parent_q) : (child_ep, child_eq)
                    to_series_p, to_series_q = parent_is_from ?
                        (child_ep, child_eq) : (parent_p, parent_q)
                    from_bus = String(transformer["bus_from"])
                    to_bus = String(transformer["bus_to"])
                    from_tm = string.(transformer["terminal_map_from"])
                    to_tm = string.(transformer["terminal_map_to"])
                    ps_to, qs_to = _l3f_derived_shunt_power(net, "noload", edge.id,
                        nothing, to_bus, to_tm, reference, w, phi)
                    to_p, to_q = _l3f_total_endpoint_power(
                        to_series_p, to_series_q, 1.0, ps_to, qs_to)
                    _l3f_add_power_circle!(model, constraints, :transformer_apparent_power,
                        (edge.subtype, edge.id, "from", phi), from_p, from_q,
                        get(transformer, "s_rating", nothing))
                    _l3f_add_power_circle!(model, constraints, :transformer_apparent_power,
                        (edge.subtype, edge.id, "to", phi), to_p, to_q,
                        get(transformer, "s_rating", nothing))
                    _l3f_add_current_cone!(model, constraints, :transformer_current,
                        (edge.subtype, edge.id, "from", phi), from_p, from_q,
                        w[(from_bus, from_tm[phi])],
                        _l3f_scalar_or_indexed(transformer, "i_max_from", phi))
                    _l3f_add_current_cone!(model, constraints, :transformer_current,
                        (edge.subtype, edge.id, "to", phi), to_p, to_q,
                        w[(to_bus, to_tm[phi])],
                        _l3f_scalar_or_indexed(transformer, "i_max_to", phi))
                end
            end
            continue
        end
        line = net["line"][edge.id]
        fallback = haskey(line, "linecode") ?
            get(get(net, "linecode", Dict()), String(line["linecode"]), nothing) : nothing
        Z = _l3f_series_matrix(net, line, length(edge.parent_map), edge.id)
        vbar = ComplexF64[reference.voltage[(edge.parent, terminal)]
                          for terminal in edge.parent_map]
        drop = line_drop_coefficients(Z, vbar)
        for phi in eachindex(edge.parent_map)
            rhs = JuMP.AffExpr(0.0)
            JuMP.add_to_expression!(rhs, 1.0, w[(edge.parent, edge.parent_map[phi])])
            for psi in eachindex(edge.parent_map)
                JuMP.add_to_expression!(rhs, -drop.active[phi, psi], p_line[(edge.family, edge.id, psi)])
                JuMP.add_to_expression!(rhs, -drop.reactive[phi, psi], q_line[(edge.family, edge.id, psi)])
            end
            _l3f_register_constraint!(constraints, :line_voltage_drop,
                (edge.id, phi), @constraint(model,
                    w[(edge.child, edge.child_map[phi])] == rhs))
            p, q = p_line[(edge.family, edge.id, phi)], q_line[(edge.family, edge.id, phi)]
            smax = _l3f_rating(line, fallback, "s_max", phi)
            imax = _l3f_rating(line, fallback, "i_max", phi)
            parent_side, child_side = edge.reversed ? ("to", "from") : ("from", "to")
            psh_parent, qsh_parent = _l3f_derived_shunt_power(net, "lineshunt",
                edge.id, parent_side, edge.parent, edge.parent_map,
                reference, w, phi)
            psh_child, qsh_child = _l3f_derived_shunt_power(net, "lineshunt",
                edge.id, child_side, edge.child, edge.child_map,
                reference, w, phi)
            parent_p, parent_q = _l3f_total_endpoint_power(
                p, q, 1.0, psh_parent, qsh_parent)
            child_p, child_q = _l3f_total_endpoint_power(
                p, q, -1.0, psh_child, qsh_child)
            # Stamp both physical ends even without a pi shunt. In exact AC a
            # shunt-free series current gives |S_from|²/w_from = |S_to|²/w_to.
            # Lossless L3F reuses one S while allowing w to drop, so the two
            # inferred currents can differ. Requiring both is a deliberate
            # conservative endpoint contract, not exact BMOPFTools parity.
            _l3f_add_power_circle!(model, constraints, :line_apparent_power,
                (edge.id, parent_side, phi), parent_p, parent_q, smax)
            _l3f_add_power_circle!(model, constraints, :line_apparent_power,
                (edge.id, child_side, phi), child_p, child_q, smax)
            _l3f_add_current_cone!(model, constraints, :line_current,
                (edge.id, parent_side, phi), parent_p, parent_q,
                w[(edge.parent, edge.parent_map[phi])], imax)
            _l3f_add_current_cone!(model, constraints, :line_current,
                (edge.id, child_side, phi), child_p, child_q,
                w[(edge.child, edge.child_map[phi])], imax)
        end
    end

    balance_p = Dict(key => JuMP.AffExpr(0.0) for key in keys(w))
    balance_q = Dict(key => JuMP.AffExpr(0.0) for key in keys(w))
    for edge in topology
        p = [p_line[(edge.family, edge.id, k)] for k in eachindex(edge.child_map)]
        q = [q_line[(edge.family, edge.id, k)] for k in eachindex(edge.child_map)]
        if edge.family == :transformer
            transformer = net["transformer"][edge.subtype][edge.id]
            transform = _l3f_transformer_oriented_map(edge, transformer)
            vbar = ComplexF64[reference.voltage[(edge.parent, terminal)]
                              for terminal in edge.parent_map]
            H = connection_power_map(transform, vbar)
            for terminal in eachindex(edge.parent_map)
                pt, qt = _l3f_mapped_terminal_power(H, p, q, terminal)
                JuMP.add_to_expression!(balance_p[(edge.parent, edge.parent_map[terminal])], -pt)
                JuMP.add_to_expression!(balance_q[(edge.parent, edge.parent_map[terminal])], -qt)
            end
        else
            for terminal in eachindex(edge.parent_map)
                JuMP.add_to_expression!(balance_p[(edge.parent, edge.parent_map[terminal])],
                                        -1.0, p[terminal])
                JuMP.add_to_expression!(balance_q[(edge.parent, edge.parent_map[terminal])],
                                        -1.0, q[terminal])
            end
        end
        for terminal in eachindex(edge.child_map)
            JuMP.add_to_expression!(balance_p[(edge.child, edge.child_map[terminal])],
                                    1.0, p[terminal])
            JuMP.add_to_expression!(balance_q[(edge.child, edge.child_map[terminal])],
                                    1.0, q[terminal])
        end
    end
    load_power = Dict{Tuple{String,Int},Tuple{JuMP.AffExpr,JuMP.AffExpr}}()
    for (lid_raw, load) in get(net, "load", Dict())
        lid = String(lid_raw)
        bus = String(load["bus"]); tm = string.(load["terminal_map"])
        p = Float64.(load["p_nom"]); q = Float64.(load["q_nom"])
        D = _l3f_connection_incidence(String(load["configuration"]), length(tm), length(p))
        vbar = ComplexF64[reference.voltage[(bus, t)] for t in tm]
        H = connection_power_map(D, vbar)
        for channel in eachindex(p)
            load_power[(lid, channel)] = _l3f_load_channel_power(
                load, D, vbar, bus, tm, w, channel)
        end
        for terminal in eachindex(tm), channel in eachindex(p)
            pc, qc = load_power[(lid, channel)]
            JuMP.add_to_expression!(balance_p[(bus, tm[terminal])],
                -H.real_part[terminal, channel], pc)
            JuMP.add_to_expression!(balance_p[(bus, tm[terminal])],
                H.imag_part[terminal, channel], qc)
            JuMP.add_to_expression!(balance_q[(bus, tm[terminal])],
                -H.imag_part[terminal, channel], pc)
            JuMP.add_to_expression!(balance_q[(bus, tm[terminal])],
                -H.real_part[terminal, channel], qc)
        end
    end
    variables[:load_power] = load_power
    shunt_power = Dict{Tuple{String,Int},Tuple{JuMP.AffExpr,JuMP.AffExpr}}()
    for (sid_raw, shunt) in get(net, "shunt", Dict())
        sid = String(sid_raw); bus = String(shunt["bus"]); tm = string.(shunt["terminal_map"])
        for phi in eachindex(tm)
            ps, qs = _l3f_shunt_power(shunt, bus, tm, reference, w, phi)
            shunt_power[(sid, phi)] = (ps, qs)
            JuMP.add_to_expression!(balance_p[(bus, tm[phi])], -ps)
            JuMP.add_to_expression!(balance_q[(bus, tm[phi])], -qs)
        end
    end
    variables[:shunt_power] = shunt_power
    for (gid_raw, gen) in get(net, "generator", Dict())
        gid = String(gid_raw); bus = String(gen["bus"]); tm = string.(gen["terminal_map"])
        nch = length(gen["p_min"])
        D = _l3f_connection_incidence(String(gen["configuration"]), length(tm), nch)
        vbar = ComplexF64[reference.voltage[(bus, t)] for t in tm]
        H = connection_power_map(D, vbar)
        for terminal in eachindex(tm), channel in 1:nch
            JuMP.add_to_expression!(balance_p[(bus, tm[terminal])],
                H.real_part[terminal, channel], p_generator[(gid, channel)])
            JuMP.add_to_expression!(balance_p[(bus, tm[terminal])],
                -H.imag_part[terminal, channel], q_generator[(gid, channel)])
            JuMP.add_to_expression!(balance_q[(bus, tm[terminal])],
                H.imag_part[terminal, channel], p_generator[(gid, channel)])
            JuMP.add_to_expression!(balance_q[(bus, tm[terminal])],
                H.real_part[terminal, channel], q_generator[(gid, channel)])
        end
    end
    for (sid_raw, source) in get(net, "voltage_source", Dict())
        sid = String(sid_raw); bus = String(source["bus"]); tm = string.(source["terminal_map"])
        for k in eachindex(tm)
            JuMP.add_to_expression!(balance_p[(bus, tm[k])], 1.0, p_source[(sid, k)])
            JuMP.add_to_expression!(balance_q[(bus, tm[k])], 1.0, q_source[(sid, k)])
        end
    end
    for key in sort!(collect(keys(w)))
        _l3f_register_constraint!(constraints, :nodal_active_balance, key,
            @constraint(model, balance_p[key] == 0.0))
        _l3f_register_constraint!(constraints, :nodal_reactive_balance, key,
            @constraint(model, balance_q[key] == 0.0))
    end

    objective = JuMP.AffExpr(0.0)
    if options.objective == :cost
        for (gid_raw, gen) in get(net, "generator", Dict())
            gid = String(gid_raw)
            for k in eachindex(gen["p_min"])
                JuMP.add_to_expression!(objective,
                    _l3f_cost_coefficient(gen, k, "generator", gid), p_generator[(gid, k)])
            end
        end
        for (sid_raw, source) in get(net, "voltage_source", Dict())
            sid = String(sid_raw)
            for k in eachindex(source["terminal_map"])
                JuMP.add_to_expression!(objective,
                    _l3f_cost_coefficient(source, k, "voltage source", sid), p_source[(sid, k)])
            end
        end
    elseif options.objective == :source_import
        for variable in values(p_source)
            JuMP.add_to_expression!(objective, 1.0, variable)
        end
    end
    @objective(model, Min, objective)
    L3FBuild(model, variables, constraints, physical_reference, reference,
             report, options, physical_network, net, topology, bases)
end

"""
    build_l3f_opf(net, optimizer=Clarabel.Optimizer;
                  options=L3FOptions(), reference=nothing)

Build the L3F-BMOPF lossless radial LP/SOCP. Explicit-neutral inputs are
Kron-reduced on a copy when enabled. Inapplicable networks raise
[`L3FInapplicableError`](@ref), whose report contains stable diagnostics.
Input and extracted results are SI. `options.per_unit` selects only the model's
working coordinates.
"""
function build_l3f_opf(net, optimizer=Clarabel.Optimizer;
                       options::L3FOptions=L3FOptions(), reference=nothing)
    prepared = _l3f_prepare(net; options, reference)
    is_l3f_applicable(prepared.applicability) ||
        throw(L3FInapplicableError(prepared.applicability))
    physical_reference = try
        _l3f_reference(prepared.network, prepared.topology, reference,
                       options.reference_policy)
    catch err
        throw(L3FInapplicableError(_l3f_report_error(prepared.applicability,
            "E.L3F.REFERENCE_MISSING", sprint(showerror, err))))
    end
    working, working_reference, bases = _l3f_working_coordinates(
        prepared.network, physical_reference, options)
    _l3f_build_model(working, prepared.topology, working_reference,
                     prepared.applicability, optimizer, options;
                     physical_network=prepared.network,
                     physical_reference=physical_reference, bases)
end

"""Classify a continuous L3F model as `:LP`, `:QP`, or `:SOCP`."""
function l3f_model_class(build::L3FBuild)
    any(v -> JuMP.is_binary(v) || JuMP.is_integer(v), JuMP.all_variables(build.model)) &&
        return :unsupported
    objective_type = JuMP.objective_function_type(build.model)
    quadratic = objective_type <: JuMP.GenericQuadExpr
    quadratic || objective_type <: JuMP.GenericAffExpr ||
        objective_type <: JuMP.VariableRef || objective_type <: Real || return :unsupported
    conic = false
    for (function_type, set_type) in JuMP.list_of_constraint_types(build.model)
        if set_type <: JuMP.MOI.SecondOrderCone ||
           set_type <: JuMP.MOI.RotatedSecondOrderCone
            conic = true
        elseif !(function_type <: JuMP.GenericAffExpr ||
                 function_type <: JuMP.VariableRef)
            return :unsupported
        end
    end
    conic ? :SOCP : quadratic ? :QP : :LP
end

function _l3f_extract(build::L3FBuild, outcome::SolveOutcome)
    publish = outcome.optimal
    value(variable) = publish ? JuMP.value(variable) : NaN
    voltage_scale(bus) = build.bases === nothing ? 1.0 : build.bases.v_base[bus]
    power_scale() = build.bases === nothing ? 1.0 : build.options.s_base
    buses = Dict{String,Any}()
    for (busid_raw, bus) in get(build.network, "bus", Dict())
        busid = String(busid_raw); terminals = Dict{String,Any}()
        for terminal in string.(bus["terminal_names"])
            wv = value(build.variables[:w][(busid, terminal)]) * voltage_scale(busid)^2
            terminals[terminal] = Dict{String,Any}(
                "w" => wv,
                "vm" => isfinite(wv) ? sqrt(max(0.0, wv)) : NaN,
                "reference_angle" => angle(build.reference.voltage[(busid, terminal)]),
            )
        end
        buses[busid] = terminals
    end
    lines = Dict{String,Any}(); transformers = Dict{String,Any}()
    for edge in build.topology
        output = Dict{String,Any}(
            "parent" => edge.parent, "child" => edge.child,
            "reversed_from_input" => edge.reversed,
            "terminal_map_parent" => edge.parent_map,
            "terminal_map_child" => edge.child_map,
            "p" => [value(build.variables[:p_line][(edge.family, edge.id, k)]) * power_scale() for k in eachindex(edge.parent_map)],
            "q" => [value(build.variables[:q_line][(edge.family, edge.id, k)]) * power_scale() for k in eachindex(edge.parent_map)],
        )
        if edge.family == :line
            lines[edge.id] = output
        else
            output["subtype"] = edge.subtype
            output["effective_ratio_from_to"] =
                _l3f_transformer_neff(edge.subtype,
                    build.network["transformer"][edge.subtype][edge.id])
            transformers[edge.id] = output
        end
    end
    generators = Dict{String,Any}()
    for (gid_raw, gen) in get(build.network, "generator", Dict())
        gid = String(gid_raw); n = length(gen["p_min"])
        generators[gid] = Dict{String,Any}(
            "terminal_map" => string.(gen["terminal_map"]),
            "pg" => [value(build.variables[:p_generator][(gid, k)]) * power_scale() for k in 1:n],
            "qg" => [value(build.variables[:q_generator][(gid, k)]) * power_scale() for k in 1:n],
        )
    end
    sources = Dict{String,Any}()
    for (sid_raw, source) in get(build.network, "voltage_source", Dict())
        sid = String(sid_raw); n = length(source["terminal_map"])
        sources[sid] = Dict{String,Any}(
            "terminal_map" => string.(source["terminal_map"]),
            "pg" => [value(build.variables[:p_source][(sid, k)]) * power_scale() for k in 1:n],
            "qg" => [value(build.variables[:q_source][(sid, k)]) * power_scale() for k in 1:n],
        )
    end
    buses, lines, transformers, generators, sources
end

function _l3f_physical_objective(build::L3FBuild, outcome::SolveOutcome)
    outcome.optimal || return NaN
    value = JuMP.objective_value(build.model)
    build.bases === nothing && return value
    build.options.objective == :source_import ? value * build.options.s_base : value
end

function _l3f_fix_dispatch!(net, result::L3FResult)
    for (gid, dispatch) in result.generators
        haskey(get(net, "generator", Dict()), gid) || continue
        gen = net["generator"][gid]
        gen["p_min"] = copy(dispatch["pg"]); gen["p_max"] = copy(dispatch["pg"])
        gen["q_min"] = copy(dispatch["qg"]); gen["q_max"] = copy(dispatch["qg"])
    end
    net
end

"""
    l3f_reference_from_powerflow(net; options=L3FOptions(), dispatch=nothing,
                                 nonlinear_optimizer=Ipopt.Optimizer,
                                 solver_options=()) -> L3FReferenceState

Build the linearization reference from a converged nonlinear power flow instead
of the flat propagated profile.

The default reference carries the source phasors outward through the topology
applying ratios but no line drop, so on a loaded feeder every coefficient is
formed at roughly nominal voltage. Linearizing at an actual operating point
tightens the fixed-angle closure and the frozen channel-to-terminal split, which
is the direct remedy for accuracy loss under heavy loading or strong unbalance.

The returned state has provenance `:power_flow` and can be passed straight back
as the `reference` argument, giving the standard successive-linearization loop:

```julia
first = solve_l3f_opf(net, Clarabel.Optimizer; options)
better = l3f_reference_from_powerflow(net; options, dispatch=first)
second = solve_l3f_opf(net, Clarabel.Optimizer; options, reference=better)
```

A power flow is a determined problem, so BMOPFTools requires every generator to
be a fixed setpoint. Pass an [`L3FResult`](@ref) as `dispatch` to pin them at a
previous solution, or fix `p_min == p_max` and `q_min == q_max` in the input.
Transformer nameplates are stripped from the private working copy first: a
reference is a linearization point, not a feasibility claim, and an overloaded
coil would otherwise make the reference solve fail rather than report the
operating point being linearized at.

Source terminals are restored to their exactly declared phasors, so the result
always satisfies the reference/source consistency check.
"""
function l3f_reference_from_powerflow(net; options::L3FOptions=L3FOptions(),
                                      dispatch=nothing,
                                      nonlinear_optimizer=Ipopt.Optimizer,
                                      solver_options=())
    # This helper produces the explicit reference and its provenance.
    # Preparation must not require neutral-reduction provenance before the PF,
    # including for an already-reduced case or for
    # `reference_policy=:source_propagated`.  An explicit policy also needs a
    # provisional source-propagated policy because the helper has no reference
    # argument until it returns its PF state.
    prepare_policy = options.reference_policy == :explicit ?
        :source_propagated : options.reference_policy
    prepare_options = _l3f_with_options(options;
        reference_policy=prepare_policy, require_neutral_provenance=false)
    prepared = _l3f_prepare(net; options=prepare_options)
    is_l3f_applicable(prepared.applicability) ||
        throw(L3FInapplicableError(prepared.applicability))
    working = deepcopy(prepared.network)
    dispatch === nothing || _l3f_fix_dispatch!(working, dispatch)
    for (_, table) in get(working, "transformer", Dict())
        table isa AbstractDict || continue
        for (_, transformer) in table
            transformer isa AbstractDict && delete!(transformer, "s_rating")
        end
    end
    for (gid, generator) in get(working, "generator", Dict())
        generator isa AbstractDict || continue
        lo, hi = Float64.(get(generator, "p_min", Float64[])), Float64.(get(generator, "p_max", Float64[]))
        qlo, qhi = Float64.(get(generator, "q_min", Float64[])), Float64.(get(generator, "q_max", Float64[]))
        lo == hi && qlo == qhi || throw(ArgumentError(
            "generator '$gid' is a P/Q range, but a reference power flow is a " *
            "determined problem. Pass `dispatch=<L3FResult>` to pin it at a " *
            "previous solution, or set p_min == p_max and q_min == q_max."))
    end

    pf = BMOPFTools.solve_pf(working; optimizer=nonlinear_optimizer, solver_options)
    termination = String(get(pf, "termination_status", "UNKNOWN"))
    termination in ("OPTIMAL", "LOCALLY_SOLVED") || throw(ErrorException(
        "reference power flow did not converge (status $termination)"))

    voltage = Dict{Tuple{String,String},ComplexF64}()
    for (busid, bus) in get(prepared.network, "bus", Dict())
        pf_bus = get(get(pf, "bus", Dict()), String(busid), nothing)
        pf_bus isa AbstractDict || throw(ErrorException(
            "reference power flow returned no voltages for bus '$busid'"))
        for terminal in string.(get(bus, "terminal_names", String[]))
            entry = get(pf_bus, terminal, nothing)
            entry isa AbstractDict || throw(ErrorException(
                "reference power flow returned no voltage for $busid.$terminal"))
            voltage[(String(busid), terminal)] =
                ComplexF64(Float64(entry["vr"]), Float64(entry["vi"]))
        end
    end
    # The source is fixed data, not a solved quantity; restoring the declared
    # phasors keeps the state exactly consistent with the network it describes.
    for (_, source) in get(prepared.network, "voltage_source", Dict())
        bus = String(source["bus"])
        tm = string.(source["terminal_map"])
        vm = Float64.(source["v_magnitude"])
        va = Float64.(source["v_angle"])
        for k in eachindex(tm)
            voltage[(bus, tm[k])] = vm[k] * cis(va[k])
        end
    end
    state = L3FReferenceState(voltage, :power_flow, _l3f_reference_hash(voltage),
                              Symbol(termination))
    _l3f_validate_reference(state, prepared.network)
end

"""
    validate_l3f_solution(result; nonlinear_optimizer=Ipopt.Optimizer,
                          voltage_tolerance=nothing)

Reapply the L3F generator dispatch to `result.network` — the SI, Kron-reduced
snapshot the model was built from — and run BMOPFTools' nonlinear power flow.

`status` is `"replayed"` when the nonlinear solve converged, `"failed"` when it
did not, and `"not_run"` when the L3F solve was not optimal. **It is a statement
about the replay, not about accuracy**: the linearization omits series losses,
so a converged replay always differs from the linear solution by some margin.
Supply `voltage_tolerance` (in volts) to have that margin judged; the result
then carries `within_tolerance`. Physical limits are never certified here.

`replayed_network` is `"as_supplied"` normally and `"projected"` when
`unsupported=:approximate` substituted a load law, tap, or limit. In that case
both the linear model and the replay describe the *substituted* network, so the
reported error measures the linearization but not the projection.
"""
function validate_l3f_solution(result::L3FResult;
                               nonlinear_optimizer=Ipopt.Optimizer,
                               voltage_tolerance=nothing)
    result.solve.optimal || return Dict{String,Any}(
        "status" => "not_run", "reason" => "L3F solve was not optimal",
        "physical_limits" => "unassessed")

    working = _l3f_fix_dispatch!(deepcopy(result.network), result)
    # Under `unsupported=:approximate` the snapshot is the projected network, so
    # the replay measures the linearization error but NOT the projection error:
    # both sides are solving the same substituted physics. The flag says which.
    replayed_network = any(f -> startswith(f.code, "A.L3F."),
                           result.applicability.findings) ? "projected" : "as_supplied"
    try
        pf = BMOPFTools.solve_pf(working; optimizer=nonlinear_optimizer)
        termination = String(get(pf, "termination_status", "UNKNOWN"))
        solved = termination in ("OPTIMAL", "LOCALLY_SOLVED")
        errors = Float64[]
        if solved
            for (bus, terminals) in result.buses
                pf_bus = get(get(pf, "bus", Dict()), bus, Dict())
                for (terminal, values) in terminals
                    pf_terminal = get(pf_bus, terminal, nothing)
                    pf_terminal isa AbstractDict || continue
                    vm = hypot(Float64(pf_terminal["vr"]), Float64(pf_terminal["vi"]))
                    push!(errors, abs(vm - Float64(values["vm"])))
                end
            end
        end
        maximum_error = isempty(errors) ? NaN : maximum(errors)
        out = Dict{String,Any}(
            "status" => solved ? "replayed" : "failed",
            "nonlinear_solve_status" => termination,
            "maximum_voltage_magnitude_error" => maximum_error,
            "rms_voltage_magnitude_error" => isempty(errors) ? NaN : sqrt(sum(abs2, errors) / length(errors)),
            "compared_terminals" => length(errors),
            "physical_limits" => "unassessed",
            "replayed_network" => replayed_network,
            "claim" => "nonlinear replay comparison, not a physical-feasibility certificate",
        )
        if voltage_tolerance !== nothing
            out["voltage_tolerance"] = Float64(voltage_tolerance)
            out["within_tolerance"] = solved && isfinite(maximum_error) &&
                maximum_error <= Float64(voltage_tolerance)
        end
        out
    catch err
        Dict{String,Any}(
            "status" => "failed",
            "nonlinear_solve_status" => "ERROR",
            "error" => sprint(showerror, err),
            "physical_limits" => "unassessed",
        )
    end
end

"""
    solve_l3f_opf(net, optimizer=Clarabel.Optimizer; options=L3FOptions(),
                  reference=nothing, nonlinear_optimizer=Ipopt.Optimizer,
                  solver_options=())

Build and solve the L3F-BMOPF LP/SOCP. The result includes applicability,
reference provenance, stable semantic outputs, and optional nonlinear replay.
"""
function solve_l3f_opf(net, optimizer=Clarabel.Optimizer;
                       options::L3FOptions=L3FOptions(), reference=nothing,
                       nonlinear_optimizer=Ipopt.Optimizer, solver_options=(),
                       voltage_tolerance=nothing)
    build = build_l3f_opf(net, optimizer; options, reference)
    _set_solver_options!(build.model, solver_options)
    JuMP.optimize!(build.model)
    outcome = _solve_outcome(build.model)
    status = SolveStatus(outcome)
    buses, lines, transformers, generators, sources = _l3f_extract(build, outcome)
    result = L3FResult(buses, lines, transformers, generators, sources,
        _l3f_physical_objective(build, outcome),
        Dict{String,Any}(
            "name" => "L3F-BMOPF", "version" => "0.1-prototype",
            "problem_class" => String(l3f_model_class(build)),
            "working_units" => options.per_unit ? "per_unit" : "SI",
            "per_unit" => options.per_unit,
            "s_base" => options.s_base,
            "result_units" => "SI",
            "series_losses" => "omitted",
            "unsupported_policy" => String(options.unsupported),
            "lowered" => build.applicability.lowered,
            "reference_provenance" => String(build.reference.provenance),
            "reference_hash" => build.reference.source_hash,
            "voltage_angles" => "fixed reference coefficients; not decision variables",
        ), build.reference, build.applicability,
        Dict{String,Any}("status" => "not_requested"),
        build.network, status)
    if options.validate_nonlinear
        validation = validate_l3f_solution(result; nonlinear_optimizer, voltage_tolerance)
        reference = L3FReferenceState(result.reference.voltage, result.reference.provenance,
            result.reference.source_hash, Symbol(validation["status"]))
        result = L3FResult(result.buses, result.lines, result.transformers, result.generators,
            result.sources, result.objective, result.formulation,
            reference, result.applicability, validation,
            result.network, result.solve)
    end
    result
end
