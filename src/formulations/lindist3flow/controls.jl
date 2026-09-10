# Restricted generator and IBR controls for the one-shot LinDist3Flow model.
# These helpers deliberately contain no power-flow, warm-start, IV, droop, or
# controller iteration. Supported IBRs are lowered to the generator vocabulary
# on the compiler's private network copy, and all controls stamp affine or SOC
# constraints in S-W coordinates.

_l3f_control_value(x, k::Int, n::Int, field::String) = x isa Real ? Float64(x) :
    x isa AbstractVector && length(x) in (1, n) ? Float64(x[length(x) == 1 ? 1 : k]) :
    throw(DimensionMismatch("$field must be scalar or have $n entries"))

function _l3f_finite_control_value(x, field)
    value = Float64(x)
    isfinite(value) || throw(ArgumentError("$field must be finite"))
    value
end

_l3f_valid_number(x) = x isa Real && isfinite(Float64(x))
_l3f_valid_number_vector(x) = x isa AbstractVector && !isempty(x) &&
    all(_l3f_valid_number, x)
_l3f_permissive_wellformed(x) = x isa Real ? _l3f_valid_number(x) :
    x isa AbstractString || x isa Bool || x === nothing ? true :
    x isa AbstractVector ? all(_l3f_permissive_wellformed, x) :
    x isa AbstractDict ? all(kv -> kv.first isa AbstractString &&
        _l3f_permissive_wellformed(kv.second), pairs(x)) : false

"""
Drop well-formed unsupported IBR physics for `unsupported=:permissive`.

This operates only on the compiler's private copy. Power boxes, phase ratings,
costs, signed fixed PF, voltage targets, and aggregate capability remain. An
isolated lossless shared-DC link also retains its static net-P bounds. Values
that are malformed are deliberately left in place so ordinary validation still
reports an error rather than converting bad data into a permissive warning.
"""
function _l3f_permissive_restricted_controls!(findings, net)
    ibrs = get(net, "ibr", nothing)
    ibrs isa AbstractDict || return net
    profiles = get(net, "control_profile", Dict())
    dropped_profiles = Set{String}()
    for (iid_raw, raw) in ibrs
        raw isa AbstractDict || continue
        iid = String(iid_raw)
        dropped = Dict{String,Any}()
        drop_if!(field, predicate) = begin
            haskey(raw, field) && predicate(raw[field]) || return
            dropped[field] = deepcopy(raw[field])
            delete!(raw, field)
        end

        drop_if!("r_filter", _l3f_valid_number_vector)
        drop_if!("x_filter", _l3f_valid_number_vector)
        drop_if!("b_filter_shunt", _l3f_valid_number)
        drop_if!("v_ref_internal", x -> _l3f_valid_number(x) || _l3f_valid_number_vector(x))
        drop_if!("voltage_aggregation", x -> x isa AbstractString)
        drop_if!("time_series", x -> x isa AbstractDict && _l3f_permissive_wellformed(x))
        positive_number(x) = _l3f_valid_number(x) && Float64(x) > 0
        drop_if!("i_neutral_max", positive_number)
        drop_if!("neutral_i_max", positive_number)
        drop_if!("neutral_current_control", x ->
            (x isa AbstractDict || x isa AbstractString) && _l3f_permissive_wellformed(x))
        # FOUR_LEG i_max may carry a final neutral-conductor rating. When the
        # bus conventions prove that terminal to be neutral, retain all phase
        # current caps and drop only the unsupported neutral entry.
        topology = uppercase(String(get(raw, "topology", "FOUR_LEG")))
        if topology in ("FOUR_LEG", "SINGLE_PHASE") &&
           get(raw, "i_max", nothing) isa AbstractVector
            tm = string.(get(raw, "terminal_map", String[]))
            imax = raw["i_max"]
            neutral_labels = _l3f_neutral_labels(net)
            neutral_positions = findall(in(neutral_labels), tm)
            if length(imax) == length(tm) && length(neutral_positions) == 1 &&
               _l3f_valid_number_vector(imax) && all(x -> Float64(x) > 0, imax)
                ni = only(neutral_positions)
                retained = Float64.(imax[setdiff(eachindex(imax), [ni])])
                # Both conductors carry the same SINGLE_PHASE current, so the
                # tighter of phase and return ratings is the preserved cap.
                if topology == "SINGLE_PHASE"
                    effective = min(minimum(retained), Float64(imax[ni]))
                    dropped["i_max_neutral_combined"] = Dict(
                        "original" => deepcopy(imax[ni]), "effective_cap" => effective)
                    raw["i_max"] = [effective]
                else
                    dropped["i_max_neutral"] = deepcopy(imax[ni])
                    raw["i_max"] = retained
                end
            end
        end
        if get(raw, "grid_forming", false) === true
            dropped["grid_forming"] = true
            delete!(raw, "grid_forming")
        end

        coupled = get(raw, "dc_link_coupled", false)
        if coupled isa Bool
            # The isolated link's p_dc bounds remain a safe affine aggregate
            # constraint. External buses and dynamic/master controls do not.
            drop_if!("dc_bus", x -> x isa AbstractString)
            drop_if!("dc_terminal_map", x -> x isa AbstractVector && all(y -> y isa AbstractString, x))
            drop_if!("dc_control", x -> x isa AbstractString)
            for field in ("dc_v_set", "dc_p_ref", "dc_droop", "dc_deadband")
                drop_if!(field, _l3f_valid_number)
            end
            if !coupled
                for field in ("p_dc_min", "p_dc_max")
                    drop_if!(field, _l3f_valid_number)
                end
            end
        end

        cp_id = get(raw, "control_profile", nothing)
        if cp_id isa AbstractString
            cp = get(profiles, String(cp_id), nothing)
            if cp isa AbstractDict
                laws = Set(String.(keys(cp)))
                pure_pf = laws == Set(["power_factor"]) &&
                    get(cp, "power_factor", nothing) isa AbstractDict
                well_formed_unsupported = !pure_pf && !isempty(laws) &&
                    all(value -> value isa AbstractDict, values(cp)) &&
                    _l3f_permissive_wellformed(cp)
                if well_formed_unsupported
                    dropped["control_profile"] = Dict(
                        "id" => String(cp_id), "value" => deepcopy(cp))
                    delete!(raw, "control_profile")
                    push!(dropped_profiles, String(cp_id))
                end
            end
        end
        # Name the fields that survived alongside those that did not: the
        # useful question after a permissive drop is what the solved device
        # still enforces, which a constant flag cannot answer.
        retained = sort!([field for field in String.(keys(raw))
                          if !startswith(field, "_")])
        isempty(dropped) || _l3f_warning!(findings,
            "A.L3F.IBR_FIELDS_DROPPED", :ibr, iid,
            "unsupported IBR fields were dropped under unsupported=:permissive",
            evidence=Dict("original_fields" => dropped,
                          "retained_fields" => retained))
    end

    if profiles isa AbstractDict
        referenced = Set(String(raw["control_profile"]) for raw in values(ibrs)
            if raw isa AbstractDict && get(raw, "control_profile", nothing) isa AbstractString)
        for id in setdiff(dropped_profiles, referenced)
            pop!(profiles, id, nothing)
        end
    end
    net
end

"""Lower the explicitly supported native IBR slice into generators in-place."""
function _l3f_lower_restricted_controls!(findings, net)
    ibrs = get(net, "ibr", nothing)
    ibrs isa AbstractDict || return net
    generators = _l3f_table!(net, "generator")
    profiles = get(net, "control_profile", Dict())
    used_profiles = Set{String}()
    for (iid_raw, raw) in sort!(collect(ibrs); by=first)
        iid = String(iid_raw)
        raw isa AbstractDict || begin
            _l3f_error!(findings, "E.L3F.IBR_INVALID", :ibr, iid, "IBR must be an object")
            continue
        end
        try
            supported = Set(("bus", "terminal_map", "topology", "prime_mover",
                "s_max", "i_max", "p_avail", "p_min", "p_max", "q_min", "q_max",
                "cost", "control_profile", "dc_link_coupled", "p_dc_min", "p_dc_max",
                "aggregate_p_min", "aggregate_p_max", "aggregate_q_min",
                "aggregate_q_max", "aggregate_s_max", "v_target", "fixed_pf",
                "grid_forming"))
            unknown = setdiff(Set(String.(keys(raw))), supported)
            isempty(unknown) || throw(ArgumentError(
                "unsupported IBR fields: $(join(sort!(collect(unknown)), ", "))"))
            forbidden = [field for field in
                ("r_filter", "x_filter", "b_filter_shunt", "v_ref_internal",
                 "voltage_aggregation", "time_series", "dc_terminal_map",
                 "dc_v_set", "dc_p_ref", "dc_droop", "dc_deadband", "dc_control") if haskey(raw, field)]
            isempty(forbidden) || throw(ArgumentError(
                "unsupported IBR fields: $(join(forbidden, ", "))"))
            haskey(raw, "grid_forming") && !(raw["grid_forming"] isa Bool) &&
                throw(ArgumentError("grid_forming must be Bool"))
            get(raw, "grid_forming", false) == true &&
                throw(ArgumentError("grid-forming IBR controls are unsupported"))
            haskey(raw, "dc_link_coupled") && !(raw["dc_link_coupled"] isa Bool) &&
                throw(ArgumentError("dc_link_coupled must be Bool"))
            topology = uppercase(String(get(raw, "topology", "FOUR_LEG")))
            tm = string.(get(raw, "terminal_map", String[]))
            configuration, gtm, nch = if topology == "FOUR_LEG"
                !isempty(tm) || throw(ArgumentError("FOUR_LEG requires retained phase terminals"))
                isempty(intersect(Set(tm), _l3f_neutral_labels(net))) ||
                    throw(ArgumentError("FOUR_LEG terminal_map still contains a declared neutral"))
                ("WYE", tm, length(tm))
            elseif topology == "SINGLE_PHASE"
                length(tm) in (1, 2) || throw(ArgumentError("SINGLE_PHASE requires one retained phase or a terminal pair"))
                isempty(intersect(Set(tm), _l3f_neutral_labels(net))) ||
                    throw(ArgumentError("SINGLE_PHASE terminal_map still contains a declared neutral"))
                (length(tm) == 1 ? "WYE" : "SINGLE_PHASE", tm, 1)
            elseif topology == "THREE_LEG"
                length(tm) == 3 || throw(ArgumentError("THREE_LEG requires three terminals"))
                ("DELTA", tm, 3)
            else
                throw(ArgumentError("unsupported IBR topology '$topology'"))
            end
            haskey(generators, iid) && throw(ArgumentError(
                "IBR id collides with generator id '$iid'"))
            gen = Dict{String,Any}(
                "bus" => String(raw["bus"]), "terminal_map" => gtm,
                "configuration" => configuration, "_l3f_original_ibr" => iid,
                "_l3f_original_ibr_topology" => topology)
            sraw = get(raw, "s_max", nothing)
            sraw isa AbstractVector || throw(ArgumentError(
                "restricted native IBR requires vector s_max"))
            length(sraw) == nch || throw(DimensionMismatch(
                "s_max must have one entry per retained power channel"))
            svals = Float64.(sraw)
            all(x -> isfinite(x) && x > 0, svals) || throw(ArgumentError(
                "s_max entries must be positive and finite"))
            gen["s_max"] = svals
            gen["cost"] = haskey(raw, "cost") ?
                [_l3f_control_value(raw["cost"], k, nch, "cost") for k in 1:nch] :
                zeros(nch)
            for (field, fallback) in (("p_min", -svals), ("p_max", svals),
                                      ("q_min", -svals), ("q_max", svals))
                source = get(raw, field, fallback)
                gen[field] = [_l3f_control_value(source, k, nch, field) for k in 1:nch]
            end
            for field in ("i_max",)
                haskey(raw, field) || continue
                source = raw[field]
                # A FOUR_LEG neutral i_max has no generator-channel equivalent.
                field == "i_max" && source isa AbstractVector && length(source) == nch + 1 &&
                    throw(ArgumentError("FOUR_LEG neutral i_max cannot be preserved by generator lowering"))
                gen[field] = [_l3f_control_value(source, k, nch, field) for k in 1:nch]
            end
            if get(raw, "dc_link_coupled", false) == true
                haskey(raw, "dc_bus") && throw(ArgumentError("external DC buses are unsupported"))
                gen["aggregate_p_min"] = _l3f_finite_control_value(
                    get(raw, "p_dc_min", 0.0), "p_dc_min")
                gen["aggregate_p_max"] = _l3f_finite_control_value(
                    get(raw, "p_dc_max", 0.0), "p_dc_max")
            elseif any(haskey(raw, key) for key in ("p_dc_min", "p_dc_max", "dc_bus", "dc_control"))
                throw(ArgumentError("DC-side fields require a supported isolated dc_link_coupled IBR"))
            end
            cp_id = get(raw, "control_profile", nothing)
            cp_id !== nothing && haskey(raw, "fixed_pf") && throw(ArgumentError(
                "IBR cannot declare both fixed_pf and control_profile"))
            if cp_id !== nothing
                cp_id = String(cp_id); cp = get(profiles, cp_id, nothing)
                cp isa AbstractDict || throw(ArgumentError("unknown control_profile '$cp_id'"))
                laws = Set(String.(keys(cp)))
                laws == Set(["power_factor"]) || throw(ArgumentError(
                    "only a pure power_factor control_profile is supported"))
                pfobj = cp["power_factor"]
                pfobj isa AbstractDict || throw(ArgumentError("power_factor must be an object"))
                pf = _l3f_finite_control_value(get(pfobj, "pf", NaN), "power_factor.pf")
                0 < abs(pf) <= 1 || throw(ArgumentError("power_factor.pf must satisfy 0 < abs(pf) <= 1"))
                gen["fixed_pf"] = pf
                push!(used_profiles, cp_id)
            end
            # Local aggregate P bounds intersect isolated DC-link limits.
            for field in ("aggregate_q_min", "aggregate_q_max", "aggregate_s_max",
                          "v_target", "fixed_pf")
                haskey(raw, field) && (gen[field] = deepcopy(raw[field]))
            end
            haskey(raw, "aggregate_p_min") && (gen["aggregate_p_min"] = max(
                Float64(get(gen, "aggregate_p_min", -Inf)), Float64(raw["aggregate_p_min"])))
            haskey(raw, "aggregate_p_max") && (gen["aggregate_p_max"] = min(
                Float64(get(gen, "aggregate_p_max", Inf)), Float64(raw["aggregate_p_max"])))
            generators[iid] = gen
            _l3f_info!(findings, "L.L3F.IBR_TO_GENERATOR", :ibr, iid,
                "restricted $topology IBR lowered to an equivalent $configuration generator",
                evidence=Dict("generator_id" => iid, "original_topology" => topology))
        catch err
            _l3f_error!(findings, "E.L3F.IBR_UNSUPPORTED", :ibr, iid, sprint(showerror, err))
        end
    end
    empty!(ibrs)
    if profiles isa AbstractDict
        for id in used_profiles
            pop!(profiles, id, nothing)
        end
    end
    net
end

"""Validate local generator controls that extend the base BMOPF schema."""
function _l3f_validate_restricted_controls!(findings, net)
    for (gid, gen) in get(net, "generator", Dict())
        n = length(get(gen, "p_min", Any[]))
        try
            for (lo, hi) in (("aggregate_p_min", "aggregate_p_max"),
                             ("aggregate_q_min", "aggregate_q_max"))
                lov = haskey(gen, lo) ? _l3f_finite_control_value(gen[lo], lo) : -Inf
                hiv = haskey(gen, hi) ? _l3f_finite_control_value(gen[hi], hi) : Inf
                lov <= hiv || throw(ArgumentError("$lo exceeds $hi"))
            end
            if haskey(gen, "aggregate_s_max")
                s = _l3f_finite_control_value(gen["aggregate_s_max"], "aggregate_s_max")
                s > 0 || throw(ArgumentError("aggregate_s_max must be positive"))
            end
            if haskey(gen, "fixed_pf")
                pf = _l3f_finite_control_value(gen["fixed_pf"], "fixed_pf")
                0 < abs(pf) <= 1 || throw(ArgumentError("fixed_pf must satisfy 0 < abs(pf) <= 1"))
            end
            if haskey(gen, "v_target")
                pmin, pmax = Float64.(gen["p_min"]), Float64.(gen["p_max"])
                qmin, qmax = Float64.(gen["q_min"]), Float64.(gen["q_max"])
                all(pmin .== pmax) || throw(ArgumentError("v_target requires fixed active power"))
                all(qmin .<= qmax) && any(qmin .< qmax) || throw(ArgumentError(
                    "v_target requires free, finitely bounded reactive power"))
                for k in 1:n
                    v = _l3f_control_value(gen["v_target"], k, n, "v_target")
                    isfinite(v) && v > 0 || throw(ArgumentError("v_target must be positive and finite"))
                end
            end
        catch err
            _l3f_error!(findings, "E.L3F.GENERATOR_CONTROL_INVALID", :generator, gid,
                        sprint(showerror, err))
        end
    end
end

"""Restore added control fields in working PU coordinates from the SI input."""
function _l3f_rescale_restricted_controls!(working, physical, bases, options)
    options.per_unit || return working
    for (gid, gen) in get(working, "generator", Dict())
        original = get(get(physical, "generator", Dict()), gid, nothing)
        original isa AbstractDict || continue
        for field in ("aggregate_p_min", "aggregate_p_max", "aggregate_q_min",
                      "aggregate_q_max", "aggregate_s_max")
            haskey(original, field) && (gen[field] = Float64(original[field]) / options.s_base)
        end
        if haskey(original, "v_target")
            n = length(gen["p_min"]); bus = String(gen["bus"])
            gen["v_target"] = [_l3f_control_value(original["v_target"], k, n, "v_target") /
                               bases.v_base[bus] for k in 1:n]
        end
    end
    working
end

"""Stamp aggregate capability, signed fixed-PF, and hard winding-voltage targets."""
function _l3f_stamp_restricted_controls!(model, constraints, net, reference, w,
                                         p_generator, q_generator)
    for (gid_raw, gen) in get(net, "generator", Dict())
        gid = String(gid_raw); n = length(gen["p_min"])
        ps = [p_generator[(gid, k)] for k in 1:n]
        qs = [q_generator[(gid, k)] for k in 1:n]
        for (field, sense, family, vars) in
            (("aggregate_p_min", :ge, :generator_aggregate_p_lower, ps),
             ("aggregate_p_max", :le, :generator_aggregate_p_upper, ps),
             ("aggregate_q_min", :ge, :generator_aggregate_q_lower, qs),
             ("aggregate_q_max", :le, :generator_aggregate_q_upper, qs))
            haskey(gen, field) || continue
            value = Float64(gen[field])
            c = sense == :ge ? @constraint(model, sum(vars) >= value) :
                               @constraint(model, sum(vars) <= value)
            _l3f_register_constraint!(constraints, family, gid, c)
        end
        # This bounds |sum_k(P_k+jQ_k)|. It is intentionally distinct from
        # sum_k|S_k| and from the concurrently retained per-channel nameplates.
        haskey(gen, "aggregate_s_max") && _l3f_register_constraint!(constraints,
            :generator_aggregate_apparent_power, gid,
            @constraint(model, [Float64(gen["aggregate_s_max"]), sum(ps), sum(qs)] in
                               JuMP.SecondOrderCone()))
        if haskey(gen, "fixed_pf")
            pf = Float64(gen["fixed_pf"]); slope = tan(acos(abs(pf)))
            for k in 1:n
                _l3f_register_constraint!(constraints, :generator_fixed_pf, (gid, k),
                    @constraint(model, sign(pf) * qs[k] + slope * ps[k] == 0))
            end
        end
        if haskey(gen, "v_target")
            tm = string.(gen["terminal_map"]); bus = String(gen["bus"])
            D = _l3f_connection_incidence(String(gen["configuration"]), length(tm), n)
            vbar = ComplexF64[reference.voltage[(bus, terminal)] for terminal in tm]
            for k in 1:n
                v2 = _l3f_winding_voltage_squared(D, k, vbar, bus, tm, w)
                target = _l3f_control_value(gen["v_target"], k, n, "v_target")
                _l3f_register_constraint!(constraints, :generator_voltage_target, (gid, k),
                    @constraint(model, v2 == target^2))
            end
        end
    end
    constraints
end

"""Attach original native-IBR identity to lowered generator result rows."""
function _l3f_trace_restricted_controls!(generators, network)
    for (gid, row) in generators
        gen = get(get(network, "generator", Dict()), gid, nothing)
        gen isa AbstractDict || continue
        haskey(gen, "_l3f_original_ibr") || continue
        row["original_component"] = "ibr"
        row["original_id"] = String(gen["_l3f_original_ibr"])
        row["original_topology"] = String(gen["_l3f_original_ibr_topology"])
    end
    generators
end

"""Whether nonlinear replay can preserve every restricted control law."""
function _l3f_restricted_replay_reason(net)
    any(haskey(gen, "v_target") for gen in values(get(net, "generator", Dict()))) &&
        return "unavailable: nonlinear dispatch replay cannot maintain hard generator voltage targets"
    any(!isempty(get(get(net, "transformer", Dict()), subtype, Dict()))
        for subtype in _L3F_LOCAL_BANK_SUBTYPES) &&
        return "unavailable: nonlinear replay has no BMOPFTools component for L3F-local fixed transformer banks"
    nothing
end
