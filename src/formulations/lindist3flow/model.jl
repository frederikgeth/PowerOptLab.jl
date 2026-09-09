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
            scale = if edge.family == :line
                1.0
            else
                data = net["transformer"][edge.subtype][edge.id]
                neff = _l3f_transformer_neff(edge.subtype, data)
                edge.reversed ? neff : inv(neff)
            end
            for k in eachindex(edge.parent_map)
                key = (edge.child, edge.child_map[k])
                value = scale * voltage[(edge.parent, edge.parent_map[k])]
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

function _l3f_reference(net, topology, reference)
    result = reference === nothing ? _l3f_source_reference(net, topology) :
             _l3f_explicit_reference(reference)
    _l3f_validate_reference(result, net)
end

function _l3f_report_error(report::L3FApplicabilityReport, code, message)
    findings = copy(report.findings)
    _l3f_error!(findings, code, :network, nothing, message)
    L3FApplicabilityReport(:inapplicable, findings, report.roots,
                           report.islands, report.kron_reduced)
end

function _l3f_add_bounds!(variable, data, index::Int, n::Int,
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

function _l3f_rating(component, fallback, field::String, k::Int)
    data = haskey(component, field) ? component : fallback
    data === nothing || !haskey(data, field) ? nothing : Float64(data[field][k])
end

function _l3f_add_power_circle!(model, constraints, family, key, p, q, radius)
    radius === nothing && return
    _l3f_register_constraint!(constraints, family, key,
        @constraint(model, [Float64(radius), p, q] in JuMP.SecondOrderCone()))
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

function _l3f_cost_coefficient(component, k::Int, kind::String, id::String)
    cost = get(component, "cost", nothing)
    cost === nothing && return 0.0
    cost isa AbstractVector || throw(ArgumentError("$kind '$id' cost must be a per-phase vector"))
    length(cost) >= k || throw(ArgumentError("$kind '$id' cost does not cover phase $k"))
    Float64(cost[k]) / 1000.0
end

function _l3f_build_model(net, topology, reference, report, optimizer, options)
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
        ubar = D * vbar
        for k in 1:nch
            p = @variable(model, base_name=_l3f_name("l3f_pg", gid, k))
            q = @variable(model, base_name=_l3f_name("l3f_qg", gid, k))
            _l3f_add_bounds!(p, gen, k, length(tm), "p_min", "p_max")
            _l3f_add_bounds!(q, gen, k, length(tm), "q_min", "q_max")
            p_generator[(gid, k)] = p; q_generator[(gid, k)] = q
            smax = _l3f_rating(gen, nothing, "s_max", k)
            imax = _l3f_rating(gen, nothing, "i_max", k)
            _l3f_add_power_circle!(model, constraints, :generator_apparent_power,
                (gid, k), p, q, smax)
            _l3f_add_power_circle!(model, constraints, :generator_reference_current,
                (gid, k), p, q,
                imax === nothing ? nothing : imax * abs(ubar[k]))
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
            _l3f_add_bounds!(p, source, k, length(tm), "p_min", "p_max")
            _l3f_add_bounds!(q, source, k, length(tm), "q_min", "q_max")
            p_source[(sid, k)] = p; q_source[(sid, k)] = q
            key = (String(source["bus"]), tm[k])
            _l3f_register_constraint!(constraints, :source_voltage, (sid, tm[k]),
                @constraint(model, w[key] == abs2(reference.voltage[key])))
            smax = _l3f_rating(source, nothing, "s_max", k)
            imax = _l3f_rating(source, nothing, "i_max", k)
            _l3f_add_power_circle!(model, constraints, :source_apparent_power,
                (sid, k), p, q, smax)
            _l3f_add_power_circle!(model, constraints, :source_reference_current,
                (sid, k), p, q,
                imax === nothing ? nothing : imax * abs(reference.voltage[key]))
        end
    end
    variables[:p_source] = p_source; variables[:q_source] = q_source

    for edge in topology
        if edge.family == :transformer
            transformer = net["transformer"][edge.subtype][edge.id]
            neff = _l3f_transformer_neff(edge.subtype, transformer)
            scale = edge.reversed ? neff : inv(neff)
            for phi in eachindex(edge.parent_map)
                _l3f_register_constraint!(constraints, :transformer_voltage_ratio,
                    (edge.subtype, edge.id, phi), @constraint(model,
                        w[(edge.child, edge.child_map[phi])] ==
                        scale^2 * w[(edge.parent, edge.parent_map[phi])]))
                p = p_line[(edge.family, edge.id, phi)]
                q = q_line[(edge.family, edge.id, phi)]
                _l3f_add_power_circle!(model, constraints, :transformer_apparent_power,
                    (edge.subtype, edge.id, phi), p, q,
                    get(transformer, "s_rating", nothing))
                from_vm = abs(reference.voltage[(String(transformer["bus_from"]),
                                                  string.(transformer["terminal_map_from"])[phi])])
                to_vm = abs(reference.voltage[(String(transformer["bus_to"]),
                                                string.(transformer["terminal_map_to"])[phi])])
                _l3f_add_power_circle!(model, constraints, :transformer_from_reference_current,
                    (edge.subtype, edge.id, phi), p, q,
                    haskey(transformer, "i_max_from") ? Float64(transformer["i_max_from"]) * from_vm : nothing)
                _l3f_add_power_circle!(model, constraints, :transformer_to_reference_current,
                    (edge.subtype, edge.id, phi), p, q,
                    haskey(transformer, "i_max_to") ? Float64(transformer["i_max_to"]) * to_vm : nothing)
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
            _l3f_add_power_circle!(model, constraints, :line_apparent_power,
                (edge.id, phi), p, q, smax)
            _l3f_add_power_circle!(model, constraints, :line_reference_current,
                (edge.id, phi), p, q, imax === nothing ? nothing :
                    imax * abs(reference.voltage[(edge.parent, edge.parent_map[phi])]))
        end
    end

    balance_p = Dict(key => JuMP.AffExpr(0.0) for key in keys(w))
    balance_q = Dict(key => JuMP.AffExpr(0.0) for key in keys(w))
    for edge in topology, k in eachindex(edge.parent_map)
        JuMP.add_to_expression!(balance_p[(edge.parent, edge.parent_map[k])], -1.0, p_line[(edge.family, edge.id, k)])
        JuMP.add_to_expression!(balance_q[(edge.parent, edge.parent_map[k])], -1.0, q_line[(edge.family, edge.id, k)])
        JuMP.add_to_expression!(balance_p[(edge.child, edge.child_map[k])], 1.0, p_line[(edge.family, edge.id, k)])
        JuMP.add_to_expression!(balance_q[(edge.child, edge.child_map[k])], 1.0, q_line[(edge.family, edge.id, k)])
    end
    for (lid, load) in get(net, "load", Dict())
        bus = String(load["bus"]); tm = string.(load["terminal_map"])
        p = Float64.(load["p_nom"]); q = Float64.(load["q_nom"])
        D = _l3f_connection_incidence(String(load["configuration"]), length(tm), length(p))
        vbar = ComplexF64[reference.voltage[(bus, t)] for t in tm]
        H = connection_power_map(D, vbar).matrix
        terminal_power = H * complex.(p, q)
        for k in eachindex(tm)
            JuMP.add_to_expression!(balance_p[(bus, tm[k])], -real(terminal_power[k]))
            JuMP.add_to_expression!(balance_q[(bus, tm[k])], -imag(terminal_power[k]))
        end
    end
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
    L3FBuild(model, variables, constraints, reference, report, options,
             net, topology)
end

"""
    build_l3f_opf(net, optimizer=Clarabel.Optimizer;
                  options=L3FOptions(), reference=nothing)

Build the L3F-BMOPF lossless radial LP/SOCP. Explicit-neutral inputs are
Kron-reduced on a copy when enabled. Inapplicable networks raise
[`L3FInapplicableError`](@ref), whose report contains stable diagnostics.
"""
function build_l3f_opf(net, optimizer=Clarabel.Optimizer;
                       options::L3FOptions=L3FOptions(), reference=nothing)
    prepared = _l3f_prepare(net; options, reference)
    is_l3f_applicable(prepared.applicability) ||
        throw(L3FInapplicableError(prepared.applicability))
    ref = try
        _l3f_reference(prepared.network, prepared.topology, reference)
    catch err
        throw(L3FInapplicableError(_l3f_report_error(prepared.applicability,
            "E.L3F.REFERENCE_MISSING", sprint(showerror, err))))
    end
    _l3f_build_model(prepared.network, prepared.topology, ref,
                     prepared.applicability, optimizer, options)
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
        if set_type <: JuMP.MOI.SecondOrderCone
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
    buses = Dict{String,Any}()
    for (busid_raw, bus) in get(build.network, "bus", Dict())
        busid = String(busid_raw); terminals = Dict{String,Any}()
        for terminal in string.(bus["terminal_names"])
            wv = value(build.variables[:w][(busid, terminal)])
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
            "p" => [value(build.variables[:p_line][(edge.family, edge.id, k)]) for k in eachindex(edge.parent_map)],
            "q" => [value(build.variables[:q_line][(edge.family, edge.id, k)]) for k in eachindex(edge.parent_map)],
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
            "pg" => [value(build.variables[:p_generator][(gid, k)]) for k in 1:n],
            "qg" => [value(build.variables[:q_generator][(gid, k)]) for k in 1:n],
        )
    end
    sources = Dict{String,Any}()
    for (sid_raw, source) in get(build.network, "voltage_source", Dict())
        sid = String(sid_raw); n = length(source["terminal_map"])
        sources[sid] = Dict{String,Any}(
            "terminal_map" => string.(source["terminal_map"]),
            "pg" => [value(build.variables[:p_source][(sid, k)]) for k in 1:n],
            "qg" => [value(build.variables[:q_source][(sid, k)]) for k in 1:n],
        )
    end
    buses, lines, transformers, generators, sources
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
    validate_l3f_solution(net, result; nonlinear_optimizer=Ipopt.Optimizer)

Reapply L3F generator dispatch to the Kron-reduced BMOPF snapshot and run the
nonlinear BMOPFTools power flow. The returned dictionary reports voltage error
and explicitly leaves full physical-limit certification unassessed.
"""
function validate_l3f_solution(net, result::L3FResult;
                               nonlinear_optimizer=Ipopt.Optimizer)
    result.solve.optimal || return Dict{String,Any}(
        "status" => "not_run", "reason" => "L3F solve was not optimal")
    working = _l3f_fix_dispatch!(deepcopy(result.network), result)
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
        Dict{String,Any}(
            "status" => solved ? "passed" : "failed",
            "nonlinear_solve_status" => termination,
            "maximum_voltage_magnitude_error" => isempty(errors) ? NaN : maximum(errors),
            "rms_voltage_magnitude_error" => isempty(errors) ? NaN : sqrt(sum(abs2, errors) / length(errors)),
            "physical_limits" => "unassessed",
            "claim" => "nonlinear replay comparison, not a physical-feasibility certificate",
        )
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
                       nonlinear_optimizer=Ipopt.Optimizer, solver_options=())
    build = build_l3f_opf(net, optimizer; options, reference)
    _set_solver_options!(build.model, solver_options)
    JuMP.optimize!(build.model)
    outcome = _solve_outcome(build.model)
    status = SolveStatus(outcome)
    buses, lines, transformers, generators, sources = _l3f_extract(build, outcome)
    result = L3FResult(buses, lines, transformers, generators, sources,
        outcome.optimal ? JuMP.objective_value(build.model) : NaN,
        Dict{String,Any}(
            "name" => "L3F-BMOPF", "version" => "0.1-prototype",
            "problem_class" => String(l3f_model_class(build)),
            "series_losses" => "omitted",
            "voltage_angles" => "fixed reference coefficients; not decision variables",
        ), build.reference, build.applicability,
        Dict{String,Any}("status" => "not_requested"),
        build.network, status)
    if options.validate_nonlinear
        validation = validate_l3f_solution(build.network, result; nonlinear_optimizer)
        result = L3FResult(result.buses, result.lines, result.transformers, result.generators,
            result.sources, result.objective, result.formulation,
            result.reference, result.applicability, validation,
            result.network, result.solve)
    end
    result
end
