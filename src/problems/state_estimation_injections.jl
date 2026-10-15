# Explicit unknown branch currents extend the voltage formulation when its
# elimination assumptions do not apply. Positive current injects into the
# network at the original customer terminal and leaves at its return.
struct _SECurrentPowerSpec
    kind::Symbol
    terminal::Int
    reference::Int
    real_state::Int
    imag_state::Int
end

struct _SEKnownPowerSpec
    kind::Symbol
    terminal::Int
    reference::Int
    branches::Vector{Tuple{Symbol,Int}}
end

function _se_requires_current_states(connections, groups, measurements, devices, index, fixed)
    isempty(groups) || return true
    isempty(devices) || isempty(connections) || return true
    returns = Dict{Int,Int}()
    for c in connections
        i = index[c.positive]
        j = c.negative === nothing ? 0 : index[c.negative]
        haskey(returns, i) && returns[i] != j && return true
        returns[i] = j
    end
    any(j -> j != 0 && haskey(returns, j), values(returns))
end

function _se_current_measurements(specs, measurements, neutral, map, devices, index, nv, nu)
    known=Dict{Tuple{TerminalID,Int},Vector{Tuple{Symbol,Int}}}()
    cursor=0
    for d in devices
        kind=d.model isa ConstantPowerDevice ? :constant_power :
             d.model isa ConstantCurrentDevice ?
             (d.model.angle_reference===:earth ? :constant_current : :constant_magnitude_current) :
             :zip
        for c in d.model.connections
            cursor+=1
            j=c.negative===nothing ? 0 : index[c.negative]
            push!(get!(known, (c.positive, j), Tuple{Symbol,Int}[]), (kind, cursor))
        end
    end
    result=copy(specs)
    for (row, m) in enumerate(measurements)
        m isa Measurement && m.kind in (:pinj, :qinj) || continue
        node=(m.bus, m.terminal)
        spec=specs[row]
        if haskey(map, node)
            k=map[node]
            result[row]=_SECurrentPowerSpec(m.kind, spec.terminal, spec.reference, nv+k, nv+nu+k)
        else
            result[row]=_SEKnownPowerSpec(
                m.kind,
                spec.terminal,
                spec.reference,
                known[(node, spec.reference)],
            )
        end
    end
    result
end

function _se_current_power(spec, vr, vi, x)
    a=vr[spec.terminal]-(spec.reference==0 ? 0.0 : vr[spec.reference])
    b=vi[spec.terminal]-(spec.reference==0 ? 0.0 : vi[spec.reference])
    c=x[spec.real_state]
    d=x[spec.imag_state]
    spec.kind===:pinj ? (a*c+b*d, (c, d, a, b)) : (b*c-a*d, (-d, c, b, -a))
end

function _se_known_power(spec, vr, vi, p)
    a=vr[spec.terminal]-(spec.reference==0 ? 0.0 : vr[spec.reference])
    b=vi[spec.terminal]-(spec.reference==0 ? 0.0 : vi[spec.reference])
    c=d=ca=cb=da=db=0.0
    for (kind, k) in spec.branches
        if kind in (:constant_power, :zip)
            r, i, ra, rb, ia, ib=_power_current_and_derivative(p.device_powers[k], a, b, p)
            c-=r
            d-=i
            ca-=ra
            cb-=rb
            da-=ia
            db-=ib
        end
        if kind===:constant_magnitude_current
            r, i, ra, rb, ia, ib=_se_voltage_following_current(p.device_currents[k], a, b)
            c-=r
            d-=i
            ca-=ra
            cb-=rb
            da-=ia
            db-=ib
        elseif kind in (:constant_current, :zip)
            c-=real(p.device_currents[k])
            d-=imag(p.device_currents[k])
        end
        if kind===:zip
            g=real(p.device_admittances[k])
            h=imag(p.device_admittances[k])
            c-=g*a-h*b
            d-=h*a+g*b
            ca-=g
            cb+=h
            da-=h
            db-=g
        end
    end
    spec.kind===:pinj ? (a*c+b*d, (c+a*ca+b*da, d+a*cb+b*db)) :
    (b*c-a*d, (-d+b*ca-a*da, c+b*cb-a*db))
end

function _se_power_state_jacobian(spec, s, p, x, vr, vi)
    n=length(s.nodes)
    ja=s.voltage_state_rows[:, spec.terminal]
    jb=s.voltage_state_rows[:, n+spec.terminal]
    if spec.reference!=0
        ja-=s.voltage_state_rows[:, spec.reference]
        jb-=s.voltage_state_rows[:, n+spec.reference]
    end
    if spec isa _SECurrentPowerSpec
        value, g=_se_current_power(spec, vr, vi, x)
        row=g[1]*ja+g[2]*jb
        row[spec.real_state]+=g[3]
        row[spec.imag_state]+=g[4]
    else
        value, g=_se_known_power(spec, vr, vi, p)
        row=g[1]*ja+g[2]*jb
    end
    value, row
end

function _se_initial_currents(s, p, x)
    isempty(s.injection_state_map) && return x
    nv=2length(s.free_state_map)
    nu=size(s.injection_incidence, 2)
    length(x)==nv && (x=vcat(x, zeros(2nu)))
    vr, vi, _, _=_se_parts(s, p, x)
    for (node, k) in s.injection_state_map
        c=only(filter(c->c.positive==node, s.injection_connections))
        i=s.node_index[node]
        j=c.negative===nothing ? 0 : s.node_index[c.negative]
        u=complex(vr[i], vi[i])-(j==0 ? 0im : complex(vr[j], vi[j]))
        power=0im
        for kind in (:pinj, :qinj)
            rows=findall(
                spec->spec isa _SECurrentPowerSpec && spec.real_state==nv+k && spec.kind===kind,
                s.measurement_pattern,
            )
            isempty(rows) && continue
            weights=1 ./ p.covariance_values[rows] .^ 2
            value=sum(weights .* p.measurement_values[rows])/sum(weights)
            power+=kind===:pinj ? value : im*value
        end
        current=abs(u)>0 ? conj(power/u) : 0im
        x[nv+k]=real(current)
        x[nv+nu+k]=imag(current)
    end
    x
end

"""
    injection_state_indices(structure)

Map original unknown customer terminals to `(real, imag)` SI current-state
indices. Empty for the eliminated voltage formulation. Voltage-state indices
remain in `free_state_map`; current entries follow all voltage entries.
"""
function injection_state_indices(s::SEStructure)
    nv=2length(s.free_state_map)
    nu=size(s.injection_incidence, 2)
    Dict(node=>(real = nv+k, imag = nv+nu+k) for (node, k) in s.injection_state_map)
end

# Native WLS keeps explicit customer currents. Eliminate their local tangent
# directions only for the VOLTAGE observability diagnostic, including the
# nuisance directions of customers aliased to one electrical node.
function _se_native_current_observability(net, measurements, neutral, bus)
    s=compile_state_estimator(net, measurements; neutral, formulation = :current)
    covered=Set(s.constraint_pattern)
    fixed=Set(s.node_index[node] for node in keys(s.reference_map))
    extra=[node for (i, node) in enumerate(s.nodes) if !(i in covered) && !(i in fixed)]
    s=compile_state_estimator(
        net,
        measurements;
        neutral,
        formulation = :current,
        zero_injection = extra,
    )
    p=SEParameters(s, measurements)
    x=copy(s.nominal_state)
    nv=2length(s.free_state_map)
    nf=nv÷2
    for (node, i) in s.node_index
        i==0 && continue
        canonical=s.nodes[i]
        haskey(s.free_state_map, canonical) || continue
        k=s.free_state_map[canonical]
        v=get(get(bus, node[1], Dict()), node[2], Dict())
        x[k]=get(v, "vr", x[k])
        x[nf+k]=get(v, "vi", x[nf+k])
    end
    # Conditional weighted fit is linear in the unknown currents at fixed V.
    C=collect(Float64, constraint_jacobian(s, p, x))
    A=C[:, (nv+1):end]
    c=evaluate_state_estimator(s, p, x).constraints
    step=_se_least_squares(A, -c)
    x[(nv+1):end]+=step
    _, Z=_se_rank_nullspace(A)
    H=collect(Float64, residual_jacobian(s, p, x))
    Hj=H[:, (nv+1):end]
    r=evaluate_state_estimator(s, p, x).residual
    x[(nv+1):end]-=Z*_se_least_squares(Hj*Z, r)
    H=collect(Float64, residual_jacobian(s, p, x))
    Hj=H[:, (nv+1):end]
    _, L=_se_rank_nullspace(transpose(A))
    T=-pinv(A)*C[:, 1:nv]
    effective=H[:, 1:nv]+Hj*T
    W=Hj*Z
    rankW, F=_se_svd_rank(W; reference_scale = norm(Hj))
    rankW==0 || (effective-=F.U[:, 1:rankW]*(F.U[:, 1:rankW]'*effective))
    _se_observability_metrics(
        effective,
        L'*C[:, 1:nv];
        nuisance_dimension = rankW,
        constraint_reference_scale = norm(C[:, 1:nv]),
    )
end

# Current at unit connection angle is a complex power-factor coefficient.
# Its magnitude remains constant while the voltage angle rotates it.
function _se_voltage_following_current(coefficient, a, b)
    iszero(coefficient) && return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    magnitude=hypot(a, b)
    magnitude>0 || throw(
        DomainError(
            magnitude,
            "voltage-following current angle is undefined at zero connection voltage",
        ),
    )
    g=real(coefficient)
    h=imag(coefficient)
    r=(g*a-h*b)/magnitude
    i=(h*a+g*b)/magnitude
    (
        r,
        i,
        g/magnitude-r*a/magnitude^2,
        -h/magnitude-r*b/magnitude^2,
        h/magnitude-i*a/magnitude^2,
        g/magnitude-i*b/magnitude^2,
    )
end

function _se_known_connections(devices, index)
    Set(
        (c.positive, c.negative===nothing ? 0 : index[c.negative]) for d in devices for
        c in d.model.connections
    )
end

function _se_check_explicit_device_unknowns(connections, devices, index, fixed)
    connections===nothing && return nothing
    free(node) = node===nothing ? 0 : index[node]
    known=Set(
        free(node) for d in devices for c in d.model.connections for
        node in (c.positive, c.negative) if free(node)!=0 && !haskey(fixed, free(node))
    )
    any(c->any(node->free(node) in known, (c.positive, c.negative)), connections) && throw(
        ArgumentError(
            "explicit unknown currents sharing an exact-device conductor require formulation=:current or :auto",
        ),
    )
    nothing
end
