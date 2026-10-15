# Conservative dependency sets. Include both rectangular components even when
# one numerical derivative happens to vanish, and every device parameter mode
# even when its current snapshot coefficient is zero.
function _se_symbolic_jacobians(s, p)
    n=length(s.nodes)
    ns=size(s.voltage_state_jacobian, 2)
    voltage(i) =
        i==0 ? Int[] :
        unique(vcat(findnz(s.voltage_state_rows[:, i])[1], findnz(s.voltage_state_rows[:, n+i])[1]))
    current=[
        unique(
            vcat(findnz(s.current_state_rows[:, i])[1], findnz(s.current_state_rows[:, n+i])[1]),
        ) for i = 1:n
    ]
    hi=Int[]
    hj=Int[]
    for (row, spec) in enumerate(s.measurement_pattern)
        dependencies=Int[]
        if spec isa _SEFunctionalSpec
            for i in spec.nodes
                append!(dependencies, voltage(i))
            end
        elseif spec isa _SEBranchMeasurementSpec
            for i in spec.nodes
                append!(dependencies, voltage(i))
            end
            append!(dependencies, voltage(spec.reference))
        else
            append!(dependencies, voltage(spec.terminal))
            append!(dependencies, voltage(spec.reference))
            if spec isa _SECurrentPowerSpec
                append!(dependencies, (spec.real_state, spec.imag_state))
            elseif !(spec isa _SEKnownPowerSpec) &&
                   spec.kind in (:pinj, :qinj, :ire, :iim, :imag) &&
                   spec.terminal!=0
                append!(dependencies, current[spec.terminal])
            end
        end
        for j in unique(dependencies)
            push!(hi, row)
            push!(hj, j)
        end
    end
    for (k, j) in enumerate(p.prior_indices)
        push!(hi, length(s.measurement_pattern)+k)
        push!(hj, j)
    end
    for d in s.device_pattern, (a, b) in zip(d.positive, d.negative)
        dependencies=union(voltage(a), voltage(b))
        append!(current[a], dependencies)
        b==0 || append!(current[b], dependencies)
    end
    nv=2length(s.free_state_map)
    nu=size(s.injection_incidence, 2)
    for (i, k, _) in zip(findnz(s.injection_incidence)...)
        append!(current[i], (nv+k, nv+nu+k))
    end
    ci=Int[]
    cj=Int[]
    for (k, i) in enumerate(s.constraint_pattern), row in (2k-1, 2k), j in unique(current[i])
        push!(ci, row)
        push!(cj, j)
    end
    offset=length(s.constraint_pattern)
    for r in axes(s.return_current_pattern, 1)
        dependencies=Int[]
        for i in findnz(s.return_current_pattern[r, :])[1]
            append!(dependencies, current[i])
        end
        for row in (2(offset+r)-1, 2(offset+r)), j in unique(dependencies)
            push!(ci, row)
            push!(cj, j)
        end
    end
    sparse(hi, hj, ones(length(hi)), length(s.measurement_pattern)+length(p.prior_indices), ns),
    sparse(ci, cj, ones(length(ci)), _se_constraint_count(s), ns)
end
