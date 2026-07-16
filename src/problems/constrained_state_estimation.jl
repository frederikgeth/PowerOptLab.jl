# Compiled four-wire voltage-state evaluation for the constrained-NLLS estimator.
#
# This is deliberately separate from the legacy JuMP/Ipopt WLS prototype in
# `state_estimation.jl`.  It owns only immutable topology/equation structure and
# mutable numerical parameters; a dense composite-step solver is layered on top
# in the next development phase.

const TerminalID = Tuple{String,String}

"""
    ExactInjectionSpecification

Marker hierarchy for information which is genuinely exact.  Meter readings and
forecasts are intentionally absent: they belong in the stochastic residual
model, not in this hierarchy.
"""
abstract type ExactInjectionSpecification end
struct NoExactInjection <: ExactInjectionSpecification end
struct ExactZeroInjection <: ExactInjectionSpecification end
struct ExactDeviceEquation{M} <: ExactInjectionSpecification
    model::M
end

"""One oriented device branch: current is positive from `positive` to `negative`."""
struct TerminalConnection
    positive::TerminalID
    negative::Union{TerminalID,Nothing}
end

abstract type ExactDeviceModel end

"""Exact constant-power branches; positive power denotes consumption."""
struct ConstantPowerDevice <: ExactDeviceModel
    connections::Vector{TerminalConnection}
    powers::Vector{ComplexF64}
    function ConstantPowerDevice(connections::AbstractVector{TerminalConnection},
                                 powers::AbstractVector{<:Complex})
        length(connections) == length(powers) ||
            throw(ArgumentError("constant-power device needs one power per connection"))
        isempty(connections) && throw(ArgumentError("constant-power device needs at least one connection"))
        new(collect(connections), ComplexF64.(powers))
    end
end

"""Exact constant-current branches; positive current flows into the device."""
struct ConstantCurrentDevice <: ExactDeviceModel
    connections::Vector{TerminalConnection}
    currents::Vector{ComplexF64}
    function ConstantCurrentDevice(connections::AbstractVector{TerminalConnection},
                                   currents::AbstractVector{<:Complex})
        length(connections) == length(currents) ||
            throw(ArgumentError("constant-current device needs one current per connection"))
        isempty(connections) && throw(ArgumentError("constant-current device needs at least one connection"))
        new(collect(connections), ComplexF64.(currents))
    end
end

"""Exact ZIP branches: `conj(S/V) + I + YV`, all in SI branch quantities."""
struct ZIPDevice <: ExactDeviceModel
    connections::Vector{TerminalConnection}
    powers::Vector{ComplexF64}
    currents::Vector{ComplexF64}
    admittances::Vector{ComplexF64}
    function ZIPDevice(connections::AbstractVector{TerminalConnection},
                       powers::AbstractVector{<:Complex}, currents::AbstractVector{<:Complex},
                       admittances::AbstractVector{<:Complex})
        n = length(connections)
        length(powers) == n && length(currents) == n && length(admittances) == n ||
            throw(ArgumentError("ZIP device needs one S, I, and Y value per connection"))
        n > 0 || throw(ArgumentError("ZIP device needs at least one connection"))
        new(collect(connections), ComplexF64.(powers), ComplexF64.(currents), ComplexF64.(admittances))
    end
end

struct _SEMeasurementSpec{Ti<:Integer}
    kind::Symbol
    terminal::Ti
    reference::Ti                  # 0 denotes the earth reference
end

struct _SEExactDeviceSpec{Ti<:Integer}
    kind::Symbol
    positive::Vector{Ti}
    negative::Vector{Ti}           # 0 denotes earth
    parameter_range::UnitRange{Int}
end

"""
    SEStructure{Ti}

Immutable compiled structure for the voltage-only, four-wire state-estimation
formulation.  It imports BMOPFTools' passive, conductor-to-earth Ybus exactly
once.  `free_state_map` maps a free conductor to its position in the rectangular
state `[real(V_free); imag(V_free)]`; ideal-source conductors are held in the
parameter vector instead.  Closed ideal switches are already represented by the
node aliases in `ybus_passive`.

The evaluator supports terminal voltage components/magnitudes and terminal
active/reactive injection measurements, exact zero-injection equations, and
connection-aware exact constant-power, constant-current, and ZIP devices.
Branch telemetry and sparse linear algebra belong to subsequent phases.
"""
struct SEStructure{Ti<:Integer}
    nodes::Vector{TerminalID}
    node_index::Dict{TerminalID,Ti}
    passive_pattern::SparseMatrixCSC{ComplexF64,Ti}
    measurement_pattern::Vector{_SEMeasurementSpec{Ti}}
    constraint_pattern::Vector{Ti}
    device_pattern::Vector{_SEExactDeviceSpec{Ti}}
    free_state_map::Dict{TerminalID,Ti}
    reference_map::Dict{TerminalID,ComplexF64}
    voltage_state_jacobian::SparseMatrixCSC{Float64,Ti}
    current_state_jacobian::SparseMatrixCSC{Float64,Ti}
    fixed_voltage_state::Vector{Float64}
    fixed_current_state::Vector{Float64}
end

"""
    SEParameters

Mutable numerical data paired with an [`SEStructure`](@ref).  Updating the
measurement values, standard deviations, or fixed source phasors does not alter
terminal ordering or symbolic sparsity.  Standard deviations are the diagonal
whitening factors for this first implementation.
"""
mutable struct SEParameters
    fixed_voltages::Vector{ComplexF64}
    measurement_values::Vector{Float64}
    covariance_values::Vector{Float64}
    magnitude_epsilon::Float64
    device_powers::Vector{ComplexF64}
    device_currents::Vector{ComplexF64}
    device_admittances::Vector{ComplexF64}
    voltage_min_model::Float64
    regularization_voltage::Float64
    continuation_alpha::Float64
end

"""Values of the compiled residual/constraint model at one voltage state."""
struct SEEvaluation
    voltage::Vector{ComplexF64}
    current::Vector{ComplexF64}
    device_current::Vector{ComplexF64}
    predicted::Vector{Float64}
    residual::Vector{Float64}
    constraints::Vector{Float64}
end

"""
    ConstrainedStateEstimationResult

Result from the dense composite-step reference solver.  `status` distinguishes a
numerically converged but underobserved estimate from failure to establish the
exact equations.  `history` records the scaled trust-region radius, merit value,
measurement objective, and exact-constraint norm at accepted iterates.
"""
struct ConstrainedStateEstimationResult
    status::Symbol
    state::Vector{Float64}
    evaluation::SEEvaluation
    iterations::Int
    constraint_rank::Int
    tangent_dimension::Int
    observable_dimension::Int
    history::Vector{NamedTuple}
end

"""Result and per-stage diagnostics from constant-power continuation."""
struct ContinuationStateEstimationResult
    status::Symbol
    result::ConstrainedStateEstimationResult
    alphas::Vector{Float64}
    stages::Vector{ConstrainedStateEstimationResult}
end

_se_node_index(s::SEStructure, node::TerminalID) = get(s.node_index, node, 0)

function _source_phasors(net::Dict{String,Any}, node_index)
    fixed = Dict{Int,ComplexF64}()
    for (_, source) in get(net, "voltage_source", Dict())
        bus = String(get(source, "bus", ""))
        terminals = String.(get(source, "terminal_map", String[]))
        magnitudes = Float64.(get(source, "v_magnitude", Float64[]))
        angles = Float64.(get(source, "v_angle", zeros(length(terminals))))
        length(magnitudes) == length(terminals) ||
            throw(ArgumentError("voltage source at bus '$bus' must provide one v_magnitude per terminal"))
        for k in eachindex(terminals)
            i = get(node_index, (bus, terminals[k]), 0)
            i == 0 && throw(ArgumentError("voltage source terminal ($bus, $(terminals[k])) is earth-referenced or absent from Ybus"))
            v = magnitudes[k] * cis(k <= length(angles) ? angles[k] : 0.0)
            if haskey(fixed, i) && !isapprox(fixed[i], v; rtol=1e-10, atol=1e-10)
                throw(ArgumentError("conflicting ideal source phasors at aliased Ybus node ($bus, $(terminals[k]))"))
            end
            fixed[i] = v
        end
    end
    fixed
end

function _compile_measurements(measurements, node_index, neutral)
    specs = _SEMeasurementSpec{Int}[]
    for m in measurements
        m isa Measurement || throw(ArgumentError("measurements must contain `Measurement` values"))
        i = get(node_index, (m.bus, m.terminal), 0)
        i == 0 && throw(ArgumentError("measurement terminal ($(m.bus), $(m.terminal)) is not a free Ybus conductor"))
        ref = _resolve_ref(m, neutral)
        j = ref === nothing ? 0 : get(node_index, (m.bus, ref), 0)
        ref !== nothing && !haskey(node_index, (m.bus, ref)) &&
            throw(ArgumentError("measurement reference terminal ($(m.bus), $ref) is absent from Ybus"))
        m.kind in (:vr, :vi, :vmag, :pinj, :qinj) ||
            throw(ArgumentError("compiled state estimator does not support measurement kind :$(m.kind)"))
        push!(specs, _SEMeasurementSpec(m.kind, i, j))
    end
    specs
end

function _compile_exact_devices(exact_devices, node_index)
    specs = _SEExactDeviceSpec{Int}[]
    constrained = Int[]
    cursor = 1
    for wrapped in exact_devices
        wrapped isa ExactDeviceEquation ||
            throw(ArgumentError("exact_devices must contain ExactDeviceEquation(model) values"))
        model = wrapped.model
        model isa ExactDeviceModel ||
            throw(ArgumentError("unsupported exact device model $(typeof(model))"))
        kind = model isa ConstantPowerDevice ? :constant_power :
               model isa ConstantCurrentDevice ? :constant_current : :zip
        positive = Int[]; negative = Int[]
        for connection in model.connections
            i = get(node_index, connection.positive, 0)
            i == 0 && throw(ArgumentError("exact device positive terminal $(connection.positive) is earth-referenced or absent from Ybus"))
            j = connection.negative === nothing ? 0 : get(node_index, connection.negative, 0)
            connection.negative !== nothing && j == 0 &&
                throw(ArgumentError("exact device negative terminal $(connection.negative) is earth-referenced or absent from Ybus"))
            i == j && throw(ArgumentError("exact device connection $(connection.positive) has identical endpoints"))
            push!(positive, i); push!(negative, j); push!(constrained, i)
            j != 0 && push!(constrained, j)
        end
        rng = cursor:(cursor + length(positive) - 1)
        cursor += length(positive)
        push!(specs, _SEExactDeviceSpec(kind, positive, negative, rng))
    end
    specs, constrained
end

function _flatten_device_parameters(exact_devices)
    powers = ComplexF64[]; currents = ComplexF64[]; admittances = ComplexF64[]
    for wrapped in exact_devices
        model = wrapped.model
        if model isa ConstantPowerDevice
            append!(powers, model.powers)
            append!(currents, zeros(ComplexF64, length(model.connections)))
            append!(admittances, zeros(ComplexF64, length(model.connections)))
        elseif model isa ConstantCurrentDevice
            append!(powers, zeros(ComplexF64, length(model.connections)))
            append!(currents, model.currents)
            append!(admittances, zeros(ComplexF64, length(model.connections)))
        elseif model isa ZIPDevice
            append!(powers, model.powers)
            append!(currents, model.currents)
            append!(admittances, model.admittances)
        else
            throw(ArgumentError("unsupported exact device model $(typeof(model))"))
        end
    end
    powers, currents, admittances
end

"""
    compile_state_estimator(net, measurements=Measurement[];
                            neutral="n", zero_injection=String[], exact_devices=[]) -> SEStructure

Compile the immutable voltage-state evaluator.  The network is represented by
BMOPFTools' passive `I = YV` relation in SI units.  Every source terminal with a
specified phasor is eliminated from the state; ungrounded neutrals and floating
conductors remain explicit states, so gauge/reference deficiencies are visible
to the later rank diagnostics rather than silently grounded.
"""
function compile_state_estimator(net::Dict{String,Any}, measurements::AbstractVector=Measurement[];
                                 neutral::Union{String,Nothing}="n",
                                 zero_injection=String[], exact_devices=Any[])
    ybus = ybus_passive(net)
    nodes = TerminalID.(ybus.nodes)
    node_index = Dict{TerminalID,Int}(TerminalID(k) => Int(v) for (k, v) in ybus.index)
    fixed = _source_phasors(net, node_index)
    n = length(nodes)
    free_indices = [i for i in 1:n if !haskey(fixed, i)]
    free_state_map = Dict{TerminalID,Int}(nodes[i] => k for (k, i) in enumerate(free_indices))

    # Map a rectangular free-voltage state into all Ybus node voltages.
    nf = length(free_indices)
    I = Int[]; J = Int[]; V = Float64[]
    for (k, i) in enumerate(free_indices)
        push!(I, i);     push!(J, k);      push!(V, 1.0)
        push!(I, n + i); push!(J, nf + k); push!(V, 1.0)
    end
    E = sparse(I, J, V, 2n, 2nf)
    Yr = sparse(real(ybus.Y)); Yi = sparse(imag(ybus.Y))
    K = [Yr -Yi; Yi Yr]
    M = sparse(K * E)

    fixed_voltage = zeros(Float64, 2n)
    for (i, v) in fixed
        fixed_voltage[i] = real(v)
        fixed_voltage[n + i] = imag(v)
    end
    fixed_current = Vector{Float64}(K * fixed_voltage)

    specs = _compile_measurements(measurements, node_index, neutral)
    zi = _zero_injection_set(net, zero_injection, neutral)
    constraint_nodes = Int[]
    for node in zi
        i = get(node_index, node, 0)
        i == 0 && throw(ArgumentError("zero-injection terminal $node is earth-referenced or absent from Ybus"))
        haskey(fixed, i) && throw(ArgumentError("zero-injection terminal $node is fixed by an ideal source"))
        push!(constraint_nodes, i)
    end
    unique!(constraint_nodes)
    devices, device_nodes = _compile_exact_devices(exact_devices, node_index)
    append!(constraint_nodes, device_nodes)
    unique!(constraint_nodes)
    sort!(constraint_nodes)

    SEStructure(nodes, node_index, sparse(ybus.Y), specs, constraint_nodes, devices,
                free_state_map, Dict(nodes[i] => v for (i, v) in fixed),
                E, M, fixed_voltage, fixed_current)
end

function SEParameters(s::SEStructure, measurements::AbstractVector=Measurement[];
                      exact_devices=Any[], magnitude_epsilon::Real=0.0,
                      voltage_min_model::Real=1e-3, regularization_voltage::Real=1e-3,
                      continuation_alpha::Real=1.0)
    length(measurements) == length(s.measurement_pattern) ||
        throw(ArgumentError("SEParameters needs the same number of measurements used for compilation"))
    epsmag = Float64(magnitude_epsilon)
    epsmag >= 0 && isfinite(epsmag) ||
        throw(ArgumentError("magnitude_epsilon must be finite and non-negative"))
    vmin = Float64(voltage_min_model); vreg = Float64(regularization_voltage)
    α = Float64(continuation_alpha)
    vmin > 0 && isfinite(vmin) || throw(ArgumentError("voltage_min_model must be finite and > 0"))
    vreg > 0 && isfinite(vreg) || throw(ArgumentError("regularization_voltage must be finite and > 0"))
    0 <= α <= 1 && isfinite(α) || throw(ArgumentError("continuation_alpha must lie in [0, 1]"))
    length(exact_devices) == length(s.device_pattern) ||
        throw(ArgumentError("SEParameters needs the exact_devices used for compilation"))
    # Store one phasor per fixed Ybus node in the deterministic node order.
    fixed_by_node = fill(ComplexF64(NaN, NaN), length(s.nodes))
    for (node, v) in s.reference_map
        fixed_by_node[s.node_index[node]] = v
    end
    values = Float64[m.value for m in measurements]
    sigmas = Float64[m.sigma for m in measurements]
    isempty(measurements) && (values = Float64[]; sigmas = Float64[])
    powers, currents, admittances = _flatten_device_parameters(exact_devices)
    SEParameters(fixed_by_node, values, sigmas, epsmag, powers, currents, admittances,
                 vmin, vreg, α)
end

function _validate_parameters(s::SEStructure, p::SEParameters)
    length(p.fixed_voltages) == length(s.nodes) ||
        throw(ArgumentError("fixed_voltages must contain one entry per compiled Ybus node"))
    length(p.measurement_values) == length(s.measurement_pattern) ||
        throw(ArgumentError("measurement_values length does not match compiled measurement pattern"))
    length(p.covariance_values) == length(s.measurement_pattern) ||
        throw(ArgumentError("covariance_values length does not match compiled measurement pattern"))
    all(x -> isfinite(x) && x > 0, p.covariance_values) ||
        throw(ArgumentError("all measurement standard deviations must be finite and > 0"))
    p.magnitude_epsilon >= 0 && isfinite(p.magnitude_epsilon) ||
        throw(ArgumentError("magnitude_epsilon must be finite and non-negative"))
    ndevice = sum((length(d.positive) for d in s.device_pattern); init=0)
    length(p.device_powers) == ndevice || throw(ArgumentError("device_powers length does not match compiled device pattern"))
    length(p.device_currents) == ndevice || throw(ArgumentError("device_currents length does not match compiled device pattern"))
    length(p.device_admittances) == ndevice || throw(ArgumentError("device_admittances length does not match compiled device pattern"))
    p.voltage_min_model > 0 && isfinite(p.voltage_min_model) ||
        throw(ArgumentError("voltage_min_model must be finite and > 0"))
    p.regularization_voltage > 0 && isfinite(p.regularization_voltage) ||
        throw(ArgumentError("regularization_voltage must be finite and > 0"))
    0 <= p.continuation_alpha <= 1 && isfinite(p.continuation_alpha) ||
        throw(ArgumentError("continuation_alpha must lie in [0, 1]"))
    return nothing
end

function _se_parts(s::SEStructure, p::SEParameters, x::AbstractVector{<:Real})
    length(x) == size(s.voltage_state_jacobian, 2) ||
        throw(DimensionMismatch("state has length $(length(x)); expected $(size(s.voltage_state_jacobian, 2))"))
    _validate_parameters(s, p)
    # Fixed source values are deliberately read from parameters: time-series
    # source updates retain the compilation and all Jacobian sparsity.
    b = copy(s.fixed_voltage_state)
    for node in keys(s.reference_map)
        i = s.node_index[node]
        v = p.fixed_voltages[i]
        (isfinite(real(v)) && isfinite(imag(v))) ||
            throw(ArgumentError("fixed source phasor at node $node must be finite"))
        b[i] = real(v); b[length(s.nodes) + i] = imag(v)
    end
    u = b + s.voltage_state_jacobian * x
    q = Vector{Float64}(s.fixed_current_state + s.current_state_jacobian * x)
    n = length(s.nodes)
    u[1:n], u[n+1:end], q[1:n], q[n+1:end]
end

function _measurement_value(spec, vr, vi, ir, ii, epsmag)
    i, j = spec.terminal, spec.reference
    dvr = vr[i] - (j == 0 ? 0.0 : vr[j])
    dvi = vi[i] - (j == 0 ? 0.0 : vi[j])
    spec.kind === :vr   && return dvr
    spec.kind === :vi   && return dvi
    spec.kind === :vmag && return sqrt(dvr^2 + dvi^2 + epsmag^2)
    spec.kind === :pinj && return dvr * ir[i] + dvi * ii[i]
    spec.kind === :qinj && return dvi * ir[i] - dvr * ii[i]
    error("unsupported compiled measurement kind $(spec.kind)")
end

function _power_current_and_derivative(S, vr, vi, p::SEParameters)
    d = vr^2 + vi^2
    if p.continuation_alpha > 0 && sqrt(d) < p.voltage_min_model
        throw(DomainError(sqrt(d), "constant-power device voltage is below voltage_min_model"))
    end
    P, Q = real(S), imag(S)
    nr = P * vr + Q * vi
    ni = P * vi - Q * vr
    dreg = d + p.regularization_voltage^2
    α = p.continuation_alpha
    denom = α == 1 ? d : dreg
    jr = nr / denom; ji = ni / denom
    djr_dvr = (P * denom - 2vr * nr) / denom^2
    djr_dvi = (Q * denom - 2vi * nr) / denom^2
    dji_dvr = (-Q * denom - 2vr * ni) / denom^2
    dji_dvi = (P * denom - 2vi * ni) / denom^2
    if 0 < α < 1
        # Continuation blends a smooth internal law with the physical law.  A
        # caller must finish at α=1 before accepting the final estimate.
        er = nr / d; ei = ni / d
        der_vr = (P * d - 2vr * nr) / d^2
        der_vi = (Q * d - 2vi * nr) / d^2
        dei_vr = (-Q * d - 2vr * ni) / d^2
        dei_vi = (P * d - 2vi * ni) / d^2
        jr = (1 - α) * jr + α * er; ji = (1 - α) * ji + α * ei
        djr_dvr = (1 - α) * djr_dvr + α * der_vr
        djr_dvi = (1 - α) * djr_dvi + α * der_vi
        dji_dvr = (1 - α) * dji_dvr + α * dei_vr
        dji_dvi = (1 - α) * dji_dvi + α * dei_vi
    end
    jr, ji, djr_dvr, djr_dvi, dji_dvr, dji_dvi
end

function _device_parts(s::SEStructure, p::SEParameters, vr, vi)
    n = length(s.nodes); ns = size(s.voltage_state_jacobian, 2)
    dr = zeros(Float64, n); di = zeros(Float64, n)
    Jdr = zeros(Float64, n, ns); Jdi = zeros(Float64, n, ns)
    E = s.voltage_state_jacobian
    for device in s.device_pattern
        for (branch_index, k) in enumerate(device.parameter_range)
            i, j = device.positive[branch_index], device.negative[branch_index]
            dvr = vr[i] - (j == 0 ? 0.0 : vr[j])
            dvi = vi[i] - (j == 0 ? 0.0 : vi[j])
            jdvr = Vector(E[i, :]); jdvi = Vector(E[n + i, :])
            if j != 0
                jdvr .-= Vector(E[j, :]); jdvi .-= Vector(E[n + j, :])
            end
            jr = ji = djr_dvr = djr_dvi = dji_dvr = dji_dvi = 0.0
            if device.kind === :constant_power || device.kind === :zip
                jr, ji, djr_dvr, djr_dvi, dji_dvr, dji_dvi =
                    _power_current_and_derivative(p.device_powers[k], dvr, dvi, p)
            end
            if device.kind === :constant_current || device.kind === :zip
                I = p.device_currents[k]
                jr += real(I); ji += imag(I)
            end
            if device.kind === :zip
                Y = p.device_admittances[k]
                G, B = real(Y), imag(Y)
                jr += G * dvr - B * dvi; ji += B * dvr + G * dvi
                djr_dvr += G; djr_dvi -= B; dji_dvr += B; dji_dvi += G
            end
            gjr = djr_dvr .* jdvr .+ djr_dvi .* jdvi
            gji = dji_dvr .* jdvr .+ dji_dvi .* jdvi
            dr[i] += jr; di[i] += ji; Jdr[i, :] .+= gjr; Jdi[i, :] .+= gji
            if j != 0
                dr[j] -= jr; di[j] -= ji; Jdr[j, :] .-= gjr; Jdi[j, :] .-= gji
            end
        end
    end
    dr, di, Jdr, Jdi
end

"""
    evaluate_state_estimator(structure, parameters, x) -> SEEvaluation

Evaluate SI phasors, whitened stochastic residuals `(h(x)-z)/σ`, and exact
zero-injection residuals.  The constraint vector is never whitened or otherwise
softened.
"""
function evaluate_state_estimator(s::SEStructure, p::SEParameters, x::AbstractVector{<:Real})
    vr, vi, ir, ii = _se_parts(s, p, x)
    dir, dii, _, _ = _device_parts(s, p, vr, vi)
    predicted = [_measurement_value(spec, vr, vi, ir, ii, p.magnitude_epsilon)
                 for spec in s.measurement_pattern]
    residual = (predicted .- p.measurement_values) ./ p.covariance_values
    constraints = Vector{Float64}(undef, 2length(s.constraint_pattern))
    for (k, i) in enumerate(s.constraint_pattern)
        constraints[2k - 1] = ir[i] + dir[i]
        constraints[2k] = ii[i] + dii[i]
    end
    SEEvaluation(ComplexF64.(vr .+ im .* vi), ComplexF64.(ir .+ im .* ii),
                 ComplexF64.(dir .+ im .* dii),
                 predicted, residual, constraints)
end

"""Analytic Jacobian of the whitened stochastic residual vector."""
function residual_jacobian(s::SEStructure, p::SEParameters, x::AbstractVector{<:Real})
    vr, vi, ir, ii = _se_parts(s, p, x)
    n = length(s.nodes)
    E = s.voltage_state_jacobian
    M = s.current_state_jacobian
    nstate = size(E, 2)
    J = zeros(Float64, length(s.measurement_pattern), nstate)
    for (row, spec) in enumerate(s.measurement_pattern)
        i, j = spec.terminal, spec.reference
        jdvr = Vector(E[i, :]); jdvi = Vector(E[n + i, :])
        if j != 0
            jdvr .-= Vector(E[j, :]); jdvi .-= Vector(E[n + j, :])
        end
        dvr = vr[i] - (j == 0 ? 0.0 : vr[j])
        dvi = vi[i] - (j == 0 ? 0.0 : vi[j])
        if spec.kind === :vr
            J[row, :] .= jdvr
        elseif spec.kind === :vi
            J[row, :] .= jdvi
        elseif spec.kind === :vmag
            mag = sqrt(dvr^2 + dvi^2 + p.magnitude_epsilon^2)
            mag > 0 || throw(DomainError(mag, "voltage-magnitude derivative is undefined at zero; set magnitude_epsilon > 0"))
            J[row, :] .= (dvr .* jdvr .+ dvi .* jdvi) ./ mag
        elseif spec.kind === :pinj
            J[row, :] .= ir[i] .* jdvr .+ ii[i] .* jdvi .+
                         dvr .* Vector(M[i, :]) .+ dvi .* Vector(M[n + i, :])
        elseif spec.kind === :qinj
            J[row, :] .= ii[i] .* jdvr .- ir[i] .* jdvi .+
                         dvi .* Vector(M[i, :]) .- dvr .* Vector(M[n + i, :])
        end
        J[row, :] ./= p.covariance_values[row]
    end
    J
end

"""Analytic Jacobian of exact zero-injection constraints."""
function constraint_jacobian(s::SEStructure, p::SEParameters, x::AbstractVector{<:Real})
    vr, vi, _, _ = _se_parts(s, p, x)
    _, _, Jdr, Jdi = _device_parts(s, p, vr, vi)
    n = length(s.nodes)
    C = zeros(Float64, 2length(s.constraint_pattern), size(s.current_state_jacobian, 2))
    for (k, i) in enumerate(s.constraint_pattern)
        C[2k - 1, :] .= s.current_state_jacobian[i, :] .+ Jdr[i, :]
        C[2k, :] .= s.current_state_jacobian[n + i, :] .+ Jdi[i, :]
    end
    C
end

# ── dense composite-step reference solver ───────────────────────────────────

function _se_rank_nullspace(A::AbstractMatrix{<:Real}; rtol::Real=sqrt(eps(Float64)))
    n = size(A, 2)
    size(A, 1) == 0 && return (0, Matrix{Float64}(I, n, n))
    F = svd(Matrix{Float64}(A); full=true)
    smax = isempty(F.S) ? 0.0 : maximum(F.S)
    tol = max(Float64(rtol) * max(size(A)...) * smax, eps(Float64))
    rank = count(>(tol), F.S)
    rank, Matrix(F.V[:, rank+1:end])
end

_se_merit(e::SEEvaluation, μ) = 0.5 * sum(abs2, e.residual) + μ * norm(e.constraints)

function _se_scaled_normal_step(C, c, scale, radius)
    isempty(c) && return zeros(Float64, size(C, 2))
    Cs = C * Diagonal(1.0 ./ scale)
    y = -(Cs \ c)
    ny = norm(y)
    ny > radius && ny > 0 && (y .*= radius / ny)
    y ./ scale
end

function _se_soc_step(C, defect, scale, radius)
    isempty(defect) && return zeros(Float64, size(C, 2))
    y = (C * Diagonal(1.0 ./ scale)) \ defect
    ny = norm(y)
    ny > radius && ny > 0 && (y .*= radius / ny)
    y ./ scale
end

function _se_history_entry(iteration, radius, merit, e)
    (iteration=iteration, radius=radius, merit=merit,
     measurement_objective=0.5 * sum(abs2, e.residual),
     constraint_norm=norm(e.constraints))
end

"""
    solve_compiled_state_estimator(structure, parameters, x0; kwargs...)
        -> ConstrainedStateEstimationResult

Dense reference implementation of the plan's equality-constrained
Gauss--Newton method.  It uses a scaled Byrd--Omojokun composite step: a normal
least-squares step for exact-equation violation, followed by a null-space
tangential measurement step.  An exact-penalty merit function globalises both
quantities; rejected nonlinear constraint steps receive one trust-region-limited
second-order correction before the radius contracts.

This intentionally transparent solver is for small-system verification.  It
uses dense SVD rank/null-space diagnostics; the compiled evaluator preserves the
sparsity needed by the planned sparse QR/Hachtel implementation.
"""
function solve_compiled_state_estimator(s::SEStructure, p::SEParameters,
                                        x0::AbstractVector{<:Real};
                                        max_iterations::Integer=100,
                                        initial_radius::Real=0.25,
                                        max_radius::Real=4.0,
                                        normal_fraction::Real=0.8,
                                        penalty::Real=10.0,
                                        acceptance_threshold::Real=0.1,
                                        constraint_tolerance::Real=1e-8,
                                        optimality_tolerance::Real=1e-8,
                                        min_radius::Real=1e-10,
                                        rank_rtol::Real=sqrt(eps(Float64)))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    0 < initial_radius <= max_radius || throw(ArgumentError("require 0 < initial_radius ≤ max_radius"))
    0 < normal_fraction < 1 || throw(ArgumentError("normal_fraction must lie in (0, 1)"))
    penalty > 0 || throw(ArgumentError("penalty must be positive"))
    0 < acceptance_threshold < 1 || throw(ArgumentError("acceptance_threshold must lie in (0, 1)"))

    x = Float64.(x0)
    e = evaluate_state_estimator(s, p, x)
    # A uniform nominal-voltage scale gives the rectangular state a physically
    # meaningful trust-region norm even when an initial imaginary component is 0.
    nominal = maximum(abs.(p.fixed_voltages[isfinite.(real.(p.fixed_voltages))]); init=1.0)
    # `scale` is the diagonal D in ||D*s|| ≤ Δ, hence it has inverse-voltage
    # units.  A unit scaled step is one nominal-voltage state increment.
    scale = 1.0 ./ max.(abs.(x), max(nominal, 1.0))
    radius = Float64(initial_radius)
    μ = Float64(penalty)
    history = [_se_history_entry(0, radius, _se_merit(e, μ), e)]
    rejected = 0

    for iteration in 1:max_iterations
        H = residual_jacobian(s, p, x)
        C = constraint_jacobian(s, p, x)
        rankC, Z = _se_rank_nullspace(C; rtol=rank_rtol)
        c = e.constraints
        r = e.residual
        # Continuation stages and externally supplied warm starts can already
        # satisfy both tests.  Do not manufacture a zero predicted-reduction
        # failure merely because there is no step left to take.
        gred0 = size(Z, 2) == 0 ? 0.0 : norm(Z' * (H' * r), Inf)
        if norm(c) <= constraint_tolerance && gred0 <= optimality_tolerance
            observable, _ = _se_rank_nullspace(H * Z; rtol=rank_rtol)
            status = observable == size(Z, 2) ? :converged_unique : :converged_underobserved
            return ConstrainedStateEstimationResult(status, x, e, iteration - 1, rankC,
                                                    size(Z, 2), observable, history)
        end
        n = _se_scaled_normal_step(C, c, scale, normal_fraction * radius)

        # The tangential component is an exact linearised-null-space step.  Its
        # own radius leaves room for the normal component already taken.
        reduced_radius = sqrt(max(radius^2 - norm(scale .* n)^2, 0.0))
        t = zeros(Float64, length(x))
        if size(Z, 2) > 0 && reduced_radius > 0
            B = H * Z
            q = -(B \ (r + H * n))
            nq = norm(scale .* (Z * q))
            nq > reduced_radius && nq > 0 && (q .*= reduced_radius / nq)
            t .= Z * q
        end
        step = n + t
        model_r = r + H * step
        model_c = c + C * step
        predicted = _se_merit(e, μ) - (0.5 * sum(abs2, model_r) + μ * norm(model_c))

        # A zero/negative predicted reduction means the local model is no longer
        # useful; contract rather than accepting a coincidental raw decrease.
        if !(isfinite(predicted) && predicted > 0)
            radius *= 0.25
            rejected += 1
            radius < min_radius && break
            continue
        end

        trial_x = x + step
        trial = try
            evaluate_state_estimator(s, p, trial_x)
        catch err
            err isa DomainError || rethrow()
            nothing
        end
        accepted = false
        ρ = -Inf
        if trial !== nothing
            actual = _se_merit(e, μ) - _se_merit(trial, μ)
            ρ = actual / predicted
            accepted = isfinite(ρ) && ρ >= acceptance_threshold
        end

        # One second-order correction repairs nonlinear constraint curvature;
        # it is bounded and evaluated through the same merit acceptance test.
        if !accepted && trial !== nothing && !isempty(c)
            defect = -trial.constraints + c + C * step
            soc = _se_soc_step(C, defect, scale, 0.5 * reduced_radius)
            if norm(scale .* (step + soc)) <= radius * (1 + 1e-12)
                soc_trial = evaluate_state_estimator(s, p, x + step + soc)
                actual = _se_merit(e, μ) - _se_merit(soc_trial, μ)
                ρsoc = actual / predicted
                if isfinite(ρsoc) && ρsoc >= acceptance_threshold
                    step .+= soc
                    trial = soc_trial
                    ρ = ρsoc
                    accepted = true
                end
            end
        end

        if accepted
            x .+= step
            e = trial
            rejected = 0
            radius = ρ > 0.75 ? min(2radius, Float64(max_radius)) : radius
            push!(history, _se_history_entry(iteration, radius, _se_merit(e, μ), e))

            # Tangent-space first-order stationarity and exact-equation
            # feasibility are separate conditions; satisfying only one is not
            # convergence for a constrained estimator.
            Cnow = constraint_jacobian(s, p, x)
            _, Znow = _se_rank_nullspace(Cnow; rtol=rank_rtol)
            gred = size(Znow, 2) == 0 ? 0.0 : norm(Znow' * (residual_jacobian(s, p, x)' * e.residual), Inf)
            if norm(e.constraints) <= constraint_tolerance && gred <= optimality_tolerance
                Hnow = residual_jacobian(s, p, x)
                observable, _ = _se_rank_nullspace(Hnow * Znow; rtol=rank_rtol)
                rankC, _ = _se_rank_nullspace(Cnow; rtol=rank_rtol)
                status = observable == size(Znow, 2) ? :converged_unique : :converged_underobserved
                return ConstrainedStateEstimationResult(status, x, e, iteration, rankC,
                                                        size(Znow, 2), observable, history)
            end
        else
            radius *= 0.25
            rejected += 1
            radius < min_radius && break
        end
    end

    C = constraint_jacobian(s, p, x)
    rankC, Z = _se_rank_nullspace(C; rtol=rank_rtol)
    status = norm(e.constraints) > constraint_tolerance ?
        (rejected > 0 ? :constraint_restoration_failed : :infeasible_constraints) :
        (radius < min_radius ? :trust_region_stalled : :max_iterations)
    ConstrainedStateEstimationResult(status, x, e, max_iterations, rankC,
                                    size(Z, 2), 0, history)
end

"""
    solve_with_continuation(structure, parameters, x0; alphas=0:0.25:1, kwargs...)

Advance exact constant-power/ZIP constraints from the regularised internal model
(`α=0`) to their physical, unregularised equations (`α=1`).  Each stage warm
starts the dense constrained solver.  A returned `:power_flow_initialization_failed`
status means a stage could not establish feasibility; no final estimate is
claimed unless the last stage is at α=1.
"""
function solve_with_continuation(s::SEStructure, p::SEParameters,
                                 x0::AbstractVector{<:Real};
                                 alphas=collect(0.0:0.25:1.0), kwargs...)
    αs = Float64.(collect(alphas))
    !isempty(αs) || throw(ArgumentError("continuation needs at least one alpha"))
    all(α -> isfinite(α) && 0 <= α <= 1, αs) ||
        throw(ArgumentError("continuation alphas must lie in [0, 1]"))
    issorted(αs) || throw(ArgumentError("continuation alphas must be nondecreasing"))
    αs[end] == 1.0 || throw(ArgumentError("continuation must finish at α=1"))

    x = Float64.(x0)
    stages = ConstrainedStateEstimationResult[]
    for α in αs
        p.continuation_alpha = α
        result = solve_compiled_state_estimator(s, p, x; kwargs...)
        push!(stages, result)
        if !(result.status in (:converged_unique, :converged_underobserved))
            return ContinuationStateEstimationResult(:power_flow_initialization_failed,
                                                     result, αs, stages)
        end
        x = result.state
    end
    ContinuationStateEstimationResult(stages[end].status, stages[end], αs, stages)
end
