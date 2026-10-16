"""
    NetInjectionMeasurement(; kind, bus, terminal="1", reference=nothing, value, sigma)

Net current (`:ire`, `:iim`, `:imag`) or power (`:pinj`, `:qinj`) into the
passive network at one conductor. Unlike customer `Measurement` powers this
observes the **total** at a merged node and declares no unknown injection path.
Declare physical paths with `injection_connections` independently of telemetry.
Power uses terminal-to-`reference` voltage; current is reference-independent.
Ideal-ground reaction currents require element telemetry, not reduced Ybus.
"""
struct NetInjectionMeasurement <: AbstractMeasurement
    kind::Symbol
    bus::String
    terminal::String
    reference::Union{String,Nothing}
    value::Float64
    sigma::Float64
    function NetInjectionMeasurement(; kind::Symbol, bus::String, terminal::String="1",
            reference::Union{String,Nothing}=nothing, value::Real, sigma::Real)
        kind in (:ire, :iim, :imag, :pinj, :qinj) || throw(ArgumentError("unsupported net-injection kind :$kind"))
        _se_validate_meter(kind, terminal, reference, value, sigma)
        isempty(bus) && throw(ArgumentError("bus must be non-empty"))
        kind in (:ire, :iim, :imag) && reference !== nothing && throw(ArgumentError("current telemetry has no voltage reference"))
        new(kind, bus, terminal, reference, Float64(value), Float64(sigma))
    end
end

"""
    TransformerMeasurement(; kind, transformer, variant, port=1, terminal,
                             reference=nothing, value, sigma)

Current (`:ire`, `:iim`, `:imag`) or power (`:pflow`, `:qflow`) into a finite
transformer port. The id is looked up across BMOPF subtypes; optional `variant` disambiguates
a duplicated id. Ordinary transformers use
port 1 (from) or 2 (to); a center-tap's port 2 contains both secondary legs and
neutral. For `n_winding`, port is the winding-list index. Power is referenced
to earth unless a terminal on the selected port's bus is named explicitly.
Shared original nodes in n-winding ports cannot be disaggregated by Yprim.
"""
struct TransformerMeasurement <: AbstractMeasurement
    kind::Symbol
    transformer::String
    variant::String
    port::Int
    terminal::String
    reference::Union{String,Nothing}
    value::Float64
    sigma::Float64
    function TransformerMeasurement(; kind::Symbol, transformer::String, variant::Union{Nothing,String}=nothing,
            port::Integer=1, terminal::String, reference::Union{String,Nothing}=nothing,
            value::Real, sigma::Real)
        _se_validate_meter(kind, terminal, reference, value, sigma; passive=true)
        isempty(transformer) && throw(ArgumentError("transformer must be non-empty"))
        variant===nothing || !isempty(variant) || throw(ArgumentError("variant must be non-empty when provided"))
        port > 0 || throw(ArgumentError("port must be positive"))
        new(kind, transformer, something(variant,""), Int(port), terminal, reference, Float64(value), Float64(sigma))
    end
end

"""
    ShuntMeasurement(; kind, shunt, terminal, reference=nothing, value, sigma)

Current or power entering a named finite BMOPF shunt, including a neutral
grounding electrode. Current is positive from the conductor into the shunt;
positive power is absorbed by it. Uses the shunt's full coupled admittance.
The kinds and voltage-reference convention match `BranchMeasurement`.
"""
struct ShuntMeasurement <: AbstractMeasurement
    kind::Symbol
    shunt::String
    terminal::String
    reference::Union{String,Nothing}
    value::Float64
    sigma::Float64
    function ShuntMeasurement(; kind::Symbol, shunt::String, terminal::String,
            reference::Union{String,Nothing}=nothing, value::Real, sigma::Real)
        _se_validate_meter(kind, terminal, reference, value, sigma; passive=true)
        isempty(shunt) && throw(ArgumentError("shunt must be non-empty"))
        new(kind, shunt, terminal, reference, Float64(value), Float64(sigma))
    end
end

function _se_validate_meter(kind, terminal, reference, value, sigma; passive=false)
    passive && !(kind in _BRANCH_MEASUREMENT_KINDS) && throw(ArgumentError("unsupported element telemetry kind :$kind"))
    _validate_measurement_scalar(kind, value, sigma)
    isempty(terminal) && throw(ArgumentError("terminal must be non-empty"))
    reference isa String && isempty(reference) && throw(ArgumentError("reference must be non-empty"))
    nothing
end

"""
    FunctionalMeasurement(; name, nodes, prediction, gradient=nothing, value, sigma)

A pure scalar measurement function of selected SI conductor voltages. Its input
is `[real(V[nodes]); imag(V[nodes])]`, retaining duplicates and grounded zeros.
`prediction(z)` returns one finite real number; supply `gradient(z)` or support
ForwardDiff dual numbers for automatic first derivatives. Values/sigmas may
change through `SEParameters`; functions, node order and captured constants
belong to the compiled pattern. Callbacks must be pure and must not mutate
inputs or captured data. They add no injection paths or exact equations.
All three shared-model solvers use this same function and Jacobian.
"""
struct FunctionalMeasurement{F,G} <: AbstractMeasurement
    kind::Symbol
    name::Symbol
    nodes::Vector{TerminalID}
    prediction::F
    gradient::G
    value::Float64
    sigma::Float64
    function FunctionalMeasurement(; name::Symbol, nodes, prediction, gradient=nothing,
            value::Real, sigma::Real)
        _validate_measurement_scalar(:functional, value, sigma)
        ns = TerminalID.(collect(nodes))
        isempty(ns) && throw(ArgumentError("functional measurement needs at least one node"))
        all(n -> !isempty(n[1]) && !isempty(n[2]), ns) || throw(ArgumentError("node names must be non-empty"))
        new{typeof(prediction),typeof(gradient)}(:functional, name, ns, prediction, gradient, Float64(value), Float64(sigma))
    end
end

struct _SEFunctionalSpec{F,G}
    kind::Symbol
    nodes::Vector{Int}
    prediction::F
    gradient::G
end

const _SEPassiveTelemetry = Union{BranchMeasurement,TransformerMeasurement,ShuntMeasurement}
const _SESupportedMeasurement = Union{Measurement,NetInjectionMeasurement,_SEPassiveTelemetry,FunctionalMeasurement}
_resolve_ref(m::NetInjectionMeasurement, neutral) = m.reference
_se_measurement_key(m::NetInjectionMeasurement) = (:net, m.kind, m.bus, m.terminal, m.reference)
_se_measurement_key(m::TransformerMeasurement) = (:transformer, m.kind, m.variant, m.transformer, m.port, m.terminal, m.reference)
_se_measurement_key(m::ShuntMeasurement) = (:shunt, m.kind, m.shunt, m.terminal, m.reference)
_se_measurement_key(m::FunctionalMeasurement) = (:functional, m.name, Tuple(m.nodes), m.prediction, m.gradient)

_se_measurement_description(m::Union{Measurement,NetInjectionMeasurement}, neutral) =
    (kind=m.kind, location="$(m.bus)/$(m.terminal)", reference=_resolve_ref(m, neutral))
_se_measurement_description(m::BranchMeasurement, neutral) =
    (kind=m.kind, location="line/$(m.line)/$(m.side)/$(m.terminal)", reference=m.reference)
_se_measurement_description(m::TransformerMeasurement, neutral) =
    (kind=m.kind, location="transformer/$(m.variant)/$(m.transformer)/$(m.port)/$(m.terminal)", reference=m.reference)
_se_measurement_description(m::ShuntMeasurement, neutral) =
    (kind=m.kind, location="shunt/$(m.shunt)/$(m.terminal)", reference=m.reference)
_se_measurement_description(m::FunctionalMeasurement, neutral) =
    (kind=:functional, location=string(m.name), reference=Tuple(m.nodes))

function _se_meter_primitive(m::BranchMeasurement, net)
    el = get(get(net, "line", Dict()), m.line, nothing)
    el isa AbstractDict || throw(ArgumentError("line '$(m.line)' is absent"))
    side = m.side === :from ? "from" : "to"
    bus = String(el["bus_$side"]); ts = String.(el["terminal_map_$side"])
    m.terminal in ts || throw(ArgumentError("terminal '$(m.terminal)' is not on line '$(m.line)' $side side"))
    nodes, Y = line_yprim(el, get(net, "linecode", Dict()))
    k = findfirst(==(m.terminal), ts)
    row = m.side === :from ? k : length(el["terminal_map_from"]) + k
    nodes, Y, row, bus
end

function _se_meter_primitive(m::TransformerMeasurement, net)
    variants=[kind for (kind,group) in get(net,"transformer",Dict()) if haskey(group,m.transformer)]
    variant=isempty(m.variant) ? (length(variants)==1 ? only(variants) :
        throw(ArgumentError("transformer id must be present and unique; provide variant to disambiguate"))) : m.variant
    el = get(get(get(net, "transformer", Dict()), variant, Dict()), m.transformer, nothing)
    el isa AbstractDict || throw(ArgumentError("transformer $(m.variant)/$(m.transformer) is absent"))
    if variant == "n_winding"
        ws = BMOPFTools._nw_windings(el)
        m.port <= length(ws) || throw(ArgumentError("transformer port $(m.port) is absent"))
        w = ws[m.port]; bus = String(w.bus); ts = String.(w.terminal_map)
        any(k -> k != m.port && ws[k].bus == bus && m.terminal in ws[k].terminal_map, eachindex(ws)) &&
            throw(ArgumentError("n-winding telemetry cannot separate ports sharing the same original terminal"))
        nodes, Y = BMOPFTools.nwinding_yprim(el)
    else
        m.port in (1, 2) || throw(ArgumentError("ordinary transformer ports are 1 (from) or 2 (to)"))
        side = m.port == 1 ? "from" : "to"
        bus = String(el["bus_$side"]); ts = String.(el["terminal_map_$side"])
        if variant == "open_delta_regulator" && (bus,m.terminal) in _se_open_delta_bond(el)
            throw(ArgumentError("open-delta shared-phase port current includes an ideal-wire reaction not recoverable from voltage alone; meter adjacent finite branches instead"))
        end
        nodes, Y = transformer_yprim(el, variant)
    end
    m.terminal in ts || throw(ArgumentError("terminal '$(m.terminal)' is absent on transformer port $(m.port)"))
    matches = findall(==((bus, m.terminal)), nodes)
    length(matches) == 1 || throw(ArgumentError("transformer primitive cannot identify this port conductor uniquely"))
    nodes, Y, only(matches), bus
end

function _se_meter_primitive(m::ShuntMeasurement, net)
    el = get(get(net, "shunt", Dict()), m.shunt, nothing)
    el isa AbstractDict || throw(ArgumentError("shunt '$(m.shunt)' is absent"))
    nodes, Y = BMOPFTools._shunt_yprim(el)
    bus = String(el["bus"])
    row = findfirst(==((bus, m.terminal)), nodes)
    row === nothing && throw(ArgumentError("terminal '$(m.terminal)' is absent on shunt '$(m.shunt)'"))
    nodes, Y, row, bus
end

function _se_compile_passive_measurement(m, net, node_index)
    (data, warnings) = _se_capture_primitive(() -> _se_meter_primitive(m, net))
    isempty(warnings) || throw(ArgumentError("telemetry primitive warned: $warnings"))
    nodes, Y, row, bus = data
    isempty(nodes) && throw(ArgumentError("element has no finite telemetry primitive"))
    size(Y) == (length(nodes), length(nodes)) && all(isfinite, Y) || throw(ArgumentError("invalid telemetry primitive"))
    all(node -> haskey(node_index, node), nodes) || throw(ArgumentError("telemetry primitive has absent nodes"))
    ref = m.reference === nothing ? 0 : get(node_index, (bus, m.reference), -1)
    ref >= 0 || throw(ArgumentError("telemetry reference '$(m.reference)' is absent on bus '$bus'"))
    ns = Int[node_index[node] for node in nodes]
    _SEBranchMeasurementSpec(m.kind, ns, ns[row], ref, ComplexF64.(Y[row, :]))
end

function _se_functional_input(spec, vr, vi)
    vcat([i == 0 ? 0.0 : vr[i] for i in spec.nodes],
         [i == 0 ? 0.0 : vi[i] for i in spec.nodes])
end
function _se_functional_value(spec, vr, vi)
    value = spec.prediction(_se_functional_input(spec, vr, vi))
    value isa Real || throw(ArgumentError("functional prediction must be scalar and real"))
    isfinite(value) || throw(DomainError(value, "nonfinite functional prediction"))
    Float64(value)
end
function _se_functional_jacobian(spec, s, vr, vi)
    z = _se_functional_input(spec, vr, vi)
    g = spec.gradient === nothing ? ForwardDiff.gradient(spec.prediction, z) : spec.gradient(z)
    g isa AbstractVector && length(g) == length(z) || throw(DimensionMismatch("functional gradient must have $(length(z)) entries"))
    all(v -> v isa Real && isfinite(v), g) || throw(DomainError(g, "nonfinite functional gradient"))
    row = spzeros(Float64, size(s.voltage_state_jacobian, 2))
    n = length(s.nodes); k = length(spec.nodes)
    for (j, i) in enumerate(spec.nodes)
        i == 0 && continue
        row += g[j] .* s.voltage_state_rows[:, i] + g[k+j] .* s.voltage_state_rows[:, n+i]
    end
    row
end

_se_spec_value(spec::_SEFunctionalSpec, vr, vi, ir, ii, p) = _se_functional_value(spec, vr, vi)
_se_spec_value(spec::_SEBranchMeasurementSpec, vr, vi, ir, ii, p) = _branch_measurement_value(spec, vr, vi, p.current_epsilon)
_se_spec_value(spec::_SEMeasurementSpec, vr, vi, ir, ii, p) = _measurement_value(spec, vr, vi, ir, ii,
    spec.kind === :imag ? p.current_epsilon : p.magnitude_epsilon)
