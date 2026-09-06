# Fundamental-frequency generators. All public data are SI RMS quantities.
# Internal EMFs, neutral/earth returns, and powers are expressions. In particular,
# no internal voltage nodes are introduced merely to represent a series drop.
const _GeneratorBound = Union{Nothing,Float64,Vector{Float64}}

"""
    GeneratorVoltageLaw(mode=:none; phasor=[], angles=[], magnitude_min=nothing,
                        magnitude_max=nothing)

Internal EMF law. Modes are `:none`, `:fixed_phasor` (b),
`:fixed_magnitudes` (c, rotating `phasor` template), `:common_magnitude` (d),
and `:phase_magnitudes` (e). `phasor` is an SI RMS complex vector for b/c.
Angles are relative offsets in radians for d/e (default 0 for one port and
0,-2π/3,2π/3 for three). Two or other port counts require explicit offsets.
Magnitude bounds are optional SI volts, scalar or per-port vectors; common
magnitude uses scalar bounds. Mode e requires explicit strictly positive lower
bounds to exclude the undefined-angle domain. No law adds an absolute reference
except b. See the generalized-generator model guide for degrees of freedom.
"""
Base.@kwdef struct GeneratorVoltageLaw
    mode::Symbol = :none
    phasor::Vector{ComplexF64} = ComplexF64[]
    angles::Vector{Float64} = Float64[]
    magnitude_min::_GeneratorBound = nothing
    magnitude_max::_GeneratorBound = nothing
end
GeneratorVoltageLaw(mode::Symbol; kwargs...) = GeneratorVoltageLaw(; mode, kwargs...)

"""
    GeneratorControl(; p=nothing, q=nothing, power_location=:poc,
                       voltage_target=nothing, voltage_metric=:positive_sequence)

Fixed active/reactive power in W/var: a scalar specifies aggregate power, a
vector specifies oriented port powers. Omission leaves dispatch free. Scalar
PCC power for a grounded source includes its neutral conductor contribution.
`power_location` is `:poc` or `:internal`. An optional positive voltage target
regulates PCC positive-sequence magnitude (three ports), or `:phase` magnitude
(one port). This hard target requires fixed P, free Q and a common-magnitude
source law; reaching a capability limit may make it infeasible. No automatic
PV/PQ switching, slack allocation, or controller saturation is implied.
"""
Base.@kwdef struct GeneratorControl
    p::_GeneratorBound = nothing
    q::_GeneratorBound = nothing
    power_location::Symbol = :poc
    voltage_target::Union{Nothing,Float64} = nothing
    voltage_metric::Symbol = :positive_sequence
end

"""
    GeneratorSequenceLimits(; location=:poc, voltage_min=nothing,
        voltage_max=nothing, current_min=nothing, current_max=nothing)

Optional three-entry magnitude limits in order **zero, positive, negative**.
Voltage location is `:poc` or `:internal`; currents are the three oriented port
currents. Uses amplitude-invariant Fortescue (1/3 normalization). Three ports
must be explicitly ordered a,b,c and share the same return. Bounds are SI V/A.
Zero upper bounds become linear real/imaginary equalities. These limits do not
replace phase, neutral, earth, or converter hardware limits.
"""
Base.@kwdef struct GeneratorSequenceLimits
    location::Symbol = :poc
    voltage_min::Union{Nothing,Vector{Float64}} = nothing
    voltage_max::Union{Nothing,Vector{Float64}} = nothing
    current_min::Union{Nothing,Vector{Float64}} = nothing
    current_max::Union{Nothing,Vector{Float64}} = nothing
end

"""
    GeneratorCapability(; kwargs...)

Optional SI limits. `p_min/max`, `q_min/max`, `s_max` are scalar aggregate or
per-port vectors at `power_location=:poc` or `:internal`. `i_max` limits port
currents; `terminal_i_max` limits all external conductors in result terminal
order. Scalars apply to every current/voltage entry. `voltage_min/max` constrain
PCC port magnitudes; internal magnitudes belong to `GeneratorVoltageLaw`.
`earth_i_max` is a scalar ground-current limit for `SourceGenerator` only.
`sequence` optionally holds `GeneratorSequenceLimits`. Omitted bounds add no
constraints; apparent-power limits use P/Q lifts to keep every row quadratic.
"""
Base.@kwdef struct GeneratorCapability
    p_min::_GeneratorBound = nothing
    p_max::_GeneratorBound = nothing
    q_min::_GeneratorBound = nothing
    q_max::_GeneratorBound = nothing
    s_max::_GeneratorBound = nothing
    power_location::Symbol = :poc
    i_max::_GeneratorBound = nothing
    terminal_i_max::_GeneratorBound = nothing
    voltage_min::_GeneratorBound = nothing
    voltage_max::_GeneratorBound = nothing
    earth_i_max::Union{Nothing,Float64} = nothing
    sequence::Union{Nothing,GeneratorSequenceLimits} = nothing
end

"""
    GeneralizedGenerator(; id, bus, connections, impedance=0, voltage=..., control=...,
                           capability=..., voltage_scale=230, cost=0)

An arbitrary independent collection of oriented two-terminal source ports.
`connections=[("a","n"), ...]` specifies outgoing/return terminal pairs at
`bus`. Wye, single-phase line-neutral/line-line and split-phase are supported;
closed connection cycles are rejected. `impedance` is a scalar, diagonal vector,
or full port impedance matrix in ohms, including any already-reduced neutral
drop. It must be passive at the modeled frequency. This component has no earth
path: its terminal currents sum to zero. Use `SourceGenerator` for grounding.
`voltage_scale` is a positive numerical scaling/start value, never a bound.
`cost` is currency/kWh on aggregate PCC injection. Stamping uses the public
BMOPFTools staged API and returns model-unit expressions and physical bases.
"""
Base.@kwdef struct GeneralizedGenerator <: AbstractDevice
    id::String
    bus::String
    connections::Vector{Tuple{String,String}}
    impedance::Union{Number,AbstractVector,AbstractMatrix} = 0.0
    voltage::GeneratorVoltageLaw = GeneratorVoltageLaw()
    control::GeneratorControl = GeneratorControl()
    capability::GeneratorCapability = GeneratorCapability()
    voltage_scale::Float64 = 230.0
    cost::Float64 = 0.0
end

"""
    SourceGenerator(; id, bus, phase_terminals=["a","b","c"], neutral="n",
                      impedance=0, grounding=0, voltage, kwargs...)

Dedicated source with an internal star point, independent external neutral
current, and an earth return. `impedance` is the full conductor primitive in
order phases then neutral (scalar/diagonal/full, ohms). `grounding=nothing`
opens the star-earth connection; zero is an exact ideal bond; a finite complex
impedance models a star-to-remote-earth (0 V) connection. Outward currents obey
sum(phase currents)+neutral current+earth current=0. PCC-neutral grounding must
be declared separately in network data. A coincident ideal engine ground and
ideal source bond with an ideal neutral lead is rejected as duplicate grounding.
The remaining keywords have the same meaning as `GeneralizedGenerator`.
"""
Base.@kwdef struct SourceGenerator <: AbstractDevice
    id::String
    bus::String
    phase_terminals::Vector{String} = ["a", "b", "c"]
    neutral::String = "n"
    impedance::Union{Number,AbstractVector,AbstractMatrix} = 0.0
    grounding::Union{Nothing,Number} = 0.0
    voltage::GeneratorVoltageLaw = GeneratorVoltageLaw()
    control::GeneratorControl = GeneratorControl()
    capability::GeneratorCapability = GeneratorCapability()
    voltage_scale::Float64 = 230.0
    cost::Float64 = 0.0
end

const _GeneralizedDevice = Union{GeneralizedGenerator,SourceGenerator}
_gg_connections(d::GeneralizedGenerator) = d.connections
_gg_connections(d::SourceGenerator) = [(p, d.neutral) for p in d.phase_terminals]
_gg_terminals(d) = unique(vcat(first.(_gg_connections(d)), last.(_gg_connections(d))))
_gg_n(d) = length(_gg_connections(d))
_gg_vector(x, n) = x === nothing ? fill(nothing, n) : x isa Real ? fill(Float64(x), n) : Float64.(x)
function _gg_matrix(z, n)
    Z = z isa Number ? Matrix(Diagonal(fill(ComplexF64(z), n))) :
        z isa AbstractVector ? Matrix(Diagonal(ComplexF64.(z))) : Matrix{ComplexF64}(z)
    size(Z) == (n, n) || throw(ArgumentError("impedance must have size ($n,$n)"))
    all(isfinite, Z) || throw(ArgumentError("impedance must be finite"))
    eigmin(Hermitian((Z + Z') / 2)) >= -1e-12 * max(opnorm(Z), eps()) ||
        throw(ArgumentError("physical series impedance must be passive"))
    Z
end

"""
    generator_sequence_impedance(zero, positive, negative)

Return a phase-domain 3×3 impedance in ohms from diagonal sequence impedances,
ordered zero/positive/negative, with a-b-c phase order. This is a port matrix for
a three-phase `GeneralizedGenerator`, not a four-conductor source primitive.
"""
function generator_sequence_impedance(zero, positive, negative)
    T = _gg_transform()
    _gg_matrix(T \ (Diagonal(ComplexF64[zero, positive, negative]) * T), 3)
end
function _gg_transform()
    a = cis(2pi / 3)
    ComplexF64[1 1 1; 1 a a^2; 1 a^2 a] / 3
end

function _gg_check_bound(x, n, name; nonnegative=false)
    x === nothing && return
    x isa Real || x isa AbstractVector{<:Real} || throw(ArgumentError("$name must be real"))
    x isa AbstractVector && length(x) != n && throw(ArgumentError("$name needs $n entries"))
    all(v -> isfinite(v) && (!nonnegative || v >= 0), _gg_vector(x, n)) ||
        throw(ArgumentError("$name must be finite$(nonnegative ? " and nonnegative" : "")"))
end
function _gg_check_pair(lo, hi, n, name; nonnegative=false, aggregate=false)
    _gg_check_bound(lo, n, "$name lower"; nonnegative)
    _gg_check_bound(hi, n, "$name upper"; nonnegative)
    lo === nothing || hi === nothing || (aggregate && (lo isa Real)!=(hi isa Real)) || all(_gg_vector(lo,n) .<= _gg_vector(hi,n)) ||
        throw(ArgumentError("$name lower exceeds upper"))
end
function _gg_angles(law, n)
    !isempty(law.angles) && return law.angles
    n == 1 && return [0.0]
    n == 3 && return [0.0, -2pi/3, 2pi/3]
    throw(ArgumentError("explicit relative angles required for $n ports"))
end

# Reduce exact internal sequence restrictions in a rotating phase-magnitude law
# in magnitude space, before creating constraints. This avoids adding redundant
# real/imaginary sequence rows to already fixed relative-angle equations.
function _gg_shape(d)
    law=d.voltage; n=_gg_n(d)
    law.mode==:common_magnitude && return cis.(_gg_angles(law,n).-first(_gg_angles(law,n)))
    law.mode==:phase_magnitudes || return nothing
    seq=d.capability.sequence
    seq===nothing || seq.location!=:internal || seq.voltage_max===nothing || begin
        zeros_idx=findall(iszero,seq.voltage_max)
        isempty(zeros_idx) && return nothing
        phi=_gg_angles(law,n)
        M=_gg_transform()[zeros_idx,:]*Diagonal(cis.(phi))
        N=nullspace(vcat(real(M),imag(M));atol=1e-12)
        size(N,2)==1 || throw(ArgumentError("exact internal sequence zeros require a unique positive magnitude shape for :phase_magnitudes; use :none or an explicit template"))
        amplitudes=N[:,1]/N[1,1]
        all(x -> isfinite(x) && x>1e-12,amplitudes) || throw(ArgumentError("sequence zeros conflict with positive phase magnitudes and relative angles"))
        return amplitudes.*cis.(phi.-phi[1])
    end
    nothing
end
function _gg_shape_bounds(d,ratios)
    n=length(ratios); law=d.voltage
    lo,hi=nothing,nothing
    for (coef,l,u) in zip(abs.(ratios),_gg_vector(law.magnitude_min,n),_gg_vector(law.magnitude_max,n))
        l===nothing || (lo=lo===nothing ? l/coef : max(lo,l/coef))
        u===nothing || (hi=hi===nothing ? u/coef : min(hi,u/coef))
    end
    seq=d.capability.sequence
    if seq!==nothing && seq.location==:internal
        coeff=abs.(_gg_transform()*ratios)
        coeff[coeff .< 1e-12*maximum(abs,ratios)] .= 0.0
        for (coef,l,u) in zip(coeff,_gg_vector(seq.voltage_min,3),_gg_vector(seq.voltage_max,3))
            if iszero(coef)
                l===nothing || iszero(l) || throw(ArgumentError("source shape has zero sequence voltage below requested minimum"))
            else
                l===nothing || (lo=lo===nothing ? l/coef : max(lo,l/coef))
                u===nothing || (hi=hi===nothing ? u/coef : min(hi,u/coef))
            end
        end
    end
    lo===nothing || hi===nothing || lo<=hi || throw(ArgumentError("source shape has conflicting magnitude bounds"))
    lo,hi
end

function validate_device(d::_GeneralizedDevice, nets; periods::Integer=length(nets))
    periods == length(nets) || throw(ArgumentError("period count mismatch"))
    isempty(strip(d.id)) && throw(ArgumentError("generator id cannot be empty"))
    pairs = _gg_connections(d); n = length(pairs); ts = _gg_terminals(d)
    n > 0 || throw(ArgumentError("generator needs at least one port"))
    all(p -> p[1] != p[2], pairs) || throw(ArgumentError("a port cannot connect to itself"))
    C = [Float64(t == p[1]) - Float64(t == p[2]) for t in ts, p in pairs]
    rank(C) == n || throw(ArgumentError("closed/duplicate winding connections are not supported"))
    Z = _gg_matrix(d.impedance, d isa SourceGenerator ? n+1 : n)
    isfinite(d.voltage_scale) && d.voltage_scale > 0 || throw(ArgumentError("voltage_scale must be positive"))
    isfinite(d.cost) || throw(ArgumentError("cost must be finite"))
    for net in nets
        bus = get(get(net, "bus", Dict()), d.bus, nothing)
        bus === nothing && throw(ArgumentError("generator bus '$(d.bus)' not found"))
        all(t -> t in bus["terminal_names"], ts) || throw(ArgumentError("generator terminal not found"))
        if d isa SourceGenerator
            if d.grounding !== nothing
                isfinite(d.grounding) && real(d.grounding) >= 0 || throw(ArgumentError("invalid grounding impedance"))
                if iszero(d.grounding) && all(iszero, Z[end,:]) &&
                        d.neutral in get(bus, "perfectly_grounded_terminals", [])
                    throw(ArgumentError("duplicate ideal source/engine grounding at $(d.bus).$(d.neutral); declare this bond once"))
                end
            end
        end
    end
    law = d.voltage
    law.mode in (:none, :fixed_phasor, :fixed_magnitudes, :common_magnitude, :phase_magnitudes) ||
        throw(ArgumentError("unknown generator voltage mode $(law.mode)"))
    if law.mode in (:fixed_phasor, :fixed_magnitudes)
        length(law.phasor) == n && all(z -> isfinite(z) && abs(z) > 0, law.phasor) ||
            throw(ArgumentError("fixed modes require $n finite nonzero phasors"))
        isempty(law.angles) || throw(ArgumentError("fixed templates already specify their angles"))
    else
        isempty(law.phasor) || throw(ArgumentError("phasor is only used by fixed modes"))
        if law.mode != :none
            ang = _gg_angles(law, n)
            length(ang) == n && all(isfinite, ang) || throw(ArgumentError("invalid relative angles"))
        else
            isempty(law.angles) || throw(ArgumentError("angles require an angle-constrained law"))
        end
    end
    law.mode == :common_magnitude &&
        (law.magnitude_min isa Vector || law.magnitude_max isa Vector) &&
        throw(ArgumentError("common magnitude needs scalar magnitude bounds"))
    _gg_check_pair(law.magnitude_min,law.magnitude_max,n,"EMF magnitude";nonnegative=true)
    if law.mode == :phase_magnitudes
        law.magnitude_min !== nothing && all(>(0), _gg_vector(law.magnitude_min,n)) ||
            throw(ArgumentError("phase magnitudes require explicit positive minima for the angle domain"))
    end
    if law.mode in (:fixed_phasor,:fixed_magnitudes)
        lo = _gg_vector(law.magnitude_min,n); hi = _gg_vector(law.magnitude_max,n)
        all(k -> (lo[k] === nothing || lo[k] <= abs(law.phasor[k])) &&
                 (hi[k] === nothing || abs(law.phasor[k]) <= hi[k]), 1:n) ||
            throw(ArgumentError("fixed EMF violates its declared magnitude limits"))
    end
    c, cap = d.control, d.capability
    c.power_location in (:poc,:internal) && cap.power_location in (:poc,:internal) ||
        throw(ArgumentError("power location must be :poc or :internal"))
    _gg_check_bound(c.p,n,"P target"); _gg_check_bound(c.q,n,"Q target")
    _gg_check_pair(cap.p_min,cap.p_max,n,"P";aggregate=true)
    _gg_check_pair(cap.q_min,cap.q_max,n,"Q";aggregate=true)
    for (x,nx,name) in ((cap.s_max,n,"S"),(cap.i_max,n,"port current"),
                       (cap.terminal_i_max,length(ts),"terminal current"))
        _gg_check_bound(x,nx,name;nonnegative=true)
    end
    _gg_check_pair(cap.voltage_min,cap.voltage_max,n,"PCC voltage";nonnegative=true)
    _gg_check_bound(cap.earth_i_max,1,"earth current";nonnegative=true)
    cap.earth_i_max !== nothing && !(d isa SourceGenerator) &&
        throw(ArgumentError("earth current capability requires SourceGenerator"))
    seq = cap.sequence
    needs_sequence = seq !== nothing || (c.voltage_target !== nothing && c.voltage_metric == :positive_sequence)
    if needs_sequence
        n == 3 && length(unique(last.(pairs))) == 1 ||
            throw(ArgumentError("sequence quantities require three ordered ports with a common return"))
    end
    if seq !== nothing
        seq.location in (:poc,:internal) || throw(ArgumentError("invalid sequence voltage location"))
        _gg_check_pair(seq.voltage_min,seq.voltage_max,3,"sequence voltage";nonnegative=true)
        _gg_check_pair(seq.current_min,seq.current_max,3,"sequence current";nonnegative=true)
        if seq.location==:internal && law.mode in (:fixed_phasor,:fixed_magnitudes)
            mag=abs.(_gg_transform()*law.phasor)
            all(k -> (seq.voltage_min===nothing || seq.voltage_min[k]<=mag[k]+1e-10) &&
                     (seq.voltage_max===nothing || seq.voltage_max[k]>=mag[k]-1e-10),1:3) ||
                throw(ArgumentError("fixed source template violates internal sequence bounds"))
        end
    end
    shape=_gg_shape(d)
    shape===nothing || _gg_shape_bounds(d,shape)
    if c.voltage_target !== nothing
        isfinite(c.voltage_target) && c.voltage_target > 0 || throw(ArgumentError("voltage target must be positive"))
        c.p !== nothing && c.q === nothing && law.mode == :common_magnitude ||
            throw(ArgumentError("PV requires fixed P, free Q, and :common_magnitude EMF"))
        c.voltage_metric == :positive_sequence || (c.voltage_metric == :phase && n == 1) ||
            throw(ArgumentError("voltage metric must be :positive_sequence, or :phase for one port"))
    end
    nothing
end

# Each row is registered on the public engine ledger. Linear simplification
# handles fixed terms and exact-zero limits without zero-gradient norm rows.
struct _GGBuilder
    ctx
    id::String
    rows::Dict{String,Any}
    magnitudes::Dict{Any,Any}
    affine_equalities::Set{Any}
end
_GGBuilder(ctx,id,rows)=_GGBuilder(ctx,id,rows,Dict{Any,Any}(),Set{Any}())
_gg_model(b) = _opf_model(b.ctx)
function _gg_record!(b, name, row)
    b.rows[name] = row
    BMOPFTools.register_opf_constraint!(b.ctx, :generalized_generator, (b.id,name), row)
    row
end
function _gg_simplify(x)
    x isa Number && return x
    if x isa JuMP.VariableRef
        return JuMP.is_fixed(x) ? JuMP.fix_value(x) : x
    elseif x isa JuMP.AffExpr
        y = JuMP.AffExpr(JuMP.constant(x))
        for (coef,v) in JuMP.linear_terms(x)
            iszero(coef) && continue
            JuMP.is_fixed(v) ? JuMP.add_to_expression!(y,coef*JuMP.fix_value(v)) :
                              JuMP.add_to_expression!(y,coef,v)
        end
        isempty(y.terms) && return y.constant
        if iszero(y.constant) && length(y.terms)==1
            v,coef=first(y.terms)
            coef==1 && return v
        end
        return y
    end
    x
end
function _gg_eq!(b,name,x,y=0.0;scale=1.0)
    f = _gg_simplify(x-y)
    if f isa Number
        abs(f/scale) <= 1e-12 || throw(ArgumentError("inconsistent generator equality $name"))
        return nothing
    end
    f isa JuMP.VariableRef && (f=JuMP.AffExpr(0.0,f=>1.0))
    if f isa JuMP.AffExpr
        terms=sort!([(JuMP.index(v).value,c) for (c,v) in JuMP.linear_terms(f) if !iszero(c)])
        factor=last(first(terms))
        key=(JuMP.constant(f)/factor,Tuple((idx,c/factor) for (idx,c) in terms))
        key in b.affine_equalities && return nothing
        push!(b.affine_equalities,key)
    end
    if f isa JuMP.AffExpr && length(f.terms) == 1
        v, coef = first(f.terms)
        target = -f.constant/coef
        JuMP.has_lower_bound(v) && target < JuMP.lower_bound(v) - 1e-12*scale &&
            throw(ArgumentError("$name conflicts with a lower bound"))
        JuMP.has_upper_bound(v) && target > JuMP.upper_bound(v) + 1e-12*scale &&
            throw(ArgumentError("$name conflicts with an upper bound"))
        JuMP.fix(v,target;force=true)
        return _gg_record!(b,name,JuMP.FixRef(v))
    end
    _gg_record!(b,name,JuMP.@constraint(_gg_model(b), f/scale == 0))
end
function _gg_bound!(b,name,x,lo,hi;scale=1.0)
    x = _gg_simplify(x)
    if x isa Number
        (lo === nothing || x >= lo-1e-12*scale) && (hi === nothing || x <= hi+1e-12*scale) ||
            throw(ArgumentError("fixed quantity violates $name"))
    elseif lo !== nothing && hi !== nothing && lo == hi
        _gg_eq!(b,name,x,lo;scale)
    elseif x isa JuMP.VariableRef
        lo === nothing || JuMP.set_lower_bound(x, max(lo,JuMP.has_lower_bound(x) ? JuMP.lower_bound(x) : -Inf))
        hi === nothing || JuMP.set_upper_bound(x, min(hi,JuMP.has_upper_bound(x) ? JuMP.upper_bound(x) : Inf))
        lo === nothing || _gg_record!(b,name*"_lo",JuMP.LowerBoundRef(x))
        hi === nothing || _gg_record!(b,name*"_hi",JuMP.UpperBoundRef(x))
    else
        lo === nothing || _gg_record!(b,name*"_lo",JuMP.@constraint(_gg_model(b),x/scale >= lo/scale))
        hi === nothing || _gg_record!(b,name*"_hi",JuMP.@constraint(_gg_model(b),x/scale <= hi/scale))
    end
    nothing
end
function _gg_emit_mag!(b,name,r,i,lo,hi)
    lo === nothing && hi === nothing && return
    r,i = _gg_simplify(r),_gg_simplify(i)
    if r isa Number && i isa Number
        _gg_bound!(b,name,hypot(r,i),lo,hi); return
    end
    if hi !== nothing && iszero(hi)
        _gg_eq!(b,name*"_r",r); _gg_eq!(b,name*"_i",i); return
    end
    if hi !== nothing
        r isa JuMP.VariableRef && _gg_bound!(b,name*"_box_r",r,-hi,hi)
        i isa JuMP.VariableRef && _gg_bound!(b,name*"_box_i",i,-hi,hi)
    end
    if lo !== nothing && hi !== nothing && lo == hi
        _gg_eq!(b,name,(r/hi)^2+(i/hi)^2,1.0)
    else
        hi === nothing || _gg_record!(b,name*"_hi",JuMP.@constraint(_gg_model(b),(r/hi)^2+(i/hi)^2 <= 1))
        lo === nothing || iszero(lo) || _gg_record!(b,name*"_lo",JuMP.@constraint(_gg_model(b),(r/lo)^2+(i/lo)^2 >= 1))
    end
end
function _gg_mag!(b,name,r,i,lo,hi)
    lo===nothing && hi===nothing && return
    # Canonical squared affine norm identifies identical port/conductor circles,
    # including sign-reversed single-phase returns. Intersect their declarations
    # before stamping so the tight limit has exactly one physical circle.
    r,i=_gg_simplify(r),_gg_simplify(i)
    if r isa Number && i isa Number
        return _gg_emit_mag!(b,name,r,i,lo,hi)
    end
    q=r^2+i^2
    lin=sort!([(JuMP.index(v).value,c) for (c,v) in JuMP.linear_terms(q) if !iszero(c)])
    quad=sort!([(min(JuMP.index(v).value,JuMP.index(w).value),
                 max(JuMP.index(v).value,JuMP.index(w).value),c)
                for (c,v,w) in JuMP.quad_terms(q) if !iszero(c)])
    key=(JuMP.constant(q),Tuple(lin),Tuple(quad))
    if haskey(b.magnitudes,key)
        old=b.magnitudes[key]
        lo=lo===nothing ? old.lo : old.lo===nothing ? lo : max(lo,old.lo)
        hi=hi===nothing ? old.hi : old.hi===nothing ? hi : min(hi,old.hi)
        name=old.name
        r,i=old.r,old.i # retain the original variable orientation for box bounds
    end
    lo===nothing || hi===nothing || lo<=hi || throw(ArgumentError("conflicting magnitude limits for $name"))
    b.magnitudes[key]=(name=name,r=r,i=i,lo=lo,hi=hi)
    nothing
end
function _gg_flush_magnitudes!(b)
    # Stable names make row order deterministic even if Dict hashing changes.
    for item in sort!(collect(values(b.magnitudes));by=x->x.name)
        _gg_emit_mag!(b,item.name,item.r,item.i,item.lo,item.hi)
    end
end
function _gg_complex_map(A,r,i)
    ([sum(real(A[k,l])*r[l]-imag(A[k,l])*i[l] for l in axes(A,2)) for k in axes(A,1)],
     [sum(imag(A[k,l])*r[l]+real(A[k,l])*i[l] for l in axes(A,2)) for k in axes(A,1)])
end
_gg_power(r,i,cr,ci) = (r.*cr .+ i.*ci, i.*cr .- r.*ci)
_gg_start(x) = x isa Number ? Float64(x) : JuMP.value(v -> something(JuMP.start_value(v),0.0),x)

function _gg_law!(b,d,er,ei,vb)
    law = d.voltage; n=length(er); scale=d.voltage_scale/vb
    shape=_gg_shape(d)
    lo,hi = _gg_vector(law.magnitude_min,n),_gg_vector(law.magnitude_max,n)
    if law.mode == :fixed_phasor
        for k in 1:n
            _gg_eq!(b,"emf_r_$k",er[k],real(law.phasor[k])/vb;scale)
            _gg_eq!(b,"emf_i_$k",ei[k],imag(law.phasor[k])/vb;scale)
        end
    elseif law.mode==:fixed_magnitudes || shape!==nothing
        ratios = law.mode == :fixed_magnitudes ? law.phasor ./ law.phasor[1] :
            shape
        for k in 2:n
            z=ratios[k]
            _gg_eq!(b,"relative_r_$k",er[k],real(z)*er[1]-imag(z)*ei[1];scale)
            _gg_eq!(b,"relative_i_$k",ei[k],imag(z)*er[1]+real(z)*ei[1];scale)
        end
        if law.mode == :fixed_magnitudes
            mag=abs(law.phasor[1])/vb
            _gg_mag!(b,"emf_magnitude",er[1],ei[1],mag,mag)
        else
            lower,upper=_gg_shape_bounds(d,ratios)
            _gg_mag!(b,"emf_magnitude",er[1],ei[1],lower===nothing ? nothing : lower/vb,upper===nothing ? nothing : upper/vb)
        end
    else
        if law.mode == :phase_magnitudes
            phi=_gg_angles(law,n)
            for k in 2:n
                z=cis(phi[k]-phi[1])
                # Rotate e_k back to the anchor ray; cross=0, dot>=0 rejects π flips.
                rr=(real(z)*er[k]+imag(z)*ei[k])/scale
                ri=(real(z)*ei[k]-imag(z)*er[k])/scale
                _gg_eq!(b,"relative_angle_$k",ri*(er[1]/scale)-rr*(ei[1]/scale))
                _gg_bound!(b,"orientation_$k",rr*(er[1]/scale)+ri*(ei[1]/scale),0.0,nothing)
            end
        end
        for k in 1:n
            _gg_mag!(b,"emf_magnitude_$k",er[k],ei[k],lo[k]===nothing ? nothing : lo[k]/vb,hi[k]===nothing ? nothing : hi[k]/vb)
        end
    end
end

function _gg_power_limits!(b,name,pp,qq,p,q,c,base,scale_si;
        known_p=nothing,known_q=nothing,aggregate_is_port_sum=true)
    # Substitute declared controls only at the same measurement location. This
    # avoids duplicate power equalities/active inequalities at a fixed target.
    # A grounded source's complete PCC power need not equal its port-power sum.
    pp=known_p isa Vector ? known_p./base : pp
    qq=known_q isa Vector ? known_q./base : qq
    p=known_p isa Real ? known_p/base : known_p isa Vector && aggregate_is_port_sum ? sum(pp) : p
    q=known_q isa Real ? known_q/base : known_q isa Vector && aggregate_is_port_sum ? sum(qq) : q
    for (tag,vec,total,lo,hi) in (("p",pp,p,c.p_min,c.p_max),("q",qq,q,c.q_min,c.q_max))
        # Scalars are aggregate; vectors are per port, independently for each side.
        if !(lo isa Vector) && !(hi isa Vector)
            sc=max(abs(something(lo,0.0)),abs(something(hi,0.0)),scale_si)/base
            _gg_bound!(b,name*tag,total,lo===nothing ? nothing : lo/base,hi===nothing ? nothing : hi/base;scale=sc)
        elseif lo isa Vector && hi isa Vector
            for k in eachindex(vec)
                _gg_bound!(b,name*tag*"_$k",vec[k],lo[k]/base,hi[k]/base;
                    scale=max(abs(lo[k]),abs(hi[k]),scale_si)/base)
            end
        else
            for (which,val) in ((:lo,lo),(:hi,hi))
                val===nothing && continue
                if val isa Real
                    _gg_bound!(b,name*tag*"_$(which)_total",total,which==:lo ? val/base : nothing,which==:hi ? val/base : nothing;scale=max(abs(val),scale_si)/base)
                else
                    for k in eachindex(vec)
                        _gg_bound!(b,name*tag*"_$(which)_$k",vec[k],which==:lo ? val[k]/base : nothing,which==:hi ? val[k]/base : nothing;scale=max(abs(val[k]),scale_si)/base)
                    end
                end
            end
        end
    end
    c.s_max===nothing && return
    powers = c.s_max isa Real ? [(p,q,c.s_max)] : collect(zip(pp,qq,c.s_max))
    for (k,(pr,qi,lim)) in enumerate(powers)
        if iszero(lim)
            _gg_eq!(b,name*"S_zero_p_$k",pr); _gg_eq!(b,name*"S_zero_q_$k",qi)
        elseif pr isa Number && qi isa Number
            hypot(pr,qi)*base <= lim+1e-10*max(lim,scale_si) ||
                throw(ArgumentError("fixed power violates apparent-power capability"))
        else
            # Normalize the lifts themselves: their box and circle are O(1).
            m=_gg_model(b)
            pl=JuMP.@variable(m,base_name="$(b.id)_S_p_$k",lower_bound=-1,upper_bound=1)
            ql=JuMP.@variable(m,base_name="$(b.id)_S_q_$k",lower_bound=-1,upper_bound=1)
            JuMP.set_start_value(pl,clamp(_gg_start(pr)*base/lim,-1,1))
            JuMP.set_start_value(ql,clamp(_gg_start(qi)*base/lim,-1,1))
            _gg_eq!(b,name*"S_p_$k",pl,pr*(base/lim))
            _gg_eq!(b,name*"S_q_$k",ql,qi*(base/lim))
            _gg_record!(b,name*"S_circle_$k",JuMP.@constraint(m,pl^2+ql^2 <= 1))
        end
    end
end

function stamp_device!(ctx,d::_GeneralizedDevice; period::Integer=1)
    validate_device(d,(_opf_network(ctx),))
    # Data validation is structural/SI only; net voltages are never used as SI data.
    m=_opf_model(ctx); n=_gg_n(d); ts=_gg_terminals(d); pairs=_gg_connections(d)
    bs=BMOPFTools.opf_coordinate_bases(ctx,d.bus); vb,ib,sb=bs.voltage,bs.current,bs.power
    b=_GGBuilder(ctx,d.id,Dict{String,Any}())
    key=BMOPFTools.OpfModelKey(:expression,:generalized_generator,(d.id,:p))
    key in BMOPFTools.opf_object_keys(ctx) && throw(ArgumentError("generator $(d.id) already stamped"))
    vr=[_opf_voltage(ctx,d.bus,t) for t in ts]
    vi=[_opf_voltage(ctx,d.bus,t;component=:imag) for t in ts]
    C=[Float64(t==p[1])-Float64(t==p[2]) for t in ts,p in pairs]
    ur,ui=C'*vr,C'*vi
    nc=d isa SourceGenerator ? n+1 : n
    ir=JuMP.@variable(m,[1:nc],base_name="$(d.id)_ir")
    ii=JuMP.@variable(m,[1:nc],base_name="$(d.id)_ii")
    for k in 1:n
        v=complex(_gg_start(ur[k]),_gg_start(ui[k]))*vb
        abs(v)>1e-8 || (v=d.voltage_scale*cis(n==3 ? [0,-2pi/3,2pi/3][k] : 0))
        p=d.control.p===nothing ? d.voltage_scale : d.control.p isa Real ? d.control.p/n : d.control.p[k]
        q=d.control.q===nothing ? 0.0 : d.control.q isa Real ? d.control.q/n : d.control.q[k]
        seed=conj(complex(p,q)/v)/ib
        JuMP.set_start_value(ir[k],real(seed)); JuMP.set_start_value(ii[k],imag(seed))
    end
    if d isa SourceGenerator
        JuMP.set_start_value(ir[end],-sum(JuMP.start_value,ir[1:n]))
        JuMP.set_start_value(ii[end],-sum(JuMP.start_value,ii[1:n]))
    end
    Z=_gg_matrix(d.impedance,nc)/bs.impedance
    dr,di=_gg_complex_map(Z,ir,ii)
    if d isa SourceGenerator
        jr,ji=ir,ii
        star_r,star_i=vr[end]+dr[end],vi[end]+di[end]
        er,ei=vr[1:n]+dr[1:n].-star_r,vi[1:n]+di[1:n].-star_i
        gr,gi=-sum(ir),-sum(ii)
        if d.grounding===nothing
            _gg_eq!(b,"earth_open_r",gr); _gg_eq!(b,"earth_open_i",gi)
        else
            zg=ComplexF64(d.grounding)/bs.impedance
            _gg_eq!(b,"ground_drop_r",star_r,real(zg)*gr-imag(zg)*gi;scale=d.voltage_scale/vb)
            _gg_eq!(b,"ground_drop_i",star_i,imag(zg)*gr+real(zg)*gi;scale=d.voltage_scale/vb)
        end
        ground_loss=star_r*gr+star_i*gi
    else
        jr,ji=C*ir,C*ii
        er,ei=ur+dr,ui+di
        gr,gi=0.0,0.0; star_r,star_i=0.0,0.0; ground_loss=0.0
    end
    for k in eachindex(ts)
        BMOPFTools.add_terminal_injection!(ctx,d.bus,ts[k],jr[k],ji[k])
    end
    _gg_law!(b,d,er,ei,vb)
    pp,qq=_gg_power(ur,ui,ir[1:n],ii[1:n])
    pe,qe=_gg_power(er,ei,ir[1:n],ii[1:n])
    pt,qt=_gg_power(vr,vi,jr,ji)
    p,q=sum(pt),sum(qt); pint,qint=sum(pe),sum(qe)
    series_loss=sum(dr.*ir+di.*ii)
    c=d.control
    cp,cq=c.power_location==:poc ? (pp,qq) : (pe,qe)
    cpt,cqt=c.power_location==:poc ? (p,q) : (pint,qint)
    for (name,target,vec,total) in (("P",c.p,cp,cpt),("Q",c.q,cq,cqt))
        target===nothing && continue
        if target isa Real
            _gg_eq!(b,"control_"*name,total,target/sb;scale=max(abs(target),d.voltage_scale)/sb)
        else
            for k in 1:n
                _gg_eq!(b,"control_$(name)_$k",vec[k],target[k]/sb;scale=max(abs(target[k]),d.voltage_scale)/sb)
            end
        end
    end
    cap=d.capability
    xp,xq=cap.power_location==:poc ? (pp,qq) : (pe,qe)
    xt,yt=cap.power_location==:poc ? (p,q) : (pint,qint)
    same_location=c.power_location==cap.power_location
    _gg_power_limits!(b,"cap_",xp,xq,xt,yt,cap,sb,d.voltage_scale;
        known_p=same_location ? c.p : nothing,known_q=same_location ? c.q : nothing,
        aggregate_is_port_sum=d isa GeneralizedGenerator || cap.power_location==:internal)
    for (name,rr,ri,limit) in (("port_i",ir[1:n],ii[1:n],cap.i_max),("terminal_i",jr,ji,cap.terminal_i_max))
        lims=_gg_vector(limit,length(rr))
        for k in eachindex(rr)
            _gg_mag!(b,"$(name)_$k",rr[k],ri[k],nothing,lims[k]===nothing ? nothing : lims[k]/ib)
        end
    end
    _gg_mag!(b,"earth_i",gr,gi,nothing,cap.earth_i_max===nothing ? nothing : cap.earth_i_max/ib)
    vlo,vhi=_gg_vector(cap.voltage_min,n),_gg_vector(cap.voltage_max,n)
    for k in 1:n
        _gg_mag!(b,"poc_voltage_$k",ur[k],ui[k],vlo[k]===nothing ? nothing : vlo[k]/vb,vhi[k]===nothing ? nothing : vhi[k]/vb)
    end
    seqvalid=n==3 && length(unique(last.(pairs)))==1
    seqv=seqvalid ? _gg_complex_map(_gg_transform(),ur,ui) : nothing
    seqe=seqvalid ? _gg_complex_map(_gg_transform(),er,ei) : nothing
    seqi=seqvalid ? _gg_complex_map(_gg_transform(),ir[1:n],ii[1:n]) : nothing
    if cap.sequence!==nothing
        s=cap.sequence; sv=s.location==:poc ? seqv : seqe
        for (name,rr,ri,lo,hi,base) in (("sequence_v",sv[1],sv[2],s.voltage_min,s.voltage_max,vb),
                                      ("sequence_i",seqi[1],seqi[2],s.current_min,s.current_max,ib))
            if name=="sequence_v" && s.location==:internal &&
                    (d.voltage.mode in (:fixed_phasor,:fixed_magnitudes) || _gg_shape(d)!==nothing)
                continue # validated fixed values or merged into the shape magnitude
            end
            ls,hs=_gg_vector(lo,3),_gg_vector(hi,3)
            for k in 1:3
                _gg_mag!(b,"$(name)_$k",rr[k],ri[k],ls[k]===nothing ? nothing : ls[k]/base,hs[k]===nothing ? nothing : hs[k]/base)
            end
        end
    end
    if c.voltage_target!==nothing
        r,i=c.voltage_metric==:positive_sequence ? (seqv[1][2],seqv[2][2]) : (ur[1],ui[1])
        _gg_mag!(b,"voltage_target",r,i,c.voltage_target/vb,c.voltage_target/vb)
    end
    _gg_flush_magnitudes!(b)
    for (tag,expr) in ((:p,p),(:q,q),(:p_internal,pint),(:q_internal,qint),(:loss,series_loss+ground_loss))
        BMOPFTools.register_opf_object!(ctx,BMOPFTools.OpfModelKey(:expression,:generalized_generator,(d.id,tag)),expr)
    end
    (id=d.id,terminals=ts,vr=vr,vi=vi,ur=ur,ui=ui,er=er,ei=ei,
     ir=ir[1:n],ii=ii[1:n],jr=jr,ji=ji,gr=gr,gi=gi,star_r=star_r,star_i=star_i,
     p=p,q=q,p_internal=pint,q_internal=qint,p_phase=pp,q_phase=qq,
     series_loss=series_loss,ground_loss=ground_loss,sequence_voltage=seqv,
     sequence_emf=seqe,sequence_current=seqi,bases=bs,constraints=b.rows,
     cost=d.cost*(sb/1000)*p,internal_voltage_variables=0)
end

"""
    GeneratorResult

SI RMS generator measurements. Complex phasors retain signs and angles;
`p/q` are complete PCC-conductor power, `p_internal/q_internal` are EMF power.
`series_loss/ground_loss` are physical watts. Sequence vectors are ordered 0,1,2
or `nothing` when inapplicable. Unpublished solves return NaN measurements.
`power_balance_error` independently compares power difference and dissipation.
"""
struct GeneratorResult <: AbstractSolveResult
    id::String
    terminals::Vector{String}
    terminal_voltage::Vector{ComplexF64}
    port_voltage::Vector{ComplexF64}
    emf::Vector{ComplexF64}
    port_current::Vector{ComplexF64}
    terminal_current::Vector{ComplexF64}
    earth_current::ComplexF64
    star_voltage::ComplexF64
    p::Float64
    q::Float64
    p_internal::Float64
    q_internal::Float64
    series_loss::Float64
    ground_loss::Float64
    voltage_sequence::Union{Nothing,Vector{ComplexF64}}
    emf_sequence::Union{Nothing,Vector{ComplexF64}}
    current_sequence::Union{Nothing,Vector{ComplexF64}}
    power_balance_error::Float64
    solve::SolveStatus
end
solve_status(r::GeneratorResult)=r.solve
solve_diagnostics(r::GeneratorResult)=(power_balance_error=r.power_balance_error,
    conductor_current_sum=sum(r.terminal_current),earth_current=r.earth_current)

function extract_device(d::_GeneralizedDevice,h,status::SolveStatus)
    val(x)=status.publishable ? (x isa Number ? Float64(x) : JuMP.value(x)) : NaN
    cv(r,i,base)=ComplexF64[complex(val(x),val(y))*base for (x,y) in zip(r,i)]
    seq(s,base)=s===nothing ? nothing : cv(s[1],s[2],base)
    vb,ib,sb=h.bases.voltage,h.bases.current,h.bases.power
    p,q,pe,qe=val(h.p)*sb,val(h.q)*sb,val(h.p_internal)*sb,val(h.q_internal)*sb
    loss,gloss=val(h.series_loss)*sb,val(h.ground_loss)*sb
    GeneratorResult(d.id,h.terminals,cv(h.vr,h.vi,vb),cv(h.ur,h.ui,vb),cv(h.er,h.ei,vb),
        cv(h.ir,h.ii,ib),cv(h.jr,h.ji,ib),complex(val(h.gr),val(h.gi))*ib,
        complex(val(h.star_r),val(h.star_i))*vb,p,q,pe,qe,loss,gloss,
        seq(h.sequence_voltage,vb),seq(h.sequence_emf,vb),seq(h.sequence_current,ib),pe-p-loss-gloss,status)
end

# Stateless snapshots also compose with the existing multi-period device protocol.
link_device!(model,d::_GeneralizedDevice,handles::AbstractVector,sb,grid::TimeGrid)=nothing
_snapshot_device_cost(d::_GeneralizedDevice,h)=h.cost
function _snapshot_device_injection(d::_GeneralizedDevice,h,status::SolveStatus)
    val(x)=status.publishable ? JuMP.value(x)*h.bases.power : NaN
    (p=val(h.p),q=val(h.q))
end
function extract_device(d::_GeneralizedDevice,handles::AbstractVector,::Nothing,sb,status::SolveStatus)
    (snapshots=[extract_device(d,h,status) for h in handles],)
end
