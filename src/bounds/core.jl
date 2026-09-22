"""
    BoundScope(; snapshot=1)

Scope for solver-free BMOPF propagation: a static snapshot, declared operational
limits, native BMOPFTools device semantics, and no objective cutoff. Referenced
time series are resolved for this snapshot. Numerical input is interpreted
as exact binary64 values after BMOPFTools normalization. Unsupported equations
are explicitly relaxed, never replaced with guessed limits.
"""
struct BoundScope
    snapshot::Int
    function BoundScope(; snapshot::Integer=1)
        snapshot >= 1 || throw(ArgumentError("snapshot must be positive"))
        new(Int(snapshot))
    end
end

"""
    QuantityKey(family, id, quantity, index=(); component=:magnitude)

Stable physical quantity identifier. `index` is a tuple of terminal/port labels
(a single string is accepted). Voltage references and current orientations are
specified by the quantity catalogue in the bound-propagation documentation.
"""
struct QuantityKey
    family::Symbol
    id::String
    quantity::Symbol
    index::Tuple{Vararg{String}}
    component::Symbol
end
QuantityKey(family::Symbol,id::AbstractString,quantity::Symbol,index=();component::Symbol=:magnitude) =
    QuantityKey(family,String(id),quantity,index isa AbstractString ? (String(index),) : Tuple(string.(index)),component)
Base.:(==)(a::QuantityKey,b::QuantityKey) =
    a.family==b.family && a.id==b.id && a.quantity==b.quantity && a.index==b.index && a.component==b.component
Base.isequal(a::QuantityKey,b::QuantityKey) = a == b
Base.hash(k::QuantityKey,h::UInt) = hash((k.family,k.id,k.quantity,k.index,k.component),h)

"""An immutable scalar enclosure in SI units, with separate endpoint evidence IDs (0 means unresolved)."""
struct BoundDomain
    lower::Float64
    upper::Float64
    unit::Symbol
    lower_evidence::Int
    upper_evidence::Int
end

"""One outward-valid endpoint derivation. Parent IDs refer to earlier evidence records."""
struct BoundEvidence
    key::QuantityKey
    endpoint::Symbol
    value::Float64
    rule::Symbol
    parents::Vector{Int}
    sources::Vector{String}
end

"""
Result of `propagate_bounds`. Contains SI domains, endpoint evidence, normalized
input, object coverage, conflicts and a stopping status. A nonempty enclosure or
`:fixed_point` is not a feasibility certificate. Treat reports as read-only; use
`bound`, `coverage`, `explain_bound`, and `bounds_sidecar` to inspect them.
"""
struct BoundReport
    scope::BoundScope
    fingerprint::String
    network::Dict{String,Any}
    domains::Dict{QuantityKey,BoundDomain}
    evidence::Vector{BoundEvidence}
    objects::Vector{NamedTuple{(:object,:status,:reason),Tuple{String,Symbol,String}}}
    conflicts::Vector{Dict{String,Any}}
    status::Symbol
    steps::Int
    relation_count::Int
    validation::Dict{String,Any}
end

"""`bound(report, key)` returns the scalar SI enclosure for a registered quantity (unknown keys throw `KeyError`)."""
bound(r::BoundReport,k::QuantityKey) = r.domains[k]

"""Return a copy of the per-object support/relaxation inventory. This is distinct from finite-bound coverage."""
coverage(r::BoundReport) = copy(r.objects)

"""
    explain_bound(report, key; endpoint=:upper)

Return a topologically ordered copy of the endpoint's evidence DAG. An empty
vector means that endpoint remains unresolved; it is not a proof of unboundedness.
"""
function explain_bound(r::BoundReport,k::QuantityKey;endpoint::Symbol=:upper)
    endpoint in (:lower,:upper) || throw(ArgumentError("endpoint must be :lower or :upper"))
    d = bound(r,k)
    id = endpoint == :lower ? d.lower_evidence : d.upper_evidence
    seen = Set{Int}(); todo = [id]
    while !isempty(todo)
        n = pop!(todo)
        (n == 0 || n in seen) && continue
        push!(seen,n); append!(todo,r.evidence[n].parents)
    end
    deepcopy(r.evidence[sort!(collect(seen))])
end

struct _BRelation
    kind::Symbol
    ids::Vector{Int}
    data::Any
    sources::Vector{String}
end
mutable struct _BEngine
    keys::Vector{QuantityKey}
    lookup::Dict{QuantityKey,Int}
    domains::Vector{_BInterval}
    units::Vector{Symbol}
    lower::Vector{Int}
    upper::Vector{Int}
    evidence::Vector{BoundEvidence}
    relations::Vector{_BRelation}
    watchers::Vector{Vector{Int}}
    queue::Vector{Int}
    queued::BitVector
    conflicts::Vector{Dict{String,Any}}
    objects::Vector{NamedTuple{(:object,:status,:reason),Tuple{String,Symbol,String}}}
    atol::Float64
    rtol::Float64
    skipped::Bool
end
_BEngine(atol,rtol) = _BEngine(QuantityKey[],Dict{QuantityKey,Int}(),_BInterval[],Symbol[],Int[],Int[],
    BoundEvidence[],_BRelation[],Vector{Int}[],Int[],BitVector(),Dict{String,Any}[],
    NamedTuple{(:object,:status,:reason),Tuple{String,Symbol,String}}[],atol,rtol,false)

function _bevidence!(e,i,endpoint,value,rule,parents,sources)
    push!(e.evidence,BoundEvidence(e.keys[i],endpoint,value,rule,sort!(unique(filter(!iszero,parents))),copy(sources)))
    length(e.evidence)
end
function _bscalar!(e,k,unit;nonnegative=false)
    haskey(e.lookup,k) && return e.lookup[k]
    push!(e.keys,k); i=length(e.keys); e.lookup[k]=i
    push!(e.domains,nonnegative ? _BNONNEG : _BENTIRE)
    push!(e.units,unit); push!(e.lower,0); push!(e.upper,0); push!(e.watchers,Int[])
    nonnegative && (e.lower[i]=_bevidence!(e,i,:lower,0.0,:nonnegative,Int[],String[]))
    i
end
function _brelation!(e,kind,ids,data,sources)
    push!(e.relations,_BRelation(kind,collect(ids),data,collect(sources)))
    n=length(e.relations); push!(e.queued,false)
    for i in unique(ids); push!(e.watchers[i],n); end
end
function _bcomplex!(e,family,id,quantity,index,unit)
    power = quantity == :power
    r = _bscalar!(e,QuantityKey(family,id,quantity,index;component=power ? :active : :real),power ? :W : unit)
    i = _bscalar!(e,QuantityKey(family,id,quantity,index;component=power ? :reactive : :imag),power ? :var : unit)
    m = _bscalar!(e,QuantityKey(family,id,quantity,index),power ? :VA : unit;nonnegative=true)
    _brelation!(e,:norm,[r,i,m],nothing,String[])
    (r,i,m)
end
function _bupdate!(e,i,c::_BInterval,rule,parents,sources;force=false)
    isempty(e.conflicts) || return
    (isnan(c.lo) || isnan(c.hi)) && error("NaN from bound rule $rule")
    old=e.domains[i]; lo=max(old.lo,c.lo); hi=min(old.hi,c.hi)
    if lo > hi
        push!(e.conflicts,Dict("quantity"=>i,"rule"=>string(rule),"existing"=>[old.lo,old.hi],
            "candidate"=>[c.lo,c.hi],"parents"=>sort!(unique(filter(!iszero,vcat(parents,e.lower[i],e.upper[i])))),
            "sources"=>copy(sources)))
        return
    end
    if !force && isfinite(old.lo) && lo > old.lo && lo-old.lo <= e.atol+e.rtol*max(abs(lo),abs(old.lo))
        lo=old.lo; e.skipped=true
    end
    if !force && isfinite(old.hi) && hi < old.hi && old.hi-hi <= e.atol+e.rtol*max(abs(hi),abs(old.hi))
        hi=old.hi; e.skipped=true
    end
    lo == old.lo && hi == old.hi && return
    lo > old.lo && (e.lower[i]=_bevidence!(e,i,:lower,lo,rule,parents,sources))
    hi < old.hi && (e.upper[i]=_bevidence!(e,i,:upper,hi,rule,parents,sources))
    e.domains[i]=_BInterval(lo,hi)
    for j in e.watchers[i]
        if !e.queued[j]; push!(e.queue,j); e.queued[j]=true; end
    end
end

# A complex row sum(c[k]*z[k]) == 0 is lowered into two real rows plus a
# magnitude triangle contractor. Combine duplicate coordinates before lowering.
function _brow!(e,terms,sources)
    cs=Dict{NTuple{3,Int},Complex{_BRational}}()
    for (z,c) in terms
        cc=complex(_brat(real(c)),_brat(imag(c)))
        cs[z]=get(cs,z,zero(cc))+cc
    end
    zs=sort!([z for (z,c) in cs if !iszero(c)])
    isempty(zs) && return
    rr=Dict{Int,_BRational}(); ii=Dict{Int,_BRational}()
    add!(d,k,c) = (d[k]=get(d,k,_brat(0))+c)
    for z in zs
        c=cs[z]
        add!(rr,z[1],real(c)); add!(rr,z[2],-imag(c))
        add!(ii,z[1],imag(c)); add!(ii,z[2],real(c))
    end
    for d in (rr,ii)
        ids=sort!([k for (k,c) in d if !iszero(c)])
        isempty(ids) || _brelation!(e,:affine,ids,[d[k] for k in ids],sources)
    end
    coeffs=[_bsqrt(_BInterval(_bround(abs2(cs[z]),false),_bround(abs2(cs[z]),true))) for z in zs]
    _brelation!(e,:triangle,[z[3] for z in zs],coeffs,sources)
end

function _bcontract!(e,r::_BRelation)
    d=e.domains
    # Snapshot parents before updates: all premises precede every emitted node.
    parents=sort!(unique(filter(!iszero,vcat(e.lower[r.ids],e.upper[r.ids]))))
    # Snapshot intervals too; evidence must describe exactly the intervals used.
    x=copy(d[r.ids])
    update(j,c) = _bupdate!(e,r.ids[j],c,r.kind,parents,r.sources)
    if r.kind == :affine
        c=r.data
        for k in eachindex(x)
            rest=_bsum(_bscale(c[j],x[j]) for j in eachindex(x) if j!=k)
            update(k,_bscale(-inv(c[k]),rest))
        end
    elseif r.kind == :norm
        a,b,m=x
        update(3,_bhypot(a,b))
        for (k,other) in ((1,b),(2,a))
            residual=_bsub(_bsquare(_bpoint(m.hi)),_bsquare(_bpoint(_bminabs(other))))
            if residual.hi < 0
                update(3,_BInterval(Inf,Inf)); return
            end
            u=_bsqrt(residual).hi
            update(k,_BInterval(-u,u))
        end
    elseif r.kind == :triangle
        a=r.data; n=length(x)
        weighted=[_bmul(a[j],x[j]) for j in 1:n]
        for k in 1:n
            others=[j for j in 1:n if j!=k]
            u=_bsum(weighted[j] for j in others).hi
            l=0.0
            for j in others
                rest=_bsum(weighted[t] for t in others if t!=j)
                l=max(l,_bsub(_bpoint(weighted[j].lo),rest).lo)
            end
            a[k].lo > 0 && update(k,_bdiv(_BInterval(l,u),a[k]))
        end
    elseif r.kind == :power
        vr,vi,vm,ir,ii,im,p,q,sm=x
        update(7,_badd(_bmul(vr,ir),_bmul(vi,ii)))
        update(8,_bsub(_bmul(vi,ir),_bmul(vr,ii)))
        update(9,_bmul(vm,im))
        vm.lo > 0 && update(6,_bdiv(sm,vm))
        im.lo > 0 && update(3,_bdiv(sm,im))
        if vm.lo > 0
            denom=_bsquare(vm)
            update(4,_bdiv(_badd(_bmul(p,vr),_bmul(q,vi)),denom))
            update(5,_bdiv(_bsub(_bmul(p,vi),_bmul(q,vr)),denom))
        end
        if im.lo > 0
            denom=_bsquare(im)
            update(1,_bdiv(_bsub(_bmul(p,ir),_bmul(q,ii)),denom))
            update(2,_bdiv(_badd(_bmul(p,ii),_bmul(q,ir)),denom))
        end
    else
        error("Unknown bound contractor $(r.kind)")
    end
end

function _brun!(e;max_steps,time_limit)
    empty!(e.queue); fill!(e.queued,true); append!(e.queue,eachindex(e.relations))
    head=1; steps=0; started=time_ns()
    while head<=length(e.queue) && isempty(e.conflicts)
        (steps>=max_steps || (time_ns()-started)/1e9>=time_limit) && return (:budget_exhausted,steps)
        j=e.queue[head]; head+=1; e.queued[j]=false
        _bcontract!(e,e.relations[j]); steps+=1
    end
    status=isempty(e.conflicts) ? (e.skipped ? :tolerance_stopped : :fixed_point) : :conflict
    (status,steps)
end
