_bkeydict(k::QuantityKey) = Dict{String,Any}("family"=>string(k.family),"id"=>k.id,
    "quantity"=>string(k.quantity),"index"=>collect(k.index),"component"=>string(k.component))
_bkeyorder(k::QuantityKey) = (string(k.family),k.id,string(k.quantity),k.index,string(k.component))
function _bjsonfinite(x)
    x isa AbstractFloat && !isfinite(x) && return nothing
    x isa AbstractDict && return Dict(string(k)=>_bjsonfinite(v) for (k,v) in x)
    (x isa AbstractVector || x isa Tuple) && return [_bjsonfinite(v) for v in x]
    x
end

"""
    bounds_sidecar(report)

Return a versioned, JSON-compatible dictionary of domains, endpoint evidence,
coverage, conflicts and scope. Infinite endpoints are `nothing` (JSON `null`)
with status `unresolved`, never JSON infinities. This is an outer enclosure, not
a certificate that the network is feasible or that unresolved quantities are
unbounded. The result is independently mutable.
"""
function bounds_sidecar(r::BoundReport)
    endpoint(value,id) = Dict("value"=>isfinite(value) ? value : nothing,
        "status"=>id==0 ? "unresolved" : r.evidence[id].rule==:declared ? "declared" : "derived",
        "evidence"=>id==0 ? nothing : id)
    domains=[Dict("key"=>_bkeydict(k),"unit"=>string(r.domains[k].unit),
        "lower"=>endpoint(r.domains[k].lower,r.domains[k].lower_evidence),
        "upper"=>endpoint(r.domains[k].upper,r.domains[k].upper_evidence))
        for k in sort!(collect(keys(r.domains));by=_bkeyorder)]
    proofs=[Dict("id"=>i,"key"=>_bkeydict(p.key),"endpoint"=>string(p.endpoint),
        "value"=>p.value,"rule"=>string(p.rule),"parents"=>copy(p.parents),"sources"=>copy(p.sources))
        for (i,p) in enumerate(r.evidence)]
    Dict{String,Any}("format"=>"PowerOptLab.bounds","version"=>1,
        "scope"=>Dict("snapshot"=>r.scope.snapshot,"problem"=>"native_declared_limit_feasibility",
            "fingerprint"=>r.fingerprint,"adapter"=>_BOUND_ADAPTER_VERSION,"objective_cutoff"=>nothing),
        "status"=>string(r.status),"arithmetic"=>"outward_binary64",
        "domains"=>domains,"evidence"=>proofs,
        "coverage"=>[Dict("object"=>x.object,"status"=>string(x.status),"reason"=>x.reason) for x in r.objects],
        "conflicts"=>_bjsonfinite(r.conflicts),"steps"=>r.steps,"relations"=>r.relation_count,
        "validation"=>deepcopy(r.validation))
end

"""
    materialize_bounds(input, report; mode=:fill_missing)

Return `(network, patch)` without mutating either argument. The network is the
normalized, static snapshot used for propagation; time-series inputs therefore
produce a snapshot, not bounds valid over every time step. A stale input or a
conflict report is rejected. `mode=:intersect` also tightens existing limits.

Currently exports bus terminal/neutral/connection voltage limits and line
terminal current/apparent-power limits. A vector field is emitted only if every
entry has a finite valid bound. Other quantities remain in `bounds_sidecar`.
Original linecode ratings are included when materializing per-line overrides.
The patch records prior values, quantities and endpoint evidence. Output schema
validation is rerun when available.
"""
function materialize_bounds(input,r::BoundReport;mode::Symbol=:fill_missing)
    mode in (:fill_missing,:intersect) || throw(ArgumentError("mode must be :fill_missing or :intersect"))
    r.status==:conflict && throw(ArgumentError("Cannot export a conflicting bound report"))
    net=_binput(input)
    _bfingerprint(net,r.scope)==r.fingerprint || throw(ArgumentError("Stale bound report: the input or scope changed"))
    out=deepcopy(r.network); patch=Dict{String,Any}[]
    inactive=Set{Tuple{String,String}}()
    for (id,b) in get(out,"bus",Dict()), t in get(b,"perfectly_grounded_terminals",String[])
        push!(inactive,(id,string(t)))
    end
    for s in values(get(out,"voltage_source",Dict()))
        bus=string(s["bus"])
        for t in s["terminal_map"]; push!(inactive,(bus,string(t))); end
        nt=_neutral_terminal(out["bus"][bus])
        nt===nothing || push!(inactive,(bus,nt))
    end
    function put!(obj,field,keys,path;lower=false,scalar=false,inherited=nothing)
        mode==:fill_missing && haskey(obj,field) && return
        isempty(keys) && return
        all(k->haskey(r.domains,k),Iterators.flatten(keys)) || return
        # Each entry may cover both line ends: take their hull, not intersection.
        ds=[[r.domains[k] for k in group] for group in keys]
        values=[lower ? minimum(d.lower for d in group) : maximum(d.upper for d in group) for group in ds]
        all(isfinite,values) || return
        # An all-zero lower vector adds no information beyond magnitude >= 0.
        lower && all(iszero,values) && !haskey(obj,field) && return
        old=get(obj,field,nothing)
        effective=old===nothing ? inherited : old
        if effective!==nothing
            previous=scalar ? [effective] : effective
            previous isa AbstractVector && length(previous)==length(values) || return
            values=lower ? max.(Float64.(previous),values) : min.(Float64.(previous),values)
        end
        new=scalar ? only(values) : values
        isequal(old,new) && return
        obj[field]=new
        push!(patch,Dict("path"=>path*"/"*_bp(field),"old"=>deepcopy(old),"new"=>deepcopy(new),
            "inherited"=>deepcopy(inherited),"quantities"=>[[_bkeydict(k) for k in group] for group in keys],
            "evidence"=>[[lower ? d.lower_evidence : d.upper_evidence for d in group] for group in ds]))
    end
    for id in _bsorted(get(out,"bus",Dict()))
        b=out["bus"][id]; path=_bpath("bus",id); nt=_neutral_terminal(b)
        phases=[string(t) for t in b["terminal_names"] if t!=nt]
        volts=[[QuantityKey(:bus,id,:voltage,(t,))] for t in phases]
        # Do not invent combinations of fields around bounds the native model
        # exempts. Such fields may be contradictory in source data yet inactive.
        active=all(t->!((id,t) in inactive),phases)
        if active
            for (field,lower) in (("v_min",true),("v_max",false)); put!(b,field,volts,path;lower); end
        end
        if nt!==nothing
            pairs=[[QuantityKey(:bus,id,:connection_voltage,(t,nt))] for t in phases]
            if active
                for (field,lower) in (("vpn_min",true),("vpn_max",false)); put!(b,field,pairs,path;lower); end
            end
            (id,nt) in inactive || put!(b,"vn_max",[[QuantityKey(:bus,id,:voltage,(nt,))]],path;scalar=true)
        end
        pairs=[[QuantityKey(:bus,id,:connection_voltage,(phases[i],phases[j]))] for i in eachindex(phases) for j in i+1:length(phases)]
        if active
            for (field,lower) in (("vpp_min",true),("vpp_max",false)); put!(b,field,pairs,path;lower); end
        end
    end
    for id in _bsorted(get(out,"line",Dict()))
        l=out["line"][id]; path=_bpath("line",id)
        lc=get(get(out,"linecode",Dict()),get(l,"linecode",""),Dict())
        fm=string.(l["terminal_map_from"]); tm=string.(l["terminal_map_to"])
        for (field,quantity) in (("i_max",:current),("s_max",:power))
            keys=[[QuantityKey(:line,id,quantity,("from",fm[k])),QuantityKey(:line,id,quantity,("to",tm[k]))] for k in eachindex(fm)]
            put!(l,field,keys,path;inherited=get(lc,field,nothing))
        end
    end
    findings=BMOPFTools.Finding[]
    checked=BMOPFTools.schema_check(out,findings)
    get(checked,"jsonschema_valid",true) || throw(ArgumentError("Exported bounds do not satisfy the BMOPF schema"))
    return out,patch
end
